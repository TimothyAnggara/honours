# ============================================================
# Helper Functions
# ============================================================

run_spicy_contrast = function(spe, cond_A, cond_B,
                              level = c("sub", "major", "lineage"),
                              cores = 4) {
  level    = match.arg(level)
  cellCol  = switch(level,
                    sub     = "cellType_sub",
                    major   = "cellType_major",
                    lineage = "lineage"
  )
  
  spe2 = spe[, spe$condition %in% c(cond_A, cond_B)]
  
  spe2$condition    = droplevels(factor(spe2$condition))
  spe2$imageID      = droplevels(factor(spe2$imageID))
  spe2[[cellCol]]   = droplevels(factor(spe2[[cellCol]]))
  
  spicy(
    spe2,
    condition = "condition",
    subject   = "subject",
    cellType  = cellCol,
    imageID   = "imageID",
    cores     = cores,
    weights   = FALSE
  )
}

canonical_pair = function(x) {
  ab = strsplit(x, "__", fixed = TRUE)
  vapply(ab, function(z) paste(sort(z), collapse = "__"), character(1))
}

extract_leaf_stats = function(spicy_sub, conditionB = "signal", eps = 1e-12) {
  colB = paste0("condition", conditionB)
  tibble(
    interaction = rownames(spicy_sub$p.value),
    leaf        = canonical_pair(rownames(spicy_sub$p.value)),
    pval        = as.numeric(spicy_sub$p.value[, colB]),
    coef        = as.numeric(spicy_sub$coefficient[, colB])
  ) %>%
    # spicyR returns both A__B and B__A; keep only the alphabetically sorted
    # form so each unordered pair appears exactly once.
    filter(interaction == leaf) %>%
    mutate(pval = pmin(pmax(pval, eps), 1 - eps))
}

extract_major_pvals = function(spicy_major, conditionB = "signal", eps = 1e-12) {
  colB = paste0("condition", conditionB)
  tibble(
    fam     = canonical_pair(rownames(spicy_major$p.value)),
    p_major = as.numeric(spicy_major$p.value[, colB])
  ) %>%
    distinct(fam, .keep_all = TRUE) %>%
    mutate(p_major = pmin(pmax(p_major, eps), 1 - eps))
}

# ── Flat BH ─────────────────────
bh_reject = function(p, q = 0.05) {
  p   = as.numeric(p)
  m   = length(p)
  o   = order(p)
  ps  = p[o]
  thr = (seq_len(m) / m) * q
  k   = max(c(0L, which(ps <= thr)))
  out = rep(FALSE, m)
  if (k > 0) out[o[seq_len(k)]] = TRUE
  out
}

# ============================================================
# TreeBH
# ============================================================
build_treebh_inputs = function(spicy_sub, conditionB = "signal",
                               major_map = MAJOR_MAP, eps = 1e-12) {
  dfLeaf = extract_leaf_stats(spicy_sub, conditionB, eps) %>%
    tidyr::separate(leaf, into = c("a", "b"), sep = "__", remove = FALSE) %>%
    mutate(
      ma  = unname(major_map[a]),
      mb  = unname(major_map[b]),
      fam = canonical_pair(paste(ma, mb, sep = "__"))
    ) %>%
    select(leaf, fam, pval, coef)
  
  groups = cbind(
    level1 = as.integer(factor(dfLeaf$fam)),
    level2 = as.integer(factor(paste(dfLeaf$fam, dfLeaf$leaf, sep = "/")))
  )
  
  list(df_leaf = dfLeaf, pvals = dfLeaf$pval, groups = groups)
}

run_treebh = function(spicy_sub, conditionB = "signal",
                      q = c(0.05, 0.05),
                      major_map = MAJOR_MAP, eps = 1e-12) {
  inp = build_treebh_inputs(spicy_sub, conditionB, major_map, eps)
  
  res = TreeBH::get_TreeBH_selections(
    pvals  = inp$pvals,
    groups = inp$groups,
    q      = q
  )
  
  inp$df_leaf %>%
    mutate(method = "TreeBH", reject = (res[, 2] == 1))
}

# ============================================================
# DART
# ============================================================
extract_leaf_for_dart = function(spicy_sub, conditionB = "signal", eps = 1e-12) {
  df = extract_leaf_stats(spicy_sub, conditionB, eps) %>%
    arrange(leaf) %>%
    filter(!is.na(pval), !is.na(coef))  # drop pairs with NA coefs/pvals
  
  list(
    df    = df,
    pvals = setNames(df$pval, df$leaf),
    coefs = setNames(df$coef, df$leaf)
  )
}

make_hypothesis_distance = function(leaf_names, major_map = MAJOR_MAP) {
  cell_dist = function(x, y) {
    if (x == y)                           return(0L)
    if (major_map[[x]] == major_map[[y]]) return(1L)
    2L
  }
  parse_pair = function(s) {
    ab = sort(strsplit(s, "__", fixed = TRUE)[[1]])
    list(a = ab[1], b = ab[2])
  }
  
  pairs = lapply(leaf_names, parse_pair)
  m     = length(leaf_names)
  D     = matrix(0L, m, m, dimnames = list(leaf_names, leaf_names))
  
  for (i in seq_len(m)) {
    ai = pairs[[i]]$a; bi = pairs[[i]]$b
    for (j in seq_len(m)) {
      aj = pairs[[j]]$a; bj = pairs[[j]]$b
      D[i, j] = min(
        cell_dist(ai, aj) + cell_dist(bi, bj),
        cell_dist(ai, bj) + cell_dist(bi, aj)
      )
    }
  }
  D
}

run_dart_one = function(spicy_sub, spe, conditionB = "signal",
                        alpha = 0.05, eps = 1e-12, Mgroup = 2,
                        subject_col = "subject") {
  
  leaf      = extract_leaf_for_dart(spicy_sub, conditionB, eps)
  pvals     = leaf$pvals
  coefs     = leaf$coefs
  leafNames = names(pvals)
  
  cat(">>> [DART] n_hypotheses (m):", length(pvals), "\n")
  cat(">>> [DART] any NA pvals:", any(is.na(pvals)), "\n")
  cat(">>> [DART] any NA coefs:", any(is.na(coefs)), "\n")
  cat(">>> [DART] pval range: [", min(pvals), ",", max(pvals), "]\n")
  
  D    = make_hypothesis_distance(leafNames)
  maxD = max(D)
  D2   = if (maxD > 0) D / maxD else D
  
  cat(">>> [DART] maxD:", maxD, "\n")
  cat(">>> [DART] any NA in D2:", any(is.na(D2)), "\n")
  
  z  = qnorm(1 - as.numeric(pvals) / 2) * sign(as.numeric(coefs))
  T1 = matrix(z, nrow = 1)
  m  = length(pvals)
  
  cat(">>> [DART] any NA/Inf in z:", any(is.na(z) | is.infinite(z)), "\n")
  cat(">>> [DART] z range: [", min(z, na.rm = TRUE), ",", max(z, na.rm = TRUE), "]\n")
  
  n = length(unique(spe[[subject_col]]))
  cat(">>> [DART] n_subjects (n):", n, "\n")
  cat(">>> [DART] subject_col used:", subject_col, "\n")
  
  loglogm = log(log(m))
  grids   = 16 / sqrt(n * log(m) * loglogm)
  cat(">>> [DART] log(m):", log(m), "\n")
  cat(">>> [DART] log(log(m)):", loglogm, "\n")
  cat(">>> [DART] grids (before clamp):", grids, "\n")
  grids = max(grids, 0.5)
  cat(">>> [DART] grids (after clamp):", grids, "\n")
  
  cat(">>> [DART] calling A.tree.mult...\n")
  Atree = DART::A.tree.mult(grids = grids, Dist0 = D2, Mgroup = Mgroup)
  cat(">>> [DART] A.tree.mult done\n")
  
  res = DART::test.mult(alpha = alpha, Llist = Atree$Llist,
                        Dist0 = D2, T1 = T1)
  cat(">>> [DART] test.mult done\n")
  
  dartIdx    = res[[length(res)]]
  dartReject = rep(FALSE, m)
  if (length(dartIdx) > 0) dartReject[dartIdx] = TRUE
  
  cat(">>> [DART] n_rejected:", sum(dartReject), "\n")
  
  tibble(leaf = leafNames, pval = as.numeric(pvals),
         coef = as.numeric(coefs), method = "DART", reject = dartReject)
}

# ============================================================
# p-filter
# ============================================================
extract_pfilter_inputs = function(spicy_sub, spicy_major,
                                  conditionB = "signal", eps = 1e-12) {
  colB = paste0("condition", conditionB)
  stopifnot(colB %in% colnames(spicy_sub$p.value))
  stopifnot(colB %in% colnames(spicy_major$p.value))
  
  majorDf = tibble(
    fam     = canonical_pair(rownames(spicy_major$p.value)),
    p_major = as.numeric(spicy_major$p.value[, colB])
  ) %>%
    group_by(fam) %>% dplyr::slice(1) %>% ungroup() %>%
    filter(!is.na(p_major)) %>%                              # drop NA pvals
    mutate(p_major = pmin(pmax(p_major, eps), 1 - eps)) %>%
    arrange(fam)
  
  subDf = tibble(
    leaf  = canonical_pair(rownames(spicy_sub$p.value)),
    p_sub = as.numeric(spicy_sub$p.value[, colB])
  ) %>%
    group_by(leaf) %>% dplyr::slice(1) %>% ungroup() %>%
    filter(!is.na(p_sub)) %>%                                # drop NA pvals
    tidyr::separate(leaf, into = c("a", "b"), sep = "__", remove = FALSE) %>%
    mutate(
      fam   = canonical_pair(paste0(MAJOR_MAP[a], "__", MAJOR_MAP[b])),
      p_sub = pmin(pmax(p_sub, eps), 1 - eps)
    ) %>%
    select(leaf, fam, p_sub) %>%
    filter(fam %in% majorDf$fam) %>%                        # only keep subtypes whose family survived
    arrange(fam, leaf)
  
  list(major_df = majorDf, sub_df = subDf)
}

run_pfilter_one = function(spicy_sub, spicy_major,
                           conditionB  = "signal",
                           alpha_major = 0.05,
                           alpha_leaf  = 0.05,
                           eps         = 1e-12) {
  inp     = extract_pfilter_inputs(spicy_sub, spicy_major, conditionB, eps)
  majorDf = inp$major_df
  subDf   = inp$sub_df
  
  # Drop subtypes whose family couldn't be matched — avoids NA group indices
  parentIdx = match(subDf$fam, majorDf$fam)
  subDf     = subDf[!is.na(parentIdx), ]
  parentIdx = parentIdx[!is.na(parentIdx)]
  
  nMajor = nrow(majorDf)
  nSub   = nrow(subDf)
  nTotal = nMajor + nSub
  
  pAll   = c(majorDf$p_major, subDf$p_sub)
  groups = matrix(0L, nrow = nTotal, ncol = 2)
  
  groups[seq_len(nMajor), 1]      = seq_len(nMajor)
  groups[(nMajor + 1):nTotal, 1] = parentIdx
  groups[, 2]                     = seq_len(nTotal)
  
  res        = pfilter(P = pAll, alphas = c(alpha_major, alpha_leaf),
                       groups = groups)
  leafReject = res[(nMajor + 1):nTotal] == 1
  
  tibble(leaf = subDf$leaf, fam = subDf$fam,
         pval = subDf$p_sub, method = "pfilter", reject = leafReject)
}

# ============================================================
# Our Method (three-level hierarchy: lineage -> major -> subtype)
# ============================================================

# Build a three-level tree data frame from spicyR results at all three levels.
#
# Level 1: lineage pairs  (e.g. "Lymphoid__Myeloid") — root nodes, no parent
# Level 2: major pairs    (e.g. "NK__T")             — parent is a lineage pair
# Level 3: subtype pairs  (e.g. "CD4_T__NK")         — parent is a major pair
#
# Each node stores its p-value from the corresponding spicyR run.
make_tree_df_from_spicy = function(spicy_lineage, spicy_major, spicy_sub,
                                   conditionB = "signal", eps = 1e-12) {
  colB = paste0("condition", conditionB)
  stopifnot(colB %in% colnames(spicy_lineage$p.value))
  stopifnot(colB %in% colnames(spicy_major$p.value))
  stopifnot(colB %in% colnames(spicy_sub$p.value))
  
  # ── Level 1: lineage pairs ───────────────────────────────────────────
  lineageDf = tibble(
    node   = canonical_pair(rownames(spicy_lineage$p.value)),
    level  = 1L,
    parent = NA_character_,
    pvalue = as.numeric(spicy_lineage$p.value[, colB])
  ) %>%
    group_by(node) %>% dplyr::slice(1) %>% ungroup() %>%
    mutate(pvalue = pmin(pmax(pvalue, eps), 1 - eps))
  
  # ── Level 2: major pairs ─────────────────────────────────────────────
  # Parent of a major pair is determined by mapping each major type to its
  # lineage via MAJOR_LINEAGE_MAP, then taking the canonical lineage pair.
  majorDf = tibble(
    node   = canonical_pair(rownames(spicy_major$p.value)),
    level  = 2L,
    pvalue = as.numeric(spicy_major$p.value[, colB])
  ) %>%
    group_by(node) %>% dplyr::slice(1) %>% ungroup() %>%
    tidyr::separate(node, into = c("a", "b"), sep = "__", remove = FALSE) %>%
    mutate(
      parent = canonical_pair(paste0(
        MAJOR_LINEAGE_MAP[a], "__", MAJOR_LINEAGE_MAP[b]
      )),
      pvalue = pmin(pmax(pvalue, eps), 1 - eps)
    ) %>%
    select(node, level, parent, pvalue) %>%
    filter(parent %in% lineageDf$node)
  
  # ── Level 3: subtype pairs ───────────────────────────────────────────
  # Parent of a subtype pair is determined by mapping each subtype to its
  # major type via MAJOR_MAP, then taking the canonical major pair.
  subtypeDf = tibble(
    node   = canonical_pair(rownames(spicy_sub$p.value)),
    level  = 3L,
    pvalue = as.numeric(spicy_sub$p.value[, colB])
  ) %>%
    group_by(node) %>% dplyr::slice(1) %>% ungroup() %>%
    tidyr::separate(node, into = c("a", "b"), sep = "__", remove = FALSE) %>%
    mutate(
      parent = canonical_pair(paste0(MAJOR_MAP[a], "__", MAJOR_MAP[b])),
      pvalue = pmin(pmax(pvalue, eps), 1 - eps)
    ) %>%
    select(node, level, parent, pvalue) %>%
    filter(parent %in% majorDf$node)
  
  bind_rows(lineageDf, majorDf, subtypeDf) %>%
    arrange(level, parent, node)
}

# Gate through three levels: L1 (lineage) -> L2 (major) -> L3 (subtype).
# At each level, BH is applied only over nodes whose parent was rejected.
# L1 nodes have no parent so all are tested.
run_levelwise_bh = function(tree_df, q = 0.05, l1_method = c("BY", "BH")) {
  l1_method = match.arg(l1_method)
  
  # ── Level 1: lineage pairs — gate with BY (arbitrary dependence) ─────
  L1 = tree_df %>%
    filter(level == 1) %>%
    mutate(
      tested = TRUE,
      reject = if (l1_method == "BY") {
        p.adjust(pvalue, method = "BY") <= q
      } else {
        bh_reject(pvalue, q)
      }
    )
  
  sigLineage = L1$node[L1$reject]
  
  # ── Level 2: major pairs — test only under rejected lineage parents ──
  L2 = tree_df %>%
    filter(level == 2) %>%
    mutate(tested = parent %in% sigLineage, reject = FALSE)
  
  if (any(L2$tested)) {
    idx = which(L2$tested)
    L2$reject[idx] = bh_reject(L2$pvalue[idx], q)
  }
  
  sigMajor = L2$node[L2$reject]
  
  # ── Level 3: subtype pairs — test only under rejected major parents ──
  L3 = tree_df %>%
    filter(level == 3) %>%
    mutate(tested = parent %in% sigMajor, reject = FALSE)
  
  if (any(L3$tested)) {
    idx = which(L3$tested)
    L3$reject[idx] = bh_reject(L3$pvalue[idx], q)
  }
  
  bind_rows(L1, L2, L3) %>% arrange(level, parent, node)
}

run_levelwise_bh_twolevel = function(major_df, sub_df, q = 0.05, l1_method = c("BY", "BH")) {
  l1_method = match.arg(l1_method)
  
  # Level 1: major pairs — test all
  L1 = major_df %>%
    mutate(
      tested = TRUE,
      reject = if (l1_method == "BY") {
        p.adjust(pvalue, method = "BY") <= q
      } else {
        bh_reject(pvalue, q)
      }
    )
  
  sigMajor = L1$node[L1$reject]
  
  # Level 2: subtype pairs — test only under rejected major parents
  L2 = sub_df %>%
    mutate(tested = parent %in% sigMajor, reject = FALSE)
  
  if (any(L2$tested)) {
    idx = which(L2$tested)
    L2$reject[idx] = bh_reject(L2$pvalue[idx], q)
  }
  
  bind_rows(L1, L2) %>% arrange(level, parent, node)
}

run_ours_one = function(spicy_lineage, spicy_major, spicy_sub,
                        conditionB = "signal", q = 0.05,
                        l1_method = c("BY", "BH")) {
  l1_method = match.arg(l1_method)
  treeDF    = make_tree_df_from_spicy(spicy_lineage, spicy_major, spicy_sub,
                                      conditionB)
  
  majorDf  = treeDF %>% filter(level == 2)  # major pairs
  subtypeDf = treeDF %>% filter(level == 3) # subtype pairs, already has parent column
  
  resOurs = run_levelwise_bh_twolevel(majorDf, subtypeDf, q = q, l1_method = l1_method)
  list(tree_df = treeDF, res_ours = resOurs)
}

# ============================================================
# ClusterTrim (BY) — corrected implementation v2
# Reference: Benjamini & Heller (2007), Procedure 4.1
# ============================================================

# --------------------------------------------------
# cluster_test_major()
#
# Runs BH or BY at the major-family level and returns:
#   sig_fams : names of significant families
#   u1       : largest rejected raw p-value (the
#              "cut-off point" from Procedure 4.1 step 1d)
#   m0_m     : estimated proportion of null families,
#              from Procedure 4.1 step 2a:
#              m0_hat = (m - k) / (1 - q), so
#              m0_m   = m0_hat / m = (m - k) / (m * (1 - q))
# --------------------------------------------------
cluster_test_major = function(p_major, q = 0.05, method = c("BY", "BH")) {
  method = match.arg(method)
  adj    = p.adjust(p_major, method = method)
  sig    = adj <= q
  m      = length(p_major)
  k      = sum(sig)
  
  u1   = if (k > 0) max(p_major[sig]) else NA_real_  # restore this
  m0_m = if (k > 0) min((m - k) / (m * (1 - q)), 1) else 1
  
  list(sig_fams = names(p_major)[sig], u1 = u1, m0_m = m0_m, adj = adj)
}

# --------------------------------------------------
# calculate_trimming_p()
#
# Implements equation (6) from Benjamini & Heller (2007)
# under the conservative µi -> 0 approximation.
#
# Full equation (6):
#
#   pli = [m0_m * Φ̃((Φ̃⁻¹(u1) - ρ·zli) / √(1-ρ²))]
#         / [m0_m·u1 + (1-m0_m)·Φ̃(Φ̃⁻¹(u1) - µi)]
#
# Conservative approximation (µi -> 0):
#   Φ̃(Φ̃⁻¹(u1) - 0) = Φ̃(Φ̃⁻¹(u1)) = u1
#   => denominator = m0_m·u1 + (1-m0_m)·u1 = u1
#
# So the conservative trimming p-value is:
#
#   pli = m0_m * Φ̃((Φ̃⁻¹(u1) - ρ·zli) / √(1-ρ²)) / u1
#
# where Φ̃ = pnorm(..., lower.tail = FALSE)
#       Φ̃⁻¹ = qnorm(..., lower.tail = FALSE)
#
# Note: m0_m is now passed in from cluster_test_major()
# rather than hardcoded, so it is estimated from data.
# --------------------------------------------------
calculate_trimming_p = function(z_subtype, u1_cutoff, rho, m0_m) {
  # Φ̃⁻¹(u1): right-tail quantile
  threshold = qnorm(u1_cutoff, lower.tail = FALSE)
  
  # Numerator: m0_m * Φ̃((threshold - ρ·z) / √(1-ρ²))
  num = m0_m * pnorm(
    (threshold - rho * z_subtype) / sqrt(1 - rho^2),
    lower.tail = FALSE
  )
  
  # Denominator: u1 under µi -> 0 approximation
  den = u1_cutoff
  
  pmin(pmax(num / den, 0), 1)
}

# --------------------------------------------------
# run_cluster_trim_one()
#
# Adaptation of Procedure 4.1 from Benjamini & Heller
# (2007) to cell-type interaction pairs.
#
# Mapping from paper to this setting:
#   clusters    <-> major-family interaction groups
#   locations   <-> subtype interaction pairs within families
#   FDRCLUST    <-> FDR on major families
#   FDRLOC      <-> FDR on subtype pairs
#
# Departures from the original paper (documented):
#   1. Two-sided tests: z_subtype uses sign(coef) to
#      handle bidirectional interactions (co-localisation
#      vs avoidance). The paper assumes one-sided tests.
#   2. rho is fixed rather than estimated per family from
#      a variogram. A single pooled value is used as an
#      approximation since we lack the spatial coordinates
#      needed for variogram estimation at this level.
#   3. BH (not two-stage BH) is used in the trimming stage
#      as a tractable approximation.
# --------------------------------------------------
run_cluster_trim_one = function(spicy_sub, spicy_major,
                                conditionB   = "signal",
                                q            = 0.05,
                                major_method = c("BY", "BH"),
                                rho          = 0.5,
                                eps          = 1e-12) {
  major_method = match.arg(major_method)
  inp     = extract_pfilter_inputs(spicy_sub, spicy_major, conditionB, eps)
  majorDf = inp$major_df
  subDf   = inp$sub_df
  
  # Pull coef for sign information and join onto subDf
  leafStats = extract_leaf_stats(spicy_sub, conditionB, eps) %>%
    select(leaf, coef)
  subDf = subDf %>% left_join(leafStats, by = "leaf")
  
  pMajor = setNames(majorDf$p_major, majorDf$fam)
  maj    = cluster_test_major(pMajor, q = q, method = major_method)
  
  # Early exit: no families passed the major-level correction
  if (length(maj$sig_fams) == 0) {
    return(subDf %>% transmute(leaf, pval = p_sub,
                               method = "ClusterTrim_BY", reject = FALSE))
  }
  
  subDf = subDf %>%
    mutate(
      in_sig_family = fam %in% maj$sig_fams,
      # Signed z-score: converts two-sided p-value to directional z.
      # zli = Φ̃⁻¹(p̃li/2) * sign(coef) to distinguish co-localisation
      # from avoidance. Adaptation from the paper's one-sided setting.
      z_subtype = qnorm(1 - p_sub / 2) * sign(coef),
      p_trimmed = ifelse(
        in_sig_family,
        vapply(z_subtype, calculate_trimming_p, numeric(1),
               u1_cutoff = maj$u1,
               rho       = rho,
               m0_m      = maj$m0_m),
        NA_real_
      )
    )
  
  # Trimming stage: apply BH to trimming p-values within significant
  # families. This approximates the two-stage procedure from step 2c
  # of Procedure 4.1 in the paper.
  inSig     = which(subDf$in_sig_family & !is.na(subDf$p_trimmed))
  rejectVec = rep(FALSE, nrow(subDf))
  if (length(inSig) > 0) {
    rejectVec[inSig] = bh_reject(subDf$p_trimmed[inSig], q = q)
  }
  
  subDf %>%
    mutate(reject = rejectVec) %>%
    transmute(leaf, pval = p_sub, method = "ClusterTrim_BY", reject)
}

#' Run all methods on spicy results at all three levels.
#'
#' @param spicy_sub     spicyR result at the subtype level.
#' @param spicy_major   spicyR result at the major type level.
#' @param spicy_lineage spicyR result at the lineage level.
#' @param spe           SpatialExperiment object (needed for DART).
#' @param true_pairs    Character vector of ground-truth canonical pair names.
#'
#' @return A list with two elements:
#'   \item{raw}{Tidy tibble, one row per method x leaf.}
#'   \item{metrics}{Tidy tibble, one row per method with TPR/FDR/etc.}
run_all_methods = function(spicy_sub, spicy_major, spicy_lineage, spe,
                           conditionB = "signal",
                           true_pairs = TRUE_PAIRS) {
  
  # ── Flat BH ──────────────────────────────────────────────────────────
  flatBh = extract_leaf_stats(spicy_sub, conditionB) %>%
    mutate(method = "Flat BH", reject = bh_reject(pval, q = 0.05))
  
  # ── Flat BY ──────────────────────────────────────────────────────────
  flatBy = extract_leaf_stats(spicy_sub, conditionB) %>%
    mutate(method = "Flat BY",
           reject = p.adjust(pval, method = "BY") <= 0.05)
  
  # ── TreeBH ───────────────────────────────────────────────────────────
  treebh = tryCatch(
    run_treebh(spicy_sub, conditionB = conditionB),
    error = function(e) {
      message("TreeBH failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>% mutate(method = "TreeBH", reject = FALSE)
    }
  )
  
  # ── DART ─────────────────────────────────────────────────────────────
  dart = tryCatch(
    run_dart_one(spicy_sub, spe, conditionB = conditionB),
    error = function(e) {
      message("DART failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>% mutate(method = "DART", reject = FALSE)
    }
  )
  
  # ── p-filter ─────────────────────────────────────────────────────────
  pfilterRes = tryCatch(
    run_pfilter_one(spicy_sub, spicy_major, conditionB = conditionB),
    error = function(e) {
      message("pfilter failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>%
        mutate(method = "pfilter", reject = FALSE)
    }
  )
  
  # ── Ours (three-level: lineage -> major -> subtype) ──────────────────
  oursRes = tryCatch({
    out = run_ours_one(spicy_lineage, spicy_major, spicy_sub,
                       conditionB = conditionB, l1_method = "BY")
    out$res_ours %>%
      filter(level == 3) %>%
      transmute(leaf = node, pval = pvalue, method = "Ours", reject)
  }, error = function(e) {
    message("Ours failed: ", conditionMessage(e))
    extract_leaf_stats(spicy_sub, conditionB) %>% mutate(method = "Ours", reject = FALSE)
  })
  
  # ── ClusterTrim (BY) ─────────────────────────────────────────────────
  clustertrim = tryCatch(
    run_cluster_trim_one(spicy_sub, spicy_major, conditionB = conditionB, major_method = "BH"),
    error = function(e) {
      message("ClusterTrim_BY failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>%
        mutate(method = "ClusterTrim_BY", reject = FALSE)
    }
  )
  
  # ── Combine ──────────────────────────────────────────────────────────
  allRaw = list(flatBh, flatBy, treebh, dart, pfilterRes, oursRes, clustertrim) %>%
    purrr::map_dfr(function(df) {
      df %>%
        select(any_of(c("method", "leaf", "pval", "coef", "reject"))) %>%
        mutate(
          reject  = as.logical(reject),
          is_true = leaf %in% true_pairs
        )
    })
  
  # ── Metrics ──────────────────────────────────────────────────────────
  metrics = allRaw %>%
    group_by(method) %>%
    group_modify(~ compute_metrics(.x, true_pairs)) %>%
    ungroup()
  
  list(raw = allRaw, metrics = metrics)
}

run_all_methods_keren = function(spicy_sub, spicy_major, spicy_lineage, spe,
                           conditionB = "signal",
                           true_pairs = TRUE_PAIRS) {
  
  message("\n=== run_all_methods started ===")
  message("Timestamp: ", Sys.time())
  
  # ── Flat BH ──────────────────────────────────────────────────────────
  message("\n[1/7] Running Flat BH...")
  flatBh = extract_leaf_stats(spicy_sub, conditionB) %>%
    mutate(method = "Flat BH", reject = bh_reject(pval, q = 0.05))
  message("  Done. Rejections: ", sum(flatBh$reject))
  
  # ── Flat BY ──────────────────────────────────────────────────────────
  message("\n[2/7] Running Flat BY...")
  flatBy = extract_leaf_stats(spicy_sub, conditionB) %>%
    mutate(method = "Flat BY",
           reject = !is.na(pval) & p.adjust(pval, method = "BY") <= 0.05)
  message("  Done. Rejections: ", sum(flatBy$reject))
  
  # ── TreeBH ───────────────────────────────────────────────────────────
  message("\n[3/7] Running TreeBH...")
  treebh = tryCatch(
    {
      res = run_treebh(spicy_sub, conditionB = conditionB)
      message("  Done. Rejections: ", sum(res$reject))
      res
    },
    error = function(e) {
      message("  TreeBH failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>% mutate(method = "TreeBH", reject = FALSE)
    }
  )
  
  # ── DART ─────────────────────────────────────────────────────────────
  message("\n[4/7] Running DART...")
  message("  (this may hang on large datasets — kill if no output after 2 min)")
  dart = tryCatch(
    {
      run_dart_one(spicy_sub, spe, conditionB = conditionB, subject_col = "DONOR_NO")
      message("  Done. Rejections: ", sum(res$reject))
      res
    },
    error = function(e) {
      message("  DART failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>% mutate(method = "DART", reject = FALSE)
    }
  )
  
  # ── p-filter ─────────────────────────────────────────────────────────
  message("\n[5/7] Running p-filter...")
  pfilterRes = tryCatch(
    {
      res = run_pfilter_one(spicy_sub, spicy_major, conditionB = conditionB)
      message("  Done. Rejections: ", sum(res$reject))
      res
    },
    error = function(e) {
      message("  pfilter failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>%
        mutate(method = "pfilter", reject = FALSE)
    }
  )
  
  # ── Ours ─────────────────────────────────────────────────────────────
  message("\n[6/7] Running Ours...")
  oursRes = tryCatch({
    out = run_ours_one(spicy_lineage, spicy_major, spicy_sub,
                       conditionB = conditionB, l1_method = "BY")
    res = out$res_ours %>%
      filter(level == 3) %>%
      transmute(leaf = node, pval = pvalue, method = "Ours", reject)
    message("  Done. Rejections: ", sum(res$reject))
    res
  }, error = function(e) {
    message("  Ours failed: ", conditionMessage(e))
    extract_leaf_stats(spicy_sub, conditionB) %>% mutate(method = "Ours", reject = FALSE)
  })
  
  # ── ClusterTrim (BY) ─────────────────────────────────────────────────
  message("\n[7/7] Running ClusterTrim...")
  clustertrim = tryCatch(
    {
      res = run_cluster_trim_one(spicy_sub, spicy_major, conditionB = conditionB, major_method = "BH")
      message("  Done. Rejections: ", sum(res$reject))
      res
    },
    error = function(e) {
      message("  ClusterTrim_BY failed: ", conditionMessage(e))
      extract_leaf_stats(spicy_sub, conditionB) %>%
        mutate(method = "ClusterTrim_BY", reject = FALSE)
    }
  )
  
  # ── Combine ──────────────────────────────────────────────────────────
  message("\nCombining results...")
  allRaw = list(flatBh, flatBy, treebh, dart, pfilterRes, oursRes, clustertrim) %>%
    purrr::map_dfr(function(df) {
      df %>%
        select(any_of(c("method", "leaf", "pval", "coef", "reject"))) %>%
        mutate(
          reject  = as.logical(reject),
          is_true = leaf %in% true_pairs
        )
    })
  
  # ── Metrics ──────────────────────────────────────────────────────────
  message("Computing metrics...")
  metrics = allRaw %>%
    group_by(method) %>%
    group_modify(~ compute_metrics(.x, true_pairs)) %>%
    ungroup()
  
  message("\n=== run_all_methods complete ===")
  message("Timestamp: ", Sys.time())
  
  list(raw = allRaw, metrics = metrics)
}

#' Compute TPR, FDR, FPR, Precision, and F1 from a results tibble.
compute_metrics = function(result_df, true_pairs = TRUE_PAIRS) {
  result_df = result_df %>%
    mutate(
      is_true  = leaf %in% true_pairs,
      rejected = as.logical(reject)
    )
  
  TP = sum( result_df$is_true  &  result_df$rejected)
  FP = sum(!result_df$is_true  &  result_df$rejected)
  FN = sum( result_df$is_true  & !result_df$rejected)
  TN = sum(!result_df$is_true  & !result_df$rejected)
  
  TPR       = if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  FDR       = if ((TP + FP) > 0) FP / (TP + FP) else 0
  FPR       = if ((FP + TN) > 0) FP / (FP + TN) else NA_real_
  Precision = if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  F1        = if (!is.na(Precision) && !is.na(TPR) && (Precision + TPR) > 0)
    2 * Precision * TPR / (Precision + TPR) else NA_real_
  
  tibble(TPR, FDR, FPR, Precision, F1, TP, FP, FN, TN)
}