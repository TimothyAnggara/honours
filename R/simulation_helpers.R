library(tidyverse)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(SpatialExperiment)
library(purrr)
library(dplyr)

# --------------------------------------------------
# safe_bw_diggle()
#
# Purpose:
# Compute a kernel bandwidth for a point pattern using
# bw.diggle(), with a configurable upper search limit.
#
# Input:
#   pp   = spatstat point pattern (ppp)
#   hmax = maximum bandwidth allowed
#
# Output:
#   numeric bandwidth value sigma
# --------------------------------------------------
safe_bw_diggle <- function(pp, hmax = 150) {
  sigma <- bw.diggle(pp, hmax = hmax)
  sigma
}

# --------------------------------------------------
# choose_anchor_n()
#
# Purpose:
# Choose the number of anchor points used to generate
# a spatial density for a branch.
#
# Idea:
# Anchor count is proportional to the final cell count
# but constrained between min_n and max_n.
#
# Input:
#   n_final = final cell count
#   min_n   = minimum anchor points
#   scale   = proportion of n_final
#   max_n   = maximum anchor points
#
# Output:
#   integer-like anchor count
# --------------------------------------------------
choose_anchor_n <- function(n_final, min_n = 40, scale = 0.35, max_n = 250) {
  max(min_n, min(max_n, round(n_final * scale)))
}

# --------------------------------------------------
# split_count()
#
# Purpose:
# Split a total count n into integer subgroup counts
# according to given proportions.
#
# Idea:
# Uses multinomial sampling to randomly allocate n
# observations across groups with the given proportions.
# The result is random, not deterministic — counts vary
# across calls with the same inputs.
#
# Input:
#   n     = total count
#   props = vector of proportions
#   names = optional names for output groups
#
# Output:
#   integer vector summing to n
# --------------------------------------------------
split_count <- function(n, props, group_names = NULL) {
  if (length(n) != 1 || is.na(n) || !is.finite(n) || n < 0) {
    stop("n must be a single finite non-negative number.")
  }
  
  n <- as.integer(round(unname(n)))
  
  if (any(is.na(props)) || any(props < 0)) {
    stop("props must be non-missing and non-negative.")
  }
  
  props <- props / sum(props)
  
  out <- as.vector(rmultinom(1, size = n, prob = props))
  
  if (!is.null(group_names)) {
    names(out) <- group_names
  } else if (!is.null(names(props))) {
    names(out) <- names(props)
  }
  
  out
}

# --------------------------------------------------
# normalize_im()
#
# Purpose:
# Convert a spatial image into a proper probability density
# so that it can be used safely in rpoint().
# --------------------------------------------------
normalize_im = function(im_obj, tol = 1e-12, strict = FALSE) {
  if (!inherits(im_obj, "im")) stop("im_obj must be a spatstat image of class 'im'.")
  v = im_obj$v
  v[!is.finite(v)] = 0
  v[v < 0] = 0
  im_obj$v = v
  total_mass = integral.im(im_obj)
  if (!is.finite(total_mass) || total_mass <= tol)
    stop("Image cannot be normalized because total mass is not positive after cleaning.")
  im_obj / total_mass
}

# --------------------------------------------------
# label_cells()
#
# Purpose:
# Convert a spatstat point pattern into a tibble and
# assign a leaf cell type label.
#
# Inputs:
#   pp           = spatstat point pattern (ppp object)
#   cell_subtype = most specific cell label
#
# Output:
#   tibble with x, y, cellType_sub
# --------------------------------------------------
label_cells <- function(pp, cell_subtype) {
  if (!inherits(pp, "ppp")) {
    stop("pp must be a spatstat point pattern object of class 'ppp'.")
  }
  
  tibble(
    x = pp$x,
    y = pp$y,
    cellType_sub = cell_subtype
  )
}

# --------------------------------------------------
# add_tree_labels()
#
# Purpose:
# Use TREE_MAP to recover broader hierarchy labels
# from the leaf subtype label.
#
# Inputs:
#   cells    = tibble with a column called cellType_sub
#   tree_map = lookup table mapping subtype to hierarchy
#
# Output:
#   tibble with hierarchy columns added
# --------------------------------------------------
add_tree_labels <- function(cells, tree_map = TREE_MAP) {
  if (!"cellType_sub" %in% colnames(cells)) {
    stop("cells must contain a column called 'cellType_sub'.")
  }
  
  tree_renamed <- tree_map
  colnames(tree_renamed)[colnames(tree_renamed) == "major"]   <- "cellType_major"
  colnames(tree_renamed)[colnames(tree_renamed) == "subtype"] <- "cellType_sub"
  
  out <- merge(
    x    = cells,
    y    = tree_renamed,
    by   = "cellType_sub",
    all.x = TRUE
  )
  
  missing_types <- unique(
    out$cellType_sub[is.na(out$cellType_major) | is.na(out$lineage)]
  )
  
  if (length(missing_types) > 0) {
    stop(
      paste(
        "These subtypes are missing from TREE_MAP:",
        paste(missing_types, collapse = ", ")
      )
    )
  }
  
  out
}

# --------------------------------------------------
# simulate_cells_from_density()
#
# Purpose:
# Simulate n cells from a spatial density and label them
# with a leaf subtype.
#
# Inputs:
#   n            = number of cells
#   cell_subtype = leaf cell label
#   density      = spatstat image
#   win          = simulation window
#
# Output:
#   tibble with x, y, cellType_sub
# --------------------------------------------------
simulate_cells_from_density <- function(n, cell_subtype, density, win) {
  if (n <= 0) return(NULL)
  
  d <- normalize_im(density)
  pp <- rpoint(n = n, f = d, win = win)
  
  label_cells(pp, cell_subtype)
}

# --------------------------------------------------
# simulate_mixture_cells()
#
# Purpose:
# Simulate n cells from a mixture of two densities and
# label them with a leaf subtype.
#
# Inputs:
#   n            = number of cells
#   cell_subtype = leaf cell label
#   density1     = first density
#   density2     = second density
#   weight       = mixing weight in [0, 1]
#   win          = simulation window
#
# Output:
#   tibble with x, y, cellType_sub
# --------------------------------------------------
simulate_mixture_cells <- function(n, cell_subtype,
                                   density1, density2, weight, win) {
  if (n <= 0) return(NULL)
  
  if (!is.numeric(weight) || length(weight) != 1 || weight < 0 || weight > 1) {
    stop("weight must be a single number between 0 and 1.")
  }
  
  d1 <- normalize_im(density1)
  d2 <- normalize_im(density2)
  
  mixed_density <- (1 - weight) * d1 + weight * d2
  pp <- rpoint(n = n, f = mixed_density, win = win)
  
  label_cells(pp, cell_subtype)
}

# -------------------------------------------------- 
# simulate_pp() 
# 
# Purpose: 
# Simulate a completely spatially random point pattern 
# inside a given window.
# 
# Idea: 
# This uses a homogeneous Poisson point process, meaning: 
# - points are placed independently of one another 
# - every location in the window is equally likely 
# - there is no clustering or repulsion 
# 
# The intensity is chosen so that the expected total number # of points in the window is n. 
# 
# Inputs: 
# n = desired expected number of points 
# win = spatstat window object defining the simulation region 
# 
# Output: 
# A tibble with one row per simulated point, containing: 
# x, y = simulated spatial coordinates 
# 
# Notes: 
# - The number of simulated points is random, not fixed exactly at n 
# - The actual number of points is Poisson distributed with mean n 
# - This is a useful null model for complete spatial randomness 
# -------------------------------------------------- 

simulate_pp <- function(n, win) { 
  rpoispp(n / area.owin(win), win = win) 
}

# --------------------------------------------------
# resimulate_celltype_with_mixture()
#
# Purpose:
# Post hoc resimulate one cell subtype using a mixture of:
#   1. its own realized density
#   2. another cell subtype's realized density
#
# This lets you introduce an artificial spatial signal
# after all_cells has already been constructed.
#
# Inputs:
#   all_cells     = full tibble of simulated cells
#   target_type   = subtype to remove and resimulate
#   partner_type  = one or more subtype labels whose
#                   density will be mixed in
#   weight        = mixture weight for partner density
#   win           = spatstat window
#   keep_metadata = if TRUE, preserve extra metadata
#                   columns such as sample_id, condition
#
# Output:
#   updated tibble with target_type resimulated
# --------------------------------------------------
resimulate_celltype_with_mixture <- function(
    all_cells,
    target_type,
    partner_type,
    weight,
    win,
    keep_metadata = TRUE
) {
  if (!all(c("x", "y", "cellType_sub") %in% colnames(all_cells))) {
    stop("all_cells must contain columns: x, y, cellType_sub")
  }
  
  if (!is.character(target_type) || length(target_type) != 1) {
    stop("target_type must be a single subtype name.")
  }
  
  if (!is.character(partner_type) || length(partner_type) < 1) {
    stop("partner_type must be one or more subtype names.")
  }
  
  if (!is.numeric(weight) || length(weight) != 1 || weight < 0 || weight > 1) {
    stop("weight must be a single number between 0 and 1.")
  }
  
  target_cells <- all_cells %>%
    filter(cellType_sub == target_type)
  
  partner_cells <- all_cells %>%
    filter(cellType_sub %in% partner_type)
  
  remaining_cells <- all_cells %>%
    filter(cellType_sub != target_type)
  
  if (nrow(target_cells) == 0) {
    stop(paste("No cells found for target_type:", target_type))
  }
  
  if (nrow(partner_cells) == 0) {
    stop(
      paste(
        "No cells found for partner_type:",
        paste(partner_type, collapse = ", ")
      )
    )
  }
  
  # Convert realized cells into point patterns
  target_pp <- ppp(
    x = target_cells$x,
    y = target_cells$y,
    window = win
  )
  
  partner_pp <- ppp(
    x = partner_cells$x,
    y = partner_cells$y,
    window = win
  )
  
  # Estimate realized densities using safe_bw_diggle()
  target_density <- density.ppp(
    target_pp,
    sigma = safe_bw_diggle(target_pp)
  )
  
  partner_density <- density.ppp(
    partner_pp,
    sigma = safe_bw_diggle(partner_pp)
  )
  
  # Resimulate the target cell type from the mixture
  target_cells_new <- simulate_mixture_cells(
    n = nrow(target_cells),
    cell_subtype = target_type,
    density1 = target_density,
    density2 = partner_density,
    weight = weight,
    win = win
  )
  
  # Rebuild dataset
  out <- bind_rows(
    remaining_cells %>% select(x, y, cellType_sub),
    target_cells_new
  ) %>%
    add_tree_labels()
  
  # Reattach metadata if needed
  if (keep_metadata) {
    meta_cols <- setdiff(
      colnames(all_cells),
      c("x", "y", "cellType_sub", "cellType_major", "lineage")
    )
    
    if (length(meta_cols) > 0) {
      meta_template <- all_cells %>%
        select(all_of(meta_cols))
      
      # Assert all metadata columns are constant across rows
      non_constant <- meta_cols[
        vapply(meta_template, function(col) length(unique(col)) > 1, logical(1))
      ]
      if (length(non_constant) > 0) {
        stop(paste(
          "resimulate_celltype_with_mixture() expects constant metadata",
          "across all cells, but these columns vary:",
          paste(non_constant, collapse = ", ")
        ))
      }
      
      meta_template <- meta_template[1, , drop = FALSE]
      
      out <- bind_cols(
        out,
        meta_template[rep(1, nrow(out)), , drop = FALSE]
      )
    }
  }
  
  out
}


# --------------------------------------------------
# make_realistic_counts()
#
# Purpose:
# Generate realistic cell counts across the hierarchy
# from a chosen total number of cells.
#
# Idea:
# Split the total count step by step through the tree:
# top level, major branches, then leaf subtypes.
#
# Input:
#   total_cells = total number of cells to simulate
#
# Output:
#   list containing counts for each level of the
#   hierarchy, including final leaf cell counts
# --------------------------------------------------
make_realistic_counts <- function(total_cells = 5000) {
  
  total_cells <- as.integer(round(unname(total_cells)))
  
  # -----------------------------
  # Top level
  # -----------------------------
  top <- split_count(
    total_cells,
    props = c(0.45, 0.25, 0.30),
    group_names = c("Immune", "Stromal", "Tissue")
  )
  
  # -----------------------------
  # Immune split
  # -----------------------------
  immune <- split_count(
    unname(top["Immune"]),
    props = c(0.55, 0.45),
    group_names = c("Lymphoid", "Myeloid")
  )
  
  # -----------------------------
  # Lymphoid split
  # -----------------------------
  lymphoid <- split_count(
    unname(immune["Lymphoid"]),
    props = c(0.20, 0.65, 0.15),
    group_names = c("B", "T", "NK")
  )
  
  b_branch <- split_count(
    unname(lymphoid["B"]),
    props = c(0.70, 0.30),
    group_names = c("Naive_B", "Memory_B")
  )
  
  t_branch <- split_count(
    unname(lymphoid["T"]),
    props = c(0.45, 0.35, 0.20),
    group_names = c("CD4_T", "CD8_T", "Memory_T")
  )
  
  # -----------------------------
  # Myeloid split
  # -----------------------------
  myeloid <- split_count(
    unname(immune["Myeloid"]),
    props = c(0.20, 0.25, 0.25, 0.20, 0.10),
    group_names = c("Monocyte", "Macrophage", "Granulocyte", "Dendritic", "Mast")
  )
  
  macrophage <- split_count(
    unname(myeloid["Macrophage"]),
    props = c(0.45, 0.55),
    group_names = c("M1_like", "M2_like")
  )
  
  granulocyte <- split_count(
    unname(myeloid["Granulocyte"]),
    props = c(0.70, 0.20, 0.10),
    group_names = c("Neutrophil", "Eosinophil", "Basophil")
  )
  
  dendritic <- split_count(
    unname(myeloid["Dendritic"]),
    props = c(0.30, 0.45, 0.25),
    group_names = c("cDC1", "cDC2", "pDC")
  )
  
  # -----------------------------
  # Stromal split
  # -----------------------------
  stromal <- split_count(
    unname(top["Stromal"]),
    props = c(0.50, 0.30, 0.20),
    group_names = c("Fibroblast", "Endothelial", "Pericyte")
  )
  
  # -----------------------------
  # Tissue
  # -----------------------------
  tissue <- c(Epithelial = unname(top["Tissue"]))
  
  # -----------------------------
  # Final leaf counts
  # -----------------------------
  leaf_counts <- c(
    b_branch,
    t_branch,
    NK = unname(lymphoid["NK"]),
    Monocyte = unname(myeloid["Monocyte"]),
    macrophage,
    granulocyte,
    dendritic,
    Mast = unname(myeloid["Mast"]),
    stromal,
    tissue
  )
  
  # final safety checks
  if (any(is.na(leaf_counts))) {
    stop("Some leaf counts are NA.")
  }
  
  if (sum(leaf_counts) != total_cells) {
    stop("Leaf counts do not sum to total_cells.")
  }
  
  list(
    total_cells = total_cells,
    top = top,
    immune = immune,
    lymphoid = lymphoid,
    myeloid = myeloid,
    leaf = leaf_counts
  )
}

make_spe_from_df <- function(all_sims) {
  all_sims <- all_sims %>%
    mutate(cell_id = paste0("cell_", row_number()))
  
  cell_metadata <- all_sims %>%
    select(
      cell_id,
      subject_id,
      image_id,
      condition,
      cellType_sub,
      cellType_major,
      lineage
    )
  
  sample_metadata <- all_sims %>%
    distinct(subject_id, image_id, condition)
  
  spe <- SpatialExperiment(
    colData = S4Vectors::DataFrame(cell_metadata),
    spatialCoords = as.matrix(all_sims[, c("x", "y")])
  )
  
  colnames(spe) <- all_sims$cell_id
  spe$subject   <- factor(all_sims$subject_id)
  spe$imageID   <- factor(all_sims$image_id)
  spe$condition <- factor(all_sims$condition)
  
  metadata(spe)$sample_metadata <- sample_metadata
  
  spe
}

make_subject_scaffold <- function(
    total_cells = 5000,
    win = WIN
) {
  message("  [scaffold] generating counts and densities...")
  counts <- make_realistic_counts(total_cells)
  
  main_pp <- simulate_pp(
    n = choose_anchor_n(total_cells, min_n = 150, scale = 0.12, max_n = 500),
    win = win
  )
  main_density <- density.ppp(main_pp, safe_bw_diggle(main_pp))
  
  immune_pp <- rpoint(
    n = choose_anchor_n(counts$top["Immune"]),
    f = normalize_im(main_density),
    win = win
  )
  immune_density <- density.ppp(immune_pp, safe_bw_diggle(immune_pp))
  
  stromal_pp <- rpoint(
    n = choose_anchor_n(counts$top["Stromal"]),
    f = normalize_im(main_density),
    win = win
  )
  stromal_density <- density.ppp(stromal_pp, safe_bw_diggle(stromal_pp))
  
  tissue_pp <- rpoint(
    n = choose_anchor_n(counts$top["Tissue"]),
    f = normalize_im(main_density),
    win = win
  )
  tissue_density <- density.ppp(tissue_pp, safe_bw_diggle(tissue_pp))
  
  lymphoid_pp <- rpoint(
    n = choose_anchor_n(counts$immune["Lymphoid"]),
    f = normalize_im(immune_density),
    win = win
  )
  lymphoid_density <- density.ppp(lymphoid_pp, safe_bw_diggle(lymphoid_pp))
  
  myeloid_pp <- rpoint(
    n = choose_anchor_n(counts$immune["Myeloid"]),
    f = normalize_im(immune_density),
    win = win
  )
  myeloid_density <- density.ppp(myeloid_pp, safe_bw_diggle(myeloid_pp))
  
  b_total <- counts$leaf["Naive_B"] + counts$leaf["Memory_B"]
  b_pp <- rpoint(
    n = choose_anchor_n(b_total),
    f = normalize_im(lymphoid_density),
    win = win
  )
  b_density <- density.ppp(b_pp, safe_bw_diggle(b_pp))
  
  t_total <- counts$leaf["CD4_T"] + counts$leaf["CD8_T"] + counts$leaf["Memory_T"]
  t_pp <- rpoint(
    n = choose_anchor_n(t_total),
    f = normalize_im(lymphoid_density),
    win = win
  )
  t_density <- density.ppp(t_pp, safe_bw_diggle(t_pp))
  
  nk_pp <- rpoint(
    n = choose_anchor_n(counts$leaf["NK"]),
    f = normalize_im(lymphoid_density),
    win = win
  )
  nk_density <- density.ppp(nk_pp, safe_bw_diggle(nk_pp))
  
  monocyte_pp <- rpoint(
    n = choose_anchor_n(counts$leaf["Monocyte"]),
    f = normalize_im(myeloid_density),
    win = win
  )
  monocyte_density <- density.ppp(monocyte_pp, safe_bw_diggle(monocyte_pp))
  
  macrophage_total <- counts$leaf["M1_like"] + counts$leaf["M2_like"]
  macrophage_pp <- rpoint(
    n = choose_anchor_n(macrophage_total),
    f = normalize_im(myeloid_density),
    win = win
  )
  macrophage_density <- density.ppp(macrophage_pp, safe_bw_diggle(macrophage_pp))
  
  granulocyte_total <- counts$leaf["Neutrophil"] + counts$leaf["Eosinophil"] + counts$leaf["Basophil"]
  granulocyte_pp <- rpoint(
    n = choose_anchor_n(granulocyte_total),
    f = normalize_im(myeloid_density),
    win = win
  )
  granulocyte_density <- density.ppp(granulocyte_pp, safe_bw_diggle(granulocyte_pp))
  
  dendritic_total <- counts$leaf["cDC1"] + counts$leaf["cDC2"] + counts$leaf["pDC"]
  dendritic_pp <- rpoint(
    n = choose_anchor_n(dendritic_total),
    f = normalize_im(myeloid_density),
    win = win
  )
  dendritic_density <- density.ppp(dendritic_pp, safe_bw_diggle(dendritic_pp))
  
  mast_pp <- rpoint(
    n = choose_anchor_n(counts$leaf["Mast"]),
    f = normalize_im(myeloid_density),
    win = win
  )
  mast_density <- density.ppp(mast_pp, safe_bw_diggle(mast_pp))
  
  list(
    counts = counts,
    win = win,
    b_density = b_density,
    t_density = t_density,
    nk_density = nk_density,
    monocyte_density = monocyte_density,
    macrophage_density = macrophage_density,
    granulocyte_density = granulocyte_density,
    dendritic_density = dendritic_density,
    mast_density = mast_density,
    stromal_density = stromal_density,
    tissue_density = tissue_density
  )
}

simulate_image_from_scaffold <- function(
    scaffold,
    image_id,
    subject_id,
    condition
) {
  message(sprintf("    [image] simulating %s", image_id))
  counts <- scaffold$counts
  win <- scaffold$win
  
  naive_b_cells <- simulate_cells_from_density(counts$leaf["Naive_B"], "Naive_B", scaffold$b_density, win)
  memory_b_cells <- simulate_cells_from_density(counts$leaf["Memory_B"], "Memory_B", scaffold$b_density, win)
  
  CD4_t_cells <- simulate_cells_from_density(counts$leaf["CD4_T"], "CD4_T", scaffold$t_density, win)
  CD8_t_cells <- simulate_cells_from_density(counts$leaf["CD8_T"], "CD8_T", scaffold$t_density, win)
  memory_t_cells <- simulate_cells_from_density(counts$leaf["Memory_T"], "Memory_T", scaffold$t_density, win)
  
  nk_cells <- simulate_cells_from_density(counts$leaf["NK"], "NK", scaffold$nk_density, win)
  
  monocyte_cells <- simulate_cells_from_density(counts$leaf["Monocyte"], "Monocyte", scaffold$monocyte_density, win)
  m1_like_cells <- simulate_cells_from_density(counts$leaf["M1_like"], "M1_like", scaffold$macrophage_density, win)
  m2_like_cells <- simulate_cells_from_density(counts$leaf["M2_like"], "M2_like", scaffold$macrophage_density, win)
  
  neutrophil_cells <- simulate_cells_from_density(counts$leaf["Neutrophil"], "Neutrophil", scaffold$granulocyte_density, win)
  eosinophil_cells <- simulate_cells_from_density(counts$leaf["Eosinophil"], "Eosinophil", scaffold$granulocyte_density, win)
  basophil_cells <- simulate_cells_from_density(counts$leaf["Basophil"], "Basophil", scaffold$granulocyte_density, win)
  
  cdc1_cells <- simulate_cells_from_density(counts$leaf["cDC1"], "cDC1", scaffold$dendritic_density, win)
  cdc2_cells <- simulate_cells_from_density(counts$leaf["cDC2"], "cDC2", scaffold$dendritic_density, win)
  pdc_cells <- simulate_cells_from_density(counts$leaf["pDC"], "pDC", scaffold$dendritic_density, win)
  
  mast_cells <- simulate_cells_from_density(counts$leaf["Mast"], "Mast", scaffold$mast_density, win)
  
  fibroblast_cells <- simulate_cells_from_density(counts$leaf["Fibroblast"], "Fibroblast", scaffold$stromal_density, win)
  endothelial_cells <- simulate_cells_from_density(counts$leaf["Endothelial"], "Endothelial", scaffold$stromal_density, win)
  pericyte_cells <- simulate_cells_from_density(counts$leaf["Pericyte"], "Pericyte", scaffold$stromal_density, win)
  
  epithelial_cells <- simulate_cells_from_density(counts$leaf["Epithelial"], "Epithelial", scaffold$tissue_density, win)
  
  bind_rows(
    naive_b_cells, memory_b_cells,
    CD4_t_cells, CD8_t_cells, memory_t_cells, nk_cells,
    monocyte_cells, m1_like_cells, m2_like_cells,
    neutrophil_cells, eosinophil_cells, basophil_cells,
    cdc1_cells, cdc2_cells, pdc_cells, mast_cells,
    fibroblast_cells, endothelial_cells, pericyte_cells,
    epithelial_cells
  ) %>%
    add_tree_labels() %>%
    mutate(
      subject_id = subject_id,
      image_id = image_id,
      condition = condition
    )
}

simulate_one_subject_consistent <- function(
    subject_id,
    n_images = 20,
    win = WIN,
    condition = "control",
    add_signal = FALSE,
    signal_weight = 0.7,
    total_cells = 5000
) {
  message(sprintf("  [subject] starting subject %s (%s)", subject_id, condition))
  scaffold <- make_subject_scaffold(
    total_cells = total_cells,
    win = win
  )
  
  purrr::map_dfr(seq_len(n_images), function(img_id) {
    img <- simulate_image_from_scaffold(
      scaffold  = scaffold,
      image_id  = paste0(condition, "_subj_", subject_id, "_img_", img_id),
      subject_id = paste0("subj_", subject_id),
      condition = condition
    )
    
    if (add_signal) {
      message(sprintf("    [signal] injecting signal into image %d", img_id))
      img <- resimulate_celltype_with_mixture(
        all_cells    = img,
        target_type  = "NK",
        partner_type = c("CD4_T", "CD8_T", "Memory_T"),
        weight       = signal_weight,
        win          = win
      )
    }
    
    img
  })
}

