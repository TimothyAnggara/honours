library(tidyverse)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(SpatialExperiment)
library(purrr)
library(dplyr)

TREE_MAP <- tribble(
  ~subtype,      ~major,      ~lineage,
  "Naive_B",     "B",         "Lymphoid",
  "Memory_B",    "B",         "Lymphoid",
  "CD4_T",       "T",         "Lymphoid",
  "CD8_T",       "T",         "Lymphoid",
  "Memory_T",    "T",         "Lymphoid",
  "NK",          "NK",        "Lymphoid",
  "Monocyte",    "Monocyte",  "Myeloid",
  "M1_like",     "Macrophage","Myeloid",
  "M2_like",     "Macrophage","Myeloid",
  "Neutrophil",  "Granulocyte","Myeloid",
  "Eosinophil",  "Granulocyte","Myeloid",
  "Basophil",    "Granulocyte","Myeloid",
  "cDC1",        "Dendritic", "Myeloid",
  "cDC2",        "Dendritic", "Myeloid",
  "pDC",         "Dendritic", "Myeloid",
  "Mast",        "Mast",      "Myeloid",
  "Fibroblast",  "Stromal",   "Stromal",
  "Endothelial", "Stromal",   "Stromal",
  "Pericyte",    "Stromal",   "Stromal",
  "Epithelial",  "Tissue",    "Tissue"
)

# subtype -> major (e.g. "CD4_T" -> "T")
MAJOR_MAP <- deframe(
  TREE_MAP %>% select(subtype, major)
)

# major -> lineage (e.g. "T" -> "Lymphoid", "NK" -> "Lymphoid")
MAJOR_LINEAGE_MAP <- deframe(
  TREE_MAP %>% distinct(major, lineage) %>% select(major, lineage)
)

# subtype -> lineage (e.g. "CD4_T" -> "Lymphoid")
LINEAGE_MAP <- deframe(
  TREE_MAP %>% select(subtype, lineage)
)


WIN = owin(xrange = c(0, 1000), yrange = c(0, 1000))

TRUE_PAIRS = c(
  "CD4_T__NK",
  "CD8_T__NK",
  "Memory_T__NK"
)