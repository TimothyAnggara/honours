library(tidyverse)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(SpatialExperiment)
library(purrr)
library(dplyr)

KEREN_TREE_MAP = tribble(
  ~subtype,         ~major,       ~lineage,
  
  # Lymphoid
  "CD4_T_cell",     "T",          "Immune",
  "CD8_T_cell",     "T",          "Immune",
  "dn_T_CD3",       "T",          "Immune",
  "Tregs",          "T",          "Immune",
  "B_cell",         "B",          "Immune",
  "NK",             "NK",         "Immune",
  
  # Myeloid
  "Macrophages",    "Macrophage", "Immune",
  "Mono_or_Neu",    "Monocyte",   "Immune",
  "DC",             "Dendritic",  "Immune",
  "DC_or_Mono",     "Dendritic",  "Immune",
  "Neutrophils",    "Granulocyte","Immune",
  
  # Other immune
  "Other_Immune",   "Other",      "Immune",
  
  # Stromal
  "Endothelial",    "Vascular",   "Stromal",
  "Mesenchymal",    "Fibroblast", "Stromal",
  
  # Tumour
  "Keratin_Tumour", "Epithelial", "Tumour",
  "Tumour",         "Epithelial", "Tumour"
)

# Lookup maps
KEREN_MAJOR_MAP = deframe(
  KEREN_TREE_MAP %>% select(subtype, major)
)

KEREN_MAJOR_LINEAGE_MAP = deframe(
  KEREN_TREE_MAP %>% distinct(major, lineage) %>% select(major, lineage)
)

KEREN_LINEAGE_MAP = deframe(
  KEREN_TREE_MAP %>% select(subtype, lineage)
)