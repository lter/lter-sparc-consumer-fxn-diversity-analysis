################################################################################
##
## Script to compute functional indices through time
##
## Camille Magneville
##
## 02/2026
##
################################################################################


# 0 - Set up ===================================================================


# Define the pipe symbol so I can use it:
`%>%` <- magrittr::`%>%`

# Load libraries
librarian::shelf(tidyverse, dplyr, funbiogeo, ggplot2, mFD, ape)

# Get set up
source("00_setup.R")

# Clear environment & collect garbage
rm(list = ls()); gc()


# 1 - Get species coordinates in the functional space ==========================


# Load our studied traits:
sp_tr_zscore <- readRDS(file.path("transformed_data",
                  "sp_tr_zscore_new.rds"))

# Subset it to tr and fish/birds for now:
sp_tr_df <- sp_tr_zscore %>% 
  dplyr::filter(taxon %in% c("Birds", "Fish")) %>% 
  dplyr::select(c("scientific_name",
                  "tr.mass.adult.zt",
                  "tr.trophic.level.zt",
                  "tr.reproduction.unified.zt",
                  "tr.age.zt")) %>% 
  dplyr::distinct() %>% 
  dplyr::ungroup() 

# One species has two names: Zonotrichia leucophyrs, arbitrarily chose the first
# TO FIX
sp_tr_df2 <- sp_tr_df %>% 
  dplyr::filter(scientific_name != "Zonotrichia leucophyrs") %>% 
  tibble::column_to_rownames(var = "scientific_name") %>% 
  dplyr::select(-c("taxon"))

# Build a dataframe gathering traits categories:
tr_nm <- colnames(sp_tr_df2)
tr_cat <- c("Q", "Q", "Q","Q")
tr_cat_df <- as.data.frame(matrix(ncol = 2, nrow = 4))
tr_cat_df[, 1] <- tr_nm
tr_cat_df[, 2] <- tr_cat
colnames(tr_cat_df) <- c("trait_name", "trait_type")


# Build functional distances (used gower as NAs):
sp_dist_df <- mFD::funct.dist(
  sp_tr         = sp_tr_df2,
  tr_cat        = tr_cat_df,
  metric        = "gower",
  scale_euclid  = "noscale",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = FALSE)

# Get functional space quality and sp coord:
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist_df,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Check that 3D space still ok:
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# Check correlation traits and axes:
tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sp_tr_df2, 
  sp_faxes_coord = sp_faxes_coord_df[ , c("PC1", "PC2", "PC3")], 
  stop_if_NA = FALSE,
  plot = TRUE)
tr_faxes

sp_faxes_coord_df <- fspaces_quality$"details_fspaces"$"sp_pc_coord"


# Save functional distance matrix:
saveRDS(sp_dist_df,
        file.path("transformed_data",
                  "funct_dist_matrix.rds"))

# Save species coordinates:
saveRDS(sp_faxes_coord_df, 
        file.path("transformed_data",
                  "sp_faxes_coord.rds"))



