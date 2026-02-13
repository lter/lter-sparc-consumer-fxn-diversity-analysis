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
                  "sp_tr_zscore.rds"))

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

sp_faxes_coord_df <- fspaces_quality$"details_fspaces"$"sp_pc_coord"
saveRDS(sp_faxes_coord_df, 
        file.path("transformed_data",
                  "sp_faxes_coord.rds"))

# Check correlation traits and axes:
tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sp_tr_df2, 
  sp_faxes_coord = sp_faxes_coord_df[ , c("PC1", "PC2", "PC3")], 
  stop_if_NA = FALSE,
  plot = TRUE)
tr_faxes

# 2 - Plot functional space with arrows ========================================

# Using Matthias Grenie functions: https://gist.github.com/Rekyt/ee15330639f8719d87aebdb8a5b095d4
# Adapted to our data

# Function 1 to add arrows - get covariance btw tr and faxes:
trait_df <- sp_tr_df2 
sp_coord <- sp_faxes_coord_df 
pc_eigenvalues <- fspaces_quality$details_fspaces$pc_eigenvalues

compute_arrows <- function(sp_coord, trait_df, pc_eigenvalues) {
  
  # Get nb of species
  n <- nrow(trait_df)
  points.stand <- scale(sp_faxes_coord_df)
  
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand, use = "pairwise.complete.obs")
  
  # Select only positive eigenvalues
  pos_eigen = pc_eigenvalues$Eigenvalues
  
  # Standardize value of covariance (see Legendre & Legendre 1998)
  std_cov <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(std_cov) <- colnames(sp_faxes_coord_df)
  
  return(std_cov)
}

std_cov <- compute_arrows(trait_df = sp_tr_df2,
                          sp_coord = sp_faxes_coord_df, 
                          pc_eigenvalues = fspaces_quality$details_fspaces$pc_eigenvalues)


# Basic plot with individuals:
# Add taxa colors:
sp_faxes_coord_df2 <- as.data.frame(sp_faxes_coord_df) %>% 
  tibble::rownames_to_column(var = "scientific_name") %>% 
  dplyr::left_join(sp_tr_zscore[, c("scientific_name",
                                    "taxon")],
                   by = "scientific_name") %>%
  dplyr::filter(taxon %in% c("Fish", "Birds")) %>% 
  dplyr::distinct() %>% 
  tibble::column_to_rownames(var = "scientific_name")
 
## Do PC1-PC2 plot: 
plot_pcoa_12 <- ggplot2::ggplot() +
  ggplot2::geom_point(sp_faxes_coord_df2, 
                      mapping = ggplot2::aes(x = PC1, y = PC2,
                                                       color = taxon),
                      alpha = 0.6) +
  ggplot2::scale_color_manual(values = c("#C9A227", "#1B9E77")) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme_bw() 
  
plot_pcoa_12

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable
std_cov2 <- as.data.frame(std_cov) %>% 
  tibble::rownames_to_column(var = "traits") %>%
  dplyr::mutate(
    PC1 = PC1 / 20,
    PC2 = PC2 / 20) %>% 
  dplyr::mutate(tr_nm = c("Body Mass", "Trophic Level", "Reproduction",
                          "Life Span"))

plot_pcoa_12 <- plot_pcoa_12 +
  ggplot2::geom_segment(data = as.data.frame(std_cov2),
               x = 0, y = 0, alpha = 0.7,
               mapping = ggplot2::aes(xend = PC1, yend = PC2),
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = std_cov2,
                  aes(x = PC1, y = PC2, label = tr_nm),
                  size = 4,
                  nudge_x = 0.01, nudge_y = 0.01) +
  ggplot2::theme(legend.position = "none") 
  
plot_pcoa_12

# Add species silhouettes?


## Do PC2-PC3 plot: 
plot_pcoa_23 <- ggplot2::ggplot() +
  ggplot2::geom_point(sp_faxes_coord_df2, 
                      mapping = ggplot2::aes(x = PC2, y = PC3,
                                             color = taxon),
                      alpha = 0.6) +
  ggplot2::scale_color_manual(values = c("#C9A227", "#1B9E77")) +
  ggplot2::theme_bw() 

plot_pcoa_23

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable
std_cov3 <- as.data.frame(std_cov) %>% 
  tibble::rownames_to_column(var = "traits") %>%
  dplyr::mutate(
    PC2 = PC2 / 100,
    PC3 = PC3 / 100) %>% 
  dplyr::mutate(tr_nm = c("Body Mass", "Trophic Level", "Reproduction",
                          "Life Span"))

plot_pcoa_23 <- plot_pcoa_23 +
  ggplot2::geom_segment(data = as.data.frame(std_cov3),
                        x = 0, y = 0, alpha = 0.7,
                        mapping = ggplot2::aes(xend = PC2, yend = PC3),
                        arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = std_cov3,
                           aes(x = PC2, y = PC3, label = tr_nm),
                           size = 4,
                           nudge_x = 0.01, nudge_y = 0.01)
plot_pcoa_23

# Add them together:
plot_fnal_space <- plot_pcoa_12 + plot_pcoa_23
plot_fnal_space

# 3 - Compute species uniqueness TO DO ===============================================


# Get functional distance between species:¨
sp_dist_all_df <- mFD::dist.to.df(list("dist" = sp_dist_df))

fish_sp_list <- sp_tr_df$scientific_name[which(sp_tr_df$taxon == "Fish")]
birds_sp_list <- sp_tr_df$scientific_name[which(sp_tr_df$taxon == "Birds")]

# Build df to keep the uniqueness of each species (based on fish or bird pool):
funiq_fish_df <- as.data.frame(matrix(ncol = 3, nrow = length(fish_sp_list), NA))
colnames(funiq_df) <- c("scientific_name", "FUni_score", "Taxa")
funiq_df$scientific_name <- unique(sp_tr_df$scientific_name)
funiq_df


#### STOPPED HERE


# 3 - Build the assemblage dataframe ===========================================


# Note: focus on the fish assemblages but still need all species in the asb df:

# Load community data:
comm_df <- read.csv(file.path("Data",
                              "community_tidy-data",
                              "04_harmonized_consumer_excretion_sparc_cnd_site.csv"))
sp_comm_df <- read.csv(file.path("transformed_data",
                              "Mack_data.csv")) #system
saveRDS(sp_comm_df,
        file.path("transformed_data", "Macks_data.rds"))

# Call species list:
# sp_list <- readRDS(file.path("transformed_data",
#                              "species_list_corrected_fish.rds"))

# Link comm and species list:
#sp_comm_df <- dplyr::left_join(comm_df, dplyr::distinct(sp_list[, c(4,12)]),
                               #by = "scientific_name")
sp_comm_df2 <- sp_comm_df %>% 
  dplyr::select(c("project", "habitat", "system", "year", "scientific_name")) 

# Add new column Project_site_year:
sp_comm_df2$samp_unit <- paste0(sp_comm_df2$project,
                                sep = "_",
                                sp_comm_df2$system,
                                sep = "_",
                                sp_comm_df2$year)
fish_asb_df <- sp_comm_df2 %>% 
  dplyr::select(c("scientific_name", "samp_unit")) %>% 
  dplyr::distinct() %>% # because several subsites inside each "site"
  dplyr::mutate(presence = 1) %>%
  tidyr::pivot_wider(names_from = scientific_name,
              values_from = presence,
              values_fill = 0)

# Then add all the other species used to build functional space:
# Get species used to build the functional space but not in fish asb:
other_sp <- setdiff(rownames(sp_faxes_coord_df),
                    colnames(fish_asb_df)[-1])
# Add these new species:
missing_sp_mat <- matrix(0,
                        nrow = nrow(fish_asb_df),
                        ncol = length(other_sp))
colnames(missing_sp_mat) <- other_sp

all_asb_temp_df <- cbind(fish_asb_df,
                         missing_sp_mat)

new_row <- rep(0, ncol(fish_asb_df))
names(new_row) <- colnames(fish_asb_df)

new_row[other_sp] <- 1

all_asb_df <- rbind(all_asb_temp_df, 
                    other = new_row)
all_asb_df$samp_unit[nrow(all_asb_df)] <- "Others"

all_asb_df <- all_asb_df %>% 
  tibble::column_to_rownames(var = "samp_unit")


# 4 - Compute functional diversity indices =====================================

# Remove sampling units which have less than 3 species:
all_asb_df$sp_richn <- rowSums(all_asb_df)
all_asb_df <- all_asb_df %>% 
  dplyr::filter(sp_richn > 3) %>% 
  dplyr::select(-c(sp_richn))

setdiff(colnames(all_asb_df), rownames(sp_faxes_coord_df))
setdiff(rownames(sp_faxes_coord_df), colnames(all_asb_df))

# Removing "Amia calva"               "Lepisosteus platyrhincus"
all_asb_df2 <- all_asb_df %>% 
  dplyr::select(-c( "Amia calva", "Lepisosteus platyrhincus"))
all_asb_df2$sp_richn <- rowSums(all_asb_df2)
all_asb_df2 <- all_asb_df2 %>% 
  dplyr::filter(sp_richn > 3) %>% 
  dplyr::select(-c(sp_richn))

alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_df[ , c("PC1", "PC2", "PC3")],
  asb_sp_w         = as.matrix(all_asb_df2),
  ind_vect         = c("fric", "fdis", "fspe", "fnnd", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

ind_df <- alpha_fd_indices$functional_diversity_indices %>% 
  tibble::rownames_to_column(var = "samp_unit")

# Combine the df with fd indices and the comm_df2 which has info about proj:
info_df <- sp_comm_df2 %>% 
  dplyr::select(-c("scientific_name")) %>% 
  dplyr::distinct()

# Remove the asb which are not studied (sp richn too low):
info_df <- dplyr::filter(info_df,
                         samp_unit %in% ind_df$samp_unit)

results_df <- ind_df %>% 
  dplyr::left_join(info_df, by = "samp_unit") %>% 
  dplyr::arrange(project, system, year) %>% 
  dplyr::filter(samp_unit != "Others")
saveRDS(results_df,
        file.path("transformed_data",
                  "fd_ind_time.rds"))


# 5 - Plot variation through time ==============================================


# Pisco:
results_df_Coastal <- results_df %>% 
  dplyr::filter(project %in% c("COASTAL_CEN", "COASTAL_SOUTH"))

fric_plot_Coastal <- ggplot2::ggplot(data = results_df_Coastal,
                             ggplot2::aes(x = year, y = fric,
                                          colour = project)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Richness - Coastal")
fric_plot_Coastal

# Pisco Central
results_df_CoastalCEN <- results_df %>% 
  dplyr::filter(project %in% c("COASTAL_CEN"))
fric_plot_CoastalCEN <- ggplot2::ggplot(data = results_df_CoastalCEN,
                                     ggplot2::aes(x = year, y = fric)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Richness - Coastal")
fric_plot_CoastalCEN

results_df_CoastalCEN_fide1 <- results_df %>% 
  dplyr::filter(project %in% c("COASTAL_CEN"))

fide1_plot_CoastalCEN <- ggplot2::ggplot(data = results_df_CoastalCEN,
                                        ggplot2::aes(x = year, y = fide_PC1)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC1 - Coastal CEN")
fide1_plot_CoastalCEN

fide2_plot_CoastalCEN <- ggplot2::ggplot(data = results_df_CoastalCEN,
                                         ggplot2::aes(x = year, y = fide_PC2)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC 2 - Coastal CEN")
fide2_plot_CoastalCEN

fide3_plot_CoastalCEN <- ggplot2::ggplot(data = results_df_CoastalCEN,
                                         ggplot2::aes(x = year, y = fide_PC3)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC 3 - Coastal CEN")
fide3_plot_CoastalCEN

# FCE
results_df_FCE <- results_df %>% 
  dplyr::filter(project %in% c("FCE"))

fric_plot_FCE <- ggplot2::ggplot(data = results_df_FCE,
                                        ggplot2::aes(x = year, y = fric)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Richness - FCE")
fric_plot_FCE

fide1_plot_FCE <- ggplot2::ggplot(data = results_df_FCE,
                                         ggplot2::aes(x = year, y = fide_PC1)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC1 - FCE")
fide1_plot_FCE

fide2_plot_FCE <- ggplot2::ggplot(data = results_df_FCE,
                                         ggplot2::aes(x = year, y = fide_PC2)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC2 - FCE")
fide2_plot_FCE

# MCR
results_df_MCR <- results_df %>% 
  dplyr::filter(project %in% c("MCR"))

fric_plot_MCR <- ggplot2::ggplot(data = results_df_MCR,
                                 ggplot2::aes(x = year, y = fric)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Richness - MCR")
fric_plot_MCR

fide1_plot_MCR <- ggplot2::ggplot(data = results_df_MCR,
                                  ggplot2::aes(x = year, y = fide_PC1)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC1 - MCR")
fide1_plot_MCR

fide2_plot_MCR <- ggplot2::ggplot(data = results_df_MCR,
                                  ggplot2::aes(x = year, y = fide_PC2)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC2 - MCR")
fide2_plot_MCR

# SBC
results_df_SBC <- results_df %>% 
  dplyr::filter(project %in% c("SBC"))

fric_plot_SBC <- ggplot2::ggplot(data = results_df_SBC,
                                 ggplot2::aes(x = year, y = fric)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Richness - SBC")
fric_plot_SBC

fide1_plot_SBC <- ggplot2::ggplot(data = results_df_SBC,
                                  ggplot2::aes(x = year, y = fide_PC1)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC1 - SBC")
fide1_plot_SBC

fide2_plot_SBC <- ggplot2::ggplot(data = results_df_SBC,
                                  ggplot2::aes(x = year, y = fide_PC2)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC2 - SBC")
fide2_plot_SBC

# VCR
results_df_VCR <- results_df %>% 
  dplyr::filter(project %in% c("VCR"))

fric_plot_VCR <- ggplot2::ggplot(data = results_df_VCR,
                                 ggplot2::aes(x = year, y = fric)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Richness - VCR")
fric_plot_VCR

fide1_plot_VCR <- ggplot2::ggplot(data = results_df_VCR,
                                  ggplot2::aes(x = year, y = fide_PC1)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC1 - VCR")
fide1_plot_VCR

fide2_plot_VCR <- ggplot2::ggplot(data = results_df_VCR,
                                  ggplot2::aes(x = year, y = fide_PC2)) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = TRUE) +
  ggplot2::facet_wrap(~ system, scales = "free_y") +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Functional Identity PC2 - VCR")
fide2_plot_VCR


# 6 - Do FD in change through time ? ===========================================






