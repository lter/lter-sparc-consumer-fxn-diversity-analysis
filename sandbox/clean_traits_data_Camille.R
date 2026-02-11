################################################################################
##
## Script to get traits data clean and z-transformed
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
librarian::shelf(tidyverse, dplyr, funbiogeo, ggplot2)

# Get set up
source("00_setup.R")

# Clear environment & collect garbage
rm(list = ls()); gc()


# 1 - Load data ================================================================

# Call the new imputed data after download from Drive (ask Shalanda's ok):
#imp_tr_df <- read.csv(file.path("Data", "traits_tidy-data", "consumer-trait-species-imputed-taxonmic-database.csv"))
# But does not work - don't have the last version don't know why
imp_tr_df <- read.csv(file.path("Data", "traits_tidy-data", "consumer-trait-species-imputed-taxonmic-database.csv"))


# Call the sp dataset
sp_list <- read.csv(file.path("Data", "species_tidy-data", "23_species_master-spp-list.csv"))


# 2 - Create a dataframe with species/hab/project/taxa =========================


sp_list_v2 <- sp_list %>%
  dplyr::mutate(across(everything(), na_if, y = ""))

sp_list_v2$sp.proj <- paste(sp_list_v2$scientific_name, sp_list_v2$project, sep=".")

sp_list_ready <- sp_list_v2 %>% distinct(sp.proj, .keep_all = T)



## add in a taxa column

# first make a taxon lookup table per project.
# we're gonna skip PIE (no traits)
# KBS_INS - insects
# KONZA - insects
project.taxon <- data.frame(project = unique(sp_list_ready$project))
project.taxon$taxa <- NA
project.taxon$taxa[which(project.taxon$project %in% c("NGA","Arctic","Palmer","CCE", "NorthLakes"))] <- "Zooplankton"
project.taxon$taxa[which(project.taxon$project %in% c("CoastalCA","FCE","SBC","MCR","VCR","RLS","FISHGLOB"))] <- "Fish"
project.taxon$taxa[which(project.taxon$project %in% c("KBS_MAM","SEV","MOHONK_MAM"))] <- "Mammals"
project.taxon$taxa[which(project.taxon$project %in% c("MOHONK_BIR","KBS_BIR","SBC_BEACH","HARVARD"))] <- "Birds"
project.taxon$taxa[which(project.taxon$project %in% c("KBS_AMP"))] <- "Amphibians"
project.taxon$taxa[which(project.taxon$project %in% c("KONZA","PIE","KBS_INS","MOHONK"))] <- "BAD"
# 
# project.taxon$kingdom <- NA
# project.taxon$kingdom[which(project.taxon$taxa)]

sp_list_ready$taxa <- project.taxon$taxa[match(sp_list_ready$project, project.taxon$project)]



# 3 - Cleaning non-focal taxa  =========================


# Get the classes of our Fish projects:
sp_list_fish_marine <- sp_list_ready %>% 
  dplyr::filter(project %in% c("CoastalCA","FCE","SBC","MCR","VCR","RLS","FISHGLOB"))
unique(sp_list_fish_marine$class)

# Keep: "Actinopterygii", "Teleostei", "Elasmobranchii", "Chondrichthyes"
sp_list_fish_marine_corrected <- sp_list_fish_marine %>% 
  dplyr::filter(class %in% c("Actinopterygii", "Teleostei", "Elasmobranchii", 
                               "Chondrichthyes"))
  
sp_list_ready_corrected <- sp_list_ready %>% 
  dplyr::filter(!project %in% c("CoastalCA","FCE","SBC","MCR","VCR","RLS","FISHGLOB")) %>% 
  dplyr::bind_rows(sp_list_fish_marine_corrected)

# Save it:
saveRDS(sp_list_ready_corrected,
        file.path("transformed_data", "species_list_corrected_fish.rds"))


#Join trait data with program species list 
program_sp_trt_data <- dplyr::left_join(sp_list_ready_corrected, imp_tr_df,
                                        by="scientific_name") %>%
  dplyr::mutate(
    order = coalesce(order.x, order.y),
    family = coalesce(family.x, family.y),
    genus = coalesce(genus.x, genus.y)) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y")) %>%
  dplyr::relocate(order, family, genus, .before = sex) %>%
  dplyr::relocate(source, .before = scientific_name) 
#some source info absent need to fix this later!!! 

program_sp_trt_data$age_life.span_years[which(program_sp_trt_data$age_life.span_years<=0)] <- NA



# Create a new column for Taxa:
program_sp_trt_data$taxa <- project.taxon$taxa[match(program_sp_trt_data$project, project.taxon$project)]

# program_sp_trt_data1 <- program_sp_trt_data %>% 
#   mutate(taxa = case_when(
#     project == "CoastalCA" ~ "Fish",
#     project == "FCE" ~ "Fish",
#     project == "SBC" ~ "Fish",
#     project == "MCR" ~ "Fish",
#     project == "VCR" ~ "Fish",
#     project == "RLS" ~ "Fish",
#     project == "FISHGLOB" ~ "Fish",
#     project == "KBS_MAM" ~ "Mammals",
#     project == "SEV" ~ "Mammals",
#     project == "MOHONK" ~ "Amphibians",
#     project == "KBS_AMP" ~ "Amphibians",
#     project %in% c("HARVARD", "KBS_BIR","SBC_BEACH") ~ "Birds",
#     project %in% c("NGA","Arctic","Palmer","CCE","NorthLakes") ~ "Zooplankton",
#     TRUE ~ NA
#     )) %>% filter(!is.na(taxa), scientific_name!= "Homo sapiens") %>% 
#   dplyr::mutate(taxa = factor(taxa))


# all_traits <- program_sp_trt_data %>% 
#   dplyr::select(c(1:9, 
#                 "age_life.span_years",
#                 "diet_trophic.level_num",
#                 "reproduction_reproductive.rate_num.offspring.per.year",
#                 "mass_adult_g",
#                 "active.time_category_ordinal", "taxa")) %>%
#   group_by(project) %>% mutate(tr.age.zp = scale(age_life.span_years)[,1],
#                                tr.trophic.level.zp = scale(diet_trophic.level_num)[,1],
#                                tr.reproductive.rate.zp = scale(reproduction_reproductive.rate_num.offspring.per.year)[,1],
#                                tr.mass.adult.zp = scale(log(mass_adult_g, 10))[,1]
#                                )

all_traits <- program_sp_trt_data %>% filter(scientific_name !="Homo sapiens", taxa !="BAD") %>%
  dplyr::select(c(1:9, 
                  "age_life.span_years",
                  "length_adult_cm",
                  "diet_trophic.level_num",
                  "reproduction_reproductive.rate_num.offspring.per.year",
                  "reproduction_fecundity_num",
                  "mass_adult_g",
                  "active.time_category_ordinal", "taxa")) %>%
  dplyr::mutate(age_life.span_years = if_else(age_life.span_years <= 0,
                                              NA, age_life.span_years)) %>% 
  dplyr::mutate(reproduction_fecundity_num = if_else(reproduction_fecundity_num < 0, 
                                                     NA, reproduction_fecundity_num)) %>% 
  group_by(project) %>% 
  mutate(tr.age.zp = scale(age_life.span_years)[,1],
         tr.trophic.level.zp = scale(diet_trophic.level_num)[,1],
         tr.reproductive.rate.zp = scale(reproduction_reproductive.rate_num.offspring.per.year)[,1],
         tr.fecundity.zp = scale(reproduction_fecundity_num)[,1],
         tr.mass.adult.zp = scale(log(mass_adult_g, 10))[,1],
         tr.length.adult.zp = scale(log(length_adult_cm, 10))[,1]
  ) %>%
  ungroup() %>% group_by(taxa) %>%
  mutate(tr.age.zt = scale(age_life.span_years)[,1],
         tr.trophic.level.zt = scale(diet_trophic.level_num)[,1],
         tr.reproductive.rate.zt = scale(reproduction_reproductive.rate_num.offspring.per.year)[,1],
         tr.fecundity.zt = scale(reproduction_fecundity_num)[,1],
         tr.mass.adult.zt = scale(log(mass_adult_g, 10))[,1],
         tr.length.adult.zt = scale(log(length_adult_cm, 10))[,1]
  )

# Order and create the Active time categories:


# all_traits$tr.active.time <- all_traits %>% ungroup() %>% 
#   dplyr::mutate(active.time_category_new = case_when(
#     active.time_category_ordinal == "diurnal" ~ "diurnal",
#     active.time_category_ordinal == "nocturnal" ~ "nocturnal",
#     ! active.time_category_ordinal %in% c("diurnal", "nocturnal", "") ~ "cathemeral")) %>% 
#   dplyr::select(active.time_category_new) %>% 
#   dplyr::rename(tr.active.time = active.time_category_new)

# ----


saveRDS(all_traits,
        file.path("transformed_data",
                  "sp_tr_zscore.rds"))







# replaced with all-in-one code above. only creates sp_tr_zscore.rds 
# 
# 
# # Add taxa to species list data:
# # Create a new column for Taxa:
# taxa_sp_list <- sp_list %>% 
#   mutate(taxa = case_when(
#     project == "CoastalCA" ~ "Fish",
#     project == "FCE" ~ "Fish",
#     project == "SBC" ~ "Fish",
#     project == "MCR" ~ "Fish",
#     project == "VCR" ~ "Fish",
#     project == "RLS" ~ "Fish",
#     project == "FISHGLOB" ~ "Fish",
#     project == "KBS_MAM" ~ "Mammals",
#     project == "SEV" ~ "Mammals",
#     project == "MOHONK" ~ "Amphibians",
#     project == "KBS_AMP" ~ "Amphibians",
#     project %in% c("HARVARD", "KBS_BIR","SBC_BEACH") ~ "Birds",
#     project %in% c("NGA","Arctic","Palmer","CCE","NorthLakes") ~ "Zooplankton")) %>% 
#   dplyr::mutate(taxa = factor(taxa))
# 
# # Check that the NA are for insect projects and PIE and rm them:
# unique(taxa_sp_list$project[which(is.na(taxa_sp_list$taxa))])
# taxa_sp_list_final <- taxa_sp_list %>% 
#   dplyr::filter(! is.na(taxa)) %>% 
#   dplyr::select(c("project", "habitat",
#                   "scientific_name", "taxa")) %>% 
#   distinct()
# 
# saveRDS(taxa_sp_list_final,
#         file.path("transformed_data", "proj_taxa_sp_list.rds"))
# 
# 
# 4 - Clean the traits data and get completedness ==============================
# 
# 
# # Keep studied traits, rm non possible values, Homo sapiens and duplicate values
# sp_tr1 <- imp_tr_df %>% 
#   dplyr::distinct() %>% 
#   dplyr::filter(scientific_name != "Homo sapiens") %>% 
#   dplyr::select(c("scientific_name",
#                   "age_life.span_years",
#                   "diet_trophic.level_num",
#                   "reproduction_reproductive.rate_num.offspring.per.year",
#                   "reproduction_fecundity_num",
#                   "length_adult_cm",
#                   "mass_adult_g")) %>% 
#   dplyr::mutate(age_life.span_years = if_else(age_life.span_years <= 0, 
#                                               NA, age_life.span_years)) %>% 
#   dplyr::mutate(reproduction_fecundity_num = if_else(reproduction_fecundity_num < 0, 
#                                               NA, reproduction_fecundity_num))
# 
# # Build one df per Taxa and get traits completedness:
# 
# # 
# 
# 
# 
# # 3 - Create a new reproduction and size column ================================
# 
# # For the two reproduction columns - chose the dominant trait:
# # Fish - fecundity
# # Mammals - repro rate
# # Amphibians - repro rate
# # Birds - repro rate
# # Zooplankton - fecundity
# 
# # Keep only data for the chosen traits for each taxa:
# sp_tr2 <- sp_tr1 %>% 
#   dplyr::left_join(taxa_sp_list_final, by = "scientific_name") %>% 
#   dplyr::mutate("reproduction_reproductive.rate_num.offspring.per.year" = if_else(taxa %in% c("Fish", "Zooplankton"),
#                                                                                   NA,   
#                                                                                   reproduction_reproductive.rate_num.offspring.per.year)) %>% 
#   dplyr::mutate("reproduction_fecundity_num" = if_else(taxa %in% c("Amphibians",
#                                                                    "Mammals",
#                                                                    "Birds"),
#                                                        NA,
#                                                        reproduction_fecundity_num)) %>% 
#   dplyr::select(c("scientific_name",
#                   "taxa",
#                   "age_life.span_years",
#                   "diet_trophic.level_num",
#                   "reproduction_reproductive.rate_num.offspring.per.year",
#                   "reproduction_fecundity_num",
#                   "length_adult_cm",
#                   "mass_adult_g")) %>% 
#   dplyr::distinct()
# 
# # Some species have taxa == NA, insect of species present only in PIE: rm
# sp_tr3 <- sp_tr2 %>% 
#   dplyr::filter(! is.na(taxa))
# 
# # ADD SIZE COLUMN HERE !!!!!
# 
# # Some species are in the list but not in the species traits df
# setdiff(taxa_sp_list_final$scientific_name,
#         sp_tr3$scientific_name)
# # All the genus ones: ok!
# 
# # Number of species to work on: (birds, mammals, fish, zooplankton, aphibians)
# length(unique(sp_tr3$scientific_name))
# 
# # Save it
# saveRDS(sp_tr3, file.path("transformed_data",
#                           "sp_tr.rds"))
# 
# 
# # 5 - Transform traits =========================================================
# 
# 
# # Z-transform the traits and log when necessary:
# sp_tr4 <- taxa_sp_list_final %>% 
#   dplyr::left_join(sp_tr3, by = "scientific_name") %>% 
#   dplyr::group_by(project) %>% 
#   dplyr::mutate(tr.age.zp = scale(age_life.span_years)[, 1],
#                 tr.trophic.level.zp = scale(diet_trophic.level_num)[,1],
#                 tr.reproductive.rate.zp = scale(reproduction_reproductive.rate_num.offspring.per.year)[,1],
#                 tr.fecundity.zp = scale(reproduction_fecundity_num)[,1],
#                 tr.mass.adult.zp = scale(log(mass_adult_g, 10))[,1],
#                 tr.length.adult.zp = scale(log(length_adult_cm, 10))[,1]) %>% 
#   dplyr::distinct() %>% 
#   dplyr::select(c("taxa.x", "scientific_name",
#                   "tr.age.zp", "tr.trophic.level.zp",
#                   "tr.reproductive.rate.zp",
#                   "tr.fecundity.zp",
#                   "tr.mass.adult.zp",
#                   "tr.length.adult.zp")) %>% 
#   dplyr::rename(taxa = "taxa.x")
#   
# # Add two reprod col together and save it:
# sp_tr5 <- sp_tr4 %>%
#   dplyr::mutate(tr.reproduction.zp = dplyr::coalesce(tr.reproductive.rate.zp, 
#                                                      tr.fecundity.zp)) %>%
#   dplyr::select(-c(tr.reproductive.rate.zp, tr.fecundity.zp))
# 
# # Do: tr completion
# # NOTE: I HAVE NAN values: check when they started to show
# # NOTE: can't synchronise to drive
# 
# saveRDS(sp_tr5, file.path("transformed_data",
#                           "sp_tr_zscore.rds"))
# 
