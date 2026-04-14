################################################################################
##
## Script to test traits imputation based on the phylogeny
##
## Camille Magneville and Shalanda Grier 
##
## 03/2026
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


# Traits data (not imputed - raw): @SHALANDA/ Should I use -copy or not copy? SG- @CAMILLW you can use -copy. Most updated for now
raw_sp_tr_df <- read.csv(file.path("Data", "traits_tidy-data", "12_traits_wrangled.csv"))

# Call the sp dataset
spp_master <- readr::read_csv(file.path("Data","species_tidy-data", "23_species_master-spp-list-copy.csv")) |>
  janitor::clean_names()


# 2 - Clean trait data =========================================================


# Only keep the studied traits and focus on birds/fish:


# Define the classes we keep:
bird_classes <- c("Aves")
fish_classes <- c(
  "Actinopterygii", "Actinopteri", "Teleostei",
  "Chondrichthyes", "Elasmobranchii", "Holocephali",
  "Chondrostei", "Holostei",
  "Myxini", "Petromyzonti"
)


dt <- spp_master %>%
  dplyr::mutate(
    taxon = case_when(
      class %in% bird_classes ~ "Birds",
      class %in% fish_classes ~ "Fish",
      TRUE ~ "others"))


# List the project we keep:
proj_allowed <- list(
  "Arctic"      = c("Birds"),
  "Palmer"      = c("Birds"),
  "FISHGLOB"    = c("Fish"),
  "FCE"         = c("Fish"),
  "SBC"         = c("Fish"),
  "SBC_BEACH"   = c("Birds"),
  "CoastalCA"   = c("Fish"),
  "MCR"         = c("Fish"),
  "VCR"         = c("Fish"),
  "KBS_BIR"     = c("Birds"),
  "HARVARD"     = c("Birds"),
  "MOHONK_BIR"  = c("Birds"),
  "PIE"         = c("Fish"),
  "RLS"         = c("Fish"),
  "AndrewsForest"  = c("Birds"),
  "Luquillo"  = c("Birds"),
  "HubbardBrook"  = c("Birds"),
  "CAP"  = c("Birds"),
  "SONGS" = c("Birds"),
  "Baltimore" = c("Birds"),
  "CCE_birds" = c("Birds"),
  "PLUM" = c("Birds")
)

dt1 <- dt |>
  ### looks up which taxon groups are allowed for each project
  mutate(allowed = map(project, ~ proj_allowed[[.x]])) |>
  ### filters those not allowed
  filter(!map_lgl(allowed, is.null)) |>
  ### keep only rows that meet allowed groups for each project
  filter(map2_lgl(taxon, allowed, ~ .x %in% .y)) |>
  ### clean up dataframe
  select(-allowed)

# rows removed by the class-level cleaning
removed <- anti_join(dt, dt1)
glimpse(removed)


# add back in projects that had NA classes but were clean to begin with (focused only on birds, fish when checking)
acceptable <- removed |> 
  filter(project %in% c('KBS_BIR', 'HARVARD',
                        'CoastalCA', "MOHONK_BIR")) |> 
  mutate(taxon = case_when(
    project == 'KBS_BIR' ~ "Birds",
    project == 'HARVARD' ~ "Birds",
    project == 'CoastalCA' ~ "Fish",
    project == 'MOHONK_BIR' ~ "Birds",
    TRUE ~ NA_character_
  ))


## Make a clean spp list 
sp_list <- rbind(dt1, acceptable)

# change NAs to be blanks
sp_list_v2 <- sp_list %>%
  dplyr::mutate(across(everything(), na_if, y = ""))

# make a species.project column to keep species that occur in multiple projects
sp_list_v2$sp.proj <- paste(sp_list_v2$scientific_name, sp_list_v2$project, sep=".")

# remove duplicated species per project (came from multiple source files)
sp_list_ready <- sp_list_v2 %>% distinct(sp.proj, .keep_all = T)


# In the species traits data, we now have several lines per species (different sources):
# Keep only the traits of interest and keep the source which has most data:
sp_tr_df1 <- raw_sp_tr_df %>% 
  dplyr::select(c("class",                                                                 
                  "order",                                                                 
                  "family",                                                                
                  "genus",                                                                 
                  "scientific_name",                                                       
                  "source",
                  "tr.age.years" = "age_life.span_years",
                  "tr.trophic.level.num" = "diet_trophic.level_num",
                  "reproduction_reproductive.rate_num.offspring.per.year",
                  "reproduction_reproductive.rate_num.litter.or.clutch.per.year",
                  "reproduction_reproductive.rate_num.offspring.per.clutch.or.litter",
                  "reproduction_fecundity_num",
                  "tr.mass.adult.g" = "mass_adult_g",
                  "length_max_cm"))

# Link it to the taxon and keep only birds and fish:
sp_tr_df2 <- sp_tr_df1 %>% 
  dplyr::left_join(sp_list_ready[, c("scientific_name", "taxon")],
                   by = "scientific_name") %>% 
  dplyr::filter(taxon %in% c("Birds", "Fish")) %>% 
  dplyr::distinct()


# Complete reproduction trait and create a unified one:
missing.offspring <- which(is.na(sp_tr_df2$reproduction_reproductive.rate_num.offspring.per.year))
sp_tr_df2$reproduction_reproductive.rate_num.offspring.per.year[missing.offspring] <- sp_tr_df2$reproduction_reproductive.rate_num.litter.or.clutch.per.year[missing.offspring] * sp_tr_df2$reproduction_reproductive.rate_num.offspring.per.clutch.or.litter[missing.offspring] 
sp_tr_df3 <- sp_tr_df2 %>% 
  dplyr::mutate(tr.reproduction.unified = case_when(
    taxon %in% c("Fish") ~ reproduction_fecundity_num,
    taxon %in% c("Birds") ~ reproduction_reproductive.rate_num.offspring.per.year,
    T ~ NA
  )) %>% 
  dplyr::mutate(tr.age.years = if_else(tr.age.years <= 0,
                                       NA, tr.age.years)) %>% 
  dplyr::mutate(reproduction_fecundity_num = if_else(reproduction_fecundity_num < 0, 
                                                     NA, reproduction_fecundity_num)) %>% 
  dplyr::select(-c("reproduction_reproductive.rate_num.offspring.per.year",
                   "reproduction_reproductive.rate_num.litter.or.clutch.per.year" ,
                   "reproduction_reproductive.rate_num.offspring.per.clutch.or.litter",
                   "reproduction_fecundity_num"))

# Gather info for each species - one row:
sp_tr_df4 <- sp_tr_df3 %>%
  dplyr::group_by(scientific_name, taxon) %>%
  dplyr::select(-c("source")) %>% 
  dplyr::summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop") %>% 
  dplyr::mutate(across(where(is.numeric), ~ replace(.x, is.nan(.x), NA)))


# How many fish species will have information on mass if we can compute it from length?
sp_tr_fish <- sp_tr_df4[which(sp_tr_df4$taxon == "Fish"), ]
length(unique(sp_tr_fish$scientific_name[which(is.na(sp_tr_fish$tr.mass.adult.g) & !is.na(sp_tr_fish$length_max_cm))]))


# Count the proprtion of NAs:
na_prop_fish <- sp_tr_fish %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(is.na(.))))
sp_tr_birds <- sp_tr_df4[which(sp_tr_df4$taxon == "Birds"), ]
na_prop_birds <- sp_tr_birds %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(is.na(.))))


# 3 - Check if we can get fish lifespan info through fishbase ==================


# Try for all our fish species with NA values:
all_fish_sp_na <- sp_tr_fish$scientific_name[which(is.na(sp_tr_fish$tr.age.years))]
fish_LS <- rfishbase::estimate(species_list = all_fish_sp_na,
                               server = c("fishbase", "sealifebase"),
                               version = "latest")

# A few AgeMax values available:
sp_tr_fish_updated <- sp_tr_fish %>% 
  dplyr::rename(Species = "scientific_name") %>% 
  dplyr::left_join(fish_LS[, c("Species", "AgeMax")]) %>% 
  dplyr::rename(scientific_name = "Species") %>% 
  dplyr::mutate(new.tr.age.years = dplyr::coalesce(tr.age.years, AgeMax))  %>% 
  dplyr::select(-c("AgeMax"))

# How did it change NA proportions?
na_prop_fish_updated <- sp_tr_fish_updated %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(is.na(.))))


# 4 - Remove genus level data: how much missing? ===============================


# Remove the "species" which are only Genus level 
sp_tr_fish_updated2 <- sp_tr_fish_updated %>% 
  dplyr::filter(stringr::str_detect(scientific_name, "^\\S+\\s+\\S+")) 

# How did it change NA proportions?
na_prop_fish_updated2 <- sp_tr_fish_updated2 %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(is.na(.))))

# CONCLUSION ON APRIL 9th:
# ... big problem with lifespan (65% NA), reproduction (83% NA),
# ... mass (47% NA but will be improved with length data)
# ... ok: trophic level (2% NA)


# 5 - Phylogenetic imputation Fish - Step 1: Get the phylogenetic tree =========


# Getting phylogeny for as many species as possible: 693 sp/923
fishtree <- fishtree::fishtree_phylogeny(species = unique(sp_tr_fish_updated2$scientific_name))

# Get complete phylogenies (missing taxa placed according to the highest level of tax resol) :
fishtree_complete <- fishtree::fishtree_complete_phylogeny(species = unique(sp_tr_fish_updated2$scientific_name)) 
# Still 74 species missing in these phylogenies (100 phylogenies): 849 species in them
sp_phylo <- gsub("_", " ", fishtree_complete[[1]]$tip.label)
missing_sp_tree <- setdiff(unique(sp_tr_fish_updated2$scientific_name),
                           sp_phylo)



# 6 - Phylogenetic imputation Fish - Step 2: Impute traits =====================


# Note: 74 species are not present in the phylogeny -> remove them from the trait dataset:
sp_tr_fish_phylo <- sp_tr_fish_updated2 %>% 
  dplyr::filter(! scientific_name %in% missing_sp_tree) %>% 
  dplyr::select(-c("taxon", "tr.age.years")) 
sp_tr_fish_phylo$scientific_name <- gsub(" ", "_", sp_tr_fish_phylo$scientific_name)
sp_tr_fish_phylo2 <- tibble::column_to_rownames(sp_tr_fish_phylo,
                                                "scientific_name")

# Restrict the traits data to the one with data for each trait:
trophic_level_tr <- sp_tr_fish_phylo2$tr.trophic.level.num[!is.na(sp_tr_fish_phylo2$tr.trophic.level.num)]
names(trophic_level_tr) <- rownames(sp_tr_fish_phylo2[!is.na(sp_tr_fish_phylo2$tr.trophic.level.num), ])

mass_tr <- sp_tr_fish_phylo2$tr.mass.adult.g[!is.na(sp_tr_fish_phylo2$tr.mass.adult.g)]
names(mass_tr) <- rownames(sp_tr_fish_phylo2[!is.na(sp_tr_fish_phylo2$tr.mass.adult.g), ])

reprod_tr <- sp_tr_fish_phylo2$tr.reproduction.unified[!is.na(sp_tr_fish_phylo2$tr.reproduction.unified)]
names(reprod_tr) <- rownames(sp_tr_fish_phylo2[!is.na(sp_tr_fish_phylo2$tr.reproduction.unified), ])

lifespan_tr <- sp_tr_fish_phylo2$new.tr.age.years[!is.na(sp_tr_fish_phylo2$new.tr.age.years)]
names(lifespan_tr) <- rownames(sp_tr_fish_phylo2[!is.na(sp_tr_fish_phylo2$new.tr.age.years), ])


# Try: Infer traits based on phylogeny 
try_trophic_level <- picante::phyEstimate(phy = fishtree_complete[[1]],
                     trait = trophic_level_tr,
                     method = "pic")
try_mass <- picante::phyEstimate(phy = fishtree_complete[[1]],
                                 trait = mass_tr,
                                 method = "pic")



# Compute mean and variance over the 100 trees 

## Define a fct to do so:
impute.trait <- function(fishtree,
                         traits,
                         trait_name) {
  
  # Final dataset:
  imputed_df <- as.data.frame(matrix(ncol = 7, nrow = 1, NA))
  colnames(imputed_df) <- c("species", "trait", "tree", "est", "se",
                            "mean", "var")
  
  # Loop on each tree:
  for (i in c(1:length(fishtree))) {
    
    print(i)
    
    # Get imputed values:
    imp_tr <- picante::phyEstimate(phy = fishtree[[i]],
                                   trait = traits,
                                   method = "pic")
    
    # Build a dataframe:
    imputed_tree_df <- imp_tr %>% 
      tibble::rownames_to_column(var = "species") %>% 
      dplyr::mutate(trait = rep(trait_name, nrow(imp_tr))) %>% 
      dplyr::mutate(tree = rep(paste0("tree", sep = "_", i), nrow(imp_tr))) %>% 
      dplyr::rename(est = estimate) %>% 
      dplyr::mutate(mean = rep(NA, nrow(imp_tr)),
                    var = rep(NA, nrow(imp_tr))) %>% 
      dplyr::select( c("species", "trait", "tree", "est", "se",
                       "mean", "var"))
    
    # Add it to the final one:
    imputed_df <- rbind(imputed_df,
                        imputed_tree_df)

  }
  
  # Compute the mean and var for each species:
  imputed_final_df <- imputed_df %>% 
    dplyr::group_by(species) %>% 
    dplyr::mutate(mean = mean(est, na.rm = TRUE),
                  var  = stats::sd(est, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  # Remove the first row:
  imputed_final_df <- imputed_final_df[-1, ]
  
  return(imputed_final_df)

}

## Impute trait for body mass:
imp_tr_mass_df <- impute.trait(fishtree = fishtree_complete,
                               traits = mass_tr,
                               trait_name = "imp.tr.mass.adult.g")

## Impute trait for reproduction:
imp_tr_reprod_df <- impute.trait(fishtree = fishtree_complete,
                               traits = reprod_tr,
                               trait_name = "imp.tr.reproduction.unified")

## Impute trait for trophic level:
imp_tr_TL_df <- impute.trait(fishtree = fishtree_complete,
                               traits = trophic_level_tr,
                               trait_name = "imp.tr.trophic.level.num")

## Impute trait for lifespan:
imp_tr_LS_df <- impute.trait(fishtree = fishtree_complete,
                               traits = lifespan_tr,
                               trait_name = "imp.new.tr.age.years")


## Bind all imputed traits together:
imputed_tr_all_df <- rbind(imp_tr_mass_df,
                           imp_tr_reprod_df,
                           imp_tr_TL_df,
                           imp_tr_LS_df)
saveRDS(imputed_tr_all_df,
        here::here("transformed_data",
                   "phylo_imputed_traits_df.rds"))


# 7 - Phylogenetic imputation Fish - Step 3: Check extreme values ==============


imputed_tr_all_df <- readRDS(here::here("transformed_data",
                                         "phylo_imputed_traits_df.rds"))

# Get species with extreme trait values:
extreme_sp_df <- imputed_tr_all_df %>%
  dplyr::group_by(trait) %>%
  dplyr::mutate(q01 = stats::quantile(mean, 0.01),
                q99 = stats::quantile(mean, 0.99)) %>%
  dplyr::filter(mean <= q01 | mean >= q99) %>%
  dplyr::select(species, trait, mean) %>%
  dplyr::mutate(mean = round(mean, 3)) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct()
