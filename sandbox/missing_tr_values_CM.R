################################################################################
##
## Explore missing traits: how much can we improve our dataset?
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


# Traits data (not imputed - raw):
raw_sp_tr_df <- read.csv(file.path("Data", "traits_tidy-data", "12_traits_wrangled-copy.csv"))

# Call the sp dataset
spp_master <- readr::read_csv(file.path("Data","species_tidy-data", "23_species_master-spp-list-copy.csv")) |>
  janitor::clean_names()

# Load imputed traits for comparison (cleaned version):
sp_tr_zscore <- readRDS(file.path("transformed_data",
                                  "sp_tr_zscore_new.rds"))


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


# 3 - Focus on fish and look for TL and lifespan ===============================


# Trophic Level: I have seen that some TL where existing in Fishbase but not our dataset
# Exploring this:
# https://www.fishbase.se/manual/key%20facts.htm
# https://www.fishbase.us/manual/English/FishbaseThe_FOOD_ITEMS_table.htm

fish <- c("Abudefduf sexfasciatus", "Abudefduf bengalensis")
rfishbase::fb_tbl("species") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% fish)

# Trophic level estimates based on models - example:
ex <- rfishbase::estimate(species_list = fish,
                    server = c("fishbase", "sealifebase"),
                    version = "latest")

# Try for all our fish species with NA values:
all_fish_sp_na <- sp_tr_fish$scientific_name[which(is.na(sp_tr_fish$tr.trophic.level.num))]
fish_TL <- rfishbase::estimate(species_list = all_fish_sp_na,
                               server = c("fishbase", "sealifebase"),
                               version = "latest")
saveRDS(fish_TL, 
        file.path("transformed_data",
                  "fish_TL_estimate_fishbase.rds")) # It also has a column "AgeMax" -> lifespan

sp_tr_fish_updated <- sp_tr_fish %>% 
  dplyr::rename(Species = "scientific_name") %>% 
  dplyr::left_join(fish_TL[, c("Species", "Troph", "AgeMax")]) %>% 
  dplyr::rename(scientific_name = "Species")
  
sp_tr_fish_updated2 <- sp_tr_fish_updated %>% 
  dplyr::mutate(new.tr.trophic.level.num = dplyr::coalesce(tr.trophic.level.num, Troph),
                new.tr.age.years = dplyr::coalesce(tr.age.years, AgeMax)) 

# Look at the change in proportion of NAs:
na_prop_fish_updated <- sp_tr_fish_updated2 %>%
  dplyr::summarise(across(where(is.numeric), ~ mean(is.na(.))))

# Take these new traits:
sp_tr_fish_updated3 <- sp_tr_fish_updated2 %>% 
  dplyr::select(-c("Troph", "AgeMax"))
hist(sp_tr_fish_updated3$new.tr.trophic.level.num, breaks = 5)
hist(sp_tr_zscore$tr.trophic.level.num, breaks = 5)


