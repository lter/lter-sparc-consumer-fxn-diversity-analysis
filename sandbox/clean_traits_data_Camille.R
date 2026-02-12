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
librarian::shelf(tidyverse, dplyr, funbiogeo, ggplot2, janitor)


# Get set up
source("00_setup.R")

# Clear environment & collect garbage
rm(list = ls()); gc()


# 1 - Load data ================================================================


# Call the new imputed data after download from Drive (ask Shalanda's ok):
imp_tr_df <- read.csv(file.path("Data", "traits_tidy-data", "consumer-trait-species-imputed-taxonmic-database.csv"))


# Call the sp dataset
spp_master <- readr::read_csv(file.path("Data","species_tidy-data", "23_species_master-spp-list.csv")) |>
  janitor::clean_names()


# 2 - Create a dataframe with species/hab/project/taxa =========================

# . taxonomic cleaning ------------------------------------------------------------
  # this was stolen directly from MW_taxon-filtering-list.R on 2.10.20



### set simple workflow functions ---
nacheck <- function(df) {
  na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
  print(na_count_per_column)
}


bird_classes <- c("Aves")
mammal_classes <- c("Mammalia")
amphib_classes <- c("Amphibia")

fish_classes <- c(
  "Actinopterygii", "Teleostei",
  "Chondrichthyes", "Elasmobranchii", "Holocephali",
  "Chondrostei", "Holostei",
  "Myxini", "Petromyzonti"
)

zoop_projects <- c("Arctic", "NGA", "CCE", "NorthLakes")

zooplankton_classes <- spp_master |>
  filter(project %in% zoop_projects) |>
  distinct(class) |>
  pull(class) |>
  na.omit() |>
  unique() |>
  setdiff(c(
    fish_classes, bird_classes, mammal_classes, amphib_classes,
    "Insecta", "Arachnida", "Diplopoda"
  ))

dt <- spp_master |>
  mutate(
    taxon = case_when(
      class %in% bird_classes ~ "Birds",
      class %in% mammal_classes ~ "Mammals",
      class %in% amphib_classes ~ "Amphibians",
      class %in% fish_classes ~ "Fish",
      class == "Insecta" ~ "Insects",
      class %in% zooplankton_classes ~ "Zooplankton",
      phylum == "Chordata" ~ "other_chordate",
      TRUE ~ "other_invertebrate"
    )
  )

proj_allowed <- list(
  "Arctic"      = c("Zooplankton", "Birds"),
  "Palmer"      = c("Zooplankton", "Birds"),
  "NorthLakes"  = c("Zooplankton"),
  "FISHGLOB"    = c("Fish"),
  "FCE"         = c("Fish"),
  "SBC"         = c("Fish", "Birds"),
  "CoastalCA"   = c("Fish", "other_invertebrates"),
  "MCR"         = c("Fish"),
  "NGA"         = c("Zooplankton"),
  "VCR"         = c("Fish"),
  "KBS_AMP"     = c("Amphibians"),
  "KBS_BIR"     = c("Birds"),
  "KBS_INS"     = c("Insects", 'other'),
  "KBS_MAM"     = c("Mammals"),
  "HARVARD"     = c("Birds"),
  "SEV"         = c("Mammals"),
  "MOHONK"      = c("Birds", "Mammals"),
  "KONZA"       = c("Insects", "Mammals"),
  "PIE"         = c("Fish"),
  "RLS"         = c("Fish"),
  "SBC_BEACH"   = c('Birds')
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


# add back in projects that had NA classes but were clean to begin with
acceptable <- removed |> 
  filter(project %in% c('KBS_AMP', 'KBS_BIR', 'KBS_MAM', 'HARVARD',
                        'KONZA', 'SEV')) |> 
  mutate(taxon = case_when(
    project == 'KBS_AMP' ~ "Amphibians",
    project == 'KBS_BIR' ~ "Birds",
    project == 'KBS_MAM' ~ "Mammals",
    project == 'HARVARD' ~ "Birds",
    project == 'KONZA' ~ "Insects",
    project == 'SEV' ~ "Mammals",
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



# ## add in a taxa column
# 
# # first make a taxon lookup table per project.
# # we're gonna skip PIE (no traits)
# # KBS_INS - insects
# # KONZA - insects
# project.taxon <- data.frame(project = unique(sp_list_ready$project))
# project.taxon$taxa <- NA
# project.taxon$taxa[which(project.taxon$project %in% c("NGA","Arctic","Palmer","CCE", "NorthLakes"))] <- "Zooplankton"
# project.taxon$taxa[which(project.taxon$project %in% c("CoastalCA","FCE","SBC","MCR","VCR","RLS","FISHGLOB"))] <- "Fish"
# project.taxon$taxa[which(project.taxon$project %in% c("KBS_MAM","SEV","MOHONK_MAM"))] <- "Mammals"
# project.taxon$taxa[which(project.taxon$project %in% c("MOHONK_BIR","KBS_BIR","SBC_BEACH","HARVARD"))] <- "Birds"
# project.taxon$taxa[which(project.taxon$project %in% c("KBS_AMP"))] <- "Amphibians"
# project.taxon$taxa[which(project.taxon$project %in% c("KONZA","PIE","KBS_INS","MOHONK"))] <- "BAD"
# # 
# # project.taxon$kingdom <- NA
# # project.taxon$kingdom[which(project.taxon$taxa)]
# 
# sp_list_ready$taxa <- project.taxon$taxa[match(sp_list_ready$project, project.taxon$project)]



# 3 - Cleaning non-focal taxa  =========================


# #--- Clean Fish (may be deprecated if cleaned upstream in sp list)
# # Get the classes of our Fish projects:
# sp_list_fish_marine <- sp_list_ready %>% 
#   dplyr::filter(project %in% c("CoastalCA","FCE","SBC","MCR","VCR","RLS","FISHGLOB"))
# unique(sp_list_fish_marine$class)
# 
# # Keep: "Actinopterygii", "Teleostei", "Elasmobranchii", "Chondrichthyes"
# sp_list_fish_marine_corrected <- sp_list_fish_marine %>% 
#   dplyr::filter(class %in% c("Actinopterygii", "Teleostei", "Elasmobranchii", 
#                                "Chondrichthyes"))
#   
# sp_list_ready_corrected <- sp_list_ready %>% 
#   dplyr::filter(!project %in% c("CoastalCA","FCE","SBC","MCR","VCR","RLS","FISHGLOB")) %>% 
#   dplyr::bind_rows(sp_list_fish_marine_corrected)

# Save it:
saveRDS(sp_list_ready,
        file.path("transformed_data", "species_list_corrected_fish.rds"))
#--- end cleaning fish 


#Join trait data with program species list 
program_sp_trt_data <- dplyr::left_join(sp_list_ready, imp_tr_df,
                                        by="scientific_name") %>%
  dplyr::mutate(
    order = coalesce(order.x, order.y),
    family = coalesce(family.x, family.y),
    genus = coalesce(genus.x, genus.y)) %>%
  dplyr::select(-ends_with(".x"), -ends_with(".y")) %>%
  dplyr::relocate(order, family, genus, .before = sex) %>%
  dplyr::relocate(source, .before = scientific_name) 
#some source info absent need to fix this later!!! 

# 0 or negative lifespans are actually NAs
program_sp_trt_data$age_life.span_years[which(program_sp_trt_data$age_life.span_years<=0)] <- NA



# Create a new column for Taxa:
# program_sp_trt_data$taxa <- project.taxon$taxa[match(program_sp_trt_data$project, project.taxon$project)]

### INFILL MISSING num.offspring.per.year (may be deprecated)
# calculate reproduction where needed
missing.offspring <- which(is.na(program_sp_trt_data$reproduction_reproductive.rate_num.offspring.per.year))
program_sp_trt_data$reproduction_reproductive.rate_num.offspring.per.year[missing.offspring] <- program_sp_trt_data$reproduction_reproductive.rate_num.litter.or.clutch.per.year[missing.offspring] * program_sp_trt_data$reproduction_reproductive.rate_num.offspring.per.clutch.or.litter[missing.offspring] 


# determine whether fecundity or offspring number is best reproductive variable per project
# xtabs(~project, program_sp_trt_data[which(program_sp_trt_data$reproduction_fecundity_num>0),])
# 
# xtabs(~project, program_sp_trt_data[which(program_sp_trt_data$reproduction_reproductive.rate_num.offspring.per.year>0),])

# offspring: 
# KBS_AMPH

# fecundity:
# CoastalCA, FCE, KONZA, MCR, NGA, SBC, VCR, PIE, RLS

# both: CCE, KBS_MAM, MOHONK, Palmer, SEV
# NGA (use fecundity)? KONZA?
# Palmer, no real data in either. identical.
# KBS, SEV, Palmer, identical an none of them real at the moment.

# check whether we need to consider reproduction and num.offspring
# repro_proj <- program_sp_trt_data %>% group_by(project, taxa) %>% summarise(
#   n.true.offspring = length(unique(reproduction_reproductive.rate_num.offspring.per.year[which(!is.na(reproduction_reproductive.rate_num.offspring.per.year))])),
#   n.true.fecundity = length(unique(reproduction_fecundity_num[which(!is.na(reproduction_fecundity_num))])), n.fecundity = length(which(!is.na(reproduction_fecundity_num))),
#   n.mass = length(which(!is.na(mass_adult_g)))
#     )
# for Fish -> use fecundity
# for Birds -> use num.offspring
# for Zooplankton -> us XXXXXXXXX

# 4 - z-score standardize  =========================

# remove projects no longer in analysis:
BAD <- c("KONZA","PIE","KBS_INS","MOHONK", "KBS_MAM", "KBS_BIR", "KBS_AMP")

# rename and z-score standardize
all_traits <- program_sp_trt_data %>% filter(!project %in% BAD) %>%
  dplyr::select(c(1:12, 
                  "age_life.span_years",
                  "length_adult_cm",
                  "diet_trophic.level_num",
                  "reproduction_reproductive.rate_num.offspring.per.year",
                  "reproduction_fecundity_num",
                  "mass_adult_g",
                  "tr.active.time" ="active.time_category_ordinal", "taxon")) %>%
  dplyr::mutate(age_life.span_years = if_else(age_life.span_years <= 0,
                                              NA, age_life.span_years)) %>% 
  dplyr::mutate(reproduction_fecundity_num = if_else(reproduction_fecundity_num < 0, 
                                                     NA, reproduction_fecundity_num)) %>% 
  dplyr::mutate(reproduction.unified = case_when(
    taxon %in% c("Fish", "Zooplnakton") ~ reproduction_fecundity_num,
    taxon %in% c("Birds","Mammals", "Amphibians") ~ reproduction_reproductive.rate_num.offspring.per.year,
    T ~ NA
  )) %>% 
  group_by(project) %>% 
  mutate(tr.age.zp = scale(age_life.span_years)[,1],
         tr.trophic.level.zp = scale(diet_trophic.level_num)[,1],
         tr.reproductive.rate.zp = scale(reproduction_reproductive.rate_num.offspring.per.year)[,1],
         tr.fecundity.zp = scale(reproduction_fecundity_num)[,1],
         tr.reproduction.unified.zp = scale(log(reproduction.unified,10))[,1],
         tr.mass.adult.zp = scale(log(mass_adult_g, 10))[,1],
         tr.length.adult.zp = scale(log(length_adult_cm, 10))[,1]
  ) %>%
  ungroup() %>% group_by(taxon) %>%
  mutate(tr.age.zt = scale(age_life.span_years)[,1],
         tr.trophic.level.zt = scale(diet_trophic.level_num)[,1],
         tr.reproductive.rate.zt = scale(reproduction_reproductive.rate_num.offspring.per.year)[,1],
         tr.fecundity.zt = scale(reproduction_fecundity_num)[,1],
         tr.reproduction.unified.zt = scale(log10(reproduction.unified))[,1],
         tr.mass.adult.zt = scale(log(mass_adult_g, 10))[,1],
         tr.length.adult.zt = scale(log(length_adult_cm, 10))[,1]
  )









# 5 - Trait QA/QC  =========================

#### . check for too many non-unique values ============
  # would indicate a problem with data processing or insufficiently fine-scale taxonomic imputation
trait.checks <- all_traits %>% group_by(project, taxon) %>% summarise(
  n.spp = n(),
  n.spp.kingprob = length(which(is.na(kingdom))),
  n.spp.classprob = length(which(is.na(class))),
  n.mass = length(which(!is.na(tr.mass.adult.zt))),
  n.unique.mass = length(unique(tr.mass.adult.zt[which(!is.na(tr.mass.adult.zt))])),
  n.trophic = length(which(!is.na(tr.trophic.level.zt))),
  n.unique.trophic = length(unique(tr.trophic.level.zt[which(!is.na(tr.trophic.level.zt))])),
  n.age = length(which(!is.na(tr.age.zt))),
  n.unique.age = length(unique(tr.age.zt[which(!is.na(tr.age.zt))])),
  n.repro = length(which(!is.na(tr.reproduction.unified.zt))),
  n.unique.repro = length(unique(tr.reproduction.unified.zt[which(!is.na(tr.reproduction.unified.zt))])),
  )



library(ggplot2)
library(ggExtra)
library(rlang)

plot_with_margins <- function(data, x, y, title) {
  
  p <- ggplot(data, aes(x = !!ensym(x), y = !!ensym(y))) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    labs(
      x = as_label(ensym(x)),
      y = as_label(ensym(y)), 
      title=title
    ) 
  
  ggMarginal(
    p,
    type = "histogram",
    bins = 30,
    fill = "gray",
    color = "black"
  )
}



#### . Visualize each site with z-score standardized values. ====
currentversion <- paste("SiteTraitCheck",format(Sys.Date(), "%Y%m%d"), sep=".")
dir.create(paste("PrelimResults/DataQA-QC/", currentversion, sep=""))

# will save 
for(i in 1:length(unique(all_traits$project))){
  my.project <- unique(all_traits$project)[i]
  proj.dat <- all_traits %>% filter(project == my.project)
  
  if(length(which(!is.na(proj.dat$tr.age.zt)))>2 & length(which(!is.na(proj.dat$tr.mass.adult.zt)))>2 ){
    p1 <- plot_with_margins(proj.dat, tr.mass.adult.zt, tr.age.zt, my.project)
  }
  if(length(which(!is.na(proj.dat$tr.age.zt)))<3 | length(which(!is.na(proj.dat$tr.mass.adult.zt)))<3 ){
    p1 <- ggplot() + labs(x="not enough age or mass data", title=my.project) + theme_bw() 
    #print(p1)
  }
  if(length(which(!is.na(proj.dat$tr.trophic.level.zt)))>2 & length(which(!is.na(proj.dat$tr.reproduction.unified.zt)))>2 ){
    p2<- plot_with_margins(proj.dat, tr.trophic.level.zt, tr.reproduction.unified.zt, my.project)
    #print(p2)
  }
  if(length(which(!is.na(proj.dat$tr.trophic.level.zt)))<3 | length(which(!is.na(proj.dat$tr.reproduction.unified.zt)))<3 ){
    p2 <- ggplot() + labs(x="not enough trophic or reproduction data", title=my.project) + theme_bw()
    #print(p2)
  }
  ggsave(plot=p1,filename=paste(my.project,"zMass-zAge.jpg", sep="_"),path = file.path("PrelimResults","DataQA-QC",currentversion,"MassAge_zscore"), create.dir = TRUE)
  ggsave(plot=p2,filename=paste(my.project,"zTrophic-zRepro.jpg", sep="_"),path = file.path("PrelimResults","DataQA-QC",currentversion,"TrophicRepro_zscore"), create.dir = TRUE)
  
}

#### . Visualize each site with raw values ====
for(i in 1:length(unique(all_traits$project))){
  my.project <- unique(all_traits$project)[i]
  proj.dat <- all_traits %>% filter(project == my.project)
  proj.dat$logmass <- log10(proj.dat$mass_adult_g)
  proj.dat$logage <- log10(proj.dat$age_life.span_years)
  proj.dat$logrepro <- log10(proj.dat$reproduction.unified)
  if(length(which(!is.na(proj.dat$tr.age.zt)))>2 & length(which(!is.na(proj.dat$tr.mass.adult.zt)))>2 ){
    p1 <- plot_with_margins(proj.dat, x=logmass, y=logage, my.project)
    #print(p1)
  }
  if(length(which(!is.na(proj.dat$tr.age.zt)))<3 | length(which(!is.na(proj.dat$tr.mass.adult.zt)))<3 ){
    p1 <- ggplot() + labs(x="not enough age or mass data", title=my.project) + theme_bw() + labs(title=project)
    #print(p1)
  }
  if(length(which(!is.na(proj.dat$tr.trophic.level.zt)))>2 & length(which(!is.na(proj.dat$tr.reproduction.unified.zt)))>2 ){
    p2<- plot_with_margins(proj.dat, diet_trophic.level_num, logrepro, my.project)
    #print(p2)
  }
  if(length(which(!is.na(proj.dat$tr.trophic.level.zt)))<3 | length(which(!is.na(proj.dat$tr.reproduction.unified.zt)))<3 ){
    p2 <- ggplot() + labs(x="not enough trophic or reproduction data", title=my.project) + theme_bw()
    #print(p2)
  }
  ggsave(plot=p1,filename=paste(my.project,"logMass-logAge.jpg", sep="_"),path = file.path("PrelimResults","DataQA-QC",currentversion,"MassAge_raw"), create.dir = TRUE)
  ggsave(plot=p2,filename=paste(my.project,"Trophic-logRepro.jpg", sep="_"),path = file.path("PrelimResults","DataQA-QC",currentversion,"TrophicRepro_raw"), create.dir = TRUE)
  
}














## 6 - Save cleaned trait data ===================
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
