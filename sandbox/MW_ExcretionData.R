###project: LTER Consumer Functional Diversity
###author(s): MW
###goal(s): Wrangling and summarizing raw CND data such that it is ready for CFD analysis
###date(s): Spring 2026
###note(s): 

# Housekeeping ------------------------------------------------------------
### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, readxl, splitstackshape, codyn, lavaan,
                 MuMIn, corrplot, performance, ggeffects, ggpubr, parameters, ggstats,
                 brms, mixedup, rstatix, sf, ggspatial, waldo, multcompView, tidySEM)

### set custom functions
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}

# Load and prepare data ----------------------------------------------------
dt <- read.csv(file.path("Data/community_tidy-data/", "04_harmonized_consumer_excretion_sparc_cnd_site.csv"),
               stringsAsFactors = F,
               ### all NAs were inititally transformed to '.'
               na.strings =".") |> 
      
      # tidy up column names
      janitor::clean_names() |> 
      
      # replace NAs/missing data with character data
      mutate(
            subsite_level1 = ifelse(is.na(subsite_level1) | subsite_level1 == "", "not available", subsite_level1),
            subsite_level2 = ifelse(is.na(subsite_level2) | subsite_level2 == "", "not available", subsite_level2),
            subsite_level3 = ifelse(is.na(subsite_level3) | subsite_level3 == "", "not available", subsite_level3)
      ) |> 
      
      # filter out projects we are using at this point
      filter(project %in% c('MCR', 'CoastalCA', 'SBC', 'VCR', 'FCE',
                            'RLS', 'NGA', 'CCE')) |> 
      
      # fix ocean misspelling
      mutate(
            habitat = case_when(
                  habitat == 'ocean ' ~ 'ocean',
                  TRUE ~ habitat
            )
      ) |> 
      
      # filter out beach habitat [talitrids at SBC]
      filter(habitat %in% c('estuary', 'ocean')) |> 
      
      # filter out sites not sampled consistently at FCE
      filter(
            !(project == "FCE" & site == "TB" & subsite_level1 == "5"),
            !(project == "FCE" & site == "RB" & subsite_level1 %in% c("17", "19"))
      ) |> 
      
      # remove weird phylum stuff for right now - fix later :)
      filter(phylum!="") |> 
      
      # set upper end for California Moray eel based on reported maximum weight
      mutate(
            dmperind_g_ind = case_when(
                  dmperind_g_ind > 9071 & scientific_name == "Gymnothorax mordax" ~ 9071,
                  TRUE ~ dmperind_g_ind
            )) |> 
      
      # remove 'biomass buster' sharks and rays from CoastalCA, SBC, and MCR [lost ~ 0.01% of data]
      group_by(project, habitat) |> 
      mutate(
            mean_dmperind = mean(dmperind_g_ind, na.rm = TRUE),
            sd_dmperind   = sd(dmperind_g_ind, na.rm = TRUE),  
            lower_bound   = mean_dmperind - 5 * sd_dmperind,  
            upper_bound   = mean_dmperind + 5 * sd_dmperind,
            outlier       = dmperind_g_ind < lower_bound | dmperind_g_ind > upper_bound,
            # sharkray      = grepl("\\bshark\\b|\\bray\\b", common_name, ignore.case = TRUE), # common name not carried forward in CFD
            elasmo        = class %in% c("Chondrichthyes", "Elasmobranchii")
      ) |> 
      ungroup() |> 
      filter(!(outlier & (elasmo))) |> 
      select(-mean_dmperind, -sd_dmperind, -lower_bound, -upper_bound, -outlier, -elasmo) |> 
      
      # coalesce density columns and remove unnecessary '..m3' column
      mutate(density_num = coalesce(density_num_m, density_num_m2, density_num_m3)) |> 
      select(-density_num_m, -density_num_m2, -density_num_m3)
glimpse(dt)
nacheck(dt)

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
# Summarize data ---------------------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

##################################################################################################
##################################################################################################
##################################################################################################
### Calculating CND, Biomass and Species Diversity metrics ---------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################

### MCR subsite_level3 has two levels (1 and 5) that should be summed across, 
### as they are conducted within same space and time so split them out and go
### sum across subsite_level2 
dt2_mcr <- dt1 |> 
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               scientific_name) |> 
      summarize(
            n       = sum(nind_ug_hr*density_num, na.rm = TRUE),
            bm      = sum(dmperind_g_ind*density_num, na.rm = TRUE),
            dens    = sum(density_num, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, 
               scientific_name) |> 
      summarize(
            species_n    = sum(n),
            species_bm   = sum(bm),
            species_dens = sum(dens),
            .groups = 'drop'
      ) |> 
      # subsite_level3 set to '15' to denote combined transects 1 and 5
      mutate(subsite_level3 = '15') |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             scientific_name, species_n, species_bm, species_dens)

glimpse(dt2_mcr)
head(dt2_mcr)
nacheck(dt2_mcr)

### all other sites should be summed down subsite_level3 as these are
### considered the 'transect'
dt2_other <- dt1 |> 
      filter(project != 'MCR') |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3, 
               scientific_name) |> 
      ### sum across unique taxa at the transect level
      ### coalesces multiple species observation to single n and bm value for each transect
      summarize(
            species_n    = sum(nind_ug_hr*density_num, na.rm = TRUE),
            species_bm   = sum(dmperind_g_ind*density_num, na.rm = TRUE),
            species_dens = sum(density_num, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             scientific_name, species_n, species_bm, species_dens)

glimpse(dt2_other)
head(dt2_other)
nacheck(dt2_other)

### join mcr and other data back together now that they are at same 'transect' scale
dt2 <- rbind(dt2_other, dt2_mcr)
glimpse(dt2)
head(dt2)
nacheck(dt2)
rm(dt2_mcr, dt2_other)

### sum across species at the transect scale to get community-level excretion and biomass
### per unit area
dt3_cnd <- dt2 |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            comm_n    = sum(species_n, na.rm = TRUE),
            comm_bm   = sum(species_bm, na.rm = TRUE),
            comm_dens = sum(species_dens, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL' ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH' ~ 'PCCS',
            TRUE ~ project
      )) |> 
      ### take everything to the true 'site' level in which transects should be
      ### averaged across
      mutate(
            site_site = case_when(
                  project == 'SBC' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
            ) 
      ) |> 
      group_by(project, site_site, year) |> 
      summarize(
            comm_n    = mean(comm_n, na.rm = TRUE),
            comm_bm   = mean(comm_bm, na.rm = TRUE),
            comm_dens = mean(comm_dens, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      rename(site = site_site)

glimpse(dt3_cnd)
head(dt3_cnd)
nacheck(dt3_cnd)
