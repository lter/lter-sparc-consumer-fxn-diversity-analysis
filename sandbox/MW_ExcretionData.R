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
      
      # remove extraordinary large schools of fish [lost < 0.0001 % of data]
      # FCE doesn't catch thousands of fish per transect, so fine for this given variable area measurements (i.e., electrofishing program)
      # could use new rule for CCE/NGA, and also don't know RLS area sampled
      mutate(
            area = case_when(
                  project == 'CoastalCA'                    ~ 60,
                  project == 'SBC'                          ~ 80,
                  project == 'VCR'                          ~ 25,
                  project == 'MCR' & subsite_level3 == '1'  ~ 50,
                  project == 'MCR' & subsite_level3 == '5'  ~ 250,
                  project == 'RLS' & subsite_level2 == '2'  ~ 50,
                  project == 'RLS' & subsite_level2 == '1'  ~ 250,
                  TRUE ~ NA_real_
            ),
            count = area*density_num_m2
      ) |> 
      filter(project %in% c("FCE", 'CCE', 'NGA') | count <10000) |> 
      select(-area, -count) |> 
      
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
      select(-density_num_m, -density_num_m2, -density_num_m3) |> 
            
      # remove high density of large-bodied fish observations that skew entire time series [lost < 0.0001% of data]
      filter(!(density_num > 1 & nind_ug_hr > 20000 & phylum == 'Chordata')) |> 
      
      # add flag for removing fishes from zooplankton sites and vice versa
      mutate(flag = case_when(
            project %in% c('NGA', 'CCE') & phylum == 'Chordata' ~ 'remove',
            project %in% c('FCE', 'SBC', 'MCR', 'VCR', 'RLS', 'CoastalCA') & phylum != 'Chordata' ~ 'remove',
            TRUE ~ 'keep'
      )) |> 
      filter(flag == 'keep') |> 
      dplyr::select(-flag)
            
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
dt_mcr <- dt |> 
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

glimpse(dt_mcr)
head(dt_mcr)
nacheck(dt_mcr)

### RLS subsite_level2 has two levels (1 and 2) that should be summed across, 
### as they are conducted within same space and time so split them out and go
### sum across subsite_level1
dt_rls <- dt |> 
      filter(project == 'RLS') |> 
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
               site, subsite_level1, 
               scientific_name) |> 
      summarize(
            species_n    = sum(n),
            species_bm   = sum(bm),
            species_dens = sum(dens),
            .groups = 'drop'
      ) |> 
      # subsite_level2 set to '12' to denote combined transects 1 and 2
      mutate(subsite_level2 = '12',
             subsite_level3 = 'not available') |> 
      select(project, habitat, year, month, site, subsite_level1, subsite_level2, subsite_level3,
             scientific_name, species_n, species_bm, species_dens)

glimpse(dt_rls)
head(dt_rls)
nacheck(dt_rls)

### all other sites should be summed down subsite_level3 as these are
### considered the 'transect'
dt_other <- dt |> 
      filter(!project %in% c('MCR', 'RLS')) |> 
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

glimpse(dt_other)
head(dt_other)
nacheck(dt_other)

### join mcr and other data back together now that they are at same 'transect' scale
dt1 <- rbind(dt_other, dt_mcr, dt_rls)
glimpse(dt1)
head(dt1)
nacheck(dt1)
rm(dt_mcr, dt_other, dt_rls)

### sum across species at the transect scale to get community-level excretion and biomass
### per unit area
dt_cnd <- dt1 |> 
      group_by(project, habitat, year, month, 
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            comm_n    = sum(species_n, na.rm = TRUE),
            comm_bm   = sum(species_bm, na.rm = TRUE),
            comm_dens = sum(species_dens, na.rm = TRUE),
            .groups   = 'drop'
      ) |> 
      mutate(project = case_when(
            project  == 'CoastalCA' & site == 'CENTRAL'   ~ 'PCCC',
            project  == 'CoastalCA' & site == 'SOUTH'     ~ 'PCCS',
            project  == 'RLS' & site == 'Maria Island'    ~ 'RLSM',
            project  == 'RLS' & site == 'Rottnest Island' ~ 'RLSR',
            project  == 'RLS' & site == 'Ningaloo Reef'   ~ 'RLSN',
            TRUE ~ project
      )) |> 
      ### take everything to the true 'site' level in which transects should be
      ### averaged across
      mutate(
            site_site = case_when(
                  project == 'SBC'  ~ site,
                  project == 'FCE'  ~ paste(site, subsite_level1, sep = ''),
                  project == 'VCR'  ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR'  ~ paste(subsite_level1, site, sep = ''),
                  project == 'CCE'  ~ paste(site, subsite_level1, sep = ''),
                  project == 'NGA'  ~ subsite_level1,
                  project == 'PCCC' ~ subsite_level2,
                  project == 'PCCS' ~ subsite_level2,
                  project == 'RLSM' ~ subsite_level1,
                  project == 'RLSR' ~ subsite_level1,
                  project == 'RLSN' ~ subsite_level1,
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

glimpse(dt_cnd)
head(dt_cnd)
nacheck(dt_cnd)

dt_cnd |> 
      mutate(project = as.factor(project)) |> 
      group_by(project, site, year) |> 
      summarize(
            mean = mean(comm_n + 1, na.rm = TRUE),
            median = median(comm_n + 1, na.rm = TRUE),
            .groups = 'drop'
      ) |> 
      ggplot(aes(x = year, y = mean, fill = site, color = site, group = site)) + 
      geom_jitter(aes(fill = site), shape = 21, width = 0.3,
                  color ='black', size = 2, stroke = 1, alpha = 0.6) +
      geom_smooth(method = 'loess', se = FALSE) +
      facet_wrap(~project, scale = 'free') + 
      labs(x = 'Project', y = 'Nitrogen Supply (ug/hr/area)',
           title = 'Mean Annual Community Nitrogen Supply by Site') + 
      theme(
            strip.text = element_text(size = 16, face = "bold", colour = "black"),
            strip.background = element_blank(),  
            axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none",
            plot.title = element_text(size = 16, face = "bold", colour = "black",
                                      hjust = 0.5)
      )