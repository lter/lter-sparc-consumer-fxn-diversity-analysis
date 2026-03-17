###project: LTER Consumer Functional Diversity
###author(s): MW
###goal(s): 
###date(s): Spring 2026
###note(s): 

# Housekeeping ------------------------------------------------------------
### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, vegan, readxl, splitstackshape, codyn, lavaan,
                 MuMIn, corrplot, performance, ggeffects, ggpubr, parameters, ggstats,
                 brms, mixedup, rstatix, sf, ggspatial, waldo, multcompView, tidySEM, skimr,
                 GGally, lme4, lmerTest)

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

# filter to only FCE for now and test non-random species losses -----------
test <- dt |> filter(project == 'FCE') |> 
      dplyr::select(project, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, 
                    kingdom, phylum, class, order, family, genus, scientific_name,
                    diet_cat, 
                    dmperind_g_ind, nind_ug_hr, density_num) |> 
      mutate(
            fisheries_species = case_when(
                  scientific_name == 'Lutjanus griseus'        ~ 'fisheries',
                  scientific_name == 'Micropterus salmoides'   ~ 'fisheries',
                  scientific_name == 'Centropomus undecimalis' ~ 'fisheries',
                  scientific_name == 'Megalops atlanticus'     ~ 'fisheries',
                  scientific_name == 'Sciaenops ocellatus'     ~ 'fisheries',
                  TRUE ~ 'other'
            )
      )

