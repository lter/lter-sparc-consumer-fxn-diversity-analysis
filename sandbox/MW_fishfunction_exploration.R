###project: CND Functional Synthesis
###author(s): MW
###goal(s): 
###date(s): February 2026
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries ---
# install.packages("librarian")
librarian::shelf(tidyverse, FSA, readr, forcats, multcompView, readxl, dplyr, broom, ggpubr, lme4,
                 patchwork, splitstackshape, ggimage, purrr, zoo, pracma, vegan, e1071, codyn, lubridate, forecast)

### set simple workflow functions ---
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}

# Housekeeping ------------------------------------------------------------

### load necessary libraries
# install.packages("librarian")
librarian::shelf(tidyverse, dplyr, splitstackshape)

### read in clean excretion and strata data from google drive
dt <- read_csv("../Collaborative/FnxSynthBase/04_harmonized_consumer_excretion_sparc_cnd_site.csv") |> 
      janitor::clean_names() |> 
      ### removing invertebrate dominant projects, plus PIE given conversations with DB, NL, AS
      filter(!project %in% c("NGA", "CCE", "PIE")) |> 
      ### removing these sites, but just for right now [02/10/26]
      filter(!project %in% c('Palmer', 'RLS', 'Arctic', 'NorthLakes')) |> 
      ### split pisco up 
      mutate(project = case_when(
            project == 'CoastalCA' & site == 'CENTRAL' ~ 'COASTAL_CEN',
            project == 'CoastalCA' & site == 'SOUTH' ~ 'COASTAL_SOUTH',
            TRUE ~ project
      )) |> 
      mutate(subsite_level1 = replace_na(subsite_level1, "Not Available"),
             subsite_level2 = replace_na(subsite_level2, "Not Available"),
             subsite_level3 = replace_na(subsite_level3, "Not Available"),
             transectarea = coalesce(transectarea_m, transectarea_m2),
             area = case_when(
                   project == "COASTAL_CEN" ~ 60,
                   project == "COASTAL_SOUTH" ~ 60,
                   project == "SBC" ~ 40,
                   project == "VCR" ~ 25,
                   project == "MCR" & subsite_level3 == "1" ~ 50,
                   project == "MCR" & subsite_level3 == "5" ~ 250,
                   TRUE ~ transectarea)) |> 
      ### remove zeros - don't need :)
      filter(nind_ug_hr != 0) |> 
      dplyr::select(project, habitat, site, subsite_level1, subsite_level2, subsite_level3, area,
                    year, month, dmperind_g_ind, nind_ug_hr, pind_ug_hr, density_num_m, density_num_m2, 
                    density_num_m3, temp_c, diet_cat, kingdom, phylum, class, order, family, genus, 
                    scientific_name)
glimpse(dt)
nacheck(dt)

### check to see NA fixes incorporated
na_count_per_column <- sapply(dt, function(x) sum(is.na(x)))
print(na_count_per_column)

dt1 <- dt |> 
      ### classify each individual as either being a vertebrate or invertebrate
      mutate(vert_1 = if_else(
            phylum == "Chordata", 
            "vertebrate", 
            "invertebrate")
      ) |> 
      mutate(vert2 = if_else(
            is.na(vert_1) & project == "CoastalCA", 
            "vertebrate", 
            vert_1)
      ) |> 
      mutate(vert = if_else(
            is.na(vert2), 
            "invertebrate", 
            vert2)) |> 
      mutate(vertebrate_n = if_else(
            vert == "vertebrate" & dmperind_g_ind != 0, 
            1, 
            0),
            invertebrate_n = if_else(
                  vert == "invertebrate" & dmperind_g_ind != 0, 
                  1, 
                  0)) |> 
      filter(vert == "vertebrate") |> 
      select(-vert, -vert_1, -vert2, -vertebrate_n, -invertebrate_n) |> 
      mutate(
            density = coalesce(density_num_m, density_num_m2),
            count = round(density*area)) |> 
      mutate(count = case_when(
            count == 0 ~ 1,
            TRUE ~ count
      )) 

unique(dt1$phylum)
unique(dt1$project)
unique(dt1$scientific_name)
nacheck(dt1)
glimpse(dt1)

dt2 <- dt1 |>
      filter(project != 'MCR') |> 
      group_by(project, habitat, year, month, scientific_name,
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            ### calculate total nitrogen supply at each sampling unit and then sum to get column with all totals
            total_nitrogen_m = sum(nind_ug_hr * density, na_rm = TRUE),
            total_nitrogen_m2 = sum(nind_ug_hr * density, na_rm = TRUE),
            # total_nitrogen_m3 = sum(nind_ug_hr * ind_density_num, na_rm = TRUE),
            ### create column with total_nitrogen contribution for each program, regardless of units
            total_n_area = coalesce(total_nitrogen_m, total_nitrogen_m2),
            ### calculate total phosphorus supply at each sampling unit and then sum to get column with all totals
            total_phosphorus_m = sum(pind_ug_hr * density, na_rm = TRUE),
            total_phosphorus_m2 = sum(pind_ug_hr * density, na_rm = TRUE),
            # total_phosphorus_m3 = sum(pind_ug_hr * ind_density_num, na_rm = TRUE),
            ### create column with total_phosphorus contribution for each program, regardless of units
            total_p_area = coalesce(total_phosphorus_m, total_phosphorus_m2),
            ### calculate total biomass at each sampling unit and then sum to get column with all totals
            total_bm_m = sum(dmperind_g_ind * density, na_rm = TRUE),
            total_bm_m2 = sum(dmperind_g_ind * density, na_rm = TRUE),
            # total_bm_m3 = sum(dmperind_g_ind*ind_density_num, na_rm = TRUE),
            ### create column with total_biomass for each program, regardless of units
            total_bm_area = coalesce(total_bm_m, total_bm_m2)) |>  
      ungroup() |> 
      arrange(project, year) |> 
      dplyr::select(project, habitat, scientific_name, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, total_n_area, total_p_area, total_bm_area)
glimpse(dt2)

dt2a <- dt1 |>
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month, scientific_name,
               site, subsite_level1, subsite_level2, subsite_level3) |> 
      summarize(
            ### calculate total nitrogen supply at each sampling unit and then sum to get column with all totals
            total_nitrogen_m = sum(nind_ug_hr * density, na_rm = TRUE),
            total_nitrogen_m2 = sum(nind_ug_hr * density, na_rm = TRUE),
            # total_nitrogen_m3 = sum(nind_ug_hr * ind_density_num, na_rm = TRUE),
            ### create column with total_nitrogen contribution for each program, regardless of units
            total_n_area = coalesce(total_nitrogen_m, total_nitrogen_m2),
            ### calculate total phosphorus supply at each sampling unit and then sum to get column with all totals
            total_phosphorus_m = sum(pind_ug_hr * density, na_rm = TRUE),
            total_phosphorus_m2 = sum(pind_ug_hr * density, na_rm = TRUE),
            # total_phosphorus_m3 = sum(pind_ug_hr * ind_density_num, na_rm = TRUE),
            ### create column with total_phosphorus contribution for each program, regardless of units
            total_p_area = coalesce(total_phosphorus_m, total_phosphorus_m2),
            ### calculate total biomass at each sampling unit and then sum to get column with all totals
            total_bm_m = sum(dmperind_g_ind * density, na_rm = TRUE),
            total_bm_m2 = sum(dmperind_g_ind * density, na_rm = TRUE),
            # total_bm_m3 = sum(dmperind_g_ind*ind_density_num, na_rm = TRUE),
            ### create column with total_biomass for each program, regardless of units
            total_bm_area = coalesce(total_bm_m, total_bm_m2),
            subsite_level3 = 'Not Available') |> 
      ungroup() |> 
      arrange(project, year) |> 
      dplyr::select(project, habitat, scientific_name, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, total_n_area, total_p_area, total_bm_area)
glimpse(dt2a)

dt2b <- rbind(dt2, dt2a)

dt2c <- dt2b |> 
      arrange(project, year) |> 
      mutate(
            system = case_when(
                  project == 'SBC' & habitat == 'ocean' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, subsite_level2, sep = ''),
                  project == 'COASTAL_CEN' ~ subsite_level2,
                  project == 'COASTAL_SOUTH' ~ subsite_level2)
      ) |> 
      select(project, habitat, year, month, system, scientific_name, total_n_area, total_p_area, total_bm_area)

write_csv(dt2c, '../Collaborative/FnxSynthBase/cleanfishdataforfnx.csv')
glimpse(dt2c)