###project: CND Functional Synthesis
###author(s): MW
###goal(s): 
###date(s): February 2026
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries ---
# install.packages("librarian")
librarian::shelf(tidyverse, FD)

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

dt_og <- dt |>
      group_by(project, habitat) |>
      mutate(
            mean_dmperind = mean(dmperind_g_ind, na.rm = TRUE),
            sd_dmperind   = sd(dmperind_g_ind, na.rm = TRUE),
            upper_bound   = mean_dmperind + 5 * sd_dmperind,
            outlier_hi = dmperind_g_ind > upper_bound,
            elasmo = class %in% c("Chondrichthyes", "Elasmobranchii") |
                  order %in% c("Rajiformes", "Myliobatiformes", "Torpediniformes",
                               "Pristiophoriformes", "Squaliformes", "Carcharhiniformes",
                               "Lamniformes", "Orectolobiformes", "Hexanchiformes",
                               "Heterodontiformes")
      ) |>
      ungroup() |>
      filter(!(outlier_hi & elasmo)) |>
      select(-mean_dmperind, -sd_dmperind, -upper_bound, -outlier_hi, -elasmo)

# test <- anti_join(dt, dt_og)

### check to see NA fixes incorporated
na_count_per_column <- sapply(dt, function(x) sum(is.na(x)))
print(na_count_per_column)

dt1 <- dt_og |> 
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
      )) |> 
      filter(count < 1000) 

unique(dt1$phylum)
unique(dt1$project)
unique(dt1$scientific_name)
nacheck(dt1)
glimpse(dt1)

dt2 <- dt1 |>
      filter(project != 'MCR') |> 
      group_by(project, habitat, year, month,
               site, subsite_level1, subsite_level2, subsite_level3, scientific_name) |> 
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
            density = sum(density)) |>  
      ungroup() |> 
      arrange(project, year) |> 
      dplyr::select(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, total_n_area, total_p_area, total_bm_area, density, scientific_name)
glimpse(dt2)

dt2a <- dt1 |>
      filter(project == 'MCR') |> 
      group_by(project, habitat, year, month,
               site, subsite_level1, subsite_level2, subsite_level3, scientific_name) |> 
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
            density = sum(density),
            subsite_level3 = 'Not Available') |>  
      ungroup() |> 
      arrange(project, year) |> 
      dplyr::select(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, total_n_area, total_p_area, total_bm_area, density, scientific_name)
glimpse(dt2a)
dt2b <- rbind(dt2, dt2a)

dt2c <- dt2b |>
      arrange(project, year) |>
      mutate(
            system = case_when(
                  project == 'SBC' & habitat == 'ocean' ~ site,
                  project == 'FCE' ~ paste(site, subsite_level1, sep = ''),
                  project == 'MCR' ~ paste(subsite_level1, site, sep = ''),
                  project == 'COASTAL_CEN' ~ subsite_level2,
                  project == 'COASTAL_SOUTH' ~ subsite_level2)
      ) |> 
      select(project, habitat, site, year, month, system, scientific_name, 
             total_n_area, total_p_area, total_bm_area, density) |> 
      group_by(project, habitat, site, year, month, system) |>
      summarize(
            n = sum(total_n_area),
            p = sum(total_p_area),
            bm = sum(total_bm_area)
      ) |>
      group_by(project, year, habitat, system) |> 
      summarize(
            n = mean(n, na.rm = TRUE),
            p = mean(p, na.rm = TRUE),
            bm = mean(n, na.rm = TRUE),
            .groups = 'drop'
      )

fdiv <- read_rds('transformed_data/fd_ind_time.rds') |> 
      select(project, habitat, system, year, everything()) |> 
      rename(
            species_richness = sp_richn,
            functional_richness = fric,
            functional_dispersion = fdis,
            functional_specialization = fspe
      )
glimpse(fdiv)
glimpse(dt2c)

# additional stuff --------------------------------------------------------
# traits <- read_csv('../Collaborative/FnxSynthBase/consumer-trait-species-imputed-taxonmic-database.csv') |> 
#       select(phylum, order, family, genus, scientific_name,
#              age_life.span_years, diet_trophic.level_num, mass_adult_g,
#              reproduction_reproductive.rate_num.offspring.per.clutch.or.litter) |> 
#       distinct() |> 
#       mutate(age_life.span_years = case_when(
#             age_life.span_years < 0 ~ NA_real_,
#             TRUE ~ age_life.span_years
#       ))
# 
# fish <- dt2c |> select(scientific_name) |> distinct()
# traits1 <- fish |> left_join(traits)
# glimpse(traits1)
# 
# all <- dt2c |> left_join(traits1) 
# glimpse(all)

# all1 <- all |> 
#       group_by(project, habitat, site, year, month, system) |> 
#       mutate(across)
#       summarize(
#             total_bm = sum(total_bm_area, na.rm = TRUE),
#             age_wm  = ifelse(total_bm > 0,
#                              sum(age_life.span_years * total_bm_area, na.rm = TRUE) / total_bm,
#                              NA_real_),
#             diet_wm = ifelse(total_bm > 0,
#                              sum(diet_trophic.level_num * total_bm_area, na.rm = TRUE) / total_bm,
#                              NA_real_),
#             mass_wm = ifelse(total_bm > 0,
#                              sum(mass_adult_g * total_bm_area, na.rm = TRUE) / total_bm,
#                              NA_real_),
#             repro_wm = ifelse(total_bm > 0,
#                               sum(reproduction_reproductive.rate_num.offspring.per.clutch.or.litter * total_bm_area,
#                                   na.rm = TRUE) / total_bm,
#                               NA_real_),
#             .groups = "drop"
#       )
# write_csv(all, '../Collaborative/dataforcamille.csv')