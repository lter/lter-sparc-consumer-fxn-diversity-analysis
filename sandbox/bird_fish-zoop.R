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
               site, subsite_level1, subsite_level2, scientific_name) |> 
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
                  project == 'COASTAL_SOUTH' ~ subsite_level2,
                  project == 'VCR' ~ paste(site, subsite_level1, sep = ''))) |> 
      select(project, habitat, site, year, month, system, 
             total_n_area, total_p_area, total_bm_area, density, scientific_name) |> 
      group_by(project, habitat, site, year, month, system, scientific_name) |>
      summarize(
            n = mean(total_n_area),
            p = mean(total_p_area),
            bm = mean(total_bm_area)
      ) |> 
      group_by(project, habitat, year, system, scientific_name) |>
      summarize(
            n = mean(n),
            p = mean(p),
            bm = mean(bm),
            .groups = 'drop'
      ) 
fish <- dt2c
# fdiv <- read_rds('transformed_data/fd_ind_time.rds') |> 
#       select(project, habitat, system, year, scientific_name, everything()) |> 
#       rename(
#             species_richness = sp_richn,
#             functional_richness = fric,
#             functional_dispersion = fdis,
#             functional_specialization = fspe,
#             functional_nearest_neighbor_distance = fnnd
#       )
# glimpse(fdiv)
# glimpse(dt2c)
# 
# test <- anti_join(dt2c, fdiv)
# 
# all <- dt2c |> left_join(fdiv) |> 
#       ### 18 observations didn't pull over, but can figure out later + now, a bunch of VCR sites...
#       filter(!is.na(functional_richness))
# glimpse(all)
# fish <-all
# all1 <- fish
rm(list = setdiff(ls(), c("fish", "nacheck")))

dt <- read_csv('../Collaborative/FnxSynthBase/04_harmonized_consumer_excretion_sparc_cnd_site02122026.csv') |> 
      janitor::clean_names() |> 
      ### removing invertebrate dominant projects, plus PIE given conversations with DB, NL, AS
      filter(project %in% c("NGA", "CCE", 'Palmer', 'Arctic', 'NorthLakes')) |> 
      # ### removing these sites, but just for right now [02/10/26]
      # filter(!project %in% c('Palmer', 'RLS', 'Arctic', 'NorthLakes')) |> 
      # ### split pisco up 
      # mutate(project = case_when(
      #       project == 'CoastalCA' & site == 'CENTRAL' ~ 'COASTAL_CEN',
      #       project == 'CoastalCA' & site == 'SOUTH' ~ 'COASTAL_SOUTH',
      #       TRUE ~ project
      # )) |> 
      mutate(subsite_level1 = replace_na(subsite_level1, "Not Available"),
             subsite_level2 = replace_na(subsite_level2, "Not Available"),
             subsite_level3 = replace_na(subsite_level3, "Not Available")) |> 
      ### remove zeros - don't need :)
      filter(nind_ug_hr != 0) |> 
      dplyr::select(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, dmperind_g_ind, nind_ug_hr, pind_ug_hr, density_num_m, density_num_m2, 
                    density_num_m3, temp_c, diet_cat, kingdom, phylum, class, order, family, genus, 
                    scientific_name)
glimpse(dt)
nacheck(dt)

spp_master <- readr::read_csv("../Collaborative/FnxSynthBase//23_species_master-spp-list.csv") |>
      janitor::clean_names()

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

dt_filter <- spp_master |>
      mutate(
            taxon_group = case_when(
                  class %in% bird_classes ~ "birds",
                  class %in% mammal_classes ~ "mammals",
                  class %in% amphib_classes ~ "amphibians",
                  class %in% fish_classes ~ "fishes",
                  class == "Insecta" ~ "insects",
                  class %in% zooplankton_classes ~ "zooplankton",
                  phylum == "Chordata" ~ "other_chordate",
                  TRUE ~ "other_invertebrate"
            )
      ) |> 
      select(scientific_name, taxon_group) |> 
      distinct()
glimpse(dt_filter)
glimpse(dt)

dt1 <- dt |> left_join(dt_filter) |> 
      filter(taxon_group == 'zooplankton')
glimpse(dt1)
unique(dt1$taxon_group)
nacheck(dt1)

dt2 <- dt1 |>
      group_by(project, habitat, year, month,
               site, subsite_level1, subsite_level2, subsite_level3, scientific_name) |> 
      mutate(
            density = coalesce(density_num_m, density_num_m2, density_num_m3),
            ) |> 
      summarize(
            total_n_area = sum(nind_ug_hr * density, na_rm = TRUE),
            total_p_area = sum(pind_ug_hr * density, na_rm = TRUE),
            total_bm_area = sum(dmperind_g_ind * density, na_rm = TRUE),
            density = sum(density)) |>  
      ungroup() |> 
      arrange(project, year) |> 
      dplyr::select(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                    year, month, total_n_area, total_p_area, total_bm_area, density, scientific_name)
glimpse(dt2)

systems <- dt2 |> select(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                         scientific_name) |> distinct()

zoop <- dt2 |> 
      arrange(project, year) |>
      mutate(
            system = case_when(
                  project == 'Arctic' ~ paste(site, subsite_level1, sep = ''),
                  project == 'CCE' ~ site,
                  project == 'NGA' ~ subsite_level1,
                  project == 'NorthLakes' ~ site,
                  project == 'Palmer' ~ site)
      ) |> 
      select(project, habitat, site, year, month, system, 
             total_n_area, total_p_area, total_bm_area, density, scientific_name) |> 
      group_by(project, habitat, site, year, month, system, scientific_name) |>
      summarize(
            n = mean(total_n_area),
            p = mean(total_p_area),
            bm = mean(total_bm_area)
      ) |> 
      group_by(project, habitat, year, system, scientific_name) |>
      summarize(
            n = mean(n),
            p = mean(p),
            bm = mean(bm),
            .groups = 'drop'
      ) 
glimpse(zoop)
glimpse(fish)
zoop_fish <- rbind(fish, zoop)
glimpse(zoop_fish)

traits <- read_rds('transformed_data/sp_tr_zscore.rds') |> ungroup() |> 
      select(
            scientific_name, 
            tr.age.zt,
            tr.trophic.level.zt,
            tr.fecundity.zt,
            tr.reproduction.unified.zt,
            tr.mass.adult.zt
            ) |> 
      distinct()
glimpse(traits)

zoop_fish_traits <- zoop_fish |> left_join(traits) |> distinct()
nacheck(zoop_fish_traits)

rm(list = setdiff(ls(), c("fish", "zoop", "nacheck", "fish_short", "zoop_fish", "zoop_fish_traits")))
glimpse(zoop_fish_traits)
tr_cols <- names(zoop_fish_traits) |> grep("^tr\\.", x = _, value = TRUE)

cwm <- zoop_fish_traits |>
      group_by(project, system, year) |>
      summarize(
            across(
                  all_of(tr_cols),
                  ~ weighted.mean(.x, w = bm, na.rm = TRUE),
                  .names = "cwm_{.col}"
            ),
            bm_tot = sum(bm, na.rm = TRUE),
            n_taxa = n_distinct(scientific_name[bm!=0]),
            .groups = "drop"
      ) |> 
      mutate(
            across(
                  starts_with("cwm_tr."),
                  ~ ifelse(is.nan(.x), 
                           NA_real_, 
                           .x)
            )
      ) |> 
      filter(!project %in% c('Arctic', 'Palmer', 'NorthLakes'))
glimpse(cwm)

trait_metrics <- names(cwm) |> grep("^cwm_tr\\.", x = _, value = TRUE)
trait_metrics

cwm |>
      pivot_longer(cols = all_of(trait_metrics), names_to = "metric", values_to = "value") |>
      ggplot(aes(x = value)) +
      geom_histogram(bins = 30, color = "black", fill = "grey80") +
      facet_wrap(~metric, scales = "free", ncol = 2) +
      theme_classic() +
      labs(x = NULL, y = "Count")

cwm_ts_long <- cwm |>
      # optional: convert NaN to NA across numeric cols if needed
      mutate(across(where(is.numeric), ~ replace(.x, is.nan(.x), NA_real_))) |>
      # z-score each metric across all rows (or you can do within project/system—see note below)
      mutate(across(all_of(trait_metrics), ~ scale(.x)[, 1], .names = "{.col}_z")) |>
      pivot_longer(
            cols = ends_with("_z"),
            names_to = "metric",
            values_to = "value"
      ) |>
      mutate(
            variable = str_remove(metric, "_z$"),
            year_scaled = (year - min(year, na.rm = TRUE)) / (max(year, na.rm = TRUE) - min(year, na.rm = TRUE))
      ) |>
      group_by(project, system, variable) |>
      filter(n_distinct(year) >= 4) |>
      ungroup() |>
      select(project, system, year, metric, value, variable, year_scaled)

glimpse(cwm_ts_long)

trend_results_cwm <- cwm_ts_long |>
      group_by(project, system, variable) |>
      nest() |>
      mutate(
            n_years = map_int(data, ~ n_distinct(.x$year)),
            model = map(data, ~ if (all(is.na(.x$value))) NULL else lm(value ~ year_scaled, data = .x)),
            model_type = "lm",
            slope = map_dbl(model, ~ if (is.null(.x)) NA_real_ else unname(coef(.x)["year_scaled"])),
            std_error = map_dbl(model, ~ if (is.null(.x)) NA_real_ else summary(.x)$coefficients["year_scaled","Std. Error"]),
            p_value   = map_dbl(model, ~ if (is.null(.x)) NA_real_ else summary(.x)$coefficients["year_scaled","Pr(>|t|)"]),
            conf_low  = slope - 1.96 * std_error,
            conf_high = slope + 1.96 * std_error,
            significant = case_when(
                  is.na(p_value) ~ NA_character_,
                  p_value < 0.05 & slope > 0 ~ "Increase",
                  p_value < 0.05 & slope < 0 ~ "Decline",
                  p_value < 0.1 ~ "Marginal",
                  TRUE ~ "Not significant"
            )
      ) |>
      select(project, system, variable, n_years, slope, std_error, p_value,
             conf_low, conf_high, model_type, significant) |> 
      filter(!is.na(slope)) |> 
      ungroup()

glimpse(trend_results_cwm)
nacheck(trend_results_cwm)

# visualize trends (CWM traits) --------------------------------------------

metric_order <- c(
      "cwm_tr.age.zt",
      "cwm_tr.trophic.level.zt",
      "cwm_tr.fecundity.zt",
      "cwm_tr.reproduction.unified.zt",
      "cwm_tr.mass.adult.zt"
      # add more here if you have them, e.g. "cwm_tr.length.adult.zt"
)

pretty_names <- c(
      "cwm_tr.age.zt"                  = "Age",
      "cwm_tr.trophic.level.zt"        = "Trophic level",
      "cwm_tr.fecundity.zt"            = "Fecundity",
      "cwm_tr.reproduction.unified.zt" = "Reproduction",
      "cwm_tr.mass.adult.zt"           = "Adult mass"
)

trend_plot_df <- trend_results_cwm |>
      mutate(
            variable = factor(variable, levels = metric_order),
            variable_pretty = pretty_names[as.character(variable)],
            measure = "Trait",
            zooplankton_site = project %in% c("CCE", "NGA")
      ) |>
      filter(!is.na(project), !is.na(variable_pretty)) |>
      mutate(significance = case_when(
            significant == "Not significant" ~ "None",
            significant == "Marginal" ~ "Marginal",
            significant %in% c("Increase", "Decline") ~ "Significant",
            TRUE ~ NA_character_
      ))

glimpse(trend_plot_df)
unique(trend_plot_df$significance)

sig_palette <- c(
      "None" = "white",
      "Marginal" = "black",
      "Significant" = "red"
)

trait_plot <- function(trait_name) {
      trend_plot_df |>
            filter(variable_pretty == trait_name) |>
            ggplot(aes(x = project, y = slope)) +
            geom_hline(yintercept = 0, linetype = "solid", color = "black") +
            geom_boxplot(width = 0.5, linewidth = 1, outlier.shape = NA, color = "black") +
            geom_jitter(
                  aes(fill = significance, shape = zooplankton_site),
                  color = "black", size = 4, width = 0.25, alpha = 1
            ) +
            scale_fill_manual(values = sig_palette) +
            scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +
            labs(y = NULL, x = NULL, fill = "Significance", shape = "Zooplankton site") +
            theme_classic() +
            theme(
                  axis.text.x = element_text(face = "bold", color = "black", size = 12),
                  axis.text.y = element_text(face = "bold", color = "black", size = 12),
                  strip.text = element_text(face = "bold", color = "black", size = 14),
                  legend.position = "right"
            )
}

a <- trait_plot("Age")
b <- trait_plot("Trophic level")
c <- trait_plot("Fecundity")
d <- trait_plot("Reproduction")
e <- trait_plot("Adult mass")

a; b; c; d; e
