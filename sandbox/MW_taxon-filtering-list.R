###project: CND Functional Synthesis
###author(s): MW
###goal(s): 
###date(s): February 2026
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries ---
# install.packages("librarian")
librarian::shelf(tidyverse, janitor)

### set simple workflow functions ---
nacheck <- function(df) {
      na_count_per_column <- sapply(df, function(x) sum(is.na(x)))
      print(na_count_per_column)
}

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

dt <- spp_master |>
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
      )

proj_allowed <- list(
      "Arctic"      = c("zooplankton", "birds"),
      "Palmer"      = c("zooplankton", "birds"),
      "NorthLakes"  = c("zooplankton"),
      "FISHGLOB"    = c("fishes"),
      "FCE"         = c("fishes"),
      "SBC"         = c("fishes", "birds"),
      "CoastalCA"   = c("fishes", "other_invertebrates"),
      "MCR"         = c("fishes"),
      "NGA"         = c("zooplankton"),
      "VCR"         = c("fishes"),
      "KBS_AMP"     = c("amphibians"),
      "KBS_BIR"     = c("birds"),
      "KBS_INS"     = c("insects", 'other'),
      "KBS_MAM"     = c("mammals"),
      "HARVARD"     = c("birds"),
      "SEV"         = c("mammals"),
      "MOHONK"      = c("birds", "mammals"),
      "KONZA"       = c("insects", "mammals"),
      "PIE"         = c("fishes"),
      "RLS"         = c("fishes"),
      "SBC_BEACH"   = c('birds')
)

dt1 <- dt |>
      ### looks up which taxon groups are allowed for each project
      mutate(allowed = map(project, ~ proj_allowed[[.x]])) |>
      ### filters those not allowed
      filter(!map_lgl(allowed, is.null)) |>
      ### keep only rows that meet allowed groups for each project
      filter(map2_lgl(taxon_group, allowed, ~ .x %in% .y)) |>
      ### clean up dataframe
      select(-allowed)

removed <- anti_join(dt, dt1)
glimpse(removed)

acceptable <- removed |> 
      filter(project %in% c('KBS_AMP', 'KBS_BIR', 'KBS_MAM', 'HARVARD',
                            'KONZA', 'SEV')) |> 
      mutate(taxon_group = case_when(
            project == 'KBS_AMP' ~ "amphibians",
            project == 'KBS_BIR' ~ "birds",
            project == 'KBS_MAM' ~ "mammals",
            project == 'HARVARD' ~ "birds",
            project == 'KONZA' ~ "insects",
            project == 'SEV' ~ "mammals",
            TRUE ~ NA_character_
      ))
      
dt2 <- rbind(dt1, acceptable)

test <- dt2 |> count(project, taxon_group) |> arrange(project, desc(n))