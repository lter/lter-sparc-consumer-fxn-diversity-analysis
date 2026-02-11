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

dt1 <- dt |>
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

# ---- project -> allowed taxon groups ----
proj_allowed <- list(
      "Arctic"      = c("zooplankton", "birds"),
      "Palmer"      = c("zooplankton", "birds"),
      "NorthLakes"  = c("zooplankton"),
      "FISHGLOB"    = c("fishes"),
      "FCE"         = c("fishes"),
      "SBC"         = c("fishes", "birds"),
      "CoastalCA"   = c("fishes"),
      "MCR"         = c("fishes"),
      "NGA"         = c("zooplankton"),
      "VCR"         = c("fishes"),
      "KBS_AMP"     = c("amphibians"),
      "KBS_BIR"     = c("birds"),
      "KBS_INS"     = c("insects"),
      "KBS_MAM"     = c("mammals"),
      "HARVARD"     = c("birds"),
      "SEV"         = c("mammals"),
      "MOHONK"      = c("birds", "mammals"),
      "KONZA"       = c("insects", "mammals")
)

# ---- filter ----
dt1 <- dt1 |>
      mutate(allowed = purrr::map(project, ~ proj_allowed[[.x]])) |>
      filter(!purrr::map_lgl(allowed, is.null)) |>
      filter(purrr::map2_lgl(taxon_group, allowed, ~ .x %in% .y)) |>
      select(-allowed)

# quick check
test <- dt1 |> count(project, taxon_group) |> arrange(project, desc(n))
