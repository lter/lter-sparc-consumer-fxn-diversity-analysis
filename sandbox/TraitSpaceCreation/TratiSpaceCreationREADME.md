## TraitSpaceCreation

This folder holds the scripts necessary to take the most up-to-date master species list:
"Data/species_tidy-data/23_species_master-spp-list.csv"
cleans out non-focal taxa (humans, non-fish, non-birds, etc.), 
and merges it with the up-to-date master trait list:
"Data/traits_tidy-data/consumer-trait-species-imputed-taxonmic-database.csv"

Run this whole code, and it will create:
"transformed_data/species_list_corrected_fish.rds"
and more importantly:
"sp_tr_zscore.rds"

^ this latter contains raw traits, log10-transformed traits, and taxon z-scored traits
-tr.XXX.XX = trait of interest
-tr.XXX.unit = raw trait data in units of measure
-tr.XXX.log = log10-transformed trait
-tr.XXX.zp = z-score standardized traits per PROJECT
-tr.XXX.zt = z-score standardized traits per TAXON
  - for z-scored traits, the following traits are log10-transformed before z-scoring:
    - biomass.adult
    - age
    - reproduction.unified

* note on 'reproduction.unified'
this was a unified column combining fecundity for zooplankton and fish
and reproductive.rate_number.offspring.per.year for birds
** This column requires z-score standardizing for any real anlysis because it combines
data of different units.

