
## SEV DATA ONLY HAS SPECIES CODES

sev_abundance <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "SEV008_long.csv"))
sev_speciescodes <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "SEV_Small_Mammals_Species_List.csv"))

sev_speciescodes #has 37 obs 

unique(sev_speciescodes$Code) #but on 35 distinct codes... 
#first make a column with genus + species 

sev_speciescodes <- sev_speciescodes %>%
  dplyr::mutate(scientific_name = paste(Genus, Species)) %>%
  distinct()

#cool, now 35. not sure why dups to begin with?

#combine with abundance data

sev_abundance_species <- sev_abundance %>%
  left_join(sev_speciescodes, by = c("species" = "Code")) %>%
  select(-Genus) %>%
  select(-Species) %>%
  filter(location %in% c("5pgrass", "goatdraw", "rslarrea", "5plarrea"))

#8 locations: "core_PJ"  "savanna"  "rsgrass"  "two22"    "5plarrea" "rslarrea" "5pgrass"  "blugrama"
#goatdraw (1992-2008), include
#blue grama (2002-2004), exclude
#rsgrass (1989-1998), exclude (9 years)
#rslarrea (1989-2009), include 
#two2 (1989-1998), exclude
#savanna (1999-2002). exclude 
# Only 2 sites have been continuously been sampled since 1989 (5pgrass and 5plarrea).
#wtf is goatdraw???
  
unique(sev_abundance_species$location)

glimpse(sev_abundance_species)
