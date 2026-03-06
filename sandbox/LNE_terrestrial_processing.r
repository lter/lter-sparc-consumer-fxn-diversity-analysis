

library(tidyverse)

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

## SBC DATA HAS DEAD ANIMALS IN IT

sbc_abundance <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "Shorebird_count_20231020.csv"))

colnames(sbc_abundance)
unique(sbc_abundance$TAXON_GROUP)
#42 taxon types

#kyle recommended shorebird and gull and probably need "."

sbc_shorbird <- sbc_abundance %>%
  filter(TAXON_GROUP %in% c("shorebird", "gull"))

unique(sbc_shorbird$TAXON_SPECIES)

#looks like 44 species 

taxo_edbird <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "eBird_taxonomy_v2025.csv"))

colnames(taxo_edbird)

taxo_edbird_limited <- taxo_edbird %>%
  dplyr::select(SCI_NAME, ORDER, FAMILY)

taxo_edbird_limited$FAMILY <- gsub("\\s*\\([^)]*\\)", "", taxo_edbird_limited$FAMILY)

unique(taxo_edbird_limited$FAMILY)

taxo_edbird_limited_commonname <- taxo_edbird %>%
  dplyr::select(SCI_NAME, ORDER, FAMILY, PRIMARY_COM_NAME)

taxo_edbird_limited_commonname$FAMILY <- gsub("\\s*\\([^)]*\\)", "", taxo_edbird_limited_commonname$FAMILY)

unique(taxo_edbird_limited$FAMILY)

#### Andrews Experimental ####

andrews_ex_abundance <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "SA02402_v4.csv"))
andrews_birdcodes <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "AndrewsForest_BirdDataCodes.csv"))
unique(andrews_ex_abundance$PLOT)
unique(andrews_ex_abundance$ENTITY)
unique(andrews_ex_abundance$REPLICATE)
unique(andrews_ex_abundance$SURVEY_DATE)
unique(andrews_ex_abundance$PERIOD) #3 periods
colnames(andrews_ex_abundance)
glimpse(andrews_ex_abundance)
#184 plots 
#why does meta data say 183? 
#DBCODE - unless 
#entity - 2 for everything 
#up to 6 replicates?

#connect bird codes with scientific names

andrews_ex_abundance_codes <- andrews_ex_abundance %>%
  left_join(andrews_birdcodes, by = c("SPECIES" = "species_code"))

#look for Nas
df_na_species <- andrews_ex_abundance_codes %>%
  dplyr::filter(is.na(species_name)) %>%
  distinct(SPECIES)

#okay, all good 
#need to just pull out distinct bird names

andrews_ex_birdspecieslist <- andrews_ex_abundance_codes %>%
  filter(!SPECIES %in% c("PIKA","DOSQ", "UNMA", "CHIP", "UNKN", "UNFI", "UNFL", 
                         "UNGR", "UNJA", "UNTH", "UNVI", "UNWA", "UNWO", "UNWR", "UNCH",
                         "UNCR", "UNKI", "UNSW", "UNOW", "UNHM", "UNDU", "UNRA", "UNCO", "UNHU")) %>% #get rid of mammals and unknown things
  distinct(species_name)
  
andrews_ex_birdspecieslist_formerge <- andrews_ex_birdspecieslist %>%
  mutate(project = "AndrewsForest",
         habitat = "terrestrial",
         raw_filename = "SA02402_v4.csv",
         class = "Aves")

andrews_ex_birdspecieslist_formerge_sciname_edits <- andrews_ex_birdspecieslist_formerge %>%
  mutate(
    species_name =   case_when(
      species_name == "Hylatomus pileatus" ~ "Dryocopus pileatus", #change to historical for merging, will need to change back. 
      species_name == "Carduelis pinus" ~ "Spinus pinus",#taxonomic change
      species_name == "Picoides villosus" ~ "Leuconotopicus villosus", #taxonomic change
      species_name == "Regulus calendula" ~ "Corthylio calendula", #taxonomic change
      species_name == "Accipiter cooperii" ~ "Astur cooperii",#taxonomic change
      species_name == "Vermivora celata" ~ "Leiothlypis celata", #taxonomic change
      species_name == "Coccothraustes vespertinus" ~ "Hesperiphona vespertina", #change to historical for merging, will need to change back. 
      species_name == "Picoides pubescens" ~ "Dryobates pubescens", #taxonomic change
      species_name == "Accipiter gentilis" ~ "Astur gentilis", #this change seems less legit. need to change back
      TRUE ~ species_name
    )
  )

andrews_ex_birdspecieslist_formerge_taxo <- andrews_ex_birdspecieslist_formerge_sciname_edits %>%
  left_join(taxo_edbird_limited, by = c("species_name" = "SCI_NAME")) %>%
  rename("order" = "ORDER",
         "family" = "FAMILY",
         "scientific_name" = "species_name") %>%
  mutate(genus = word(scientific_name, 1))

#reverse those two taxonomic name changes

andrews_ex_birdspecieslist_final <- andrews_ex_birdspecieslist_formerge_taxo %>%
  mutate(
    scientific_name =   case_when(
      scientific_name =="Dryocopus pileatus" ~ "Hylatomus pileatus", #changed from historical to OG. 
      scientific_name == "Hesperiphona vespertina" ~ "Coccothraustes vespertinus",
      TRUE ~ scientific_name
      )) %>% #changed from historical to OG. 
  select(project, habitat, raw_filename, class, order, family, genus, scientific_name)

df_na_taxo <- andrews_ex_birdspecieslist_formerge_taxo %>%
  dplyr::filter(is.na(order)) #cool, 0 

andrews_ex_birdspecieslist_final %>%
  dplyr::filter(is.na(order))

#write this out for Shalanda

#andrewsforest <- "andrewsforest_specieslist.csv"

# Export
write.csv(andrews_ex_birdspecieslist_formerge_taxo, row.names = F, na = '',
          file = file.path("sandbox", "birds_specieslist", "andrewsforest_specieslist.csv"))

write_csv()

#need to remove non new-reocrds:

#and also remove non-birds 
andrews_ex_abundance_codes_newrecords <- andrews_ex_abundance_codes %>%
  filter(NEW_RECORD == "1") %>%
  filter(DISTANCE %in% c(1, 2)) %>%
  filter(!SPECIES %in% c("PIKA","DOSQ", "UNMA", "CHIP"))

library(lubridate)

andrews_abundance_bymonth <- andrews_ex_abundance_codes_newrecords %>%
  mutate(month = format(as.Date(SURVEY_DATE), "%m")) %>%
  mutate(OBS = 1)

andrews_abundance_juneonly <- andrews_abundance_bymonth %>%
  filter(month == "06")


andrews_abundance_juneonly_total<- andrews_abundance_juneonly %>%
  group_by(YEAR, PLOT, SPECIES) %>%
  summarise(count_total = sum(OBS))

andrews_abundance_juneonly_summarized <- andrews_abundance_juneonly_total %>%
  group_by(YEAR, SPECIES) %>%
  summarise(count_avg = mean(count_total),
            se_avg = sd(count_total)/sqrt(n()),
            n = n()) 


andrews_abundance_juneonly

andrews_abundance_juneonly_summarized %>%
 filter(SPECIES %in% c("WETA", "HEWA", "PAWR", "CBCH", "CHSP", "RECR")) %>%
  ggplot(aes(x = YEAR, y =count_avg)) +
  geom_point() +
  geom_errorbar(aes(ymin = count_avg - se_avg, ymax = count_avg + se_avg))+
  facet_wrap(~SPECIES)

unique(andrews_abundance_bymonth$SURVEY_DATE)

#only pull out June?
# need relative abundance, not density

unique(andrews_ex_abundance$SPECIES)

andrews_ex_abundance_newrecordonly <- 


test <- andrews_ex_abundance %>%
  group_by(PLOT, YEAR, SURVEY_DATE, SPECIES) %>%
  summarise(abundance = n(), .groups = "drop")

test2 <- andrews_ex_abundance %>%
  group_by(PLOT, YEAR, REPLICATE, SURVEY_DATE) %>%
  summarise(abundance = n(), .groups = "drop")

test3 <- andrews_ex_abundance %>%
  group_by(PLOT, YEAR, REPLICATE) %>%
  summarise(abundance = n(), .groups = "drop")

test4 <- andrews_ex_abundance %>%
group_by(PLOT, YEAR) %>%
  summarise(n_reps = n_distinct(REPLICATE, na.rm = TRUE), .groups = "drop") %>%
  arrange(n_reps)

test5 <- andrews_ex_abundance %>%
  group_by(PLOT) %>%
  summarise(n_reps = n_distinct(YEAR, na.rm = TRUE), .groups = "drop") %>%
  arrange(n_reps)



#10 minutes VS sections 

andrews_ex_abundance$OBS <- as.numeric(1)


andrews_ex_abundance %>%
  filter(YEAR == 2010) %>%
  filter(PLOT == "PC002")

PA002_17 <- andrews_ex_abundance %>%
  filter(YEAR == 2017) %>%
  filter(PLOT == "PA002")

unique(andrews_ex_abundance$RECORD)

hja2 <- andrews_ex_abundance %>%
  group_by(YEAR,PLOT,SURVEY_DATE,REPLICATE,PERIOD,SPECIES) %>%
  summarise(
    COUNT = sum(OBS)
  ) %>% ungroup()

hja3 <- hja2 %>%
  group_by(PLOT, YEAR, PERIOD,SURVEY_DATE) %>%
  complete(REPLICATE, SPECIES, fill = list(COUNT = 0)) %>%
  ungroup()


hja4 <- hja3 %>%
  group_by(YEAR,PLOT,SURVEY_DATE,REPLICATE,SPECIES) %>%
  summarise(
    MEAN = mean(COUNT)
  ) %>% ungroup()


date_count <- hja4 %>%
  group_by(YEAR, PLOT, REPLICATE) %>%
  summarise(
    survey_count = n_distinct(SURVEY_DATE),
  ) %>% ungroup()


#### Arctic LTER ####

#all csvs are separate

#### California Current Ecosystem LTER###
#2 datasets, not sure how they relate?

#### Hubbard Brook LTER ####
#nightmare

hubbardbrook_abundance <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "HubbardBrook_Bird_MainPlot_1969-2018.csv"))

hubbardbrook_abundance_birdlist <- hubbardbrook_abundance %>%
  select(BirdSpecies)

hubbardbrook_abundance_commonname_edits <- hubbardbrook_abundance_birdlist %>%
  mutate(
    BirdSpecies  = case_when(
      BirdSpecies  == "Ruby-throated  Hummingbird" ~ "Ruby-throated Hummingbird", #hyphen or capitalization update for matching
      BirdSpecies  == "Eastern Wood Pewee" ~ "Eastern Wood-Pewee", #hyphen or capitalization update for matching 
      BirdSpecies  == "Black and White Warbler" ~ "Black-and-white Warbler", #hyphen or capitalization update for matching
      TRUE ~ BirdSpecies
    )
  )

hubbardbrook_abundance_birdlist_taxo <- hubbardbrook_abundance_commonname_edits %>%
  left_join(taxo_edbird_limited_commonname, by = c("BirdSpecies" = "PRIMARY_COM_NAME")) %>%
  rename("order" = "ORDER",
         "family" = "FAMILY",
         "scientific_name" = "SCI_NAME") %>%
  mutate(genus = word(scientific_name, 1),
          project = "HubbardBrook",
         habitat = "terrestrial",
         raw_filename = "HubbardBrook_Bird_MainPlot_1969-2018.csv",
         class = "Aves")


hubbardbrook_abundance_birdlist_taxo%>%
  dplyr::filter(is.na(SCI_NAME)) # 0 - stick 

hubbardbrook_abundance_birdlist_taxo_final <- hubbardbrook_abundance_birdlist_taxo %>%
  select(project, habitat, raw_filename, class, order, family, genus, scientific_name)
  
  



#### Luquillo LTER ####

luquillo_abundance_biotime <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo_birds.csv"))
unique(luquillo_abundance_biotime$YEAR)

#Measurements of bird abundance are taken in the 9 ha grid (5/89 to 6/90) and Luquillolterdb23-Bird point counts Forest Dynamics Plot (10/90 to present) at El Verde and in Watershed 1 and Whendee Silver's cut plots at Bisley.DATA SET METHODSCircular plot counts are used to measure relative numbers of birds over time and betweensites. 

luquillo_abundance_biotime_clean <- luquillo_abundance_biotime %>%
  mutate(plot = str_to_lower(str_extract(SAMPLE_DESC, "[^_]+$"))) %>%
  summarise(n_plots = n_distinct(plot), .by = YEAR)

glimpse(luquillo_abundance_biotime_clean)
unique(luquillo_abundance_biotime_clean$SAMPLE_DESC)

luquillo_abundance_biotime_clean <- luquillo_abundance_biotime %>%
  mutate(plot = str_to_lower(str_extract(SAMPLE_DESC, "[^_]+$"))) %>%
  mutate(date_str = str_to_lower(str_extract(SAMPLE_DESC, "^\\d{4}_\\d{1,2}_\\d{1,2}")))

unique(luquillo_abundance_biotime_clean$ABUNDANCE)

unique(luquillo_abundance_biotime_clean$plot)
#shold have plot 12, 3, 6, 13, 15, 9
#should be 5 plots per year? 
#seems combined

luquillo_abundance_EL_GT25 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo", "All_ElVerde_GT_25m.csv"))
luquillo_abundance_EL_LT25 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo", "ALL_ElVerde_LT_25m.csv"))
luquillo_abundance_BS_GT25 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo", "BisleySilver_Cut_GT_25m.csv"))
luquillo_abundance_BS_LT25 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo", "BisleySilver_Cut_LT_25m.csv"))
luquillo_abundance_BW_GT25 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo", "BisleyWatershed1_GT_25m.csv"))
luquillo_abundance_BW_LT25 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "luquillo", "BisleyWatershed1_LT_25m.csv"))
colnames(luquillo_abundance_EL_GT25) #2 plots
glimpse(luquillo_abundance_EL_GT25)
unique(luquillo_abundance_EL_GT25$PLOT1) #"A"  "C"  "E"  "G"  "I"  "K"  "3"  "6"  "9"  "12" "15"
unique(luquillo_abundance_EL_GT25$PLOT2) #  1  3  5  7  9 11 10  8  6  4  2 12 15 18 21 24
colnames(luquillo_abundance_EL_LT25) #2 plots

luquillo_birdcode <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "Luquillo_BirdDataCodes.csv"))

luquillo_abundance_EL_GT25 %>%
  distinct(PLOT1, PLOT2, YEAR) %>%   # one record per pair per year
  count(PLOT1, PLOT2, name = "n_years")

#okay, seems like Plot 1 is what we are interested in as a "site"
#wtf is up with the letters? 

luquillo_abundance_EL_GT25 %>%
  filter(PLOT1 == "K") %>%
  distinct(YEAR) 
# A - 1989 and 1990
# C - 1989 and 1990
# E - 1989 and 1990
# G - 1989 and 1990
# I - 1989 and 1990
# K - 1989 and 1990

#okay, going to ignore first 2 years, they were sampled differently. 

luquillo_abundance_EL_GT25 %>%
  filter(PLOT1 == "C") %>%
  distinct(YEAR)


#don't care since only two year.... 
unique(luquillo_abundance_EL_GT25$YEAR) #1989-2021
unique(luquillo_abundance_EL_LT25$YEAR) #1989-2021
colnames(luquillo_abundance_BS_GT25) #2 plots
colnames(luquillo_abundance_BS_LT25) #2plots
unique(luquillo_abundance_BS_GT25$YEAR) #1989-1990 only... 
unique(luquillo_abundance_BS_LT25$YEAR) #1989-1990 only...
colnames(luquillo_abundance_BW_GT25) #2 plots
colnames(luquillo_abundance_BW_LT25) #2 plots
unique(luquillo_abundance_BW_GT25$YEAR) #1989-1991 only... 
unique(luquillo_abundance_BW_LT25$YEAR) #1989-1991 only...

#this tracks with the metadata. 
#so only luquillo_abundance_EL_GT25 & luquillo_abundance_EL_LT25

#first go wide to long... 

#going to remove the first 2 years because sampling seems inconsistent...
#also removing 1991 because it is werid sampling dates as well 

luquillo_abundance_EL_GT25_sub <- luquillo_abundance_EL_GT25 %>%
  filter(!YEAR %in% c(1989, 1990, 1991), 
  dplyr::between(JULIAN, 121, 182)) %>%
  select(-c(WEATHER, WIND, RAIN, TIME, JULIAN, OBSERVER))

unique(luquillo_abundance_EL_GT25_sub$YEAR) #nice 

luquillo_abundance_EL_GT25_sub %>%
  filter(YEAR == "1992") %>%
  distinct(DATE, PLOT1)

luquillo_abundance_EL_LT25_sub <- luquillo_abundance_EL_LT25 %>%
  filter(!YEAR %in% c(1989, 1990)) %>%
  select(-c(WEATHER, WIND, RAIN, TIME, JULIAN, OBSERVER))


#pull out species list 
luquillo_abundance_EL_GT25_sub


unique(luquillo_abundance_EL_LT25_sub$YEAR) #nice

id_cols <- c("YEAR", "PLOT1", "PLOT2", "DATE", "PLACE")   # add others you have (SITE, DATE, TRANSECT, etc.)

luquillo_abundance_EL_GT25_sub %>%
  filter(YEAR == "2012") %>%
  filter(PLOT1 == "12")



luquillo_abundance_EL_LT25_sub %>%
  filter(YEAR == "2012") %>%
  filter(PLOT1 == "12")

luquillo_abundance_EL_GT25_LONG <- luquillo_abundance_EL_GT25_sub %>%
  pivot_longer(
    cols = -all_of(id_cols),
    names_to = "species_code",
    values_to = "abundance",
  ) %>%
  mutate(abundance = replace_na(abundance, 0)) %>%
  mutate(dataset = "GT25")

glimpse(luquillo_abundance_EL_GT25_LONG)

luquillo_abundance_EL_LT25_LONG <- luquillo_abundance_EL_LT25_sub %>%
  pivot_longer(
    cols = -all_of(id_cols),
    names_to = "species_code",
    values_to = "abundance"
  ) %>%
  mutate(abundance = replace_na(abundance, 0)) %>%
  mutate(dataset = "LT25")


luquillo_abundance_EL_all_LONG <- rbind(luquillo_abundance_EL_GT25_LONG, luquillo_abundance_EL_LT25_LONG)

uni

luquillo_abundance_EL_all_LONG_above0 <- luquillo_abundance_EL_all_LONG %>%
  dplyr::filter(abundance > 0) %>%
  distinct(species_code) 

luquillo_abundance_EL_all_LONG_all<- luquillo_abundance_EL_all_LONG %>%
  distinct(species_code) 

setdiff(luquillo_abundance_EL_all_LONG_all$species_code, luquillo_abundance_EL_all_LONG_above0$species_code)
# "BFGRQ" "KENWA" "PRAIW" "WEVIR" "YFGRQ" --> never appear? 


luquillo_abundance_EL_all_LONG_all_specieslist <- luquillo_abundance_EL_all_LONG_all %>%
  filter(!species_code %in% c("UNKNO","UNWAR", "WATTH", "TOTAL"))

luquillo_birdspecieslist_formerge <- luquillo_abundance_EL_all_LONG_all_specieslist %>%
  mutate(project = "Luquillo",
         habitat = "terrestrial",
         raw_filename = "2filesmergedtogetherneedtofigureout",
         class = "Aves")

luquillo_birdspecieslist_formerge_commonname <- luquillo_birdspecieslist_formerge %>%
  left_join(luquillo_birdcode, by = "species_code")

luquillo_birdspecieslist_formerge_commonname %>%
  dplyr::filter(is.na(common_name)) #0 :) 

#update spelling mistakes and some taxonomic changes and changing common names to match bird base... 
luquillo_birdspecieslist_formerge_commonname_edits <- luquillo_birdspecieslist_formerge_commonname %>%
  mutate(
    common_name = case_when(
      common_name == "Plain Pegeon" ~ "Plain Pigeon", #spelling error in data??
      common_name == "Blue-hooded Euphonia" ~ "Elegant Euphonia", #common name update for matching 
      common_name == "Black and white Warbler" ~ "Black-and-white Warbler", #hyphen or capitalization update for matching
      common_name == "Black-cowled oriole" ~ "Black-cowled Oriole", #hyphen or capitalization update for matching
      common_name == "Parula Warbler" ~ "Northern Parula",#common name update for matching 
      common_name == "Puerto Rico Lizard Cuckoo" ~ "Puerto Rican Lizard-Cuckoo", #hyphen or capitalization update for matching
      common_name == "Red-legged Thrush" ~ "Western Red-legged Thrush", ##common name update for matching 
      common_name == "Puerto Rican Screech Owl" ~ "Puerto Rican Owl", #common name update for matching 
      common_name == "Stripe-headed Tanager" ~ "Western Spindalis", #taxonomic change
      TRUE ~ common_name
    )
  )

luquillo_birdspecieslist_formerge_taxo <- luquillo_birdspecieslist_formerge_commonname_edits %>%
  left_join(taxo_edbird_limited_commonname, by = c("common_name" = "PRIMARY_COM_NAME")) %>%
  rename("order" = "ORDER",
         "family" = "FAMILY",
         "scientific_name" = "SCI_NAME") %>%
  mutate(genus = word(scientific_name, 1))

luquillo_birdspecieslist_final <- luquillo_birdspecieslist_formerge_taxo %>%
  select(project, habitat, raw_filename, class, order, family, genus, scientific_name)
  

#fix final one by hand, then split scientifi name into gensu



luquillo_birdspecieslist_formerge_taxo %>%
  dplyr::filter(is.na(scientific_name)) #missing 0 :) 

#

#write this out for Shalanda

#andrewsforest <- "andrewsforest_specieslist.csv"

# Export





#above 0 = 41 --> 
#all = 46 --> which species never appear? 
  
luquillo_abundance_EL_all_LONG %>%
 # filter(!SPECIES %in% c("PIKA","DOSQ", "UNMA", "CHIP", "UNKN", "UNFI", "UNFL", 
          #               "UNGR", "UNJA", "UNTH", "UNVI", "UNWA", "UNWO", "UNWR", "UNCH",
           #              "UNCR", "UNKI", "UNSW", "UNOW", "UNHM", "UNDU", "UNRA", "UNCO", "UNHU")) %>% #get rid of mammals and unknown things
  distinct(species_code)


#46 10 +36
  
  df_na_species <- andrews_ex_abundance_codes %>%
  dplyr::filter(is.na(species_name)) %>%
  distinct(SPECIES)
  
  

colnames(luquillo_abundance_EL_all_LONG)
luquillo_abundance_EL_final_LONG_abundance <-luquillo_abundance_EL_all_LONG %>%
  group_by(YEAR, PLACE, PLOT1, PLOT2, species_code, DATE) %>%
  summarise(abundance_total = sum(abundance))

luquillo_abundance_biotime_clean %>%
  filter(YEAR == 2008) %>%
 # filter(date_str == 2006) %>%
  filter(plot == "plot15") %>%
  filter(valid_name == "Spindalis zena")
 
unique(luquillo_abundance_biotime_clean$plot)
#2 sampling data

unique(luquillo_abundance_EL_final_LONG_abundance$abundance_total)

luquillo_abundance_EL_final_LONG_abundance %>%
  filter(YEAR == 2008) %>%
  filter(PLOT1 == 15) %>%
  filter(species_code == "WWDOV") %>%
  filter(abundance_total > 0)


unique(luquillo_abundance_EL_final_LONG_abundance$species_code)

#### CAP ####

cap_main <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "cap", "46_bird_observations.csv"))
cap_saltriver <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "cap", "641_bird_observations.csv"))


cap_main_commonnames <- cap_main %>%
  distinct(common_name) 
#310
#212

cap_saltriver_commonnames <-cap_saltriver %>%
  distinct(common_name)
#522

cap_allbirds_dups <- rbind(cap_main_commonnames, cap_saltriver_commonnames)

cap_allbirds <- cap_allbirds_dups %>%
  distinct(common_name) %>%
  filter(!str_detect(common_name, "Unidentified"))


#cool, seems like 313 unique birds across both
#51 w/ "unidentfied in them"
#313-51 = 262


setdiff(cap_saltriver_commonnames,cap_main_commonnames) #3 diff
setdiff(cap_main_commonnames, cap_saltriver_commonnames) #101


#let's try to merge..


caps_birdspecieslist_formerge <- cap_allbirds %>%
  mutate(project = "CAP",
         habitat = "terrestrial",
         raw_filename = "twodatasetswillneedtofigurethatout",
         class = "Aves")

caps_birdspecieslist_formerge_commonnames_edits <- caps_birdspecieslist_formerge %>%
  mutate(
    common_name =   case_when(
      common_name == "Audubon's Warbler" ~ "Yellow-rumped Warbler (Audubon's)", #common name chage to help with merging.
      common_name == "Oregon Junco" ~ "Dark-eyed Junco (Oregon)",#NEED TO TALK TO EXPERTS
      common_name == "Sage Sparrow" ~ "Sagebrush Sparrow", #NEED TO TALK TO EXPERTS
      common_name == "Common Ground-Dove" ~ "Common Ground Dove", #hyphen or capitalization update for matching
      common_name == "Black-crowned Night-Heron" ~ "Black-crowned Night Heron", #hyphen or capitalization update for matching
      common_name == "Cattle Egret" ~ "Western Cattle-Egret",#NEED TO TALK TO EXPERTS
      common_name == "Slate-colored Junco" ~ "Dark-eyed Junco (Slate-colored)", #need to check with experts 
      common_name == "Domestic Duck" ~ "Mallard", #need to check with experts 
      common_name == "American Green-winged Teal" ~ "Green-winged Teal (American)", #hyphen or capitalization update for matching
      common_name == "Cordilleran Flycatcher" ~ "Western Flycatcher (Cordilleran)", #need to check with experts
      common_name == "Gray-headed Junco" ~ "Dark-eyed Junco (Gray-headed)",  #need to check with experts
      common_name == "Yellow Warbler" ~ "Mangrove Yellow Warbler",  #need to check with experts
      common_name == "Pacific Slope Flycatcher" ~ "Western Flycatcher", #need to check with experts
      common_name == "Black-bellied Whistling Duck" ~ "Black-bellied Whistling-Duck", #hyphen or capitalization update for matching
      common_name == "House Wren" ~ "Northern House Wren", #need to check with experts
      common_name == "Pink-sided Junco" ~ "Dark-eyed Junco (Pink-sided)", #should check with experts 
      common_name == "Warbling Vireo" ~ "Eastern Warbling Vireo", ##should check with experts 
      common_name == "Chinese Goose" ~ "Swan Goose", # ##should check with experts 
      common_name == "Domestic Goose" ~ "Graylag Goose", ##should check with experts 
      common_name == "Myrtle Warbler" ~ "Yellow-rumped Warbler", ##should check with experts 
      common_name == "Common Canary" ~ "Island Canary", ##should check with experts 
      common_name == "Ridgeway's Rail" ~ "Ridgway's Rail", ##spelling update for matching
      common_name == "Barn Owl" ~ "Western Barn Owl", ##should check with experts 
      TRUE ~ common_name
    )
  )

cap_allbirds_taxo <- caps_birdspecieslist_formerge_commonnames_edits %>%
  left_join(taxo_edbird_limited_commonname, by = c("common_name" = "PRIMARY_COM_NAME")) %>%
  rename("order" = "ORDER",
         "family" = "FAMILY",
         "scientific_name" = "SCI_NAME") %>%
  mutate(genus = word(scientific_name, 1))
  
cap_allbirds_taxo %>%
  dplyr::filter(is.na(SCI_NAME))  #0
  


CAP_birdspecieslist_final <- cap_allbirds_taxo %>%
  select(project, habitat, raw_filename, class, order, family, genus, scientific_name)

df_na_taxo <- andrews_ex_birdspecieslist_formerge_taxo %>%
  dplyr::filter(is.na(order)) #cool, 0 

andrews_ex_birdspecieslist_final %>%
  dplyr::filter(is.na(order))

#### Arctic LTER ####

arctic_abundance_birdlist <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "arctic_lter_birds.csv"))

arctic_abundance_birdlist_commonname_edits <- arctic_abundance_birdlist %>%
  mutate(
    common_name  = case_when(
      common_name  == "American Golden Plover" ~ "American Golden-Plover", #hyphen or capitalization update for matching
      common_name  == "American Pipet" ~ "American Pipit", #spelling change for matching
      common_name  == "Common Redpoll" ~ "Redpoll", #hyphen or capitalization update for matching
      common_name  == "Whimbrel" ~ "Eurasian Whimbrel", #need to check with experts
      TRUE ~ common_name
    )
  )

arctic_abundance_birdlist_taxo <- arctic_abundance_birdlist_commonname_edits %>%
  left_join(taxo_edbird_limited_commonname, by = c("common_name" = "PRIMARY_COM_NAME")) %>%
  rename("order" = "ORDER",
         "family" = "FAMILY",
         "scientific_name" = "SCI_NAME") %>%
  mutate(genus = word(scientific_name, 1),
         project = "Arctic",
         habitat = "terrestrial",
         raw_filename = "onedatasetperyearneedtofigurethisout",
         class = "Aves")


arctic_abundance_birdlist_taxo%>%
  dplyr::filter(is.na(scientific_name)) # 0 - stick 

arctic_abundance_birdlist_taxo_final <- arctic_abundance_birdlist_taxo %>%
  select(project, habitat, raw_filename, class, order, family, genus, scientific_name)

#### SONGS ####

#### Arctic LTER ####

songs_birdlist <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "SONGS_birdlist.csv"))
songs_birdlist_raw <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "SONGS_birdlist_total.csv"))

songs_birdlist_sciname1 <- songs_birdlist %>%
  select(scientific_name)

colnames(songs_birdlist)

songs_birdlist_cleaned <- songs_birdlist_raw %>%
  mutate(scientific_name_new = str_extract(scientific_name, "(?<=:)\\s*\\w+\\s+\\w+")) %>%
  mutate(scientific_name_new = str_trim(scientific_name_new)) 
  
songs_birdlist_cleaned <- songs_birdlist_cleaned[!is.na(songs_birdlist_cleaned$scientific_name_new), ]

unique(songs_birdlist$scientific_name)

songs_birdlist_sciname2 <- songs_birdlist_cleaned %>%
  select(scientific_name_new) %>%
  rename(scientific_name = scientific_name_new)

songs_birdlist_sciname_all <- rbind(songs_birdlist_sciname1, songs_birdlist_sciname2) #491

songs_birdlist_sciname_all_formerge <- songs_birdlist_sciname_all %>%
  distinct(scientific_name) %>%
  filter(!str_detect(scientific_name, "Unidentified"))

songs_abundance_birdlist_taxo <- songs_birdlist_sciname_all_formerge %>%
  left_join(taxo_edbird_limited_commonname, by = c("scientific_name" = "SCI_NAME")) %>%
  rename("order" = "ORDER",
         "family" = "FAMILY") %>%
  mutate(genus = word(scientific_name, 1),
         project = "SONGS",
         habitat = "terrestrial",
         raw_filename = "willneedtofigurethisout",
         class = "Aves")

songs_abundance_birdlist_taxo %>%
  filter(!is.na(ORDER))

unique(songs_abundance_birdlist_taxo$scientific_name)

songs_abundance_birdlist_taxo_final <- songs_abundance_birdlist_taxo %>%
  select(project, habitat, raw_filename, class, order, family, genus, scientific_name)


#188 species

#### STACK TOGETHEER ####

#andrews: 77
#luqillo: 42
#hubbard: 33
#cap: 262
#art: 20 
#songs: 184 

colnames(andrews_ex_birdspecieslist_final) #8 cols
colnames(luquillo_birdspecieslist_final) #8 cols
colnames(hubbardbrook_abundance_birdlist_taxo_final) #8

birdlist_LNEv1 <- rbind(andrews_ex_birdspecieslist_final, 
                      luquillo_birdspecieslist_final, 
                      hubbardbrook_abundance_birdlist_taxo_final,
                      CAP_birdspecieslist_final,
                      arctic_abundance_birdlist_taxo_final,
                      songs_abundance_birdlist_taxo_final) 

#434?
#618
#119, nice 



write.csv(birdlist_LNEv1, row.names = F, na = '',
          file = file.path("sandbox", "birds_specieslist", "LNE_specieslist_03042026_v3.csv"))

  