
librarian::shelf(lter/ltertools)

key <- read.csv(file = file.path("data", "-keys", "community_datakey.csv"))

list_raw <- ltertools::read(raw_folder = file.path("data", "terrestrial_community_raw-data"), 
                            data_format = "csv")


#investigating if these birds are dead or alive?
#sbc_v1 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "Shorebird_count_20231020.csv"))
#colnames(sbc_v1)
#unique(sbc_v1$SURVEY)
#sbc_v1 %>%
 # filter(SURVEY == ".")

sev_v1 <- read.csv(file = file.path("data", "terrestrial_community_raw-data", "SEV008_long.csv"))
unique(sev_v1$species)

sev_v1 %>%
  

#colnames(sbc_v1)


# Make a list to store standardized outputs
list_std <- list()

namelist <- sort(intersect(x = key$source, y = names(list_raw)))

for (focal_src in sort(intersect(x = key$source, y = names(list_raw)))){
  
  #focal_src= "FISHGLOB_NorthSea_EnglishChannel.csv"
  # Progress message
  message("Standarizing file: '", focal_src, "'")
  
  
  # Standardize this file
  focal_v1 <- ltertools::standardize(focal_file = focal_src, 
                                     key = key, 
                                     df_list = list_raw)
  
  
  #  Identify valid (non-blank) column names
  nms  <- names(focal_v1)
  good <- !is.na(nms) & nzchar(str_trim(nms))
  
  # Keep only those columns (avoids passing "" anywhere)
  focal_v1_clean <- focal_v1 %>%
    select(all_of(nms[good]))
  
  # Do some bonus processing if taxa are in wide format
  if(any(stringr::str_detect(string = names(focal_v1_clean), pattern = "orig_taxa_"))){
    
    # Flip it to long format & tidy up taxa names
    focal_v2 <- focal_v1_clean %>% 
      tidyr::pivot_longer(cols = dplyr::starts_with("orig_taxa_"),
                          names_to = "species",
                          values_to = "density") %>% 
      dplyr::mutate(species = gsub(pattern = "orig_taxa_", replacement = "", x = species))
    
  } else { focal_v2 <- focal_v1_clean }
  
  # Make a final object
  focal_std <- focal_v2
  
  # Add to output list
  list_std[[focal_src]] <- focal_std
  
} # Close loop


# Unlist outputs 
combo_v1 <- purrr::list_rbind(x = list_std)

# Check structure
dplyr::glimpse(combo_v1)


## --------------------------- ##
# Add on Key Metadata ----
## --------------------------- ##

# Grab some important metadata stored in the key
key_meta <- key %>% 
  dplyr::select(project, data_type, habitat, source) %>% 
  dplyr::distinct()

# Check that out
key_meta

# Attach that to the combined data using the 'source' column
combo_v2 <- combo_v1 %>% 
  dplyr::left_join(y = key_meta, by = "source") %>% 
  dplyr::relocate(project:habitat, .before = source)

# Check structure
dplyr::glimpse(combo_v2)

#Fix SBC Bird Data
#has sea lions and separate columns for genus and species
combo_v2_SBC <- combo_v2 %>%
  filter(project == "SBC") %>%
  filter(class != "Mammalia") %>%
  dplyr::mutate(scientific_name = paste(genus, species)) %>%
  select(-species) 

#class -- . = no birds 
# 
combo_v2_SBC %>%
  filter(is.na(class))

unique(combo_v2_SBC$class)

colnames(combo_v2_SBC)
glimpse(combo_v2_SBC)
unique(combo_v2_SBC$site) 
#6 sites 
colnames(combo_v2_SBC)
unique(combo_v2_SBC$survey) 

combo_v2_SEV <- combo_v2 %>%
  filter(project == "SEV") %>%
  filter(class != "Mammalia") %>%
  dplyr::mutate(scientific_name = paste(genus, species)) %>%
  select(-species) 

combo_v2_allothers <- combo_v2 %>%
  filter(project != "SBC") %>%
  filter(project != "SEV")  %>% #taking out SEV because it's fucked
  rename(scientific_name = species) # rename for merging to match SBC
  
unique(combo_v2_allothers$project)

combo_v3 <- rbind(combo_v2_allothers, combo_v2_SBC)

combo_v3 <-combo_v3 %>%
  rename(species = scientific_name)  #rename back to match aquatic community data


#don't need this step for this data as of rn, 
# if there is a date column but the year column is NA, extract year from date, same for month and day
#combo_v3 <- combo_v2 %>%
 # dplyr::mutate(year = ifelse(is.na(year) & !is.na(date),
                   #           lubridate::year(lubridate::ymd(date)),
                      #        year),
               # month = ifelse(is.na(month) & !is.na(date),
                   #            lubridate::month(lubridate::ymd(date)),
                    #           month),
              #  day = ifelse(is.na(day) & !is.na(date),
                          #   lubridate::day(lubridate::ymd(date)),
                         #    day))

## --------------------------- ##
# Export ----
## --------------------------- ##

# Make a final object
combo_v99 <- combo_v3

# Double check its structure
dplyr::glimpse(combo_v99)

# Make a filename for it
combo_file_terrestrial <- "01_terrestrial_community_harmonized.csv"

# Export
write.csv(x = combo_v99, row.names = F, na = '',
          file = file.path("data", "community_tidy-data", combo_file_terrestrial))

# End ----
