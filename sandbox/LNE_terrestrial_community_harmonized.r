
# read in harmonized comm data 

#com_dt<- read.csv(file = file.path("Data", "community_tidy-data", "01_community_harmonized.csv"))

# for now

com_dt <- combo_v99

#filter

unique(com_dt$project)

hubbard_brook <- com_dt %>%
  filter(project == "Hubbard_Brook") %>%
  mutate(site = "Hubb_all") # site is currently the year, because only one site for this project 

unique(hubbard_brook$site)


palmer <- com_dt %>%
  filter(project == "Palmer") %>%
  mutate(site = sub("^(?:[^_]*_){3}", "", site)) #19 sites

glimpse(palmer)

unique(palmer$site)


#skipping for now, only goes through 2026
luquillo <- com_dt %>%
  filter(project == "Luquillo")

unique(luquillo$year)

sbc_birds <- com_dt %>%
  filter(project == "SBC")

unique(sbc_birds$site) #already normal, I love SBC!

CAP <- com_dt %>%
  filter(project == "CAP") %>%
  separate(site,
           into = c("a","b","c", "d", "site","t1","t2"),
           sep = "_",
           remove = FALSE) %>%
  select(-a,-b,-c,-d, -t1,-t2)

unique(CAP$site) #37 sites?

comm_terrest_data <- rbind(CAP, sbc_birds, palmer, hubbard_brook)

colnames(CAP)
colnames(sbc_birds)

glimpse(comm_terrest_data)

combo_file_terrestrial_wrang <- "02_terrestrial_community_wrangled.csv"

write.csv(x = comm_terrest_data, row.names = F, na = '',
          file = file.path("data", "terrestrial_community_tidy-data", combo_file_terrestrial_wrang))
