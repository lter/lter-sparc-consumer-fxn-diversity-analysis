###

### IMPORTANT: requires running the script LDLA_dataexploration.R


communities <- read.csv(file = "Data/community_tidy-data/04_harmonized_consumer_excretion_sparc_cnd_site.csv") 

comm.filtered.abu <- communities %>%
  select(project, year, month, habitat, temp_c, site, subsite_level1, subsite_level2, subsite_level3,
         scientific_name, diet_cat, nind_ug.hr, pind_ug.hr, count_num, density_num.m, density_num.m2,
         density_num.m3,biomass_g, dmperind_g.ind, kingdom, phylum, class, order, family, genus) %>%
  mutate(taxa = case_when(
    project == "CoastalCA" ~ "Fish",
    project == "FCE" ~ "Fish",
    project == "SBC" ~ "Fish",
    project == "MCR" ~ "Fish",
    project == "VCR" ~ "Fish",
    project == "RLS" ~ "Fish",
    project == "FISHGLOB" ~ "Fish",
    project == "KBS_MAM" ~ "Mammals",
    project == "SEV" ~ "Mammals",
    project == "MOHONK" ~ "Amphibians",
    project == "KBS_AMP" ~ "Amphibians",
    project %in% c("HARVARD", "KBS_BIR","SBC_BEACH") ~ "Birds",
    project %in% c("NGA","Arctic","Palmer","CCE","NorthLakes") ~ "Zooplankton",
    TRUE ~ NA
  )) %>%
  filter(habitat != "beach") %>%
  filter(taxa == "Fish") %>%
  mutate(project = case_when(project == 'CoastalCA' & site == 'CENTRAL' ~ 'COASTAL_CEN',
                             project == 'CoastalCA' & site == 'SOUTH' ~ 'COASTAL_SOUTH',
                             TRUE ~ project)) %>%
  mutate(density_num.m = case_when(is.na(density_num.m) ~ 0,
                                   TRUE ~ density_num.m),
         density_num.m2 = case_when(is.na(density_num.m2) ~ 0,
                                    TRUE ~ density_num.m2)) %>%
  mutate(density_coll = density_num.m+density_num.m2) %>%
  filter(density_coll > 0) %>%
  mutate(biomass_coll = dmperind_g.ind*density_coll) %>%
  unite("ID", c(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                year, month), sep = "_", remove = F) %>%
  group_by(ID, project, scientific_name) %>%
  summarize(density = sum(density_coll, na.rm = TRUE),
            biomass = sum(biomass_coll)) %>%
  group_by(project) %>%
  filter(!is.na(density), !is.na(biomass)) %>%
  mutate(z_density = scale(log(density))[,1],
         z_biomass = scale(log(biomass+1))[,1]) %>%
  ungroup()









comm.filtered.abu <- communities %>%
  select(project, year, month, habitat, temp_c, site, subsite_level1, subsite_level2, subsite_level3,
         scientific_name, diet_cat, nind_ug.hr, pind_ug.hr, count_num, density_num.m, density_num.m2,
         density_num.m3,biomass_g, dmperind_g.ind, kingdom, phylum, class, order, family, genus) %>%
  mutate(taxa = case_when(
    project == "CoastalCA" ~ "Fish",
    project == "FCE" ~ "Fish",
    project == "SBC" ~ "Fish",
    project == "MCR" ~ "Fish",
    project == "VCR" ~ "Fish",
    project == "RLS" ~ "Fish",
    project == "FISHGLOB" ~ "Fish",
    project == "KBS_MAM" ~ "Mammals",
    project == "SEV" ~ "Mammals",
    project == "MOHONK" ~ "Amphibians",
    project == "KBS_AMP" ~ "Amphibians",
    project %in% c("HARVARD", "KBS_BIR","SBC_BEACH") ~ "Birds",
    project %in% c("NGA","Arctic","Palmer","CCE","NorthLakes") ~ "Zooplankton",
    TRUE ~ NA
  )) %>%
  filter(habitat != "beach") %>%
  filter(taxa == "Fish") %>%
  mutate(project = case_when(project == 'CoastalCA' & site == 'CENTRAL' ~ 'COASTAL_CEN',
                             project == 'CoastalCA' & site == 'SOUTH' ~ 'COASTAL_SOUTH',
                             TRUE ~ project)) %>%
  mutate(density_num.m = case_when(is.na(density_num.m) ~ 0,
                                   TRUE ~ density_num.m),
         density_num.m2 = case_when(is.na(density_num.m2) ~ 0,
                                   TRUE ~ density_num.m2)) %>%
  mutate(density_coll = density_num.m+density_num.m2) %>%
  filter(density_coll > 0) %>%
  mutate(biomass_coll = dmperind_g.ind*density_coll) %>%
  unite("ID", c(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                year, month), sep = "_", remove = F) %>%
  group_by(ID, project, scientific_name) %>%
  summarize(density = sum(density_coll, na.rm = TRUE),
            biomass = sum(biomass_coll)) %>%
  group_by(project) %>%
  filter(!is.na(density), !is.na(biomass)) %>%
  mutate(z_density = scale(log(density))[,1],
         z_biomass = scale(log(biomass+1))[,1]) %>%
  ungroup()

test <- comm.filtered %>%
  filter(project == "FCE")

density.plot <- ggplot(comm.filtered, aes(x = z_biomass)) +
  geom_histogram(aes(fill = project)) +
  facet_grid(project~., scales = "free_y")
density.plot


traits.df.parsed <- taxa.df %>%
  ungroup() %>%
  select(scientific_name, PC1, PC2)

traitsxcomms <- comm.filtered %>%
  ungroup() %>%
  inner_join(traits.df.parsed, by = "scientific_name") %>%
  distinct()

anti <- comm.filtered %>%
  anti_join(taxa.df)

### function for KDE

weighted_kde2d <- function(x, y, w, n = 200, bandwidth_factor = 6) {
  gx <- seq(min(x), max(x), length.out = n)
  gy <- seq(min(y), max(y), length.out = n)
  
  # broader bandwidths ⇒ smoother, prettier blobs
  hx <- bw.nrd(x) * bandwidth_factor
  hy <- bw.nrd(y) * bandwidth_factor
  
  z <- matrix(0, nrow = n, ncol = n)
  
  for (i in seq_along(x)) {
    z <- z + w[i] * outer(
      dnorm(gx, mean = x[i], sd = hx),
      dnorm(gy, mean = y[i], sd = hy)
    )
  }
  
  # normalize to [0, 1]
  z <- z / max(z, na.rm = TRUE)
  
  list(x = gx, y = gy, z = z)
}


#### create dataframe
df_long <- traitsxcomms %>%
  filter(!is.na(PC1), !is.na(PC2)) 

### run KDE function for each project
kde_by_proj <- df_long %>%
  group_by(project) %>%
  group_modify(~ {
    x <- .x$PC1
    y <- .x$PC2
    w <- .x$density
    
    kde <- weighted_kde2d(x, y, w, n = 200)
    
    expand.grid(
      PC1 = kde$x,
      PC2 = kde$y
    ) %>%
      mutate(density = as.vector(kde$z))
  }) %>%
  ungroup()


### plot for each project
list_of_dfs_tidy <- kde_by_proj %>%
  group_by(project) %>%
  group_split()



### create contour plots

contour_plots1 <- ggplot(list_of_dfs_tidy[[1]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[1]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots1


contour_plots2 <- ggplot(list_of_dfs_tidy[[2]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[2]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots2

contour_plots3 <- ggplot(list_of_dfs_tidy[[3]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[3]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots3

contour_plots4 <- ggplot(list_of_dfs_tidy[[4]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[4]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots4
x <- list_of_dfs_tidy[[5]]
contour_plots5 <- ggplot(list_of_dfs_tidy[[5]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 10,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[5]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots5

contour_plots6 <- ggplot(list_of_dfs_tidy[[6]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 10,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[6]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots6

contour_plots7 <- ggplot(list_of_dfs_tidy[[7]], aes(PC1, PC2, z = density)) +
  geom_contour(
    aes(z = density,
        color = after_stat(level)),
    bins = 10,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[7]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots7

contour.plots.density <- (contour_plots1 + contour_plots2 + contour_plots3) /  (contour_plots4 + contour_plots5 + contour_plots6) + plot_layout(guides = "collect")

contour.plots.density








### run KDE function for each project
kde_by_proj <- df_long %>%
  group_by(project) %>%
  group_modify(~ {
    x <- .x$PC1
    y <- .x$PC2
    w <- .x$biomass
    
    kde <- weighted_kde2d(x, y, w, n = 200)
    
    expand.grid(
      PC1 = kde$x,
      PC2 = kde$y
    ) %>%
      mutate(biomass = as.vector(kde$z))
  }) %>%
  ungroup()


### plot for each project
list_of_dfs_tidy <- kde_by_proj %>%
  group_by(project) %>%
  group_split()



### create contour plots

contour_plots1 <- ggplot(list_of_dfs_tidy[[1]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[1]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots1


contour_plots2 <- ggplot(list_of_dfs_tidy[[2]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[2]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots2

contour_plots3 <- ggplot(list_of_dfs_tidy[[3]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[3]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots3

contour_plots4 <- ggplot(list_of_dfs_tidy[[4]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 8,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[4]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )

contour_plots4
x <- list_of_dfs_tidy[[5]]
contour_plots5 <- ggplot(list_of_dfs_tidy[[5]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 10,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[5]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots5

contour_plots6 <- ggplot(list_of_dfs_tidy[[6]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 10,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[6]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots6

contour_plots7 <- ggplot(list_of_dfs_tidy[[7]], aes(PC1, PC2, z = biomass)) +
  geom_contour(
    aes(z = biomass,
        color = after_stat(level)),
    bins = 10,
    size = 0.5
  ) +
  scale_color_viridis() +
  geom_point(data = taxa.df, aes(x = PC1, y = PC2), inherit.aes = FALSE, color = "grey", alpha = 0.3) +
  labs(
    title = list_of_dfs_tidy[[7]]$project,
    x = "PC1",
    y = "PC2",
    fill = "Biomass"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10)
  )
contour_plots7

contour.plots.biomass <- (contour_plots1 + contour_plots2 + contour_plots3) /  (contour_plots4 + contour_plots5 + contour_plots6) + plot_layout(guides = "collect")

contour.plots.biomass


contour.plots.all <- contour.plots.density | contour.plots.biomass

ggsave(contour.plots.all, file = "contour.plots.all.png", width = 16, height = 10)














model1 <- brm(abundance ~ PC1 + PC2 + (1|project), family = "negbinomial", data = traitsxcomms,
              cores = 4, chains = 4, backend = "cmdstanr")
summary(model1)
plot(model1)

params <- model1 %>%
  spread_draws(b_Intercept, b_PC1, b_PC2) %>%
  pivot_longer(starts_with("b"), names_to = "Parameter", values_to = "value")

param.plot <- ggplot(params, aes(y = Parameter, x = value)) +
  geom_halfeyeh(fill = "#4B0082", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()
param.plot

m1p1 <- conditional_effects(model1, resolution = 1000)
df1.1 <- as.data.frame(m1p1[[1]])

plot1.1 <- ggplot(df1.1) +
  geom_line(aes(x = effect1__, y = estimate__), color = "black", lwd = 1.2) +
  geom_ribbon(aes(x= effect1__, ymin = lower__, ymax = upper__), color = "black",fill = alpha("#4B0082", 0.3), lty = 2, lwd = 0.3) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(color="black", size = 8),
        axis.title=element_text(color="black", size = 8),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black")) 
plot1.1


model2 <- brm(biomass ~ PC1 + PC2 + (1|project), family = "gamma", data = traitsxcomms,
              cores = 4, chains = 4, backend = "cmdstanr")
summary(model2)
plot(model2)

params2 <- model2 %>%
  spread_draws(b_Intercept, b_PC1, b_PC2) %>%
  pivot_longer(starts_with("b"), names_to = "Parameter", values_to = "value")

param.plot2 <- ggplot(params2, aes(y = Parameter, x = value)) +
  geom_halfeyeh(fill = "gold3", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_bw()
param.plot2

m1p2 <- conditional_effects(model2, resolution = 1000)
df1.2 <- as.data.frame(m1p2[[2]])

plot1.2 <- ggplot(df1.2) +
  geom_line(aes(x = effect1__, y = estimate__), color = "black", lwd = 1.2) +
  geom_ribbon(aes(x= effect1__, ymin = lower__, ymax = upper__), color = "black",fill = alpha("gold3", 0.3), lty = 2, lwd = 0.3) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text=element_text(color="black", size = 8),
        axis.title=element_text(color="black", size = 8),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black")) 
plot1.2

param.plots <- param.plot + param.plot2