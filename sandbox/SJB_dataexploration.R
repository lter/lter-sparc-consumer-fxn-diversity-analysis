###

### IMPORTANT: requires running the script LDLA_dataexploration.R


communities <- read.csv(file = "Data/community_tidy-data/04_harmonized_consumer_excretion_sparc_cnd_site.csv") 

comm.filtered.pa <- communities %>%
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
  mutate(scientific_name = str_remove(scientific_name, " spp.")) %>%
  mutate(project = case_when(project == 'CoastalCA' & site == 'CENTRAL' ~ 'COASTAL_CEN',
                             project == 'CoastalCA' & site == 'SOUTH' ~ 'COASTAL_SOUTH',
                             TRUE ~ project)) %>%
  mutate(density_num.m = case_when(is.na(density_num.m) ~ 0,
                                   TRUE ~ density_num.m),
         density_num.m2 = case_when(is.na(density_num.m2) ~ 0,
                                    TRUE ~ density_num.m2)) %>%
  mutate(density_coll = density_num.m+density_num.m2) %>%
  mutate(pres.abs = case_when(density_coll > 0 ~ 1,
                              TRUE ~ 0)) %>%
  unite("ID", c(project, habitat, site, subsite_level1, subsite_level2, subsite_level3,
                year, month), sep = "_", remove = F) %>%
  filter(pres.abs > 0)


traits.df.parsed <- taxa.df.t %>%
  ungroup() %>%
  select(scientific_name, PC1, PC2)

traitsxcomms <- comm.filtered.pa %>%
  ungroup() %>%
  inner_join(traits.df.parsed, by = "scientific_name") %>%
  distinct()

summary(traitsxcomms)

centroids <- traitsxcomms %>%
  filter(pres.abs == 1) %>%
  unite("site_fine", c(project, site, subsite_level1), remove = F, sep = "_") %>%
  group_by(project, site, site_fine, year) %>%              # ID uniquely identifies site-year
  summarise(
    centroid_PC1 = mean(PC1, na.rm = TRUE),
    centroid_PC2 = mean(PC2, na.rm = TRUE),
    n_species = n(),                        # nice to keep track
    .groups = "drop"
  )

### calculate trajectores and step lengths using Pythagorean theorem
centroids.traj <- centroids %>%
  arrange(project, site, site_fine, year) %>%
  group_by(project, site, site_fine) %>%
  mutate(
    lag_PC1 = lag(centroid_PC1),
    lag_PC2 = lag(centroid_PC2),
    dPC1 = centroid_PC1 - lag_PC1,
    dPC2 = centroid_PC2 - lag_PC2,
    step_length = sqrt(dPC1^2 + dPC2^2)      # movement through space
  ) %>%
  mutate(angle = atan2(dPC2, dPC1)) %>%
  ungroup()


### plot steps
centroid.plot <- ggplot(centroids, aes(centroid_PC1, centroid_PC2, group = site_fine, color = project)) +
  geom_path(alpha = 0.95) +
  geom_point(aes(fill = year), size = 2, shape = 21, color = "black") +
  coord_equal() +
  scale_color_viridis_d() +
  theme_bw(base_size = 14) +
  facet_wrap(.~fct_reorder(site, project), ncol = 9) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Site",
    fill = "Year",
    title = "Community centroid shifts through functional space"
  )
centroid.plot


angle.plot <- ggplot(centroids.traj, aes(year, angle, group = site_fine, color = project)) +
  geom_line() +
  scale_y_continuous(
    breaks = c(-pi, -pi/2, 0, pi/2, pi),
    labels = c("−π", "−π/2", "0", "π/2", "π")
  ) +
  facet_wrap(.~site) +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Direction of movement",
    title = "Year-to-year direction of functional centroid shifts"
  )

angle.plot

centroid.tests <- centroids.traj %>%
  filter(!is.na(angle)) %>%
  group_by(site_fine) %>%
  summarise(
    p = circular::rayleigh.test(circular(angle))$p.value,
    mean_dir = as.numeric(mean(circular(angle))),
    .groups = "drop"
  ) %>%
  arrange(p)
centroid.tests

rosette.plot <- ggplot(centroids.traj, aes(angle, fill = project)) +
  geom_histogram(binwidth = pi/6) +
  coord_polar() +
  facet_wrap(~ site, scales = "free_y", ncol = 9) +
  theme_bw() +
  scale_fill_viridis_d()
  
rosette.plot

centroid.traj.dir <- centroids.traj %>%
  filter(!is.na(dPC1)) %>%
  pivot_longer(cols = c(dPC1, dPC2),
               names_to = "axis",
               values_to = "delta")

centroid.traj.plot <- ggplot(centroid.traj.dir, aes(axis, abs(delta), fill = project)) +
  geom_boxplot() +
  facet_wrap(~ site, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    y = "Magnitude of movement",
    x = "",
  )
centroid.traj.plot

axis_contrib <- centroids.traj %>%
  filter(!is.na(dPC1)) %>%
  group_by(project, site) %>%
  summarise(
    var_PC1 = var(dPC1, na.rm = TRUE),
    var_PC2 = var(dPC2, na.rm = TRUE),
    ratio = var_PC1 / (var_PC1 + var_PC2),
    .groups = "drop"
  )
axis_contrib


axis.cont.plot <- ggplot(axis_contrib, aes(x = fct_reorder(site, ratio), y = ratio, fill = project)) +
  geom_point(shape = 21, size = 5) +
  geom_segment(aes(x = fct_reorder(site, ratio), y = ratio, yend = 0, color = project)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  geom_hline(yintercept = 0.5, lty = 2) +
  theme_bw(base_size = 14)
axis.cont.plot


step_lengths <- ggplot(centroids.traj, aes(site_fine, step_length, fill = project)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4, outliers = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    x = "Site",
    y = "Step length"
  ) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_viridis_d()
step_lengths

step_lengths_time <- centroids.traj %>%
  filter(!is.na(step_length)) %>%
  ggplot(aes(year, step_length, colour = project)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  facet_wrap(.~site, scales = "free_y") +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Centroid step length",
  ) 
step_lengths_time


### Abundance and biomass
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
  mutate(scientific_name = str_remove(scientific_name, " spp.")) %>%
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
  group_by(ID, project, site, subsite_level1, subsite_level2, subsite_level3, year, scientific_name) %>%
  summarize(density = sum(density_coll, na.rm = TRUE),
            biomass = sum(biomass_coll)) %>%
  filter(!is.na(density), !is.na(biomass)) %>%
  ungroup()


traits.df.parsed <- taxa.df.t %>%
  ungroup() %>%
  select(scientific_name, PC1, PC2)

traitsxcomms <- comm.filtered.abu %>%
  ungroup() %>%
  inner_join(traits.df.parsed, by = "scientific_name") %>%
  distinct()

# safety: restrict to rows with traits & biomass/density
traits_metrics <- traitsxcomms %>%
  filter(!is.na(PC1), !is.na(PC2))

# define global grid breaks
nbins <- 40  # adjust: 20–40 is usually fine

pc1_breaks <- seq(min(traits_metrics$PC1), max(traits_metrics$PC1), length.out = nbins + 1)
pc2_breaks <- seq(min(traits_metrics$PC2), max(traits_metrics$PC2), length.out = nbins + 1)

grid.comms <- traits_metrics %>%
  # reshape metrics
  pivot_longer(
    cols = c(density, biomass),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value), value > 0) %>%
  mutate(
    pc1_bin = cut(PC1, breaks = pc1_breaks, include.lowest = TRUE),
    pc2_bin = cut(PC2, breaks = pc2_breaks, include.lowest = TRUE)
  ) %>%
  group_by(project, site, year, metric, pc1_bin, pc2_bin) %>%
  summarise(
    cell_value = sum(value, na.rm = TRUE),
    # use bin centers as coordinates
    PC1_center = mean(PC1, na.rm = TRUE),
    PC2_center = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )


hotspot.cells <- grid.comms %>%
  group_by(project, site, year, metric) %>%
  arrange(desc(cell_value), .by_group = TRUE) %>%
  mutate(
    total_value = sum(cell_value, na.rm = TRUE),
    cum_prop = cumsum(cell_value) / total_value,
    hotspot = cum_prop <= 0.8   # top 50% of biomass/density
  ) %>%
  ungroup() %>%
  filter(hotspot)


hotspot.centroids <- hotspot.cells %>%
  group_by(project, site, year, metric) %>%
  summarise(
    hotspot_PC1 = sum(PC1_center * cell_value) / sum(cell_value),
    hotspot_PC2 = sum(PC2_center * cell_value) / sum(cell_value),
    n_cells     = n(),
    .groups = "drop"
  )


hotspot.centroids.traj <- hotspot.centroids %>%
  arrange(project, site, metric, year) %>%
  group_by(project, site, metric) %>%
  mutate(
    lag_PC1 = lag(hotspot_PC1),
    lag_PC2 = lag(hotspot_PC2),
    dPC1    = hotspot_PC1 - lag_PC1,
    dPC2    = hotspot_PC2 - lag_PC2,
    step_length = sqrt(dPC1^2 + dPC2^2),
    angle       = atan2(dPC2, dPC1)
  ) %>%
  ungroup()

hotspot.centroid.plot <- ggplot(hotspot.centroids,
                                  aes(hotspot_PC1, hotspot_PC2,
                                      group = interaction(site, metric),
                                      colour = metric)
) +
  geom_path(arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            alpha = 0.9) +
  geom_point(aes(size = year), alpha = 0.8) +
  scale_color_viridis_d() +
  facet_wrap(~ site, ncol = 9) +
  coord_equal() +
  theme_bw(base_size = 13) +
  labs(
    x = "PC1 (traits)",
    y = "PC2 (traits)",
    colour = "Metric",
    size   = "Year",
  )
hotspot.centroid.plot

hotspot.traj.plot <- ggplot(
  hotspot.centroids.traj,
  aes(year, step_length, colour = metric)) +
  geom_line() +
  facet_wrap(~ site) +
  theme_minimal(base_size = 13) +
  labs(
    y = "Step length of hotspot trajectory",
  )

hotspot.traj.plot

hotspot.direction.tests <- hotspot.centroids.traj %>%
  filter(!is.na(angle)) %>%
  group_by(project, site, metric) %>%
  summarise(
    n_steps   = n(),
    rayleigh_p = rayleigh.test(circular(angle))$p.value,
    mean_dir   = as.numeric(mean(circular(angle))),
    .groups = "drop"
  ) %>%
  arrange(rayleigh_p)

hotspot.direction.tests


hotspot.angle.plot <- ggplot(hotspot.centroids.traj, aes(x = year, y = angle, group = interaction(site, metric), color = metric)) +
  geom_line() +
  scale_y_continuous(
    breaks = c(-pi, -pi/2, 0, pi/2, pi),
    labels = c("−π", "−π/2", "0", "π/2", "π")
  ) +
  facet_wrap(.~site, ncol = 9) +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Direction of movement",
    title = "Year-to-year direction of functional centroid shifts"
  )

hotspot.angle.plot

hotspot.rosette.plot.density <- ggplot(hotspot.centroids.traj %>% filter(metric == "density"), aes(angle, fill = project)) +
  geom_histogram(binwidth = pi/12) +
  coord_polar() +
  facet_wrap(~ site) +
  theme_bw() +
  scale_fill_viridis_d() +
  labs(title = "density")

hotspot.rosette.plot.density


hotspot.rosette.plot.biomass <- ggplot(hotspot.centroids.traj %>% filter(metric == "biomass"), aes(angle, fill = project)) +
  geom_histogram(binwidth = pi/12) +
  coord_polar() +
  facet_wrap(~ site) +
  theme_bw() +
  scale_fill_viridis_d()+
  labs(title = "biomass")

hotspot.rosette.plot.biomass

rosette.plots <- hotspot.rosette.plot.density + hotspot.rosette.plot.biomass


# hotspot.centroid.traj.dir <- hotspot.centroids.traj %>%
#   filter(!is.na(dPC1)) %>%
#   pivot_longer(cols = c(dPC1, dPC2),
#                names_to = "axis",
#                values_to = "delta")
# 
# hotspot.centroid.traj.plot <- ggplot(hotspot.centroid.traj.dir, aes(axis, abs(delta), fill = project)) +
#   geom_boxplot() +
#   facet_wrap(~ site, scales = "free_y") +
#   theme_bw(base_size = 14) +
#   labs(
#     y = "Magnitude of movement",
#     x = "",
#   )
# hotspot.centroid.traj.plot

hotspot.axis.cont <- hotspot.centroids.traj %>%
  filter(!is.na(dPC1)) %>%
  group_by(project, site, metric) %>%
  summarise(
    var_PC1 = var(dPC1, na.rm = TRUE),
    var_PC2 = var(dPC2, na.rm = TRUE),
    ratio = var_PC1 / (var_PC1 + var_PC2),
    .groups = "drop"
  )
hotspot.axis.cont


hotspot.axis.cont.plot.biomass <- ggplot(hotspot.axis.cont %>% filter(metric == "biomass"), aes(x = fct_reorder(site, ratio), y = ratio, fill = project)) +
  geom_point(shape = 21, size = 5) +
  geom_segment(aes(x = fct_reorder(site, ratio), y = ratio, yend = 0, color = project)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  geom_hline(yintercept = 0.5, lty = 2) +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(0,1))
hotspot.axis.cont.plot.biomass

hotspot.axis.cont.plot.density <- ggplot(hotspot.axis.cont %>% filter(metric == "density"), aes(x = fct_reorder(site, ratio), y = ratio, fill = project)) +
  geom_point(shape = 21, size = 5) +
  geom_segment(aes(x = fct_reorder(site, ratio), y = ratio, yend = 0, color = project)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  geom_hline(yintercept = 0.5, lty = 2) +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(0,1))
hotspot.axis.cont.plot.density


hotspot.axes.plots <- hotspot.axis.cont.plot.biomass + hotspot.axis.cont.plot.density
hotspot.axes.plots


### OLD CODE: KDE (fingerprint) plots for biomass/density hotspots
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