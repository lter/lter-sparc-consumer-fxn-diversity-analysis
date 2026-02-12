###

### IMPORTANT: requires running the script LDLA_dataexploration.R
### LOAD RDS FILES FOR TRAITS/COORDINATES AND COMMUNITIES

traits.cleaned <- as.data.frame(readRDS("transformed_data/sp_faxes_coord.rds")) %>%
  mutate(scientific_name = rownames(.))

community.data <- as.data.frame(readRDS("transformed_data/Macks_data.rds"))
table(community.data$project)
table(community.data$habitat)
table(community.data$system)

comm.meta <- community.data %>%
  dplyr::select(project, habitat, system)

comm.filtered <- community.data %>%
  filter(density > 0)

traitsxcomms <- comm.filtered %>%
  ungroup() %>%
  inner_join(traits.cleaned, by = "scientific_name") %>%
  distinct()

traitsanti <- comm.filtered %>%
  anti_join(traits.cleaned, by = "scientific_name") %>%
  dplyr::select(scientific_name) %>%
  distinct()

summary(traitsxcomms)

### --- calculate centroids
centroids <- traitsxcomms %>%
  group_by(project, habitat, system, year) %>%              # ID uniquely identifies site-year
  summarise(
    centroid_PC1 = mean(PC1, na.rm = TRUE),
    centroid_PC2 = mean(PC2, na.rm = TRUE),
    n_species = n(),                        # nice to keep track
    .groups = "drop"
  ) %>%
  mutate(t01_proj = (year - min(year, na.rm = TRUE)) / (max(year, na.rm = TRUE) - min(year, na.rm = TRUE)))



### calculate trajectores and step lengths using Pythagorean theorem
centroids.traj <- centroids %>%
  arrange(habitat, project, system, year) %>%
  group_by(habitat, project, system) %>%
  mutate(
    lag_PC1 = lag(centroid_PC1),
    lag_PC2 = lag(centroid_PC2),
    dPC1 = centroid_PC1 - lag_PC1,
    dPC2 = centroid_PC2 - lag_PC2,
    step_length = sqrt(dPC1^2 + dPC2^2)      # movement through space
  ) %>%
  mutate(angle = atan2(dPC2, dPC1)) %>%
  ungroup() %>%
  mutate(t01_proj = (year - min(year, na.rm = TRUE)) / (max(year, na.rm = TRUE) - min(year, na.rm = TRUE)))



### plot steps
centroid.plot <- ggplot(centroids, aes(centroid_PC1, centroid_PC2, group = system, color = project)) +
  geom_path(alpha = 0.95) +
  geom_point(aes(fill = t01_proj), size = 2, shape = 21, color = "black") +
  coord_equal() +
  scale_color_viridis_d() +
  theme_bw(base_size = 14) +
  facet_wrap(.~project, ncol = 9) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "System",
    fill = "Year",
    title = "Community centroid shifts through functional space"
  )
centroid.plot


angle.plot <- ggplot(centroids.traj, aes(year, angle, group = system, color = project)) +
  geom_line() +
  scale_y_continuous(
    breaks = c(-pi, -pi/2, 0, pi/2, pi),
    labels = c("−π", "−π/2", "0", "π/2", "π")
  ) +
  facet_wrap(.~project) +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Direction of movement",
    title = "Year-to-year direction of functional centroid shifts"
  )

angle.plot

centroids.traj.circ <- centroids.traj %>%
  mutate(angle.circ = circular(angle))


projects <- sort(unique(centroids.traj.circ$project))


pairs <- combn(projects, 2, simplify = FALSE)

pairwise_results <- purrr::map_dfr(pairs, function(p) {
  a1 <- centroids.traj.circ %>% filter(project == p[1]) %>% pull(angle.circ)
  a2 <- centroids.traj.circ %>% filter(project == p[2]) %>% pull(angle.circ)
  
  tst <- watson.two.test(a1, a2)
  
  tibble(
    project1 = p[1],
    project2 = p[2],
    statistic = unname(tst$statistic),
    p_value = unname(tst$alpha)
  )
})

pairwise_results %>% arrange(p_value)
pairwise_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)


centroids.traj.dir <- centroids.traj %>%
  filter(!is.na(dPC1), !is.na(dPC2)) %>%
  mutate(axis_ratio = abs(dPC1) / (abs(dPC1) + abs(dPC2)))

hist(centroids.traj.dir$axis_ratio)

m1 <- brm(axis_ratio ~ project + (1|system), data = centroids.traj.dir, family = "beta", cores = 4, chains = 4, backend = "cmdstanr")
summary(m1)

## --- projects

nd <- expand_grid(centroids.traj.dir, project = sort(unique(centroids.traj.dir$project))) %>%
  filter(project != "VCR")

ep <- fitted(
  m1,
  newdata = nd,
  re_formula = NA,     # population-level only
  summary = TRUE
) %>%
  as_tibble() %>%
  bind_cols(nd) %>%
  rename(mean = Estimate, lwr = Q2.5, upr = Q97.5)

ep

project.direction.plot <- ggplot(ep, aes(x = reorder(project, mean), y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Expected axis dominance (|dPC1| / (|dPC1| + |dPC2|))",
    title = "Projects differ in which trait axis dominates centroid movement"
  )
project.direction.plot

# ### --- time 
# nd_time <- expand_grid(
#   project = unique(centroids.traj.dir$project),
#   t01_proj = seq(0, 1, length.out = 100)) %>%
#   filter(project != "VCR")
# 
# pred_time <- fitted(
#   m1,
#   newdata = nd_time,
#   re_formula = NA,
#   summary = TRUE) %>%
#   as_tibble() %>%
#   bind_cols(nd_time)
# 
# ggplot(pred_time,
#        aes(x = t01_proj, y = Estimate, colour = project)) +
#   geom_line() +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = project),
#               alpha = 0.15, colour = NA) +
#   theme_minimal(base_size = 14) +
#   labs(x = "Relative time (0–1)",
#        y = "Predicted axis_ratio",
#        title = "Axis dominance through time by project")
# 

centroids.dir.abs <- centroids.traj.dir %>%
  mutate(adPC1 = abs(dPC1), adPC2 = abs(dPC2)) 

angle.data <- centroids.traj.circ %>%
  filter(!is.na(angle)) %>%
  mutate(angle_wrapped = angle)  # already in [-pi, pi] from atan2

rosette.plots.project <- ggplot(angle.data, aes(angle_wrapped, fill = project)) +
  geom_histogram(binwidth = pi/12, boundary = 0) +
  coord_polar(start = pi/2, direction = 1) +   # rotate 90° counterclockwise
  facet_wrap(~ project, scales = "free") +
  scale_fill_viridis_d() +
  theme_bw(base_size = 13) +
  labs(
    x = NULL, y = "Count",
    title = "Centroid movement directions"
  )
rosette.plots.project

centroid.tests <- centroids.traj %>%
  filter(!is.na(angle)) %>%
  group_by(system) %>%
  summarise(
    p = circular::rayleigh.test(circular(angle))$p.value,
    mean_dir = as.numeric(mean(circular(angle))),
    .groups = "drop"
  ) %>%
  arrange(p) %>%
  left_join(comm.meta)
centroid.tests


centroid.traj.dir <- centroids.traj %>%
  filter(!is.na(dPC1)) %>%
  pivot_longer(cols = c(dPC1, dPC2),
               names_to = "axis",
               values_to = "delta")

centroid.traj.plot <- ggplot(centroid.traj.dir, aes(axis, abs(delta), fill = project)) +
  geom_boxplot() +
  facet_wrap(~ project, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    y = "Magnitude of movement",
    x = "",
  )
centroid.traj.plot


step_lengths_time <- centroids.traj %>%
  filter(!is.na(step_length)) %>%
  ggplot(aes(year, step_length, colour = project)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  facet_wrap(.~project, scales = "free_y") +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Centroid step length",
  ) 
step_lengths_time


### --- Abundance and biomass


# define global grid breaks
nbins <- 40  # adjust: 20–40 is usually fine

pc1_breaks <- seq(min(traitsxcomms$PC1), max(traitsxcomms$PC1), length.out = nbins + 1)
pc2_breaks <- seq(min(traitsxcomms$PC2), max(traitsxcomms$PC2), length.out = nbins + 1)

grid.comms <- traitsxcomms %>%
  rename(biomass = "total_bm_area") %>%
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
  group_by(project, system, year, metric, pc1_bin, pc2_bin) %>%
  summarise(
    cell_value = sum(value, na.rm = TRUE),
    # use bin centers as coordinates
    PC1_center = mean(PC1, na.rm = TRUE),
    PC2_center = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) 

### --- one more test

p <- 2

hot_traj_all <- grid.comms %>%
  filter(metric == "biomass",
         !is.na(cell_value),
         cell_value > 0,
         !is.na(PC1_center),
         !is.na(PC2_center)) %>%
  group_by(project, year, PC1_center, PC2_center) %>%
  summarise(v = sum(cell_value), .groups = "drop") %>%
  group_by(project, year) %>%
  mutate(w = v^p) %>%
  summarise(
    hot_PC1 = weighted.mean(PC1_center, w),
    hot_PC2 = weighted.mean(PC2_center, w),
    .groups = "drop"
  )

hotspots.plot.biomass <- ggplot(hot_traj_all, aes(hot_PC1, hot_PC2)) +
  geom_path(linewidth = 0.2) +
  geom_label(aes(fill = year, label = year)) +
  facet_wrap(~ project) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) +
  labs(
    x = "PC1",
    y = "PC2",
    colour = "Year",
    title = paste0("Biomass hotspot trajectories across projects (p = ", p, ")")
  )
hotspots.plot.biomass


hot_traj_all <- grid.comms %>%
  filter(metric == "density",
         !is.na(cell_value),
         cell_value > 0,
         !is.na(PC1_center),
         !is.na(PC2_center)) %>%
  group_by(project, year, PC1_center, PC2_center) %>%
  summarise(v = sum(cell_value), .groups = "drop") %>%
  group_by(project, year) %>%
  mutate(w = v^p) %>%
  summarise(
    hot_PC1 = weighted.mean(PC1_center, w),
    hot_PC2 = weighted.mean(PC2_center, w),
    .groups = "drop"
  )

hotspots.plot.density <- ggplot(hot_traj_all, aes(hot_PC1, hot_PC2)) +
  geom_path(linewidth = 0.2) +
  geom_label(aes(fill = year, label = year)) +
  facet_wrap(~ project) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank()) +
  labs(
    x = "PC1",
    y = "PC2",
    colour = "Year",
    title = paste0("Density hotspot trajectories across projects (p = ", p, ")")
  )
hotspots.plot.density

### --- 

hotspot.cells <- grid.comms %>%
  group_by(project, system, year, metric) %>%
  arrange(desc(cell_value), .by_group = TRUE) %>%
  mutate(
    total_value = sum(cell_value, na.rm = TRUE),
    cum_prop = cumsum(cell_value) / total_value,
    hotspot = cum_prop <= 0.4   # top 50% of biomass/density
  ) %>%
  ungroup() %>%
  filter(hotspot)


hotspot.centroids <- hotspot.cells %>%
  group_by(project, system, year, metric) %>%
  summarise(
    hotspot_PC1 = sum(PC1_center * cell_value) / sum(cell_value),
    hotspot_PC2 = sum(PC2_center * cell_value) / sum(cell_value),
    n_cells     = n(),
    .groups = "drop"
  ) %>%
  mutate(t01_proj = (year - min(year, na.rm = TRUE)) / (max(year, na.rm = TRUE) - min(year, na.rm = TRUE)))


hotspot.centroids.traj <- hotspot.centroids %>%
  arrange(project, system, metric, year) %>%
  group_by(project, system, metric) %>%
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
                                      group = interaction(system, metric),
                                      colour = metric)
) +
  geom_path(arrow = arrow(type = "closed", length = unit(0.08, "inches")),
            alpha = 0.9) +
  geom_point(aes(size = year), alpha = 0.8) +
  scale_color_viridis_d() +
  facet_wrap(~ project, ncol = 3) +
  coord_equal() +
  theme_bw(base_size = 13) +
  labs(
    x = "PC1 (traits)",
    y = "PC2 (traits)",
    colour = "Metric",
    size   = "Year",
  )
hotspot.centroid.plot


hotspot.direction.tests <- hotspot.centroids.traj %>%
  filter(!is.na(angle)) %>%
  group_by(project, system, metric) %>%
  summarise(
    n_steps   = n(),
    rayleigh_p = rayleigh.test(circular(angle))$p.value,
    mean_dir   = as.numeric(mean(circular(angle))),
    .groups = "drop"
  ) %>%
  arrange(rayleigh_p)

hotspot.direction.tests



hotspots.centroids.traj.circ <- hotspot.centroids.traj %>%
  mutate(angle.circ = circular(angle))


projects <- sort(unique(hotspot.centroids.traj$project))


pairs <- combn(projects, 2, simplify = FALSE)

hotspot_pairwise_results_biomass <- purrr::map_dfr(pairs, function(p) {
  a1 <- hotspots.centroids.traj.circ %>% filter(project == p[1] & metric == "biomass") %>% pull(angle.circ)
  a2 <- hotspots.centroids.traj.circ %>% filter(project == p[2] & metric == "biomass") %>% pull(angle.circ)
  
  tst <- watson.two.test(a1, a2)
  
  tibble(
    project1 = p[1],
    project2 = p[2],
    statistic = unname(tst$statistic),
    p_value = unname(tst$alpha)
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
hotspot_pairwise_results_biomass

hotspot_pairwise_results_density<- purrr::map_dfr(pairs, function(p) {
  a1 <- hotspots.centroids.traj.circ %>% filter(project == p[1] & metric == "density") %>% pull(angle.circ)
  a2 <- hotspots.centroids.traj.circ %>% filter(project == p[2] & metric == "density") %>% pull(angle.circ)
  
  tst <- watson.two.test(a1, a2)
  
  tibble(
    project1 = p[1],
    project2 = p[2],
    statistic = unname(tst$statistic),
    p_value = unname(tst$alpha)
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
hotspot_pairwise_results_density


hotspot.centroids.traj.dir <- hotspot.centroids.traj %>%
  filter(!is.na(dPC1), !is.na(dPC2)) %>%
  mutate(axis_ratio = abs(dPC1) / (abs(dPC1) + abs(dPC2)))

hot.dirs.biomass <- hotspot.centroids.traj.dir %>%
  filter(metric == "biomass")%>%
  mutate(n = n()) %>%
  mutate(axis_ratio = (axis_ratio * (n - 1) + 0.5) / n)

m2biomass <- brm(axis_ratio ~ project + (1|system), data = hot.dirs.biomass, family = "beta", cores = 4, chains = 4, backend = "cmdstanr")
summary(m2biomass)


hot.dirs.density <- hotspot.centroids.traj.dir %>%
  filter(metric == "density") %>%
  mutate(n = n()) %>%
  mutate(axis_ratio = (axis_ratio * (n - 1) + 0.5) / n)

m2density <- brm(axis_ratio ~ project + (1|system), data = hot.dirs.density, family = "beta", cores = 4, chains = 4, backend = "cmdstanr")
summary(m2density)

### --- time 
nd_time <- expand_grid(
  project = unique(hot.dirs.biomass$project)) %>%
  filter(project != "VCR")

ep <- fitted(
  m2biomass,
  newdata = nd,
  re_formula = NA,     # population-level only
  summary = TRUE
) %>%
  as_tibble() %>%
  bind_cols(nd) %>%
  rename(mean = Estimate, lwr = Q2.5, upr = Q97.5)

ep

project.direction.plot <- ggplot(ep, aes(x = reorder(project, mean), y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Expected axis dominance (|dPC1| / (|dPC1| + |dPC2|))",
    title = "Projects differ in which trait axis dominates centroid movement"
  )
project.direction.plot


nd_time <- expand_grid(
  project = unique(hot.dirs.density$project)) %>%
  filter(project != "VCR")

ep <- fitted(
  m2density,
  newdata = nd,
  re_formula = NA,     # population-level only
  summary = TRUE
) %>%
  as_tibble() %>%
  bind_cols(nd) %>%
  rename(mean = Estimate, lwr = Q2.5, upr = Q97.5)

ep

project.direction.plot <- ggplot(ep, aes(x = reorder(project, mean), y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Expected axis dominance (|dPC1| / (|dPC1| + |dPC2|))",
    title = "Projects differ in which trait axis dominates centroid movement"
  )
project.direction.plot



centroids.dir.abs <- hot.dirs.biomass %>%
  mutate(adPC1 = abs(dPC1), adPC2 = abs(dPC2)) 

angle.data <- hotspots.centroids.traj.circ %>%
  filter(!is.na(angle)) %>%
  mutate(angle_wrapped = angle)  # already in [-pi, pi] from atan2

rosette.plots.project <- ggplot(angle.data %>% filter(metric == "biomass"), aes(angle_wrapped, fill = project)) +
  geom_histogram(binwidth = pi/8, boundary = 0) +
  coord_polar(start = pi/2, direction = 1) +   # rotate 90° counterclockwise
  facet_wrap(~ project, scales = "free_y") +
  scale_fill_viridis_d() +
  theme_bw(base_size = 13) +
  labs(
    x = NULL, y = "Count",
    title = "Centroid movement directions: biomass"
  )
rosette.plots.project

rosette.plots.project <- ggplot(angle.data %>% filter(metric == "density"), aes(angle_wrapped, fill = project)) +
  geom_histogram(binwidth = pi/8, boundary = 0) +
  coord_polar(start = pi/2, direction = 1) +   # rotate 90° counterclockwise
  facet_wrap(~ project, scales = "free") +
  scale_fill_viridis_d() +
  theme_bw(base_size = 13) +
  labs(
    x = NULL, y = "Count",
    title = "Centroid movement directions: density"
  )
rosette.plots.project


centroid.tests <- centroids.traj %>%
  filter(!is.na(angle)) %>%
  group_by(system) %>%
  summarise(
    p = circular::rayleigh.test(circular(angle))$p.value,
    mean_dir = as.numeric(mean(circular(angle))),
    .groups = "drop"
  ) %>%
  arrange(p) %>%
  left_join(comm.meta)
centroid.tests


centroid.traj.dir <- centroids.traj %>%
  filter(!is.na(dPC1)) %>%
  pivot_longer(cols = c(dPC1, dPC2),
               names_to = "axis",
               values_to = "delta")

centroid.traj.plot <- ggplot(centroid.traj.dir, aes(axis, abs(delta), fill = project)) +
  geom_boxplot() +
  facet_wrap(~ project, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    y = "Magnitude of movement",
    x = "",
  )
centroid.traj.plot


step_lengths_time <- centroids.traj %>%
  filter(!is.na(step_length)) %>%
  ggplot(aes(year, step_length, colour = project)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  facet_wrap(.~project, scales = "free_y") +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(
    x = "Year",
    y = "Centroid step length",
  ) 
step_lengths_time






#######################
######## CODE GRAVEYARD
#######################

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