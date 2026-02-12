################################################################################
##
## Trait Scatter & Correlation Plots (multi-taxa)
## Taxa: Fish, Birds, Zooplankton
## Traits: Age, Mass, Trophic Level, Reproduction Unified
##
## C Honeyman
## 02/12/2026
##
################################################################################

# 0 - Set up ===================================================================

# Define the pipe symbol so I can use it:
`%>%` <- magrittr::`%>%`

# Load libraries
librarian::shelf(tidyverse, dplyr, funbiogeo, ggplot2, mFD, tibble)

# Get set up
source("00_setup.R")

# Clear environment & collect garbage
rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# 1 - Load data =================================================================
sp_tr_zscore <- readRDS(
  "~/lter-sparc-consumer-fxn-diversity-analysis/transformed_data/sp_tr_zscore.rds"
)

# 2 - Helpers ===================================================================

`%||%` <- function(lhs, rhs) if (is.null(lhs)) rhs else lhs

prep_trait_data <- function(df, taxa_keep) {
  df %>%
    select(
      project, habitat, scientific_name, class, sp.proj,
      taxa,
      tr.age.zt,
      tr.mass.adult.zt,
      tr.reproduction.unified.zt,
      tr.trophic.level.zt
    ) %>%
    filter(!is.na(taxa), taxa %in% taxa_keep)
}

trait_scatter_lm <- function(df, x, y, color = "project", title = NULL) {
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) +
    geom_point(size = 1) +
    geom_smooth(method = "lm") +
    labs(
      title = title %||% paste0(x, " x ", y),
      x = x,
      y = y
    ) +
    theme_grey()
}

make_taxon_plot_set <- function(df_taxon, trait_pairs, taxon_label) {
  plots <- lapply(trait_pairs, function(pr) {
    trait_scatter_lm(
      df = df_taxon,
      x = pr[[1]],
      y = pr[[2]],
      color = "project",
      title = paste0(taxon_label, " | ", pr[[1]], " x ", pr[[2]])
    )
  })
  
  # Return both list + a standard 2x3 grid
  list(
    plots = plots,
    grid =
      (plots[[1]] | plots[[3]] | plots[[2]]) /
      (plots[[5]] | plots[[4]] | plots[[6]])
  )
}

# 3 - Config ====================================================================

taxa_keep <- c("Fish", "Birds", "Zooplankton")

trait_pairs <- list(
  c("tr.age.zt", "tr.mass.adult.zt"),
  c("tr.age.zt", "tr.trophic.level.zt"),
  c("tr.age.zt", "tr.reproduction.unified.zt"),
  c("tr.mass.adult.zt", "tr.trophic.level.zt"),
  c("tr.mass.adult.zt", "tr.reproduction.unified.zt"),
  c("tr.trophic.level.zt", "tr.reproduction.unified.zt")
)

# 4 - Run =======================================================================

d <- prep_trait_data(sp_tr_zscore, taxa_keep = taxa_keep)

# Split by taxa and build plot grids
taxa_groups <- split(d, d$taxa)

plot_sets <- lapply(names(taxa_groups), function(tx) {
  df_tx <- taxa_groups[[tx]]
  make_taxon_plot_set(df_tx, trait_pairs, taxon_label = tx)
})
names(plot_sets) <- names(taxa_groups)

# 5 - View results ==============================================================
# Show one taxa grid:
plot_sets$Fish$grid

# Or print all taxa grids sequentially:
for (tx in names(plot_sets)) {
  print(plot_sets[[tx]]$grid)
  invisible(readline(paste0("Next taxa (", tx, " done). Press Enter ")))
}

# ==============================================================================
# Correlation matrices + plots (run AFTER your trait scatterplots)
# Assumes you already created:
#   d <- prep_trait_data(sp_tr_zscore, taxa_keep = c("Fish","Birds","Zooplankton"))
# where d contains: taxa + the .zt trait columns used below
# ==============================================================================

library(corrplot)

traits_zt <- c(
  "tr.age.zt",
  "tr.mass.adult.zt",
  "tr.reproduction.unified.zt",
  "tr.trophic.level.zt"
)

# ---- helpers (reduce redundancy) ----
clean_numeric_matrix <- function(df, traits) {
  x <- df[, traits, drop = FALSE]
  
  x[] <- lapply(x, function(v) {
    v <- as.numeric(v)
    v[!is.finite(v)] <- NA_real_
    v
  })
  
  sds <- sapply(x, sd, na.rm = TRUE)
  x <- x[, is.finite(sds) & sds > 0, drop = FALSE]
  
  x
}

compute_cor <- function(x, method = "pearson", use = "pairwise.complete.obs") {
  if (ncol(x) < 2) return(NULL)
  cor(as.data.frame(x), method = method, use = use)
}

plot_cor_set <- function(x, cor_mat, title_prefix) {
  pairs(x, main = paste0(title_prefix, " | pairs()"))
  
  corrplot::corrplot(
    cor_mat,
    method = "number",
    main = paste0(title_prefix, " | corrplot()"),
    mar = c(0, 0, 2, 0)
  )
}

# ---- run by taxa (reuses your existing filtered data 'd') ----
stopifnot(exists("d"))
stopifnot("taxa" %in% names(d))
stopifnot(all(traits_zt %in% names(d)))

taxa_groups <- split(d, d$taxa)

corr_results <- lapply(names(taxa_groups), function(tx) {
  g <- taxa_groups[[tx]]
  x <- clean_numeric_matrix(g, traits_zt)
  
  if (nrow(x) < 3 || ncol(x) < 2) {
    message("Skipping taxa=", tx, " (too few rows or varying trait cols).")
    return(list(taxa = tx, cor = NULL, data = x))
  }
  
  corr_mat <- compute_cor(x, method = "pearson", use = "pairwise.complete.obs")
  if (is.null(corr_mat)) {
    message("Skipping taxa=", tx, " (cor matrix NULL).")
    return(list(taxa = tx, cor = NULL, data = x))
  }
  
  cat("\n====================\n")
  cat("Taxa:", tx, "\n")
  cat("n rows:", nrow(g), " | traits used:", paste(colnames(x), collapse = ", "), "\n")
  cat("====================\n")
  print(round(corr_mat, 3))
  
  plot_cor_set(
    x = x,
    cor_mat = corr_mat,
    title_prefix = paste0("Taxa=", tx, " | .zt | n=", nrow(g))
  )
  
  invisible(readline("Next taxa correlation plots? Press Enter "))
  
  list(taxa = tx, cor = corr_mat, data = x)
})
names(corr_results) <- names(taxa_groups)

# Inspect later if you want:
# corr_results$Fish$cor
# corr_results$Birds$cor
# corr_results$Zooplankton$cor
