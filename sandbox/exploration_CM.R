################################################################################
##
## Script to get abundance data to look at how abundance is distributed in the
## funactional space
## Taxa: Fish, Mammals, Amphibians
## Traits: Mass, Life Span, Trophic Level (num), Reproduction Rate
##
## Camille Magneville
##
## 02/2026
##
################################################################################

# 0 - Set up ===================================================================

# Define the pipe symbol so I can use it:
`%>%` <- magrittr::`%>%`

# Load libraries
librarian::shelf(tidyverse, dplyr, funbiogeo, ggplot2)
uniaue()
# Get set up
source("00_setup.R")

# Clear environment & collect garbage
rm(list = ls()); gc()




