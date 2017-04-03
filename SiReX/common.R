## common.R script ##

library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)

s1 <- readRDS("s_an_vs_gq_vs_int_summary.RDS") %>%
  ungroup() %>%
  gather(key = key, value = value, 7:234) %>%
  separate(key, c("method", "par", "stat"), sep = "_", extra = "merge") %>%
  mutate(stat = ifelse(is.na(stat), par, stat)) %>%
  spread(key = method, value = value)

s2 <- readRDS("s_normal_gq_summary.RDS") %>%
  ungroup() %>%
  gather(key = key, value = value, 7:158) %>%
  separate(key, c("method", "par", "stat"), sep = "_", extra = "merge") %>%
  mutate(stat = ifelse(is.na(stat), par, stat)) %>%
  spread(key = method, value = value)


