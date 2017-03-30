## common.R script ##

library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)

s1 <- readRDS("s_an_vs_gq_vs_int_summary.RDS") %>%
  ungroup() %>%
  gather(key = key, value = value, 8:121) %>%
  separate(key, c("method", "par", "stat"), sep = "_", extra = "merge") %>%
  mutate(stat = ifelse(is.na(stat), par, stat),
         key = paste0(method, ngl)) %>%
  select(-method, -ngl) %>%
  spread(key = key, value = value) %>%
  mutate(AF = (AF15 + AF35 + AF75 + AF105) / 4,
         IN = (IN15 + IN35 + IN75 + IN105) / 4) %>%
  select(-AF15, -AF35, -AF75, -AF105, -IN15, -IN35, -IN75, -IN105)

 s2 <- readRDS("s_normal_gq_summary.RDS") %>%
  ungroup() %>%
  gather(key = key, value = value, 8:45) %>%
  separate(key, c("method", "par", "stat"), sep = "_", extra = "merge") %>%
  mutate(stat = ifelse(is.na(stat), par, stat),
         key = paste0(method, ngh)) %>%
   select(-method, -ngh) %>%
   spread(key = key, value = value)
