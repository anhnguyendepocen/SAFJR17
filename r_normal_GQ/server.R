library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

data <- readRDS("s_an_vs_gq_summary.RDS")

df <- data %>%
  gather(key = key, value = value, 8:75) %>%
  mutate(method = str_sub(key, 1, 2),
         par = str_sub(key, 4, str_length(key))) %>%
  separate(par, c("name", "stat"), sep = "_", extra = "merge") %>%
  mutate(stat = ifelse(is.na(stat), name, stat))

function(input, output) {
  output$table <- renderTable(
    df %>%
      filter(n_individuals == ifelse(input$sample_size %in% c("25 individuals X 100 clusters", "25 individuals X 200 clusters"), 25,
                                     ifelse(input$sample_size %in% c("50 individuals X 200 clusters", "50 individuals X 100 clusters"), 50,
                                            ifelse(input$sample_size %in% c("250 individuals X 15 clusters", "250 individuals X 30 clusters"), 250,
                                                   ifelse(input$sample_size %in% c("500 individuals X 15 clusters", "500 individuals X 30 clusters"), 500, 1000)))),
             n_clusters == ifelse(input$sample_size %in% c("250 individuals X 15 clusters", "500 individuals X 15 clusters", "1000 individuals X 15 clusters"), 15,
                                  ifelse(input$sample_size %in% c("250 individuals X 30 clusters", "500 individuals X 30 clusters", "1000 individuals X 30 clusters"), 30,
                                         ifelse(input$sample_size %in% c("25 individuals X 100 clusters", "50 individuals X 100 clusters"), 100, 200))),
             frailty_theta == input$frailty_theta,
             treatment_effect == input$treatment_effect,
             ngl == input$ngl) %>%
      ungroup() %>%
      filter(grepl(ifelse(input$par == "Convergence",
                          "^conv",
                          ifelse(input$par == "Frailty variance",
                                 "^theta",
                                 ifelse(input$par == "Treatment effect",
                                        "^trt",
                                        ifelse(input$par == "Lambda",
                                               "^lambda",
                                               "^p")))), name)) %>%
      select(value, method, stat) %>%
      spread(key = method, value = value) %>%
      mutate(stat = factor(stat,
                           levels = c("convn", "convp", "bias", "mean", "se_mean", "median", "se_median", "empse", "mse", "covp"),
                           labels = c("N. of simulations converging", "P. of simulations converging", "Bias", "Mean estimate", "Mean SE", "Median estimate","Median SE", "Empirical SE", "MSE", "Coverage probability"))) %>%
      rename(` ` = stat, `Analytical formulae` = AF, `Gaussian quadrature` = GQ),
    digits = 10
  )
}

### Deploy to shinyapps.io with:
# if (!requireNamespace("rsconnect")) install.packages("rsconnect")
# rsconnect::deployApp("r_AF_vs_GQ")
