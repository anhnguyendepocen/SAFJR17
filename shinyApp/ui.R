library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
data <- readRDS("s_an_vs_gq_summary.RDS")

# Define UI
ui <- fluidPage(

  sidebarLayout(

    # Sidebar with a slider input
    sidebarPanel(
      selectInput("n_individuals",
                  "N. of individuals:",
                  unique(data$n_individuals)),
      selectInput("n_clusters",
                  "N. of clusters:",
                  unique(data$n_clusters)),
      selectInput("frailty_theta",
                  "Frailty variance:",
                  unique(data$frailty_theta)),
      selectInput("treatment_effect",
                  "Treatment effect:",
                  unique(data$treatment_effect)),
      selectInput("ngl",
                  "Gauss-Laguerre knots:",
                  unique(data$ngl)),
      selectInput("par",
                  "Parameter:",
                  c("Convergence", "Frailty variance", "Treatment effect", "Lambda", "P"))),

    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("table")
    )
  )
)
