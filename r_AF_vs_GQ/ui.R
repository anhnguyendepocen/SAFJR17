library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(tidyr)
data <- readRDS("s_an_vs_gq_summary.RDS")

# Define UI
ui <- fluidPage(
  # Title
  titlePanel(title = "Exploring simulations results: AF vs GQ"),

  sidebarLayout(
    # Sidebar
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

    # Show subset of table
    mainPanel(
      tableOutput("table")
    )
  ),
  theme = shinytheme("slate")
)

### Deploy to shinyapps.io with:
# if (!requireNamespace("rsconnect")) install.packages("rsconnect")
# rsconnect::deployApp("r_AF_vs_GQ")
