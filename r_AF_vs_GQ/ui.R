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
      selectInput("sample_size",
                  "Sample size:",
                  c("25 individuals X 100 clusters",
                    "25 individuals X 200 clusters",
                    "50 individuals X 100 clusters",
                    "50 individuals X 200 clusters",
                    "250 individuals X 15 clusters",
                    "250 individuals X 30 clusters",
                    "500 individuals X 15 clusters",
                    "500 individuals X 30 clusters",
                    "1000 individuals X 15 clusters",
                    "1000 individuals X 30 clusters")),
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
