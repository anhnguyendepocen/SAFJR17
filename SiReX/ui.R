## ui.R script ##
source("common.R")

header <- dashboardHeader(title = "SiReX")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Simulation 1", tabName = "sim1"),
    menuItem("Simulation 2", tabName = "sim2"))
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "sim1",
            h2("Simulation 1: analytical formulae vs Gaussian quadrature in parametric survival models with frailties"),
            p("The aim of this simulation consists in comparing estimated values from parametric survival models with a shared Gamma frailty term using either analytical formulae or numerical integration under a set of possible scenarios."),
            p("We simulate 1,000 data sets under a Weibull baseline hazard with shape parameter p = 1 and scale parameter lambda = 0.5. We vary:"),
            tags$ul(
              tags$li("Number of clusters: 15, 30, 100, 200;"),
              tags$li("Number of individuals per cluster: 25, 50, 100, 250, 500, 1000;"),
              tags$li("Treatment effect: -0.50, 0.00, 0.50;"),
              tags$li("Variance of the frailty term: 0.25, 0.50, 1.00.")
            ),
            fluidRow(
              box(
                title = "Summary statistics",
                status = "primary",
                solidHeader = TRUE,
                width = 7,
                height = 500,
                div(style = "width: 80%; margin: 0 auto;",
                    tableOutput("s1_table"))
                ),
              box(
                title = "Controls",
                status = "warning",
                solidHeader = TRUE,
                width = 5,
                height = 500,
                selectInput("s1_n_clusters", "Number of clusters:", sort(unique(s1$n_clusters))),
                selectInput("s1_n_individuals", "Number of individuals per cluster:", sort(unique(s1$n_individuals))),
                selectInput("s1_treatment_effect", "Treatment effect:", sort(unique(s1$treatment_effect))),
                selectInput("s1_frailty_theta", "Frailty variance:", sort(unique(s1$frailty_theta))),
                p("N.B.: only the possible combinations of clusters and individuals per clusters are displayed."),
                hr(),
                selectInput("s1_par", "Parameter:", choices = c(Lambda = "lambda", P = "p", Theta = "theta", `Tr. effect` = "trt"))
                )
              ),
            fluidRow(
              box(
                title = "Plots",
                status = "success",
                solidHeader = TRUE,
                width = 12,
                plotOutput("s1_plot")
              )
            )
    ),
    tabItem(tabName = "sim2",
            h2("Simulation 2: Gaussian quadrature in parametric survival models with random effects"),
            p("The aim of this simulation consists in assessing the accuracy of estimating values from parametric survival models with a random effect using numerical integration under a set of possible scenarios. Additionally, we compare the accuracy of Gauss-Hermite quadrature (the numerical integration technique of choice) with increasing number of nodes."),
            p("We simulate 1,000 data sets under a Weibull baseline hazard with shape parameter p = 1.5 and scale parameter lambda = 3.0. We include a random cluster-level treatment effect, and we vary:"),
            tags$ul(
              tags$li("Number of clusters: 15, 30, 100, 200;"),
              tags$li("Number of individuals per cluster: 25, 50, 100, 250, 500, 1000;"),
              tags$li("Treatment effect: -0.50, 0.00, 0.50;"),
              tags$li("Variance of the random treatment effect: 0.25, 0.50, 1.00.")
            ),
            fluidRow(
              box(
                title = "Summary statistics",
                status = "primary",
                solidHeader = TRUE,
                width = 7,
                height = 500,
                div(style = "width: 80%; margin: 0 auto;",
                    tableOutput("s2_table"))
              ),
              box(
                title = "Controls",
                status = "warning",
                solidHeader = TRUE,
                width = 5,
                height = 500,
                selectInput("s2_n_clusters", "Number of clusters:", sort(unique(s2$n_clusters))),
                selectInput("s2_n_individuals", "Number of individuals per cluster:", sort(unique(s2$n_individuals))),
                selectInput("s2_treatment_effect", "Treatment effect:", sort(unique(s2$treatment_effect))),
                selectInput("s2_frailty_sigma", "Random effect variance:", sort(unique(s2$frailty_sigma))),
                p("N.B.: only the possible combinations of clusters and individuals per clusters are displayed."),
                hr(),
                selectInput("s2_par", "Parameter:", choices = c(Lambda = "lambda", P = "p", Sigma = "sigma", `Tr. effect` = "trt"))
              )
            ),
            fluidRow(
              box(
                title = "Plots",
                status = "success",
                solidHeader = TRUE,
                width = 12,
                plotOutput("s2_plot")
              )
            )
    )
  )
)

dashboardPage(header, sidebar, body)
