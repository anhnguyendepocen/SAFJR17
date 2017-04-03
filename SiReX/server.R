## server.R script ##
source("common.R")

function(input, output, session) {

  # Serve only possible values:
  observe({
    # Adjust n. of individuals, simulation 1
    s1_sel_n_clusters = input$s1_n_clusters
    updateSelectInput(session,
                      "s1_n_individuals",
                      choices = sort(unique(s1$n_individuals[s1$n_clusters == input$s1_n_clusters])))

    # Adjust n. of individuals, simulation 2
    s2_sel_n_clusters = input$s2_n_clusters
    updateSelectInput(session,
                      "s2_n_individuals",
                      choices = sort(unique(s2$n_individuals[s2$n_clusters == input$s2_n_clusters])))
  })

  # Table of results, simulation 1:
  output$s1_table <- renderTable(
    s1 %>%
      filter(n_clusters == input$s1_n_clusters,
             n_individuals == input$s1_n_individuals,
             frailty_theta == input$s1_frailty_theta,
             treatment_effect == input$s1_treatment_effect,
             par %in% c("convp", "convn", input$s1_par)) %>%
      select(stat, AF, IN, GQ15, GQ35, GQ75, GQ105) %>%
      mutate(stat = factor(stat, levels = c("convn", "convp", "covp", "bias", "pbias", "mean", "se_mean", "median", "se_median", "empse", "mse"), labels = c("N. of simulations converging", "P. of simulations converging", "Coverage probability", "Bias", "Percentage bias", "Mean estimate", "Mean SE", "Median estimate","Median SE", "Empirical SE", "MSE"))) %>%
      arrange(stat) %>%
      rename(Statistic = stat),
    digits = 4)

  # Table of results, simulation 2:
  output$s2_table <- renderTable(
    s2 %>%
      filter(n_clusters == input$s2_n_clusters,
             n_individuals == input$s2_n_individuals,
             frailty_sigma == input$s2_frailty_sigma,
             treatment_effect == input$s2_treatment_effect,
             par %in% c("convp", "convn", input$s2_par)) %>%
      select(stat, GQ15, GQ35, GQ75, GQ105) %>%
      mutate(stat = factor(stat, levels = c("convn", "convp", "covp", "bias", "pbias", "mean", "se_mean", "median", "se_median", "empse", "mse"), labels = c("N. of simulations converging", "P. of simulations converging", "Coverage probability", "Bias", "Percentage bias", "Mean estimate", "Mean SE", "Median estimate","Median SE", "Empirical SE", "MSE"))) %>%
      arrange(stat) %>%
      rename(Statistic = stat),
    digits = 4)

  # Plot, simulation 1:
  output$s1_plot = renderPlot(
    s1 %>%
      filter(n_clusters == input$s1_n_clusters,
             n_individuals == input$s1_n_individuals,
             frailty_theta == input$s1_frailty_theta,
             treatment_effect == input$s1_treatment_effect,
             par == input$s1_par) %>%
      select(stat, AF, IN, GQ15, GQ35, GQ75, GQ105) %>%
      filter(stat %in% c("bias", "covp", "mse")) %>%
      mutate(stat = factor(stat, levels = c("bias", "covp", "mse"), labels = c("Bias", "Coverage probability", "MSE"))) %>%
      gather(key = key, value = value, 2:7) %>%
      mutate(key = factor(key, levels = c("AF", "IN", "GQ15", "GQ35", "GQ75", "GQ105"))) %>%
      ggplot(aes(x = key, y = value)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ stat, scales = "free_y") +
      theme_bw() +
      labs(x = "", y = ""))

  # Plot, simulation 2:
  output$s2_plot = renderPlot(
    s2 %>%
      filter(n_clusters == input$s2_n_clusters,
             n_individuals == input$s2_n_individuals,
             frailty_sigma == input$s2_frailty_sigma,
             treatment_effect == input$s2_treatment_effect,
             par == input$s2_par) %>%
      select(stat, GQ15, GQ35, GQ75, GQ105) %>%
      filter(stat %in% c("bias", "covp", "mse")) %>%
      mutate(stat = factor(stat, levels = c("bias", "covp", "mse"), labels = c("Bias", "Coverage probability", "MSE"))) %>%
      gather(key = key, value = value, 2:5) %>%
      mutate(key = factor(key, levels = c("GQ15", "GQ35", "GQ75", "GQ105"))) %>%
      ggplot(aes(x = key, y = value)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ stat, scales = "free_y") +
      theme_bw() +
      labs(x = "", y = ""))
  }
