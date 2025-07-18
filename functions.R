#  Statistical tools  ####

# Build mortality table with 95% CI (Wilson binomial)
# mortality_table dependency is a tibble created by create_mortality_table() can be named otherwise

create_mortality_table<- function(
    max_sample_size = 1000, 
    step_size = 50,
    lower_threshold = 94,
    desired_surv = 98,
    max_mortality_proportion = 0.10) {
  
  if (lower_threshold > desired_surv) {
    stop("Error: `lower_threshold` (", lower_threshold, 
         ") cannot be greater than `desired_surv` (", desired_surv, ")")
  }
  
  sample_sizes <- seq(0, max_sample_size, by = step_size)
  
  do.call(rbind, lapply(sample_sizes, function(n) {
    n_mortality <- 0:floor(max_mortality_proportion * n)
    
    if (n == 0) {
      percent <- rep(0, length(n_mortality))
      lower <- rep(0, length(n_mortality))
      upper <- rep(0, length(n_mortality))
    } else {
      conf <- binom.confint(n_mortality, n, methods = "wilson")
      percent <- round(conf$mean * 100, 2)
      lower <- round(conf$lower * 100, 2)
      upper <- round(conf$upper * 100, 2)
    }
    
    survival <- round(100 - percent, 2)
    ci_lower_surv <- round(100 - upper, 2)
    ci_upper_surv <- round(100 - lower, 2)
    
    
    
    decision <- case_when(
      survival >= desired_surv & ci_lower_surv >= lower_threshold ~ "Accept H0",
      ci_upper_surv < desired_surv + 0.1 ~ "Reject H0",
      TRUE ~ "Continue"
    )
    
    survival_formatted <- paste0(survival, "% [CI: ", ci_lower_surv, "% - ", ci_upper_surv, "%]")
    
    data.frame(
      sample_size = n,
      n_mortality = n_mortality,
      survival = survival,
      ci_lower_surv = ci_lower_surv,
      ci_upper_surv = ci_upper_surv,
      desired_surv = desired_surv,
      lower_threshold = lower_threshold,
      decision = factor(decision, levels = c("Reject H0", "Accept H0", "Continue")),
      survival_formatted = survival_formatted
    )
  }))
}

# Lookup mortality data from table

lookup_mortality_data <- function(mortality_table, sample_size, n_deaths) {
  result <- mortality_table %>%
    filter(sample_size == !!sample_size, n_mortality == !!n_deaths) %>%
    select(sample_size, n_mortality, survival, ci_lower_surv, ci_upper_surv, 
           decision, survival_formatted)
  
  if(nrow(result) == 0) {
    cat("No match found for sample_size =", sample_size, ", n_mortality =", n_deaths, "\n")
    return(NULL)
  }
  
  return(result)
}

# Simulate sequential trial design

simulate_sequential_trial <- function(simulations = 1:50, 
                                      true_survival_rates = seq(0.9, 1, by = 0.01), 
                                      max_sample = 500, 
                                      check_intervals = c(100, 200, 300, 400, 500),
                                      mortality_table,
                                      recapture_rates = c(1.0)) {  # Accept vector of rates
  
  # Initialize results
  all_results <- tibble()
  
  # Loop through each simulation
  for(sim_id in simulations) {
    # Loop through each recapture rate  
    for(recapture_rate in recapture_rates) {
      # Loop through each survival rate
      for(survival_rate in true_survival_rates) {
        
        # Initialize for this specific simulation
        results <- tibble(
          interval = integer(),
          cumulative_n = integer(),
          cumulative_survivors = integer(),
          cumulative_deaths = integer(),
          decision = character(),
          trial_stopped = logical()
        )
        
        cumulative_survivors <- 0
        
        # Run through check intervals
        for(interval in check_intervals) {
          if(interval == check_intervals[1]) {
            interval_size <- interval
          } else {
            prev_n <- check_intervals[which(check_intervals == interval) - 1]
            interval_size <- interval - prev_n
          }
          
          # Generate true outcomes
          true_new_survivors <- rbinom(1, interval_size, survival_rate)
          true_new_deaths <- interval_size - true_new_survivors
          
          # Apply recapture rate
          recaptured_alive <- rbinom(1, true_new_survivors, recapture_rate)
          recaptured_dead <- rbinom(1, true_new_deaths, recapture_rate)
          total_recaptured <- recaptured_alive + recaptured_dead
          
          # Calculate observed survival rate from recaptured fish
          if(total_recaptured > 0) {
            observed_survival_rate <- recaptured_alive / total_recaptured
          } else {
            observed_survival_rate <- survival_rate  # Fallback
          }
          
          # Estimate total outcomes using proportional assumption
          estimated_new_deaths <- interval_size * (1 - observed_survival_rate)
          cumulative_survivors <- cumulative_survivors + (interval_size - estimated_new_deaths)
          cumulative_deaths <- interval - cumulative_survivors
          
          # Round estimated deaths for lookup
          cumulative_deaths_rounded <- round(cumulative_deaths)
          
          # Look up decision using estimated totals
          decision_row <- mortality_table %>% 
            filter(sample_size == interval, n_mortality == cumulative_deaths_rounded)
          
          current_decision <- if(nrow(decision_row) > 0) decision_row$decision[1] else "Continue"
          
          # Add to results
          results <- results %>% 
            add_row(
              interval = interval,
              cumulative_n = interval,
              cumulative_survivors = cumulative_survivors,
              cumulative_deaths = cumulative_deaths,
              decision = current_decision,
              trial_stopped = current_decision != "Continue"
            )
          
          # Stop if decision reached
          if(current_decision != "Continue") break
        }
        
        # Add simulation metadata and combine
        results <- results %>%
          mutate(sim_id = sim_id, 
                 true_survival_rate = survival_rate,
                 recapture_rate = recapture_rate)  # Add this column
        
        all_results <- bind_rows(all_results, results)
      }
    }
  }
  
  return(all_results)
}

# Determine when trials stop
create_stopping_summary <- function(simulation_data) {
  
  # Get final results for each simulation
  final_results <- simulation_data %>%
    group_by(sim_id, true_survival_rate, recapture_rate) %>%
    slice_tail(n = 1) %>%  # Get last row for each simulation
    ungroup()
  
  # Stopping time analysis for trials that stopped
  stopped_trials <- final_results %>%
    filter(trial_stopped == TRUE) %>%
    group_by(true_survival_rate, recapture_rate, decision) %>%
    summarise(
      n_stopped = n(),
      median_fish = median(cumulative_n),
      mean_fish = round(mean(cumulative_n), 1),
      q25_fish = quantile(cumulative_n, 0.25),
      q75_fish = quantile(cumulative_n, 0.75),
      min_fish = min(cumulative_n),
      max_fish = max(cumulative_n),
      .groups = "drop"
    )
  
  # Overall trial summary
  trial_summary <- final_results %>%
    group_by(true_survival_rate, recapture_rate) %>%
    summarise(
      total_trials = n(),
      n_accept = sum(decision == "Accept H0", na.rm = TRUE),
      n_reject = sum(decision == "Reject H0", na.rm = TRUE), 
      n_continue = sum(decision == "Continue", na.rm = TRUE),
      prop_accept = round(n_accept / total_trials, 3),
      prop_reject = round(n_reject / total_trials, 3),
      prop_continue = round(n_continue / total_trials, 3),
      .groups = "drop"
    )
  
  # Combine results
  list(
    stopping_details = stopped_trials,
    trial_overview = trial_summary
  )
}


# Determine acceptance and rejection rates for given sample size
# e.g., when do we stop?
calculate_cumulative_stopping <- function(simulation_data) {
  
  # For each trial, find when it stopped (if it did)
  stopping_points <- simulation_data %>%
    filter(trial_stopped == TRUE) %>%
    select(sim_id, true_survival_rate, recapture_rate, interval, decision)
  
  # Create all combinations of intervals and survival/recapture rates
  all_intervals <- simulation_data %>%
    distinct(true_survival_rate, recapture_rate) %>%
    cross_join(tibble(interval = unique(simulation_data$interval))) %>%
    arrange(true_survival_rate, recapture_rate, interval)
  
  # Calculate cumulative stopping proportions
  cumulative_stats <- all_intervals %>%
    left_join(stopping_points, by = c("true_survival_rate", "recapture_rate")) %>%
    group_by(true_survival_rate, recapture_rate, interval.x) %>%
    summarise(
      total_trials = n_distinct(simulation_data$sim_id[
        simulation_data$true_survival_rate == first(true_survival_rate) & 
          simulation_data$recapture_rate == first(recapture_rate)
      ]),
      
      # Count trials that stopped by this interval
      n_accept_by_now = sum(decision == "Accept H0" & interval.y <= interval.x, na.rm = TRUE),
      n_reject_by_now = sum(decision == "Reject H0" & interval.y <= interval.x, na.rm = TRUE),
      n_stopped_by_now = sum(!is.na(decision) & interval.y <= interval.x, na.rm = TRUE),
      
      # Calculate proportions
      prop_accept_by_now = n_accept_by_now / total_trials,
      prop_reject_by_now = n_reject_by_now / total_trials,
      prop_stopped_by_now = n_stopped_by_now / total_trials,
      prop_continue = 1 - prop_stopped_by_now,
      
      .groups = "drop"
    ) %>%
    rename(interval = interval.x)
  
  return(cumulative_stats)
}

#Determine how many trials were inconclusive
create_inconclusiveness_wide <- function(cumulative_data) {
  
  max_interval <- max(cumulative_data$interval)
  
  wide_table <- cumulative_data %>%
    filter(interval == max_interval) %>%
    select(true_survival_rate, recapture_rate, prop_continue) %>%
    mutate(
      pct_inconclusive = round(prop_continue * 100, 1),
      survival_label = paste0(true_survival_rate * 100, "%")
    ) %>%
    select(survival_label, recapture_rate, pct_inconclusive) %>%
    pivot_wider(names_from = recapture_rate, 
                values_from = pct_inconclusive,
                names_prefix = "Recapture_") %>%
    arrange(survival_label)
  
  return(wide_table)
}

#  Plot functions ####

theme_JN <- function(base_size=10){ 
  theme_grey() %+replace%
    theme(
      axis.text = element_text(colour="black"),
      axis.title = element_text(colour="black"),
      axis.ticks = element_line(colour="black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(colour = "black",fill = NA),
      panel.spacing.x = unit(12, "pt")
    ) 
}

create_decision_plot <- function(mortality_table) {
  # Extract parameters from table
  desired_surv <- unique(mortality_table$desired_surv)[1]
  lower_threshold <- unique(mortality_table$lower_threshold)[1]
  
  # Calculate plot limits
  x_limits <- c(
    min(mortality_table$sample_size, na.rm = TRUE),
    max(mortality_table$sample_size, na.rm = TRUE)
  )
  
  # Define y-axis visual crop limit based on "Continue" zone
  y_crop_max <- mortality_table %>%
    filter(decision == "Continue") %>%
    pull(n_mortality) %>%
    max(na.rm = TRUE)
  
  y_crop_limits <- c(0, y_crop_max * 1.05)
  
  # Create boundary data for regions
  continue_points <- mortality_table %>%
    filter(sample_size > 0, decision == "Continue")
  
  # Build plot with regions if continue points exist
  p <- ggplot() + xlim(x_limits)  # do not include ylim here
  
  if(nrow(continue_points) > 0) {
    # Calculate smooth boundaries
    boundaries <- continue_points %>%
      group_by(sample_size) %>%
      summarise(
        upper_boundary = max(n_mortality) + 0.2,
        lower_boundary = pmax(0, min(n_mortality) - 0.2),
        .groups = "drop"
      )
    
    # Extend boundaries to exact plot limits
    x_seq <- seq(x_limits[1], x_limits[2], length.out = 150)
    
    extended_boundaries <- boundaries %>%
      bind_rows(
        data.frame(
          sample_size = c(x_limits[1], x_limits[2]),
          upper_boundary = c(first(boundaries$upper_boundary), last(boundaries$upper_boundary)),
          lower_boundary = c(first(boundaries$lower_boundary), last(boundaries$lower_boundary))
        )
      ) %>%
      arrange(sample_size)
    
    # Create smooth predictions using extended boundaries
    upper_pred <- predict(loess(upper_boundary ~ sample_size, data = extended_boundaries, span = 1), 
                          newdata = data.frame(sample_size = x_seq))
    lower_pred <- pmax(0, predict(loess(lower_boundary ~ sample_size, data = extended_boundaries, span = 1), 
                                  newdata = data.frame(sample_size = x_seq)))
    
    ribbon_data <- data.frame(x = x_seq, upper = upper_pred, lower = lower_pred)
    
    # Add all regions and boundaries at once
    p <- p + 
      geom_ribbon(data = ribbon_data, aes(x = x, ymin = upper, ymax = y_crop_limits[2]), 
                  fill = "red", alpha = 0.2) +
      geom_ribbon(data = ribbon_data, aes(x = x, ymin = lower, ymax = upper), 
                  fill = "grey", alpha = 0.2) +
      geom_ribbon(data = ribbon_data, aes(x = x, ymin = 0, ymax = lower), 
                  fill = "green", alpha = 0.2) +
      geom_line(data = ribbon_data, aes(x = x, y = upper), color = "black") +
      geom_line(data = ribbon_data, aes(x = x, y = lower), color = "black")
  }
  
  y_breaks <- function(lims) {
    y_min <- floor(lims[1])
    y_max <- ceiling(lims[2])
    range <- y_max - y_min
    
    interval <- case_when(
      y_max < 50 ~ 1,
      y_max < 100 ~ 10,
      TRUE ~ 50
    )
    
    breaks <- seq(y_min, y_max, by = interval)
    
    # Ensure min and max are included
    if (!(y_min %in% breaks)) breaks <- c(y_min, breaks)
    if (!(y_max %in% breaks)) breaks <- c(breaks, y_max)
    
    sort(unique(breaks))
  }
  
  x_breaks <- function(lims) {
    x_min <- floor(lims[1])
    x_max <- ceiling(lims[2])
    
    interval <- if (x_max < 200) 10 else 100
    
    seq(x_min, x_max, by = interval)
  }
  
  
  # Add annotations and styling
  p + 
    annotate("text", x = x_limits[2] * 0.25, y = y_crop_limits[2] * 0.85, 
             label = "Reject H0: Not fish-friendly", size = 4) +
    annotate("text", x = mean(x_limits), y = mean(y_crop_limits), 
             label = "Precision not achieved: Continue testing", size = 4) +
    annotate("text", x = x_limits[2] * 0.75, y = y_crop_limits[2] * 0.15, 
             label = "Accept H0: Fish-friendly", size = 4) +
    annotate("text", x = x_limits[2] * 0.98, y = y_crop_limits[2] * 0.98, 
             label = paste0("CI = 95% (Wilson Binomial)\nDesired survival ≥ ", desired_surv, 
                            "%\nLower bound threshold ≥ ", lower_threshold, "%"), 
             size = 4, hjust = 1, vjust = 1, fontface = "italic") +
    scale_y_continuous(
      breaks = y_breaks(y_crop_limits),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks(x_limits),
      expand = c(0, 0)
    )+
    labs(x = "Sample Size", y = "Mortality (n)") +
    theme_JN() +
    theme(
      panel.grid.major = element_line(color = alpha("black", 0.1), size = 0.5),
      panel.ontop = FALSE 
    ) +
    coord_cartesian(ylim = y_crop_limits, clip = "off")  # visual crop only
}


create_journey_plot <- function(mortality_table, 
                                simulation_data, 
                                survival_rates_to_show =  c(0.98),
                                recapture_rate_to_show = 1.0, 
                                n_journeys = 10,
                                journey_ids = NULL) {
  
  # Filter for specific survival rate AND recapture rate
  filtered_data <- simulation_data %>%
    filter(abs(true_survival_rate - survival_rates_to_show) < 1e-10,
           abs(recapture_rate - recapture_rate_to_show) < 1e-10) 
  
  # Select journeys to show
  if(!is.null(journey_ids)) {
    # Use specific journey IDs if provided
    journey_data <- filtered_data %>%
      filter(sim_id %in% journey_ids)
  } else {
    # Use first n_journeys
    journey_data <- filtered_data %>%
      filter(sim_id %in% 1:n_journeys)
  }
  
  # Start with base decision plot
  base_plot <- create_decision_plot(mortality_table)
  
  # Add journey lines and points
  journey_plot <- base_plot +
    # Grey journey lines (unchanged)
    geom_line(
      data = journey_data, 
      aes(x = cumulative_n, y = cumulative_deaths, group = sim_id),
      color = "grey", alpha = 0.5, size = 0.8
    ) +
    # Density-aware points (size = overlap count)
    geom_count(
      data = journey_data,
      aes(
        x = cumulative_n, 
        y = cumulative_deaths, 
        color = decision,
        size = after_stat(n)  # n is computed by geom_count
      ),
      alpha = 0.6,
    ) +
    scale_size_continuous(
      range = c(2, 8),
    ) +
    # Color scale (unchanged)
    scale_color_manual(
      values = c(
        "Accept H0" = "darkgreen", 
        "Reject H0" = "darkred", 
        "Continue" = "grey50"
      )
    ) +
    {if(length(survival_rates_to_show) > 1) 
      facet_wrap(~paste("True Survival:", true_survival_rate * 100, "%"), 
                 scales = "free")} +
    labs(title = paste0("Sequential Trial simulation (Recapture = ", 
                        recapture_rate_to_show * 100, "%, n = ", 
                        length(unique(journey_data$sim_id)), " simulations per scenario)"))
  return(journey_plot)
}


plot_stopping_distributions <- function(simulation_data) {
  
  stopped_data <- simulation_data %>%
    group_by(sim_id, true_survival_rate, recapture_rate) %>%
    slice_tail(n = 1) %>%
    filter(trial_stopped == TRUE) %>%
    ungroup()
  
  ggplot(stopped_data, aes(x = factor(true_survival_rate), y = cumulative_n, fill = decision)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~paste("Recapture Rate:", recapture_rate)) +
    scale_fill_manual(values = c("Accept H0" = "darkgreen", "Reject H0" = "darkred")) +
    labs(
      x = "True Survival Rate",
      y = "Sample Size at Stopping",
      fill = "Decision",
      title = "Distribution of Sample Sizes at Trial Stopping"
    ) +
    theme_JN()
}

plot_decision_proportions <- function(simulation_data) {
  
  decision_props <- simulation_data %>%
    group_by(sim_id, true_survival_rate, recapture_rate) %>%
    slice_tail(n = 1) %>%
    group_by(true_survival_rate, recapture_rate, decision) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(true_survival_rate, recapture_rate) %>%
    mutate(proportion = count / sum(count))
  
  ggplot(decision_props, aes(x = true_survival_rate, y = proportion, fill = decision)) +
    geom_col(position = "stack") +
    facet_wrap(~paste("Recapture Rate:", recapture_rate)) +
    scale_fill_manual(values = c("Accept H0" = "darkgreen", 
                                 "Reject H0" = "darkred", 
                                 "Continue" = "grey")) +
    scale_x_continuous(labels = scales::percent_format()) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      x = "True Survival Rate",
      y = "Proportion of Trials",
      fill = "Decision"
    ) +
    theme_JN()
}

# Modified plot_cumulative_stopping function to handle decision_method facet

plot_cumulative_stopping <- function(cumulative_data, 
                                     survival_filter = NULL, 
                                     recapture_filter = NULL,
                                     decision_method_filter = NULL) {
  
  # Apply filters if specified
  filtered_data <- cumulative_data
  
  if (!is.null(survival_filter)) {
    filtered_data <- filtered_data %>% 
      filter(abs(true_survival_rate - survival_filter) < 1e-10)
  }
  
  if (!is.null(recapture_filter)) {
    filtered_data <- filtered_data %>% 
      filter(abs(recapture_rate - recapture_filter) < 1e-10)
  }
  
  if (!is.null(decision_method_filter)) {
    filtered_data <- filtered_data %>% 
      filter(decision_method == decision_method_filter)
  }
  
  # Prepare plot data
  plot_data <- filtered_data %>%
    select(true_survival_rate, recapture_rate, decision_method, interval, 
           prop_accept_by_now, prop_reject_by_now, prop_continue) %>%
    pivot_longer(cols = c(prop_accept_by_now, prop_reject_by_now, prop_continue),
                 names_to = "outcome", 
                 values_to = "proportion") %>%
    mutate(
      outcome = case_when(
        outcome == "prop_accept_by_now" ~ "Accept H0",
        outcome == "prop_reject_by_now" ~ "Reject H0", 
        outcome == "prop_continue" ~ "Continue"
      ),
      percentage = proportion * 100
    )
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = interval, y = percentage, color = outcome)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE, span = 0.7) +
    scale_color_manual(values = c("Accept H0" = "darkgreen", 
                                  "Reject H0" = "darkred", 
                                  "Continue" = "grey60")) +
    scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0,100), breaks = seq(0,100,10)) +
    scale_x_continuous(limits = c(100,1000), breaks=seq(100,1000,100)) +
    labs(
      x = "Sample Interval",
      y = "Cumulative Proportion (%)",
      color = "Decision",
      title = "Cumulative Trial Outcomes by Sample Size"
    ) +
    theme_JN() +
    theme(legend.position = "bottom") +
    coord_cartesian(clip = "off")
  
  # Determine faceting based on remaining variables
  n_survival <- length(unique(filtered_data$true_survival_rate))
  n_recapture <- length(unique(filtered_data$recapture_rate))
  n_method <- length(unique(filtered_data$decision_method))
  
  # Add faceting if multiple scenarios remain
  if (n_survival > 1 && n_recapture > 1 && n_method > 1) {
    # All three variables vary - use complex faceting
    p <- p + facet_grid(
      recapture_rate ~ true_survival_rate + decision_method, 
      labeller = labeller(
        recapture_rate = function(x) paste("Recapture:", x),
        true_survival_rate = function(x) paste("Survival:", x),
        decision_method = function(x) paste("Method:", str_to_title(x))
      )
    )
  } else if (n_survival > 1 && n_recapture > 1) {
    # Survival and recapture vary
    p <- p + facet_grid(recapture_rate ~ true_survival_rate, labeller = label_both)
  } else if (n_survival > 1 && n_method > 1) {
    # Survival and method vary
    p <- p + facet_grid(decision_method ~ true_survival_rate, 
                        labeller = labeller(
                          decision_method = function(x) paste("Method:", str_to_title(x)),
                          true_survival_rate = function(x) paste("Survival:", x)
                        ))
  } else if (n_recapture > 1 && n_method > 1) {
    # Recapture and method vary
    p <- p + facet_grid(decision_method ~ recapture_rate,
                        labeller = labeller(
                          decision_method = function(x) paste("Method:", str_to_title(x)),
                          recapture_rate = function(x) paste("Recapture:", x)
                        ))
  } else if (n_survival > 1) {
    # Only survival varies
    p <- p + facet_wrap(~true_survival_rate, labeller = label_both)
  } else if (n_recapture > 1) {
    # Only recapture varies
    p <- p + facet_wrap(~recapture_rate, labeller = label_both)
  } else if (n_method > 1) {
    # Only method varies
    p <- p + facet_wrap(~decision_method, 
                        labeller = labeller(decision_method = function(x) paste("Method:", str_to_title(x))))
  }
  
  return(p)
}

# Rudimentary attempt to animate sequential trials 
create_animated_journey_plot <- function(mortality_table, 
                                         simulation_data, 
                                         survival_rates_to_show = c(0.98),
                                         recapture_rate_to_show = 1.0,
                                         n_journeys = 5,
                                         journey_ids = NULL) {
  
  # Filter and select journeys (same as before)
  filtered_data <- simulation_data %>%
    filter(true_survival_rate %in% survival_rates_to_show,
           abs(recapture_rate - recapture_rate_to_show) < 1e-10)
  
  if(!is.null(journey_ids)) {
    journey_data <- filtered_data %>%
      filter(sim_id %in% journey_ids)
  } else {
    journey_data <- filtered_data %>%
      group_by(true_survival_rate) %>%
      filter(sim_id %in% 1:n_journeys) %>%
      ungroup()
  }
  
  # Create base decision plot
  base_plot <- create_decision_plot(mortality_table)
  
  # Create animated plot
  animated_plot <- base_plot +
    geom_line(
      data = journey_data, 
      aes(x = cumulative_n, y = cumulative_deaths, group = sim_id),
      color = "grey", alpha = 0.5, size = 0.8
    ) +
    # Density-aware points (size = overlap count)
    geom_count(
      data = journey_data,
      aes(
        x = cumulative_n, 
        y = cumulative_deaths, 
        color = decision,
        size = after_stat(n)  # n is computed by geom_count
      ),
      alpha = 0.6,
    ) +
    scale_size_continuous(
      range = c(2, 8),
    ) +
    # Color scale (unchanged)
    scale_color_manual(
      values = c(
        "Accept H0" = "darkgreen", 
        "Reject H0" = "darkred", 
        "Continue" = "grey50"
      )
    )+
    {if(length(survival_rates_to_show) > 1) 
      facet_wrap(~paste("True Survival:", true_survival_rate * 100, "%"))} +
    # Animation component - reveal along cumulative_n
    transition_reveal(cumulative_n) +
    labs(title = "Sequential Trial Journeys Unfolding",
         subtitle = "Sample Size: {closest_state}") +
    ease_aes('linear')
  
  return(animated_plot)
}





calculate_ht_estimates <- function(exposed_fish, recaptured_fish, observed_deaths) {
  
  if (recaptured_fish == 0) {
    # Handle edge case of no recaptures
    return(list(
      deaths_point = exposed_fish * 0.5,  
      deaths_se = Inf,
      deaths_ci_lower = 0,
      deaths_ci_upper = exposed_fish,
      method = "No recaptures"
    ))
  }
  
  recapture_rate <- recaptured_fish / exposed_fish
  observed_survivors <- recaptured_fish - observed_deaths
  
  if (recapture_rate == 1.0) {
    # Perfect recapture - no uncertainty
    return(list(
      deaths_point = observed_deaths,
      deaths_se = 0,
      deaths_ci_lower = observed_deaths,
      deaths_ci_upper = observed_deaths,
      method = "Perfect recapture"
    ))
  }
  
  # Create data frame for survey design
  recaptured_data <- data.frame(
    death = c(rep(1, observed_deaths), rep(0, observed_survivors)),
    weight = 1 / recapture_rate,
    fpc = exposed_fish
  )
  
  # Survey design with finite population correction
  design <- svydesign(
    ids = ~1,           
    weights = ~weight,  
    fpc = ~fpc,        
    data = recaptured_data
  )
  
  # Calculate HT estimate
  ht_total <- svytotal(~death, design)
  ht_ci <- confint(ht_total)
  
  return(list(
    deaths_point = as.numeric(ht_total),
    deaths_se = as.numeric(SE(ht_total)),
    deaths_ci_lower = ht_ci[1],
    deaths_ci_upper = ht_ci[2],
    method = "HT with FPC"
  ))
}

# Updated simulation function using HT estimates
simulate_sequential_trial_ht <- function(simulations = 1:50, 
                                         true_survival_rates = seq(0.9, 1, by = 0.01), 
                                         max_sample = 500, 
                                         check_intervals = c(100, 200, 300, 400, 500),
                                         mortality_table,
                                         recapture_rates = c(1.0),
                                         decision_method = "neutral") {  # "aggressive", "neutral", "conservative"
  
  # Initialize results
  all_results <- tibble()
  
  # Loop through each simulation
  for(sim_id in simulations) {
    # Loop through each recapture rate  
    for(recapture_rate in recapture_rates) {
      # Loop through each survival rate
      for(survival_rate in true_survival_rates) {
        
        # Initialize for this specific simulation
        results <- tibble(
          interval = integer(),
          cumulative_n = integer(),
          cumulative_survivors = integer(),
          cumulative_deaths = integer(),
          decision = character(),
          trial_stopped = logical(),
          ht_deaths_point = numeric(),
          ht_deaths_ci_lower = numeric(),
          ht_deaths_ci_upper = numeric(),
          decision_method = character()
        )
        
        # Track cumulative observed data (what we actually see)
        cumulative_recaptured <- 0
        cumulative_deaths_observed <- 0
        
        # Run through check intervals
        for(interval in check_intervals) {
          if(interval == check_intervals[1]) {
            interval_size <- interval
          } else {
            prev_n <- check_intervals[which(check_intervals == interval) - 1]
            interval_size <- interval - prev_n
          }
          
          # Generate true outcomes for this interval
          true_new_survivors <- rbinom(1, interval_size, survival_rate)
          true_new_deaths <- interval_size - true_new_survivors
          
          # Apply recapture rate
          recaptured_alive <- rbinom(1, true_new_survivors, recapture_rate)
          recaptured_dead <- rbinom(1, true_new_deaths, recapture_rate)
          
          # Update cumulative observed data
          cumulative_recaptured <- cumulative_recaptured + recaptured_alive + recaptured_dead
          cumulative_deaths_observed <- cumulative_deaths_observed + recaptured_dead
          
          # Calculate HT estimates for cumulative data
          ht_result <- calculate_ht_estimates(interval, cumulative_recaptured, cumulative_deaths_observed)
          
          # Choose deaths estimate based on decision method
          if (decision_method == "aggressive") {
            deaths_for_decision <- ht_result$deaths_ci_lower
          } else if (decision_method == "conservative") {
            deaths_for_decision <- ht_result$deaths_ci_upper
          } else {  # neutral
            deaths_for_decision <- ht_result$deaths_point
          }
          
          # Round for mortality table lookup
          deaths_rounded <- round(deaths_for_decision)
          
          # Look up decision using chosen estimate
          decision_row <- mortality_table %>% 
            filter(sample_size == interval, n_mortality == deaths_rounded)
          
          current_decision <- if(nrow(decision_row) > 0) decision_row$decision[1] else "Continue"
          
          # Calculate display values using point estimates
          cumulative_survivors <- interval - ht_result$deaths_point
          
          # Add to results
          results <- results %>% 
            add_row(
              interval = interval,
              cumulative_n = interval,
              cumulative_survivors = cumulative_survivors,
              cumulative_deaths = ht_result$deaths_point,
              decision = current_decision,
              trial_stopped = current_decision != "Continue",
              ht_deaths_point = ht_result$deaths_point,
              ht_deaths_ci_lower = ht_result$deaths_ci_lower,
              ht_deaths_ci_upper = ht_result$deaths_ci_upper,
              decision_method = decision_method
            )
          
          # Stop if decision reached
          if(current_decision != "Continue") break
        }
        
        # Add simulation metadata and combine
        results <- results %>%
          mutate(sim_id = sim_id, 
                 true_survival_rate = survival_rate,
                 recapture_rate = recapture_rate)
        
        all_results <- bind_rows(all_results, results)
      }
    }
  }
  
  return(all_results)
}

# Updated journey plot function to handle HT results
create_journey_plot_ht <- function(mortality_table, 
                                   simulation_data, 
                                   survival_rates_to_show = c(0.98),
                                   recapture_rate_to_show = 1.0, 
                                   n_journeys = 10,
                                   journey_ids = NULL,
                                   decision_method = "neutral") {
  
  # Filter for specific survival rate, recapture rate, and decision method
  filtered_data <- simulation_data %>%
    filter(abs(true_survival_rate - survival_rates_to_show) < 1e-10,
           abs(recapture_rate - recapture_rate_to_show) < 1e-10,
           decision_method == !!decision_method) 
  
  # Select journeys to show
  if(!is.null(journey_ids)) {
    journey_data <- filtered_data %>%
      filter(sim_id %in% journey_ids)
  } else {
    journey_data <- filtered_data %>%
      filter(sim_id %in% 1:n_journeys)
  }
  
  # Start with base decision plot
  base_plot <- create_decision_plot(mortality_table)
  
  # Add journey lines and points
  journey_plot <- base_plot +
    # Journey lines
    geom_line(
      data = journey_data, 
      aes(x = cumulative_n, y = cumulative_deaths, group = sim_id),
      color = "grey", alpha = 0.5, size = 0.8
    ) +
    # Decision points
    geom_count(
      data = journey_data,
      aes(
        x = cumulative_n, 
        y = cumulative_deaths, 
        color = decision,
        size = after_stat(n)
      ),
      alpha = 0.6,
    ) +
    scale_size_continuous(
      range = c(2, 8),
    ) +
    scale_color_manual(
      values = c(
        "Accept H0" = "darkgreen", 
        "Reject H0" = "darkred", 
        "Continue" = "grey50"
      )
    ) +
    {if(length(survival_rates_to_show) > 1) 
      facet_wrap(~paste("True Survival:", true_survival_rate * 100, "%"), 
                 scales = "free")} +
    labs(title = paste0("Sequential Trial Simulation - HT Method (", 
                        str_to_title(decision_method), 
                        " | Recapture = ", recapture_rate_to_show * 100, 
                        "%, n = ", length(unique(journey_data$sim_id)), 
                        " simulations per scenario)"))
  
  return(journey_plot)
}

# Function to run HT simulations with all three decision scenarios
run_ht_all_scenarios <- function(mortality_table, 
                                 simulations = 1:100,
                                 true_survival_rates = seq(0.95, 1, 0.01),
                                 max_sample = 500,
                                 recapture_rates = c(0.9, 0.95, 1.0),
                                 check_intervals = c(100, 200, 300, 400, 500)) {
  
  # Run simulations for all three decision methods
  methods <- c("aggressive", "neutral", "conservative")
  all_results <- tibble()
  
  for (method in methods) {
    cat("Running simulations for", method, "method...\n")
    
    method_results <- simulate_sequential_trial_ht(
      simulations = simulations,
      true_survival_rates = true_survival_rates,
      max_sample = max_sample,
      check_intervals = check_intervals,
      mortality_table = mortality_table,
      recapture_rates = recapture_rates,
      decision_method = method
    )
    
    all_results <- bind_rows(all_results, method_results)
  }
  
  return(all_results)
}

# Updated stopping summary for HT results
create_stopping_summary_ht <- function(simulation_data) {
  
  # Get final results for each simulation
  final_results <- simulation_data %>%
    group_by(sim_id, true_survival_rate, recapture_rate, decision_method) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  # Stopping time analysis for trials that stopped
  stopped_trials <- final_results %>%
    filter(trial_stopped == TRUE) %>%
    group_by(true_survival_rate, recapture_rate, decision_method, decision) %>%
    summarise(
      n_stopped = n(),
      median_fish = median(cumulative_n),
      mean_fish = round(mean(cumulative_n), 1),
      q25_fish = quantile(cumulative_n, 0.25),
      q75_fish = quantile(cumulative_n, 0.75),
      min_fish = min(cumulative_n),
      max_fish = max(cumulative_n),
      .groups = "drop"
    )
  
  # Overall trial summary
  trial_summary <- final_results %>%
    group_by(true_survival_rate, recapture_rate, decision_method) %>%
    summarise(
      total_trials = n(),
      n_accept = sum(decision == "Accept H0", na.rm = TRUE),
      n_reject = sum(decision == "Reject H0", na.rm = TRUE), 
      n_continue = sum(decision == "Continue", na.rm = TRUE),
      prop_accept = round(n_accept / total_trials, 3),
      prop_reject = round(n_reject / total_trials, 3),
      prop_continue = round(n_continue / total_trials, 3),
      .groups = "drop"
    )
  
  return(list(
    stopping_details = stopped_trials,
    trial_overview = trial_summary
  ))
}







# Updated journey plot function with fixed sizing and comprehensive faceting

create_journey_plot_ht_comprehensive <- function(mortality_table, 
                                                 simulation_data, 
                                                 survival_rates_to_show = c(0.98),
                                                 recapture_rates_to_show = c(1.0), 
                                                 decision_methods_to_show = c("neutral"),
                                                 n_journeys = 10,
                                                 journey_ids = NULL) {
  
  # Filter for specified parameters and order factors
  filtered_data <- simulation_data %>%
    filter(true_survival_rate %in% survival_rates_to_show,
           recapture_rate %in% recapture_rates_to_show,
           decision_method %in% decision_methods_to_show) %>%
    mutate(
      # Order survival rates using the input order
      true_survival_rate = factor(true_survival_rate, 
                                  levels = survival_rates_to_show),
      # Order decision methods  
      decision_method = factor(decision_method, 
                               levels = c("aggressive", "neutral", "conservative"))
    ) 
  
  # Select journeys to show
  if(!is.null(journey_ids)) {
    journey_data <- filtered_data %>%
      filter(sim_id %in% journey_ids)
  } else {
    journey_data <- filtered_data %>%
      group_by(true_survival_rate, recapture_rate, decision_method) %>%
      filter(sim_id %in% 1:n_journeys) %>%
      ungroup()
  }
  
  # Start with base decision plot
  base_plot <- create_decision_plot(mortality_table)
  
  # Add journey lines and points with fixed sizing
  journey_plot <- base_plot +
    # Journey lines
    geom_line(
      data = journey_data, 
      aes(x = cumulative_n, y = cumulative_deaths, group = sim_id),
      color = "grey", alpha = 0.5, size = 0.8
    ) +
    # Decision points with fixed size scale
    geom_count(
      data = journey_data,
      aes(
        x = cumulative_n, 
        y = cumulative_deaths, 
        color = decision,
        size = after_stat(pmin(n, 50))  # Cap at 50, but use actual count
      ),
      alpha = 0.6,
    ) +
    # Fixed size scale 1-50, with fallback for higher values
    scale_size_continuous(
      range = c(1, 8),  # Size range from 1 to 8 pixels
      breaks = c(1, 10, 25, 50),
      labels = c("1", "10", "25", "50+"),
      limits = c(1, 50),
      guide = "none"  # Remove size legend
    ) +
    # Color scale
    scale_color_manual(
      values = c(
        "Accept H0" = "darkgreen", 
        "Reject H0" = "darkred", 
        "Continue" = "grey50"
      ),
      guide = "none"  # Remove color legend
    ) +
    # Comprehensive faceting with ordered methods
    facet_grid(
      rows = vars(paste("Recapture:", recapture_rate * 100, "%")),
      cols = vars(
        paste("Survival:", as.numeric(as.character(true_survival_rate)) * 100, "%"),
        paste("Method:", str_to_title(decision_method))
      ),
      scales = "free",
      space = "free"
    ) +
    # Updated title
    labs(
      title = paste0("Sequential Trial Simulation - HT Methods Comparison"),
      subtitle = paste0("n = ", n_journeys, " simulations per scenario")
    ) +
    # Remove all legends
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),  # Smaller facet labels
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10)
    )
  
  return(journey_plot)
}


# Alternative version with single facet_wrap if grid gets too complex
create_journey_plot_ht_wrapped <- function(mortality_table, 
                                           simulation_data, 
                                           survival_rates_to_show = c(0.98),
                                           recapture_rates_to_show = c(1.0), 
                                           decision_methods_to_show = c("neutral"),
                                           n_journeys = 10,
                                           journey_ids = NULL) {
  
  # Filter for specified parameters
  filtered_data <- simulation_data %>%
    filter(true_survival_rate %in% survival_rates_to_show,
           recapture_rate %in% recapture_rates_to_show,
           decision_method %in% decision_methods_to_show) 
  
  # Select journeys to show
  if(!is.null(journey_ids)) {
    journey_data <- filtered_data %>%
      filter(sim_id %in% journey_ids)
  } else {
    journey_data <- filtered_data %>%
      group_by(true_survival_rate, recapture_rate, decision_method) %>%
      filter(sim_id %in% 1:n_journeys) %>%
      ungroup()
  }
  
  # Create combined facet label
  journey_data <- journey_data %>%
    mutate(
      facet_label = paste0(
        "Survival: ", true_survival_rate * 100, "% | ",
        "Recapture: ", recapture_rate * 100, "% | ",
        "Method: ", str_to_title(decision_method)
      )
    )
  
  # Start with base decision plot
  base_plot <- create_decision_plot(mortality_table)
  
  # Add journey lines and points
  journey_plot <- base_plot +
    # Journey lines
    geom_line(
      data = journey_data, 
      aes(x = cumulative_n, y = cumulative_deaths, group = sim_id),
      color = "grey", alpha = 0.5, size = 0.8
    ) +
    # Decision points with capped sizing
    geom_count(
      data = journey_data,
      aes(
        x = cumulative_n, 
        y = cumulative_deaths, 
        color = decision,
        size = after_stat(pmin(n, 50))  # Cap at 50
      ),
      alpha = 0.6,
    ) +
    # Fixed size scale
    scale_size_continuous(
      range = c(1, 8),
      breaks = c(1, 10, 25, 50),
      labels = c("1", "10", "25", "50+"),
      limits = c(1, 50),
      guide = "none"
    ) +
    # Color scale
    scale_color_manual(
      values = c(
        "Accept H0" = "darkgreen", 
        "Reject H0" = "darkred", 
        "Continue" = "grey50"
      ),
      guide = "none"
    ) +
    # Wrapped faceting
    facet_wrap(
      ~facet_label, 
      scales = "free",
      ncol = 3  # Adjust as needed
    ) +
    # Remove legends and adjust text
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(
      title = "Sequential Trial Simulation - HT Methods Comparison",
      subtitle = paste0("n = ", n_journeys, " simulations per scenario")
    )
  
  return(journey_plot)
}

# Usage examples:

# Example 1: Compare all three methods for one scenario
# journey_plot_all_methods <- create_journey_plot_ht_comprehensive(
#   mortality_data,
#   ht_simulation_results_fixed,
#   survival_rates_to_show = c(0.98),
#   recapture_rates_to_show = c(1.0),
#   decision_methods_to_show = c("aggressive", "neutral", "conservative"),
#   n_journeys = 20
# )

# Example 2: Compare across recapture rates for one method
# journey_plot_recapture <- create_journey_plot_ht_comprehensive(
#   mortality_data,
#   ht_simulation_results_fixed,
#   survival_rates_to_show = c(0.98),
#   recapture_rates_to_show = c(0.9, 0.95, 1.0),
#   decision_methods_to_show = c("neutral"),
#   n_journeys = 20
# )

# Example 3: Full comparison (might be very busy)
# journey_plot_comprehensive <- create_journey_plot_ht_wrapped(
#   mortality_data,
#   ht_simulation_results_fixed,
#   survival_rates_to_show = c(0.96, 0.98, 1.0),
#   recapture_rates_to_show = c(0.9, 1.0),
#   decision_methods_to_show = c("neutral", "conservative"),
#   n_journeys = 10
# )





# Function to analyze real study data
analyze_real_study_data_enhanced <- function(exposed_fish, recaptured_fish, observed_deaths, 
                                             mortality_table, cumulative_data, 
                                             assumed_true_survival = 0.98,
                                             recapture_rate_scenario = 1.0) {
  
  # Calculate basic rates
  recapture_rate <- recaptured_fish / exposed_fish
  observed_survivors <- recaptured_fish - observed_deaths
  survival_rate_recaptured <- observed_survivors / recaptured_fish
  
  # Calculate HT estimates using your existing function
  ht_result <- calculate_ht_estimates(exposed_fish, recaptured_fish, observed_deaths)
  
  # Three decision scenarios
  scenarios <- list(
    aggressive = ht_result$deaths_ci_lower,
    neutral = ht_result$deaths_point,
    conservative = ht_result$deaths_ci_upper
  )
  
  # Look up decisions for each scenario
  decisions <- list()
  
  cat("=== REAL STUDY DATA ANALYSIS ===\n")
  cat("Exposed fish:", exposed_fish, "\n")
  cat("Recaptured fish:", recaptured_fish, "\n")
  cat("Observed deaths:", observed_deaths, "\n")
  cat("Survival rate (recaptured):", round(survival_rate_recaptured * 100, 2), "%\n")
  cat("Recapture rate:", round(recapture_rate * 100, 1), "%\n\n")
  
  cat("HT ESTIMATES:\n")
  cat("Point estimate:", round(ht_result$deaths_point, 2), "deaths\n")
  cat("95% CI: [", round(ht_result$deaths_ci_lower, 2), " - ", round(ht_result$deaths_ci_upper, 2), "] deaths\n\n")
  
  cat("DECISION SCENARIOS:\n")
  for (scenario_name in names(scenarios)) {
    deaths_estimate <- scenarios[[scenario_name]]
    deaths_rounded <- round(deaths_estimate)
    
    # Look up in mortality table
    decision_row <- mortality_table %>% 
      filter(sample_size == exposed_fish, n_mortality == deaths_rounded)
    
    if (nrow(decision_row) > 0) {
      decision <- decision_row$decision[1]
      survival_formatted <- decision_row$survival_formatted[1]
    } else {
      decision <- "Not in table"
      survival_formatted <- "Not available"
    }
    
    cat(paste(toupper(scenario_name), ":\n"))
    cat("  Deaths estimate:", round(deaths_estimate, 2), "\n")
    cat("  Rounded for lookup:", deaths_rounded, "\n")
    cat("  Decision:", decision, "\n")
    cat("  Survival estimate:", survival_formatted, "\n\n")
    
    decisions[[scenario_name]] <- list(
      deaths_estimate = deaths_estimate,
      deaths_rounded = deaths_rounded,
      decision = decision,
      survival_formatted = survival_formatted
    )
  }
  
  # PROBABILITY PREDICTIONS FROM SIMULATION DATA - ALL THREE METHODS
  cat("=== PROBABILITY PREDICTIONS ===\n")
  cat("Assuming true survival =", assumed_true_survival * 100, "%\n")
  cat("Assuming recapture rate =", recapture_rate_scenario * 100, "%\n\n")
  
  # Handle edge cases gracefully
  available_survival_rates <- unique(cumulative_data$true_survival_rate)
  available_recapture_rates <- unique(cumulative_data$recapture_rate)
  
  # Find closest available rates if exact match not found
  if (!assumed_true_survival %in% available_survival_rates) {
    closest_survival <- available_survival_rates[which.min(abs(available_survival_rates - assumed_true_survival))]
    cat("⚠️  Assumed survival rate", assumed_true_survival, "not in simulation data.\n")
    cat("   Using closest available rate:", closest_survival, "\n\n")
    assumed_true_survival <- closest_survival
  }
  
  if (!recapture_rate_scenario %in% available_recapture_rates) {
    closest_recapture <- available_recapture_rates[which.min(abs(available_recapture_rates - recapture_rate_scenario))]
    cat("⚠️  Assumed recapture rate", recapture_rate_scenario, "not in simulation data.\n")
    cat("   Using closest available rate:", closest_recapture, "\n\n")
    recapture_rate_scenario <- closest_recapture
  }
  
  # Loop through all three decision methods
  decision_methods <- c("aggressive", "neutral", "conservative")
  
  for (method in decision_methods) {
    cat(paste("=== ", toupper(method), " METHOD ===\n"))
    
    # Find current interval probabilities
    current_interval_data <- cumulative_data %>%
      filter(abs(true_survival_rate - assumed_true_survival) < 1e-10,
             abs(recapture_rate - recapture_rate_scenario) < 1e-10,
             decision_method == method,
             interval == exposed_fish)
    
    if (nrow(current_interval_data) > 0) {
      cat("CURRENT INTERVAL (", exposed_fish, " fish):\n")
      cat("  Probability of Accept H0:", round(current_interval_data$prop_accept_by_now * 100, 1), "%\n")
      cat("  Probability of Reject H0:", round(current_interval_data$prop_reject_by_now * 100, 1), "%\n")
      cat("  Probability of Continue:", round(current_interval_data$prop_continue * 100, 1), "%\n\n")
    } else {
      cat("CURRENT INTERVAL (", exposed_fish, " fish): No simulation data available\n\n")
    }
    
    # Find next interval probabilities
    next_intervals <- cumulative_data %>%
      filter(abs(true_survival_rate - assumed_true_survival) < 1e-10,
             abs(recapture_rate - recapture_rate_scenario) < 1e-10,
             decision_method == method,
             interval > exposed_fish) %>%
      arrange(interval)
    
    if (nrow(next_intervals) > 0) {
      next_interval <- next_intervals$interval[1]
      next_data <- next_intervals[1, ]
      
      cat("NEXT INTERVAL (", next_interval, " fish):\n")
      cat("  Probability of Accept H0:", round(next_data$prop_accept_by_now * 100, 1), "%\n")
      cat("  Probability of Reject H0:", round(next_data$prop_reject_by_now * 100, 1), "%\n")
      cat("  Probability of Continue:", round(next_data$prop_continue * 100, 1), "%\n\n")
      
      # Show marginal benefit of additional sampling
      if (nrow(current_interval_data) > 0) {
        marginal_accept <- (next_data$prop_accept_by_now - current_interval_data$prop_accept_by_now) * 100
        marginal_reject <- (next_data$prop_reject_by_now - current_interval_data$prop_reject_by_now) * 100
        marginal_decision <- marginal_accept + marginal_reject
        
        cat("MARGINAL BENEFIT OF ADDITIONAL", next_interval - exposed_fish, "FISH:\n")
        cat("  Additional Accept probability:", round(marginal_accept, 1), "%\n")
        cat("  Additional Reject probability:", round(marginal_reject, 1), "%\n")
        cat("  Total additional decision probability:", round(marginal_decision, 1), "%\n")
        
        if (marginal_decision < 1) {
          cat("⚠️  Marginal benefit is <1% - additional sampling may not be worthwhile\n")
        }
        cat("\n")
      }
    } else {
      cat("NEXT INTERVAL: No further intervals in simulation data\n\n")
    }
  }
  
  return(invisible(decisions))
}
