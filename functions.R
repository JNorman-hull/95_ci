

#  Statistical tools  ####

# Build mortality table with 95% CI (Wilson binomial)

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


simulate_sequential_trial <- function(simulations = 1:50, 
                                      true_survival_rates = seq(0.9, 1, by = 0.01), 
                                      max_sample = 500, 
                                      check_intervals = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
                                      mortality_table) {
  
  # Initialize results
  all_results <- tibble()
  
  # Loop through each simulation
  for(sim_id in simulations) {
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
          new_survivors <- rbinom(1, interval, survival_rate)
        } else {
          prev_n <- check_intervals[which(check_intervals == interval) - 1]
          interval_size <- interval - prev_n
          new_survivors <- rbinom(1, interval_size, survival_rate)
        }
        
        cumulative_survivors <- cumulative_survivors + new_survivors
        cumulative_deaths <- interval - cumulative_survivors
        
        # Look up decision
        decision_row <- mortality_table %>% 
          filter(sample_size == interval, n_mortality == cumulative_deaths)
        
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
        mutate(sim_id = sim_id, true_survival_rate = survival_rate)
      
      all_results <- bind_rows(all_results, results)
    }
  }
  
  return(all_results)
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


#mortality_table dependency is a tibble cereated by create_mortality_table() can be named otherwise

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
                                survival_rate_to_show = 0.98,
                                n_journeys = 10,
                                journey_ids = NULL) {
  
  # Filter for specific survival rate
  filtered_data <- simulation_data %>%
    filter(true_survival_rate == survival_rate_to_show)
  
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
    geom_line(data = journey_data, 
              aes(x = cumulative_n, y = cumulative_deaths, group = sim_id),
              color = "blue", alpha = 0.7, size = 0.8) +
    geom_point(data = journey_data, 
               aes(x = cumulative_n, y = cumulative_deaths, 
                   color = decision, shape = trial_stopped),
               size = 2) +
    scale_color_manual(values = c("Accept H0" = "darkgreen", 
                                  "Reject H0" = "darkred", 
                                  "Continue" = "grey50")) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
    labs(title = paste0("Sequential Trial Journeys (True Survival = ", 
                        survival_rate_to_show * 100, "%, n = ", 
                        length(unique(journey_data$sim_id)), " journeys)"))
  
  return(journey_plot)
}