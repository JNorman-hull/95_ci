#  Statistical tools  ####

## 95% CI tools  ####

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

#Horvitz–Thompson estimator

horvitz_thompson_est <- function(exposed_fish, recaptured_fish, observed_deaths) {
  
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
simulate_sequential_trial <- function(simulations = 1:50, 
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
          ht_result <- horvitz_thompson_est(interval, cumulative_recaptured, cumulative_deaths_observed)
          
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


# Function to run HT simulations with all three decision scenarios
simulate_all_HT <- function(mortality_table, 
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
    
    method_results <- simulate_sequential_trial(
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

## Hypotheses tools  ####

#Power analysis - Exact binomial test for compliance
find_compliance_sample_size <- function(
    compliance_threshold = 0.98,  # Target rate (e.g., 98%)
    true_rate = 0.96,             # True rate to detect (e.g., 99%)
    alpha = 0.05,                 # Significance level
    target_power = 0.8,           # Desired power
    test_direction = "less",   # "greater" or "less"
    max_n = 50000                 # Max sample size to check
) {
  for (n in 1:max_n) {
    if (test_direction == "greater") {
      # Reject H0 (compliance threshold) if survivors ≥ critical_value
      critical_value <- qbinom(1 - alpha, n, compliance_threshold) + 1
      power <- 1 - pbinom(critical_value - 1, n, true_rate)
    } else if (test_direction == "less") {
      # Reject H0 if survivors ≤ critical_value
      critical_value <- qbinom(alpha, n, compliance_threshold) - 1
      power <- pbinom(critical_value, n, true_rate)
    } else {
      stop("test_direction must be 'greater' or 'less'")
    }
    
    if (power >= target_power) {
      return(list(
        n = n,
        power = power,
        compliance_threshold = compliance_threshold,
        rejection_rule = ifelse(
          test_direction == "greater",
          paste("Reject (compliance ≤ ", compliance_threshold, ") if ≥ ", critical_value, "/", n, " survivors"),
          paste("Reject (compliance ≥ ", compliance_threshold, ") if ≤ ", critical_value, "/", n, " survivors")
        ),
        true_rate_for_power = true_rate
      ))
    }
  }
  warning("No sample size achieved target power within max_n. Try increasing max_n or relaxing constraints.")
  return(NULL)
}


# Modified Fisher's Exact Power Analysis Functions
# Uses survival-first matrix structure consistently

simulate_fisher_power <- function(n_impact, n_control, p_impact, p_control, alpha = 0.05, alternative = "less") {
  reps <- 1000  # Hardcoded
  sig_count <- 0L
  
  for (i in seq_len(reps)) {
    # Simulate survivors
    surv_impact <- rbinom(1, n_impact, p_impact)
    surv_control <- rbinom(1, n_control, p_control)
    
    # Calculate deaths
    dead_impact <- n_impact - surv_impact
    dead_control <- n_control - surv_control
    
    # Survival-first matrix: [survived, died]
    # Row 1 = impact, Row 2 = control
    survival_matrix <- matrix(c(surv_impact, dead_impact, surv_control, dead_control), 
                              nrow = 2, byrow = TRUE)
    
    # Test with specified alternative
    p_val <- fisher.test(survival_matrix, alternative = alternative)$p.value
    if (p_val < alpha) sig_count <- sig_count + 1L
  }
  
  return(sig_count / reps)
}

find_min_sample_size_fisher <- function(p_impact, p_control, alpha = 0.05, power_target = 0.8,
                                        group_ratio = 3, step_size = 10, alternative = "less") {
  
  max_n_control <- 2000  # Hardcoded
  
  for (n_control in seq(100, max_n_control, by = step_size)) {
    n_impact <- round(n_control * group_ratio)
    
    power_est <- simulate_fisher_power(n_impact, n_control, p_impact, p_control, alpha, alternative)
    
    cat("control n =", n_control, ", impact n =", n_impact, ", power =", round(power_est, 3), "\n")
    
    if (!is.na(power_est) && power_est >= power_target) {
      # Calculate final p-value for output
      final_p <- simulate_fisher_power(n_impact, n_control, p_impact, p_control, alpha, alternative)
      
      cat("\n=== FINAL RESULT ===\n")
      cat("control n =", n_control, ", impact n =", n_impact, ", p =", round(final_p, 3), ", power =", round(power_est, 3), "\n")
      
      return(list(
        n_control = n_control,
        n_impact = n_impact, 
        p_value = final_p,
        power = power_est,
        alternative = alternative
      ))
    }
  }
  
  warning("Could not achieve target power within maximum sample size. Try increasing power_target tolerance or check effect size.")
  return(NULL)
}

# function may be a replicate, check
calculate_fisher_exact <- function(total_deaths, total_n) {
  # 3:1 allocation (impact:control) 
  impact_n <- round(total_n * 0.75)
  control_n <- total_n - impact_n
  
  # All observed deaths assumed to occur in impact group (worst case)
  impact_deaths <- total_deaths
  control_deaths <- 0  # Assume perfect control survival for boundary calculation
  
  # Calculate survivors
  surv_impact <- impact_n - impact_deaths
  surv_control <- control_n - control_deaths
  
  # Create survival matrix using your existing structure
  # Row 1 = impact, Row 2 = control, Cols = [survived, died]  
  survival_matrix <- matrix(c(surv_impact, impact_deaths, surv_control, control_deaths), 
                            nrow = 2, byrow = TRUE)
  
  # One-sided test: H1 = impact survival < control survival (pump is harmful)
  # alternative = "less" tests if impact survival is significantly lower
  if (sum(survival_matrix) > 0 && all(rowSums(survival_matrix) > 0) && all(colSums(survival_matrix) > 0)) {
    test_result <- fisher.test(survival_matrix, alternative = "less")
    return(test_result$p.value)
  } else {
    return(1)  # No evidence if invalid matrix
  }
}

# Find death thresholds for each analysis using p_to_z_calculator 
find_death_thresholds <- function(sample_size, z_upper, z_lower) {
  max_deaths <- ceiling(sample_size * 0.15)  
  
  # Find efficacy threshold (stop when >= this many deaths)
  efficacy_threshold <- NA
  for (deaths in 1:max_deaths) {
    p_val <- calculate_fisher_exact(deaths, sample_size)
    if (p_val > 0 && p_val < 1) {  
      z_score <- qnorm(p_val, lower.tail = FALSE)
      if (z_score >= z_upper) {
        efficacy_threshold <- deaths
        break
      }
    }
  }
  
  # Find futility threshold (stop when <= this many deaths)
  # We want the MAXIMUM number of deaths that still allows futility stopping
  futility_threshold <- NA
  for (deaths in max_deaths:0) {  # Search from high deaths down to 0
    p_val <- calculate_fisher_exact(deaths, sample_size)
    if (p_val > 0 && p_val < 1) {  
      z_score <- qnorm(p_val, lower.tail = FALSE)  
      if (z_score <= z_lower) {  # If this death count meets futility boundary
        futility_threshold <- deaths  # This is the maximum deaths for futility
        break  # Stop at first (highest) death count that meets criteria
      }
    }
  }
  
  return(list(efficacy = efficacy_threshold, futility = futility_threshold))
}

#Enhanced function to show death estimates across multiple recapture rates
show_recapture_death_spectrum <- function(comprehensive_boundaries, 
                                          recapture_rates = seq(80, 100, by = 1)) {
  
  # Create expanded table showing death estimates at different recapture rates
  results_list <- list()
  
  for(i in 1:nrow(comprehensive_boundaries)) {
    row_data <- comprehensive_boundaries[i, ]
    
    # For futility boundary
    if(!is.na(row_data$deaths_for_futility)) {
      observed_futility <- row_data$deaths_for_futility
      
      futility_estimates <- data.frame(
        analysis = row_data$analysis,
        sample_size = row_data$sample_size,
        boundary_type = "Futility",
        original_deaths = observed_futility,
        recapture_rate = recapture_rates,
        estimated_deaths = sapply(recapture_rates, function(rate) {
          if(rate == 100) return(observed_futility)
          recaptured <- round(row_data$sample_size * rate / 100)
          ht_estimate <- (observed_futility / recaptured) * row_data$sample_size
          round(ht_estimate)
        })
      )
      results_list[[paste0("futility_", i)]] <- futility_estimates
    }
    
    # For efficacy boundary  
    if(!is.na(row_data$deaths_for_efficacy)) {
      # Use 1 less than efficacy threshold as starting observed deaths
      observed_efficacy <- row_data$deaths_for_efficacy - 1
      
      efficacy_estimates <- data.frame(
        analysis = row_data$analysis,
        sample_size = row_data$sample_size,
        boundary_type = "Efficacy",
        original_deaths = observed_efficacy,
        recapture_rate = recapture_rates,
        estimated_deaths = sapply(recapture_rates, function(rate) {
          if(rate == 100) return(observed_efficacy)
          recaptured <- round(row_data$sample_size * rate / 100)
          ht_estimate <- (observed_efficacy / recaptured) * row_data$sample_size
          round(ht_estimate)
        })
      )
      results_list[[paste0("efficacy_", i)]] <- efficacy_estimates
    }
  }
  
  # Combine all results
  full_results <- do.call(rbind, results_list)
  
  # Add interpretation
  full_results <- full_results %>%
    mutate(
      death_increase = estimated_deaths - original_deaths,
      interpretation = case_when(
        boundary_type == "Futility" & death_increase > 0 ~ paste0("Prevents futility stopping (+", death_increase, " deaths)"),
        boundary_type == "Efficacy" & death_increase > 0 ~ paste0("Enables efficacy stopping (+", death_increase, " deaths)"),
        TRUE ~ "No change"
      )
    )
  
  return(full_results)
}



find_recapture_boundary_effects <- function(comprehensive_boundaries, fixed_recapture = 90) {
  
  results <- comprehensive_boundaries %>%
    select(analysis, sample_size, deaths_for_futility, deaths_for_efficacy) %>%
    rowwise() %>%
    mutate(
      # FUTILITY PREVENTION: When does H-T push deaths UP and prevent futility stopping?
      futility_prevention_recapture = if(!is.na(deaths_for_futility)) {
        observed_deaths <- deaths_for_futility
        # Find recapture rate where observed deaths round UP to next integer
        critical_recaptured <- (observed_deaths * sample_size) / (observed_deaths + 0.5)
        critical_rate <- critical_recaptured / sample_size
        round(critical_rate * 100, 1)
      } else NA,
      
      # EFFICACY FACILITATION: When does H-T push deaths UP to reach efficacy boundary?
      efficacy_facilitation_recapture = if(!is.na(deaths_for_efficacy)) {
        observed_deaths <- deaths_for_efficacy - 1  # Start with 1 less death
        # Find recapture rate where (observed-1) deaths round UP to efficacy threshold
        critical_recaptured <- (observed_deaths * sample_size) / (observed_deaths + 0.5)
        critical_rate <- critical_recaptured / sample_size
        if(critical_rate <= 1.0) round(critical_rate * 100, 1) else NA
      } else NA,
      
      # FIXED RECAPTURE RATE ANALYSIS
      # Futility at fixed recapture rate
      futility_fixed_effect = if(!is.na(deaths_for_futility)) {
        observed_deaths <- deaths_for_futility
        recaptured <- round(sample_size * fixed_recapture / 100)
        ht_estimate <- (observed_deaths / recaptured) * sample_size
        estimated_deaths <- round(ht_estimate)
        estimated_deaths - observed_deaths
      } else NA,
      
      # Efficacy at fixed recapture rate
      efficacy_fixed_effect = if(!is.na(deaths_for_efficacy)) {
        observed_deaths <- deaths_for_efficacy - 1
        recaptured <- round(sample_size * fixed_recapture / 100)
        ht_estimate <- (observed_deaths / recaptured) * sample_size
        estimated_deaths <- round(ht_estimate)
        estimated_deaths - observed_deaths
      } else NA,
      
      futility_prevention_text = if(!is.na(futility_prevention_recapture)) {
        paste0("Recapture <", futility_prevention_recapture, "% prevents futility stopping (", 
               deaths_for_futility, "→", deaths_for_futility + 1, " deaths)")
      } else "No futility boundary",
      
      efficacy_facilitation_text = if(!is.na(efficacy_facilitation_recapture)) {
        paste0("Recapture <", efficacy_facilitation_recapture, "% enables efficacy stopping (", 
               deaths_for_efficacy - 1, "→", deaths_for_efficacy, " deaths)")
      } else "High recapture needed",
      
      # Fixed recapture rate text
      futility_prevention_text_fixed = if(!is.na(deaths_for_futility)) {
        if(futility_fixed_effect > 0) {
          paste0("At ", fixed_recapture, "% recapture: prevents futility stopping (", 
                 deaths_for_futility, "→", deaths_for_futility + futility_fixed_effect, " deaths)")
        } else {
          paste0("At ", fixed_recapture, "% recapture: no effect on futility")
        }
      } else "No futility boundary",
      
      efficacy_facilitation_text_fixed = if(!is.na(deaths_for_efficacy)) {
        if(efficacy_fixed_effect > 0) {
          paste0("At ", fixed_recapture, "% recapture: enables efficacy stopping (", 
                 deaths_for_efficacy - 1, "→", deaths_for_efficacy - 1 + efficacy_fixed_effect, " deaths)")
        } else {
          paste0("At ", fixed_recapture, "% recapture: no effect on efficacy")
        }
      } else "No efficacy boundary"
    ) %>%
    select(analysis, sample_size, deaths_for_futility, deaths_for_efficacy,
           futility_prevention_recapture, efficacy_facilitation_recapture,
           futility_prevention_text, efficacy_facilitation_text,
           futility_prevention_text_fixed, efficacy_facilitation_text_fixed)
  
  return(results)
}

extract_gsdesign_boundaries <- function(gsdesign_obj) {
  # Extract cumulative sample sizes
  sample_sizes <- ceiling(gsdesign_obj$n.I)
  
  # Extract boundaries
  upper_bounds <- gsdesign_obj$upper$bound
  
  # Check if lower bounds exist (for two-sided tests)
  if (!is.null(gsdesign_obj$lower)) {
    lower_bounds <- gsdesign_obj$lower$bound
  } else {
    lower_bounds <- rep(-Inf, length(upper_bounds))  # No futility stopping
  }
  
  # Create boundary data frame
  boundaries_df <- data.frame(
    analysis = 1:length(sample_sizes),
    sample_size = sample_sizes,
    upper_bound = upper_bounds,
    lower_bound = lower_bounds
  )
  
  return(boundaries_df)
}

p_to_z_calculator <- function(p_value) {
  if(p_value <= 0 | p_value >= 1) {
    stop("P-value must be between 0 and 1")
  }
  z_score <- qnorm(p_value, lower.tail = FALSE)
  
  cat("P-value:", p_value, "\n")
  cat("Z-score:", round(z_score, 3), "\n")
  
  return(z_score)
}

#  Summary functions ####

# Determine acceptance and rejection rates for given sample size
# e.g., when do we stop?
cumulative_stopping <- function(simulation_data) {
  
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


cumulative_stopping_sum <- function(simulation_data) {
  
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

decision_plot <- function(mortality_table) {
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
             label = "Reject: Non-compliant", size = 4) +
    annotate("text", x = mean(x_limits), y = mean(y_crop_limits), 
             label = "Precision not achieved: Continue testing", size = 4) +
    annotate("text", x = x_limits[2] * 0.75, y = y_crop_limits[2] * 0.15, 
             label = "Accept: Compliant", size = 4) +
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
        outcome == "prop_accept_by_now" ~ "Accept",
        outcome == "prop_reject_by_now" ~ "Reject", 
        outcome == "prop_continue" ~ "Continue"
      ),
      percentage = proportion * 100,
      decision_method = factor(
        decision_method,
        levels = c("aggressive", "neutral", "conservative")
      ),
      true_survival_rate = factor(
        true_survival_rate,
        levels = sort(unique(filtered_data$true_survival_rate))
      ),
      recapture_rate = factor(
        recapture_rate,
        levels = sort(unique(filtered_data$recapture_rate))
      )
    )
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = interval, y = percentage, color = outcome)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = FALSE, span = 0.7) +
    scale_color_manual(values = c("Accept" = "darkgreen", 
                                  "Reject" = "darkred", 
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
        decision_method = function(x) paste("H-T Method:", str_to_title(x))
      )
    )
  } else if (n_survival > 1 && n_recapture > 1) {
    # Survival and recapture vary
    p <- p + facet_grid(recapture_rate ~ true_survival_rate, labeller = label_both)
  } else if (n_survival > 1 && n_method > 1) {
    # Survival and method vary
    p <- p + facet_grid(decision_method ~ true_survival_rate, 
                        labeller = labeller(
                          decision_method = function(x) paste("H-T Method:", str_to_title(x)),
                          true_survival_rate = function(x) paste("Survival:", x)
                        ))
  } else if (n_recapture > 1 && n_method > 1) {
    # Recapture and method vary
    p <- p + facet_grid(decision_method ~ recapture_rate,
                        labeller = labeller(
                          decision_method = function(x) paste("H-T Method:", str_to_title(x)),
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
                        labeller = labeller(decision_method = function(x) paste("H-T Method:", str_to_title(x))))
  }
  
  return(p)
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
      true_survival_rate = as.numeric(true_survival_rate),
      facet_survival = factor(
        paste0("Survival: ", true_survival_rate * 100, "%"),
        levels = paste0("Survival: ", sort(unique(true_survival_rate)) * 100, "%")
      ),
      facet_recapture = factor(
        paste0("Recapture: ", recapture_rate * 100, "%"),
        levels = paste0("Recapture: ", sort(unique(recapture_rate)) * 100, "%")
      ),
      facet_method = factor(
        paste0("H-T Method: ", str_to_title(decision_method)),
        levels = paste0("H-T Method: ", str_to_title(c("aggressive", "neutral", "conservative")))
      )
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
  base_plot <- decision_plot(mortality_table)
  
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
      rows = vars(facet_recapture),
      cols = vars(facet_survival, facet_method),
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

plot_custom_gsdesign <- function(gsdesign_obj, title = NULL) {
  # Extract boundary data
  boundaries <- extract_gsdesign_boundaries(gsdesign_obj)
  
  # Set plot limits
  x_max <- max(boundaries$sample_size) * 1.1
  y_min <- min(boundaries$lower_bound) - 1
  y_max <- max(boundaries$upper_bound) + 0.5
  
  # Create extended sample size sequence for smooth lines starting from first boundary
  x_start <- min(boundaries$sample_size)
  x_seq <- seq(x_start, x_max, length.out = 200)
  
  # Interpolate boundaries for smooth lines
  upper_smooth <- approx(boundaries$sample_size, boundaries$upper_bound, 
                         xout = x_seq, rule = 2)$y
  lower_smooth <- approx(boundaries$sample_size, boundaries$lower_bound, 
                         xout = x_seq, rule = 2)$y
  
  # Create extended line data that includes extreme points at first x-value
  line_data_upper <- data.frame(
    x = c(x_start, x_seq),
    y = c(y_max, upper_smooth)  # Add extreme high point at first x
  )
  
  line_data_lower <- data.frame(
    x = c(x_start, x_seq), 
    y = c(y_min, lower_smooth)  # Add extreme low point at first x
  )
  
  ribbon_data <- data.frame(
    x = x_seq,
    upper = upper_smooth,
    lower = lower_smooth
  )
  
  p <- ggplot() +
    # Green continue area extends from x=0 to first boundary
    geom_ribbon(aes(x = c(0, x_start), ymin = y_min, ymax = y_max), 
                fill = "lightgreen", alpha = 0.6) +
    # Boundary-defined regions start from first analysis point
    geom_ribbon(data = ribbon_data, 
                aes(x = x, ymin = upper, ymax = y_max), 
                fill = "lightblue", alpha = 0.6) +
    geom_ribbon(data = ribbon_data, 
                aes(x = x, ymin = lower, ymax = upper), 
                fill = "lightgreen", alpha = 0.6) +
    geom_ribbon(data = ribbon_data, 
                aes(x = x, ymin = y_min, ymax = lower), 
                fill = "lightcoral", alpha = 0.6) +
    geom_vline(data = boundaries, 
               aes(xintercept = sample_size), 
               linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_line(data = line_data_upper, aes(x = x, y = y), 
              color = "darkblue", size = 0.8) +
    geom_line(data = line_data_lower, aes(x = x, y = y), 
              color = "darkred", size = 0.8) +
    geom_point(data = boundaries, 
               aes(x = sample_size, y = upper_bound), 
               color = "darkblue", size = 2) +
    geom_point(data = boundaries, 
               aes(x = sample_size, y = lower_bound), 
               color = "darkred", size = 2) +
    # Add z-score labels above upper boundary points
    geom_text(data = boundaries,
              aes(x = sample_size, y = upper_bound, 
                  label = paste(round(upper_bound, 2))),
              vjust = -0.5, hjust = -0.6, size = 3) +
    # Add z-score labels below lower boundary points  
    geom_text(data = boundaries,
              aes(x = sample_size, y = lower_bound,
                  label = paste(round(lower_bound, 2))),
              vjust = 0.6, hjust = -0.6, size = 3) +
    # Add sample size labels below the plot area
    geom_text(data = boundaries,
              aes(x = sample_size, y = y_min,
                  label = paste("n =", sample_size)),
              vjust = -0.4, size = 3, hjust = -0.3, color = "black", angle = 45) +
    geom_point(aes(x = gsdesign_obj$n.fix, y = 1.96), 
               shape = 18, size = 3, color = "black") +
    scale_x_continuous(limits = c(0, x_max),
                       breaks = seq(0, x_max, by = 50),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = seq(floor(y_min), ceiling(y_max), by = 0.5),
                       expand = c(0, 0)) +
    labs(
      x = "Total Sample Size (3:1 impact:control)",
      y = "Normal Critical Value (Z)"
    ) +
    annotate(
      "text",
      x = x_max * 0.98, 
      y = y_max * 0.98, 
      label = "Parameters:\nLan-DeMets O'Brien-Fleming\nlower bounds non-binding\nPower (1-β) = 80%\nα = 0.025 (one sided)",
      hjust = 1,
      vjust = 1,
      size = 4,
      color = "black"
    ) +
    
    # Apply your custom theme
    theme_JN() +
    theme(
      panel.grid.major = element_line(colour = "grey90", size = 0.5),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

p_to_z_visualizer <- function(target_z = NULL, p_range = c(0.0000001, 1)) {
  
  # Generate a sequence of p-values
  p_values <- seq(p_range[1], p_range[2], length.out = 1000)
  
  # Calculate corresponding z-scores
  z_scores <- qnorm(p_values, lower.tail = FALSE)
  
  # Create the plot
  plot(p_values, z_scores, type = "l", lwd = 2, col = "blue",
       xlab = "P-value", ylab = "Z-score",
       main = "Relationship between P-values and Z-scores",
       xlim = c(0, max(p_values)), ylim = c(-2, max(z_scores)))
  
  grid()
  
  # Add reference lines for common significance levels
  abline(h = c(1.645, 1.96, 2.576), col = "gray", lty = 2)
  text(x = max(p_values) * 0.9, y = 1.645, "90% (1.645)", pos = 3, cex = 0.8)
  text(x = max(p_values) * 0.9, y = 1.96, "95% (1.96)", pos = 3, cex = 0.8)
  text(x = max(p_values) * 0.9, y = 2.576, "99% (2.576)", pos = 3, cex = 0.8)
  
  # If a target z-score is specified, find and mark the corresponding p-value
  if (!is.null(target_z)) {
    # Calculate p-value for target z-score
    target_p <- pnorm(target_z, lower.tail = FALSE)
    
    # Add point and lines
    points(target_p, target_z, col = "red", pch = 19, cex = 1.5)
    abline(h = target_z, col = "red", lty = 2)
    abline(v = target_p, col = "red", lty = 2)
    
    # Add text annotation
    text(x = target_p, y = target_z, 
         labels = paste0("Z = ", target_z, "\nP = ", round(target_p, 6)),
         pos = 4, col = "red", cex = 0.9)
    
    cat("For a Z-score of", target_z, "the corresponding p-value is:", target_p, "\n")
    
    return(target_p)
  }
}

# Lookup functions ####

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


# Core lookup function with exact original logic
lookup_mortality_scenario <- function(mortality_table, 
                                      sample_size = NULL, 
                                      n_deaths = NULL, 
                                      pct_difference = NULL,
                                      decision_type = NULL,
                                      current_sample = NULL) {
  
  result <- mortality_table
  
  # Handle different lookup types with original logic
  if (!is.null(pct_difference) && is.null(sample_size) && is.null(n_deaths)) {
    # Replacement for lookup_min_sample_pct_diff
    target_survival <- 100 - pct_difference
    result <- result %>%
      filter(decision != "Continue") %>%  # Original filters out "Continue"
      mutate(survival_diff = abs(survival - target_survival)) %>%
      arrange(survival_diff, sample_size)
  }
  else if (!is.null(n_deaths) && !is.null(decision_type) && is.null(current_sample)) {
    # Replacement for lookup_min_sample_fixed_deaths
    result <- result %>%
      filter(n_mortality == n_deaths, decision == decision_type) %>%
      arrange(sample_size)
  }
  else if (!is.null(n_deaths) && !is.null(current_sample)) {
    # Replacement for lookup_next_stopping_point
    result <- result %>%
      filter(sample_size > current_sample, n_mortality == n_deaths, decision == decision_type) %>%
      arrange(sample_size)
  }
  else if (!is.null(sample_size) && !is.null(pct_difference)) {
    # Replacement for lookup_ci_at_sample_size
    mortality_rate <- pct_difference / 100
    expected_deaths <- round(sample_size * mortality_rate)
    result <- result %>%
      filter(sample_size == !!sample_size, n_mortality == !!expected_deaths)
  }
  else if (!is.null(sample_size) && !is.null(decision_type)) {
    # For rejection thresholds
    result <- result %>%
      filter(sample_size == sample_size, decision == decision_type) %>%
      arrange(n_mortality)
  }
  
  # Return appropriate result
  result <- result %>%
    select(sample_size, n_mortality, survival, ci_lower_surv, ci_upper_surv, 
           decision, survival_formatted)
  
  if(nrow(result) == 0) {
    return(NULL)
  }
  
  return(result %>% slice_head(n = 1))
}

# Exact replacement for lookup_min_sample_pct_diff
lookup_min_sample_pct_diff <- function(mortality_table, pct_difference) {
  lookup_mortality_scenario(mortality_table, pct_difference = pct_difference)
}

# Exact replacement for lookup_min_sample_fixed_deaths  
lookup_min_sample_fixed_deaths <- function(mortality_table, current_deaths) {
  lookup_mortality_scenario(mortality_table, n_deaths = current_deaths, decision_type = "Accept H0")
}

# Exact replacement for lookup_next_stopping_point
lookup_next_stopping_point <- function(mortality_table, current_sample, current_deaths) {
  lookup_mortality_scenario(mortality_table, n_deaths = current_deaths, decision_type = "Accept H0", current_sample = current_sample)
}

# Exact replacement for lookup_ci_at_sample_size
lookup_ci_at_sample_size <- function(mortality_table, sample_size, pct_difference) {
  lookup_mortality_scenario(mortality_table, sample_size = sample_size, pct_difference = pct_difference)
}

# Exact replacement for lookup_rejection_threshold
lookup_rejection_threshold <- function(mortality_table, current_sample, current_deaths, next_interval) {
  reject_threshold <- lookup_mortality_scenario(mortality_table, sample_size = next_interval, decision_type = "Reject H0")
  
  if(is.null(reject_threshold)) {
    cat("No rejection threshold found at", next_interval, "fish\n")
    return(NULL)
  }
  
  additional_deaths_needed <- reject_threshold$n_mortality - current_deaths
  
  tibble(
    current_sample = current_sample,
    current_deaths = current_deaths,
    next_interval = next_interval,
    rejection_threshold_deaths = reject_threshold$n_mortality,
    additional_deaths_to_reject = additional_deaths_needed,
    rejection_survival = reject_threshold$survival,
    rejection_ci_upper = reject_threshold$ci_upper_surv
  )
}

# Keep original show_continuation_scenarios as it's complex
show_continuation_scenarios <- function(mortality_table, intervals = c(100, 200, 300, 400)) {
  
  scenarios <- tibble()
  
  for(current_interval in intervals[-length(intervals)]) {
    next_interval <- intervals[which(intervals == current_interval) + 1]
    
    continue_cases <- mortality_table %>%
      filter(sample_size == current_interval, decision == "Continue") %>%
      filter(n_mortality <= 15) %>%
      select(sample_size, n_mortality, survival, decision)
    
    for(i in seq_len(nrow(continue_cases))) {
      current_deaths <- continue_cases$n_mortality[i]
      current_survival <- continue_cases$survival[i]
      
      next_accept <- lookup_next_stopping_point(mortality_table, current_interval, current_deaths)
      
      next_reject_threshold <- mortality_table %>%
        filter(sample_size == next_interval, decision == "Reject H0") %>%
        slice_head(n = 1)
      
      scenario <- tibble(
        current_interval = current_interval,
        current_deaths = current_deaths,
        current_survival = current_survival,
        accept_at_n = if(!is.null(next_accept)) next_accept$sample_size[1] else NA,
        accept_survival = if(!is.null(next_accept)) next_accept$survival[1] else NA,
        next_interval = next_interval,
        reject_threshold_deaths = if(nrow(next_reject_threshold) > 0) next_reject_threshold$n_mortality[1] else NA,
        additional_deaths_to_reject = if(nrow(next_reject_threshold) > 0) 
          next_reject_threshold$n_mortality[1] - current_deaths else NA
      )
      
      scenarios <- scenarios %>% bind_rows(scenario)
    }
  }
  
  return(scenarios)
}

# Keep find_decision_thresholds as is (it works well)
find_decision_thresholds <- function(mortality_table, intervals = c(100, 200, 300, 400)) {
  
  thresholds <- map_dfr(intervals, function(n) {
    
    interval_data <- mortality_table %>%
      filter(sample_size == n) %>%
      mutate(mortality_pct = round((n_mortality / sample_size) * 100, 1))
    
    last_accept <- interval_data %>% filter(decision == "Accept H0") %>% slice_tail(n = 1)
    first_reject <- interval_data %>% filter(decision == "Reject H0") %>% slice_head(n = 1)
    continue_range <- interval_data %>% filter(decision == "Continue") %>% 
      summarise(min_continue_pct = min(mortality_pct), max_continue_pct = max(mortality_pct))
    
    tibble(
      interval = n,
      accept_threshold = if(nrow(last_accept) > 0) last_accept$mortality_pct else NA,
      accept_deaths = if(nrow(last_accept) > 0) last_accept$n_mortality else NA,
      reject_threshold = if(nrow(first_reject) > 0) first_reject$mortality_pct else NA,
      reject_deaths = if(nrow(first_reject) > 0) first_reject$n_mortality else NA,
      continue_min = if(!is.infinite(continue_range$min_continue_pct)) continue_range$min_continue_pct else NA,
      continue_max = if(!is.infinite(continue_range$max_continue_pct)) continue_range$max_continue_pct else NA
    )
  })
  
  return(thresholds)
}

# Exact replacement for build_stopping_design_table
build_stopping_design_table <- function(mortality_table, intervals = c(100, 200, 300, 400)) {
  
  map_dfr(intervals, function(n) {
    
    mortality_table %>%
      filter(sample_size == n) %>%
      mutate(
        mortality_pct = round((n_mortality / sample_size) * 100, 1),
        interval = n
      ) %>%
      select(interval, n_mortality, mortality_pct, survival, ci_lower_surv, 
             ci_upper_surv, decision, survival_formatted)
  })
}

# Exact replacement for build_comprehensive_stopping_table
build_comprehensive_stopping_table <- function(mortality_table, intervals = c(100, 200, 300, 400)) {
  
  comprehensive_table <- tibble()
  
  for(i in seq_along(intervals)) {
    current_interval <- intervals[i]
    next_interval <- if(i < length(intervals)) intervals[i + 1] else NA
    
    interval_data <- mortality_table %>%
      filter(sample_size == current_interval) %>%
      mutate(mortality_pct = round((n_mortality / sample_size) * 100, 1))
    
    last_accept <- interval_data %>% filter(decision == "Accept H0") %>% slice_tail(n = 1)
    first_reject <- interval_data %>% filter(decision == "Reject H0") %>% slice_head(n = 1)
    continue_cases <- interval_data %>% filter(decision == "Continue") %>% filter(n_mortality <= 15)
    
    base_row <- tibble(
      interval = current_interval,
      accept_threshold_pct = if(nrow(last_accept) > 0) last_accept$mortality_pct else NA,
      accept_max_deaths = if(nrow(last_accept) > 0) last_accept$n_mortality else NA,
      reject_threshold_pct = if(nrow(first_reject) > 0) first_reject$mortality_pct else NA,
      reject_min_deaths = if(nrow(first_reject) > 0) first_reject$n_mortality else NA,
      continue_min_deaths = if(nrow(continue_cases) > 0) min(continue_cases$n_mortality) else NA,
      continue_max_deaths = if(nrow(continue_cases) > 0) max(continue_cases$n_mortality) else NA,
      next_interval = next_interval
    )
    
    if(nrow(continue_cases) > 0 && !is.na(next_interval)) {
      
      if(nrow(continue_cases) <= 3) {
        sample_cases <- continue_cases
      } else {
        sample_cases <- continue_cases %>%
          slice(c(1, ceiling(nrow(continue_cases)/2), nrow(continue_cases)))
      }
      
      for(j in seq_len(nrow(sample_cases))) {
        current_deaths <- sample_cases$n_mortality[j]
        current_survival <- sample_cases$survival[j]
        
        accept_point <- lookup_next_stopping_point(mortality_table, current_interval, current_deaths)
        reject_at_next <- mortality_table %>% filter(sample_size == next_interval, decision == "Reject H0") %>% slice_head(n = 1)
        
        scenario_row <- base_row %>%
          mutate(
            scenario_current_deaths = current_deaths,
            scenario_current_survival = current_survival,
            accept_if_no_more_deaths_at_n = if(!is.null(accept_point)) accept_point$sample_size[1] else NA,
            reject_if_deaths_reach = if(nrow(reject_at_next) > 0) reject_at_next$n_mortality[1] else NA,
            additional_deaths_to_reject = if(nrow(reject_at_next) > 0) reject_at_next$n_mortality[1] - current_deaths else NA
          )
        
        comprehensive_table <- bind_rows(comprehensive_table, scenario_row)
      }
    } else {
      comprehensive_table <- bind_rows(comprehensive_table, base_row)
    }
  }
  
  return(comprehensive_table)
}

# Function to analyze real study data
lookup_seq_95 <- function(exposed_fish, recaptured_fish, observed_deaths, 
                                             mortality_table, cumulative_data, 
                                             assumed_true_survival = 0.98,
                                             recapture_rate_scenario = 1.0) {
  
  # Calculate basic rates
  recapture_rate <- recaptured_fish / exposed_fish
  observed_survivors <- recaptured_fish - observed_deaths
  survival_rate_recaptured <- observed_survivors / recaptured_fish
  
  # Calculate HT estimates using your existing function
  ht_result <- horvitz_thompson_est(exposed_fish, recaptured_fish, observed_deaths)
  
  # Three decision scenarios
  scenarios <- list(
    aggressive = ht_result$deaths_ci_lower,
    neutral = ht_result$deaths_point,
    conservative = ht_result$deaths_ci_upper
  )
  
  # Look up decisions for each scenario
  decisions <- list()
  
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
      cat("  Probability of Accept:", round(current_interval_data$prop_accept_by_now * 100, 1), "%\n")
      cat("  Probability of Reject:", round(current_interval_data$prop_reject_by_now * 100, 1), "%\n")
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






###non-inferiority testing###
#
# not sure about this approach as it only considers the variability in the impact group, and not the control group
# This would be suitable if we had a known survival probability of one pump, and we wanted to ensure a new pump was not inferior to the existing one.
# Reserve code for future potential pump comparisons 
# 
# control_survivors <- 64
# control_deaths <- 17
# control_total <- 81
# control_prop <- control_survivors / control_total
# 
# # Non-inferiority margin: 2% (pump acceptable if ≥97%)
# ni_margin <- 0.02
# ni_threshold <- control_prop - ni_margin  
# 
# # Function to perform non-inferiority test
# perform_ni_test <- function(pump_survivors, pump_total, ni_threshold) {
#   # Test H0: pump_prop < ni_threshold vs H1: pump_prop >= ni_threshold
#   result <- binom.test(pump_survivors, pump_total, p = ni_threshold, 
#                        alternative = "greater")
#   return(result)
# }
# 
# # Scenario 1: Poor pump performance (95% survival)
# pump_survivors_1 <- 82
# pump_total_1 <- 105
# result_1 <- perform_ni_test(pump_survivors_1, pump_total_1, ni_threshold)
# 
# result_1
