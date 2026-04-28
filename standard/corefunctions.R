# Simple function to test Hardy-Weinberg Equilibrium
test_hwe_simple <- function(genotype_array, alpha = 0.05) {
  # genotype_array: 3D array with dimensions [individuals, alleles(2), loci]
  
  if (length(dim(genotype_array)) != 3) {
    stop("Input must be a 3D array with dimensions [individuals, alleles, loci]")
  }
  
  num_individuals <- dim(genotype_array)[1]
  num_loci <- dim(genotype_array)[3]
  
  # Storage for results
  hwe_results <- data.frame(
    locus = 1:num_loci,
    num_alleles = NA,
    observed_het = NA,
    expected_het = NA,
    chi_square = NA,
    p_value = NA,
    in_hwe = NA,
    hwe_testable = FALSE  # Can we actually test HWE?
  )
  
  for (locus in 1:num_loci) {
    # Extract genotypes for this locus
    locus_genotypes <- genotype_array[, , locus]
    
    # Get unique alleles
    all_alleles <- unique(c(locus_genotypes[, 1], locus_genotypes[, 2]))
    num_alleles <- length(all_alleles)
    hwe_results$num_alleles[locus] <- num_alleles
    
    # Skip monomorphic loci (can't test HWE)
    if (num_alleles <= 1) {
      hwe_results$observed_het[locus] <- 0
      hwe_results$expected_het[locus] <- 0
      hwe_results$in_hwe[locus] <- TRUE  # Trivially in HWE
      next
    }
    
    # Calculate observed heterozygosity
    het_count <- sum(locus_genotypes[, 1] != locus_genotypes[, 2])
    obs_het <- het_count / num_individuals
    hwe_results$observed_het[locus] <- obs_het
    
    # Calculate allele frequencies
    total_alleles <- 2 * num_individuals
    allele_freqs <- sapply(all_alleles, function(a) {
      sum(locus_genotypes == a) / total_alleles
    })
    names(allele_freqs) <- all_alleles
    
    # Expected heterozygosity under HWE
    exp_het <- 1 - sum(allele_freqs^2)
    hwe_results$expected_het[locus] <- exp_het
    
    # Count observed genotypes
    genotype_counts <- count_genotypes(locus_genotypes, all_alleles)
    
    # Calculate expected genotype counts under HWE
    expected_counts <- calculate_expected_genotype_counts(allele_freqs, num_individuals)
    
    # Chi-square test (only if expected counts are sufficient)
    if (all(expected_counts >= 5)) {
      hwe_results$hwe_testable[locus] <- TRUE
      
      # Chi-square statistic
      chi_sq <- sum((genotype_counts - expected_counts)^2 / expected_counts)
      df <- length(genotype_counts) - length(allele_freqs)  # degrees of freedom
      
      if (df > 0) {
        p_val <- 1 - pchisq(chi_sq, df)
        hwe_results$chi_square[locus] <- chi_sq
        hwe_results$p_value[locus] <- p_val
        hwe_results$in_hwe[locus] <- p_val > alpha
      }
    }
  }
  
  # Summary statistics
  testable_loci <- sum(hwe_results$hwe_testable, na.rm = TRUE)
  loci_in_hwe <- sum(hwe_results$in_hwe, na.rm = TRUE)
  polymorphic_loci <- sum(hwe_results$num_alleles > 1, na.rm = TRUE)
  
  # Mean heterozygosity (only polymorphic loci)
  polymorphic_mask <- hwe_results$num_alleles > 1
  mean_obs_het <- ifelse(sum(polymorphic_mask) > 0, 
                         mean(hwe_results$observed_het[polymorphic_mask], na.rm = TRUE), 
                         0)
  mean_exp_het <- ifelse(sum(polymorphic_mask) > 0, 
                         mean(hwe_results$expected_het[polymorphic_mask], na.rm = TRUE), 
                         0)
  
  # Inbreeding coefficient (F_IS)
  f_is <- ifelse(mean_exp_het > 0, (mean_exp_het - mean_obs_het) / mean_exp_het, 0)
  
  summary <- list(
    total_loci = num_loci,
    polymorphic_loci = polymorphic_loci,
    testable_loci = testable_loci,
    loci_in_hwe = loci_in_hwe,
    proportion_in_hwe = ifelse(testable_loci > 0, loci_in_hwe / testable_loci, NA),
    mean_observed_het = mean_obs_het,
    mean_expected_het = mean_exp_het,
    inbreeding_coefficient = f_is
  )
  
  return(list(
    summary = summary,
    by_locus = hwe_results
  ))
}

# Helper function to count genotypes
count_genotypes <- function(locus_genotypes, all_alleles) {
  num_alleles <- length(all_alleles)
  num_individuals <- nrow(locus_genotypes)
  
  # Count each genotype combination
  counts <- numeric()
  names_vec <- character()
  
  # Count homozygotes
  for (i in 1:num_alleles) {
    allele <- all_alleles[i]
    count <- sum(locus_genotypes[, 1] == allele & locus_genotypes[, 2] == allele)
    counts <- c(counts, count)
    names_vec <- c(names_vec, paste0(allele, allele))
  }
  
  # Count heterozygotes (only unique combinations)
  if (num_alleles > 1) {
    for (i in 1:(num_alleles-1)) {
      for (j in (i+1):num_alleles) {
        allele1 <- all_alleles[i]
        allele2 <- all_alleles[j]
        count <- sum((locus_genotypes[, 1] == allele1 & locus_genotypes[, 2] == allele2) |
                       (locus_genotypes[, 1] == allele2 & locus_genotypes[, 2] == allele1))
        counts <- c(counts, count)
        names_vec <- c(names_vec, paste0(allele1, allele2))
      }
    }
  }
  
  names(counts) <- names_vec
  return(counts)
}

# Helper function to calculate expected genotype counts under HWE
calculate_expected_genotype_counts <- function(allele_freqs, num_individuals) {
  allele_names <- names(allele_freqs)
  num_alleles <- length(allele_freqs)
  
  expected <- numeric()
  names_vec <- character()
  
  # Expected homozygote counts: N * p^2
  for (i in 1:num_alleles) {
    expected <- c(expected, num_individuals * allele_freqs[i]^2)
    names_vec <- c(names_vec, paste0(allele_names[i], allele_names[i]))
  }
  
  # Expected heterozygote counts: N * 2pq
  if (num_alleles > 1) {
    for (i in 1:(num_alleles-1)) {
      for (j in (i+1):num_alleles) {
        expected <- c(expected, num_individuals * 2 * allele_freqs[i] * allele_freqs[j])
        names_vec <- c(names_vec, paste0(allele_names[i], allele_names[j]))
      }
    }
  }
  
  names(expected) <- names_vec
  return(expected)
}

# Function to create a summary plot
plot_hwe_simple <- function(hwe_result) {
  library(ggplot2)
  
  # Only plot testable loci
  plot_data <- hwe_result$by_locus[hwe_result$by_locus$hwe_testable == TRUE, ]
  
  if (nrow(plot_data) == 0) {
    message("No testable loci to plot")
    return(NULL)
  }
  
  # Create plots
  p1 <- ggplot(plot_data, aes(x = factor(locus), y = -log10(p_value))) +
    geom_point(aes(color = in_hwe), size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    labs(title = "Hardy-Weinberg P-values by Locus",
         x = "Locus", y = "-log10(p-value)",
         color = "In HWE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(plot_data, aes(x = expected_het, y = observed_het)) +
    geom_point(aes(color = in_hwe), size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Observed vs Expected Heterozygosity",
         x = "Expected Heterozygosity", y = "Observed Heterozygosity",
         color = "In HWE") +
    theme_minimal()
  
  return(list(p_values = p1, heterozygosity = p2))
}

# Example usage and testing
test_hwe_functions <- function() {
  # Load your data
  temp <- read_rds("./results/rep_1_simulation_iteration_1.RDS")
  gen_data <- temp$X0_Y0.9_EnvCont0.01_HostSelInf_MicrobeSelInf$GenData[[1]]
  
  cat("Testing HWE functions:\n")
  cat("======================\n\n")
  
  # Test the simple HWE function
  hwe_result <- test_hwe_simple(gen_data$host_microbe_optima)
  
  # Test the individual heterozygosity functions
  obs_het <- calculate_observed_heterozygosity(gen_data$host_microbe_optima)
  exp_het <- calculate_expected_heterozygosity(gen_data$host_microbe_optima)
  
  # Compare with your built-in calculation
  builtin_het <- mean(gen_data$observed_heterozygosity_optima, na.rm = TRUE)
  
  cat("Results comparison:\n")
  cat("Built-in heterozygosity:    ", sprintf("%.6f", builtin_het), "\n")
  cat("Simple obs function:        ", sprintf("%.6f", obs_het$mean_heterozygosity), "\n")
  cat("HWE function (observed):    ", sprintf("%.6f", hwe_result$summary$mean_observed_het), "\n")
  cat("HWE function (expected):    ", sprintf("%.6f", hwe_result$summary$mean_expected_het), "\n\n")
  
  cat("HWE Summary:\n")
  cat("Total loci:         ", hwe_result$summary$total_loci, "\n")
  cat("Polymorphic loci:   ", hwe_result$summary$polymorphic_loci, "\n")
  cat("Testable loci:      ", hwe_result$summary$testable_loci, "\n")
  cat("Loci in HWE:       ", hwe_result$summary$loci_in_hwe, "\n")
  cat("Proportion in HWE:  ", sprintf("%.3f", hwe_result$summary$proportion_in_hwe), "\n")
  cat("Inbreeding coeff:   ", sprintf("%.3f", hwe_result$summary$inbreeding_coefficient), "\n")
  
  return(hwe_result)
}


# Create a data structure to store lineage information for diploid organisms
create_lineage_tracker <- function(initial_population_size, max_generations) {
  # Create a matrix to store parent-child relationships
  # Each row represents [generation, child_id, parent1_id, parent2_id, allele1_from_parent1, allele2_from_parent2]
  lineage_matrix <- matrix(NA, 
                           nrow = initial_population_size * max_generations,
                           ncol = 4,
                           dimnames = list(NULL, c("generation", "child_id", "parent1_id", "parent2_id")))
  
  # Initialize first generation with no parents
  lineage_matrix[1:initial_population_size, ] <- cbind(
    generation = 1,
    child_id = 1:initial_population_size,
    parent1_id = NA,
    parent2_id = NA
  )
  
  return(lineage_matrix)
}

# Update lineage information each generation for diploid organisms
update_lineage_tracker <- function(lineage_matrix, generation, parent_pairs) {
  n_offspring <- nrow(parent_pairs)
  
  # Find where to start adding new rows
  current_row <- which(is.na(lineage_matrix[,"generation"]))[1]
  if(is.na(current_row)) {
    # If matrix is full, extend it
    old_matrix <- lineage_matrix
    lineage_matrix <- matrix(NA, 
                             nrow = nrow(old_matrix) + n_offspring,
                             ncol = ncol(old_matrix),
                             dimnames = list(NULL, colnames(old_matrix)))
    lineage_matrix[1:nrow(old_matrix),] <- old_matrix
    current_row <- nrow(old_matrix) + 1
  }
  
  new_rows <- current_row:(current_row + n_offspring - 1)
  
  # Add new generation's lineage information
  new_data <- cbind(
    generation = generation,
    child_id = 1:n_offspring,
    # Extract parent IDs from the parent_pairs matrix
    parent1_id = as.numeric(gsub("Host_", "", parent_pairs[,1])),
    parent2_id = as.numeric(gsub("Host_", "", parent_pairs[,2]))
  )
  
  lineage_matrix[new_rows, ] <- new_data
  
  return(lineage_matrix)
}

# Function to convert lineage matrix to plotting format
prepare_lineage_plot_data <- function(lineage_matrix, target_lineages = NULL) {
  library(tidyverse)
  
  # Convert matrix to dataframe
  lineage_df <- as.data.frame(lineage_matrix)
  
  # If specific lineages aren't specified, track all final generation individuals
  if(is.null(target_lineages)) {
    max_gen <- max(lineage_df$generation, na.rm = TRUE)
    target_lineages <- lineage_df$child_id[lineage_df$generation == max_gen]
  }
  
  # Track ancestry recursively
  get_ancestry <- function(id, gen) {
    if(gen == 1 || is.na(id)) return(NULL)
    
    current <- lineage_df[lineage_df$generation == gen & lineage_df$child_id == id,]
    if(nrow(current) == 0) return(NULL)
    
    # Get parent connections
    p1_connections <- get_ancestry(current$parent1_id, gen - 1)
    p2_connections <- get_ancestry(current$parent2_id, gen - 1)
    
    # Create current connections
    current_connections <- data.frame(
      from_gen = gen - 1,
      to_gen = gen,
      from_id = c(current$parent1_id, current$parent2_id),
      to_id = rep(current$child_id, 2)
    )
    
    return(bind_rows(current_connections, p1_connections, p2_connections))
  }
  
  # Get all connections for target lineages
  all_connections <- map_df(target_lineages, function(id) {
    get_ancestry(id, max(lineage_df$generation, na.rm = TRUE))
  }) %>% distinct()
  
  return(all_connections)
}

# Function to create the actual plot using ggplot2
plot_lineages <- function(connections_df, highlight_lineage = NULL) {
  library(ggplot2)
  
  # Create unique IDs for each individual at each generation
  plot_data <- connections_df %>%
    mutate(
      from_unique = paste(from_gen, from_id),
      to_unique = paste(to_gen, to_id)
    )
  
  # If a specific lineage is highlighted, add this information
  if(!is.null(highlight_lineage)) {
    plot_data <- plot_data %>%
      mutate(highlight = to_id %in% highlight_lineage)
  }
  
  # Create the plot
  p <- ggplot(plot_data) +
    geom_segment(aes(x = from_id, xend = to_id,
                     y = -from_gen, yend = -to_gen,
                     alpha = if(!is.null(highlight_lineage)) highlight else NULL),
                 color = "blue") +
    geom_point(aes(x = from_id, y = -from_gen), color = "grey") +
    geom_point(aes(x = to_id, y = -to_gen), color = "grey") +
    theme_minimal() +
    labs(x = "Individual ID", y = "Generation") +
    theme(panel.grid.minor = element_blank())
  
  if(!is.null(highlight_lineage)) {
    p <- p + scale_alpha_manual(values = c(0.2, 1))
  }
  
  return(p)
}

calc_rich <- function(population) {
  population[population > 0] <- TRUE
  colSums(population)
}
fast_mean <- function(x) {
  x <- x[!is.na(x)]
  sum(x) / length(x)
}

calc_div <- function(population) {
  relabund <- population / rowSums(population,na.rm = T)
  -rowSums(relabund * log(relabund),na.rm = T)
}

fitness_func <- function(selection_parameter, optima, trait) {
  q <- exp(((trait - optima)^2) / -selection_parameter)
  q
}
sample_two_parents <- function(names, probs, N) {
  parents <- matrix(nrow = N, ncol = 2)
  
  for (i in 1:N) {
    parent1 <- sample(names, 1, prob = probs)
    names2 <- names[names != parent1]
    probs2 <- probs[names != parent1]
    probs2 <- probs2 / sum(probs2) # Normalize the probabilities
    parent2 <- sample(names2, 1, prob = probs2)
    parents[i, ] <- c(parent1, parent2)
  }
  return(parents)
}


fitness_func_bacgen <- function(selection_parameter, optima1, optima2, trait) {
  q <- exp((((trait - mean(c(optima1, optima2)))^2)) / -selection_parameter)
  q
}

mutate_trait <- function(trait = NULL, mutation_rate = NULL, 
                          mutation_range_low = NULL, mutation_range_high = NULL, 
                          TraitIsOptima = NULL, BinaryMutate = FALSE,DrawFromRange=TRUE) {
  
  # Handle array vs vector
  is_array <- !is.null(dim(trait))
  
  if (is_array) {
    # Save dimensions and dimnames
    orig_dim <- dim(trait)
    orig_dimnames <- dimnames(trait)
    
    # Flatten, process, and reconstruct
    flat_trait <- as.vector(trait)
    flat_result <- mutate_trait_cpp(
      trait = flat_trait,
      mutation_rate = mutation_rate,
      mutation_range_low = mutation_range_low, 
      mutation_range_high = mutation_range_high,
      TraitIsOptima = TraitIsOptima,
      BinaryMutate = BinaryMutate,
      DrawFromRange = DrawFromRange
    )
    
    # Reshape the result
    result <- array(flat_result, dim = orig_dim)
    dimnames(result) <- orig_dimnames
  } else {
    # Process directly for vectors
    result <- mutate_trait_cpp(
      trait = trait,
      mutation_rate = mutation_rate,
      mutation_range_low = mutation_range_low, 
      mutation_range_high = mutation_range_high,
      TraitIsOptima = TraitIsOptima,
      BinaryMutate = BinaryMutate,
      DrawFromRange = DrawFromRange
    )
  }
  
  return(result)
}

calculate_observed_heterozygosity <- function(genotype_array) {
  # genotype_array: 3D array with dimensions [individuals, alleles(2), loci]
  
  if (length(dim(genotype_array)) != 3) {
    stop("Input must be a 3D array with dimensions [individuals, alleles, loci]")
  }
  
  num_individuals <- dim(genotype_array)[1]
  num_loci <- dim(genotype_array)[3]
  
  # Calculate heterozygosity for each locus
  heterozygosity_per_locus <- numeric(num_loci)
  
  for (locus in 1:num_loci) {
    # Extract genotypes for this locus
    locus_genotypes <- genotype_array[, , locus]
    
    # Check if locus is monomorphic (all alleles the same)
    all_alleles <- unique(c(locus_genotypes[, 1], locus_genotypes[, 2]))
    
    if (length(all_alleles) <= 1) {
      # Monomorphic locus - set to NA (will be excluded from mean)
      heterozygosity_per_locus[locus] <- NA
    } else {
      # Count heterozygotes (individuals where allele1 != allele2)
      heterozygote_count <- sum(locus_genotypes[, 1] != locus_genotypes[, 2])
      heterozygosity_per_locus[locus] <- heterozygote_count / num_individuals
    }
  }
  
  # Return mean heterozygosity across polymorphic loci only
  mean_heterozygosity <- mean(heterozygosity_per_locus, na.rm = TRUE)
  
  return(list(
    mean_heterozygosity = mean_heterozygosity,
    per_locus_heterozygosity = heterozygosity_per_locus,
    polymorphic_loci = sum(!is.na(heterozygosity_per_locus)),
    total_loci = num_loci
  ))
}

# Simple function to calculate mean expected heterozygosity under Hardy-Weinberg equilibrium
calculate_expected_heterozygosity <- function(genotype_array) {
  # genotype_array: 3D array with dimensions [individuals, alleles(2), loci]
  
  if (length(dim(genotype_array)) != 3) {
    stop("Input must be a 3D array with dimensions [individuals, alleles, loci]")
  }
  
  num_individuals <- dim(genotype_array)[1]
  num_loci <- dim(genotype_array)[3]
  
  # Calculate expected heterozygosity for each locus
  expected_heterozygosity_per_locus <- numeric(num_loci)
  
  for (locus in 1:num_loci) {
    # Extract genotypes for this locus
    locus_genotypes <- genotype_array[, , locus]
    
    # Get all unique alleles and their frequencies
    all_alleles <- unique(c(locus_genotypes[, 1], locus_genotypes[, 2]))
    
    if (length(all_alleles) <= 1) {
      # Monomorphic locus - expected heterozygosity is 0, but set to NA for consistency
      expected_heterozygosity_per_locus[locus] <- NA
    } else {
      # Calculate allele frequencies
      total_alleles <- 2 * num_individuals
      allele_freqs <- sapply(all_alleles, function(a) {
        sum(locus_genotypes == a) / total_alleles
      })
      
      # Expected heterozygosity under HWE = 1 - sum(p^2)
      # This is the same as 2 * sum(p_i * p_j) for all i != j
      expected_heterozygosity_per_locus[locus] <- 1 - sum(allele_freqs^2)
    }
  }
  
  # Return mean expected heterozygosity across polymorphic loci only
  mean_expected_heterozygosity <- mean(expected_heterozygosity_per_locus, na.rm = TRUE)
  
  return(list(
    mean_expected_heterozygosity = mean_expected_heterozygosity,
    per_locus_expected_heterozygosity = expected_heterozygosity_per_locus,
    polymorphic_loci = sum(!is.na(expected_heterozygosity_per_locus)),
    total_loci = num_loci
  ))
}

Create_OffSpring_Pop_cpp <- function(host_pop, n_micro, env_pool, envpoolsize, X, fixed_envpool,
                                     selection_parameter_microbes, microbe_trait_list, host_microbe_optima,
                                     N_Species, env_condition,weightings) {
  microbe_names <- paste("Microbe", 1:N_Species, sep = "_")
  # offspring_population<-list()
  #  HostFitness_WithMicrobes<-matrix(data=NA,nrow=N_Species,ncol=ncol(host_pop))
  ENV_sampling_probability <- (env_pool) / envpoolsize # Calculate initial sampling probability based on environmental relative abundance
  #  HostFitness_WithMicrobes<-vector()
  #  weighted_samplingprob<-vector()
  names(ENV_sampling_probability) <- names(env_pool)
  offspring_population <- offspring_loopfunc(
    host_pop = host_pop, ENV_sampling_probability = ENV_sampling_probability,
    host_microbe_optima = host_microbe_optima, X = X,
    selection_parameter_microbes = selection_parameter_microbes,
    env_condition = env_condition, microbe_trait_list = microbe_trait_list,
    n_micro = n_micro, microbe_names = microbe_names,host_weighting = weightings
  ) # ,
  #                                          HostFitness_WithMicrobes = HostFitness_WithMicrobes,weighted_samplingprob = weighted_samplingprob)
  offspring_population <- list(
    offspring_population$host_pop, env_pool, microbe_trait_list, offspring_population$fitness_microbes,
    offspring_population$weighted_samplingprob
  )
  names(offspring_population) <- c("Child", "Env", "microbe_trait_list", "microbefitness", "microbe_samplingprob")
  offspring_population
}


lapply_wrapper_CPP <- function(XY,
                               HostPopulation, N_Microbes, envpoolsize, env_pool, fixed_envpool, generations, per_host_bac_gens, self_seed_prop,
                               selection_parameter_hosts, selection_parameter_microbes, selection_parameter_env,N_Species, traitpool_microbes,
                               generation_data_file, env_cond_val,
                               microbiome_importances, host_microbe_optima, mutation_rate_optima,mutation_rate_importance,
                               print_currentgen, selection_parameter_on_hosts,
                               importances_mutation, optima_mutate, nloci,nloci_importance,host_env_weightings,
                               lineage_track, Test_HWE
                               ) {
  if(lineage_track){
    lineage_tracker <- create_lineage_tracker(ncol(HostPopulation), generations)
  }  
  if (any(c(dim(microbiome_importances)[3] != nloci_importance, dim(host_microbe_optima)[3] != nloci))) {
    stop("One of the loci matrices has an incorrect number of loci")
  }
  temp_list <- list()
  print(paste("Current X Value is", XY[1], "Current EnvCon value is", XY[2]))
  gen_data <- list()
  microbe_names <- paste("Microbe", 1:N_Species, sep = "_")
  
  env_used <- NULL
  if (!nrow(env_cond_val) == generations) {
    print("env_cond_val should be a matrix of the same length as the number of generations you are simulating")
    stop()
  }
  
  
  if (!ncol(env_cond_val) == per_host_bac_gens) {
    print("Number of columns in env_conditions should be the number of bacterial generations")
    stop()
  }
  if (is.na(self_seed_prop)) {
    print("Self seeding proportion is NA, please correct")
  }
  coalescence_checker_matrix <- matrix(data = seq_len(ncol(HostPopulation)), ncol = 2, nrow = ncol(HostPopulation))
  rownames(coalescence_checker_matrix) <- paste("Host", seq_len(ncol(HostPopulation)), sep = "_")
  
  for (generation in 1:generations) {
    if (generation == 1) {
      env_used <- fixed_envpool
    }
    
    env_cond_val_used <- as.vector(env_cond_val[generation, ])
    
    # Calculate number of alleles per locus for optima traits
    optima_alleles_per_locus <- numeric(dim(host_microbe_optima)[3])
    for (locus in 1:dim(host_microbe_optima)[3]) {
      locus_genotypes <- host_microbe_optima[, , locus]
      all_alleles <- unique(c(locus_genotypes[, 1], locus_genotypes[, 2]))
      optima_alleles_per_locus[locus] <- length(all_alleles)
    }
    names(optima_alleles_per_locus) <- paste("Locus", 1:dim(host_microbe_optima)[3], sep = "_")
    
    # Calculate number of alleles per locus for importance traits
    importance_alleles_per_locus <- numeric(dim(microbiome_importances)[3])
    for (locus in 1:dim(microbiome_importances)[3]) {
      locus_genotypes <- microbiome_importances[, , locus]
      all_alleles <- unique(c(locus_genotypes[, 1], locus_genotypes[, 2]))
      importance_alleles_per_locus[locus] <- length(all_alleles)
    }
    names(importance_alleles_per_locus) <- paste("Locus", 1:dim(microbiome_importances)[3], sep = "_")

    if(Test_HWE == TRUE){
    hwe_result <- test_hwe_simple(host_microbe_optima)
    optima_Observed_heterozygosity<-hwe_result$summary$mean_observed_het
    optima_Expected_heterozygosity<-hwe_result$summary$mean_expected_het
    } else {
      hwe_result <- "Not Tested"
      optima_Observed_heterozygosity<-calculate_observed_heterozygosity(host_microbe_optima)
      optima_Expected_heterozygosity<-calculate_expected_heterozygosity(host_microbe_optima)
    }
    
    if (is.array(microbiome_importances)) {
      microbiome_importances_used <- rowSums(apply(microbiome_importances, c(1,2), mean))
      microbiome_importances_used[microbiome_importances_used<0]<-0
      microbiome_importances_used[microbiome_importances_used>1]<-1
    } else {
      microbiome_importances_used <- microbiome_importances
    }
    

    importances_Observed_heterozygosity<-calculate_observed_heterozygosity(microbiome_importances)
    importances_Expected_heterozygosity<-calculate_expected_heterozygosity(microbiome_importances)
    

    if (print_currentgen == TRUE) {
      print(paste("Current Generation is", generation))
    }
    
    if (is.array(host_env_weightings)) {
      host_env_weightings_used <-  rowSums(apply(host_env_weightings, c(1,2), mean))
      
    } else {
      host_env_weightings_used <- host_env_weightings
    }
    if (is.array(host_microbe_optima)) {
      host_microbe_optima_used <- rowSums(apply(host_microbe_optima, c(1,2), mean))
    } else {
      host_microbe_optima_used <- host_microbe_optima
    }

    for (bacgen in 1:per_host_bac_gens) {
      if (XY[2] + XY[3] > 1) {
        print("The contribution of the fixed environment (Y) and the the autocthonous environment (var_env_con) is greater than 1, please correct this ")
        stop()
      }

      gen_env_cond <- env_cond_val_used[bacgen]
      env_used <- process_microbe_probs(
        env_cond_val = gen_env_cond,
        fixed_envpool = fixed_envpool,
        HostPopulation = HostPopulation,
        N_Microbes = N_Microbes,
        envpoolsize = envpoolsize,
        selection_parameter_env = selection_parameter_env,
        XY = XY,
        traitpool_microbes = traitpool_microbes, env_used = env_used
      )

      env_fits <- env_used#weighted.mean(x = env_used[, 3], w = env_used[, 8])
      env_used <- env_used[, 8]
      names(env_used) <- names(fixed_envpool) # rownames(microbe_probs)
      # print(table(is.na(env_used)))
      if(bacgen > 1 ){
        HostPopulation <- process_host_cpp(
          HostPopulation = HostPopulation,
          selection_parameter_microbes = selection_parameter_microbes,
          host_microbe_optima = host_microbe_optima_used,
          env_condition = env_cond_val_used[bacgen],
          traitpool_microbes = traitpool_microbes,
          N_Microbes = N_Microbes,
          self_seed_prop = self_seed_prop, env_used = env_used,
          host_weighting=host_env_weightings_used,
          envpoolsize = envpoolsize
        )
      }
    }
    #  print(selection_parameter_hosts)
    # Now we've created our new environments, we now need to choose which members of a population reproduce
    # We can do this either neutrally (i.e., random chance) OR we can do this based on the fitness of a host based on what is provided by its microbiome
    host_fitnessvector <- numeric()
    
 #   print("Issue here")
    host_fitnessvector<- calculate_host_fitness_cpp(
      as.matrix(HostPopulation),
      traitpool_microbes, microbiome_importances_used,
      host_microbe_optima_used,
      selection_parameter_hosts, env_cond_val_used[per_host_bac_gens]#,cost_factor
    )
    mean_microbial_trait_vals<-host_fitnessvector[[3]]
    host_cost_factors<-host_fitnessvector[[2]]
    host_fitnessvector<-host_fitnessvector[[1]]
    host_phenotype<-generate_phenotype(as.matrix(HostPopulation),
                                       traitpool_microbes, microbiome_importances_used,
                                       host_microbe_optima_used)

    host_phenotype_nomicrobiome<-generate_phenotype(as.matrix(HostPopulation),
                                                    traitpool_microbes, rep(0,length(microbiome_importances_used)),
                                                    host_microbe_optima_used)
    
    nomicrobiomehost_fitnessvector <- numeric()
    
    nomicrobiomehost_fitnessvector <- calculate_host_fitness_cpp(
      as.matrix(HostPopulation),
      traitpool_microbes, rep(0, length(microbiome_importances)),
      host_microbe_optima_used,
      selection_parameter_hosts, env_cond_val_used[per_host_bac_gens]#,cost_factor
    )
    nomicrobiomehost_fitnessvector<-nomicrobiomehost_fitnessvector[[1]]
    
    #    names(host_fitnessvector)<-colnames(HostPopulation)
    #
    hostfitness_abs <- host_fitnessvector
    host_fitnessvector <- host_fitnessvector
    
      HostPopulationInt <- sample_two_parents(
        names = colnames(HostPopulation),
        N = ncol(HostPopulation), probs = host_fitnessvector
      )
      Ne_old<-4*ncol(HostPopulation)/(2+var(table(HostPopulationInt)))

      all_parent_ids <- colnames(HostPopulation) # or however you track the previous generation
      offspring_counts <- table(factor(HostPopulationInt, levels = all_parent_ids))
      Vk <- var(as.numeric(offspring_counts))
      Ne_new <- 4*ncol(HostPopulation) / (2 + Vk)

      
      if(lineage_track){lineage_tracker <- update_lineage_tracker(lineage_tracker, generation, HostPopulationInt)}
      num_unique_parents<-(unique(as.vector(HostPopulationInt)))
      # print(max(apply(host_microbe_optima,3,function(x){dim(table(x))})))
      new_host_microbe_optima <- simulate_mating_with_pairs_cpp(
        genotypes = host_microbe_optima, parent_pairs = HostPopulationInt,
        individual_names = colnames(HostPopulation), num_individuals = Host_PopSize,
        num_loci = dim(host_microbe_optima)[3]
      )
      
      #  print(max(apply(new_host_microbe_optima,3,function(x){dim(table(x))})))
      new_microbiome_importances <- simulate_mating_with_pairs_cpp(
        genotypes = microbiome_importances, parent_pairs = HostPopulationInt,
        individual_names = colnames(HostPopulation), num_individuals = Host_PopSize,
        num_loci = dim(microbiome_importances)[3]
      )
      
      
      new_coalescence_checker_matrix <- matrix(data = NA, ncol = 2, nrow = ncol(HostPopulation))
      rownames(coalescence_checker_matrix) <- paste("Host", seq_len(ncol(HostPopulation)), sep = "_")
      
      for (i in seq_len(nrow(HostPopulationInt))) {
        new_coalescence_checker_matrix[i, 1] <- sample(coalescence_checker_matrix[HostPopulationInt[i, 1], ], 1)
        new_coalescence_checker_matrix[i, 2] <- sample(coalescence_checker_matrix[HostPopulationInt[i, 2], ], 1)
      }
      
      new_hostpop <- matrix(
        data = NA, nrow = nrow(HostPopulation), ncol = ncol(HostPopulation),
        dimnames = dimnames(HostPopulation)
      )
      for (i in seq_len(nrow(HostPopulationInt))) {
        new_hostpop[, i] <- HostPopulation[, HostPopulationInt[i, 1]] + HostPopulation[, HostPopulationInt[i, 2]]
      }
      
      HostPopulation <- new_hostpop
      host_microbe_optima <- new_host_microbe_optima
      microbiome_importances <- new_microbiome_importances
      coalescence_checker_matrix <- new_coalescence_checker_matrix
      #host_env_weightings <- new_host_env_weightings
      # print(coalescence_checker_matrix)
      if (optima_mutate == TRUE) {
        if (mutation_rate_optima > 0) {
          host_microbe_optima <- mutate_trait(trait= host_microbe_optima,mutation_rate =  mutation_rate_optima,mutation_range_low =  -2.5,mutation_range_high = 2.5,TraitIsOptima=FALSE,BinaryMutate = F,DrawFromRange=TRUE)
        }
      }
      if (importances_mutation == TRUE) {
        if (mutation_rate_importance > 0) {
          microbiome_importances <- mutate_trait(trait = microbiome_importances, mutation_rate = mutation_rate_importance,mutation_range_low = 0,mutation_range_high = 0.5,TraitIsOptima=FALSE,BinaryMutate = F,DrawFromRange=TRUE)
        }
      }
      

      individual_names <- paste("Host", 1:Host_PopSize, sep = "_")
      loci_names <- paste("Locus", 1:nloci, sep = "_")
      importance_loci_names <- paste("Locus", 1:nloci_importance, sep = "_")  
      allele_names <- c("Allele_1", "Allele_2")
      # print(dim(microbiome_importances))
      dimnames(host_microbe_optima) <- list(individual_names, allele_names, loci_names)
      dimnames(microbiome_importances) <- list(individual_names, allele_names, importance_loci_names)
      #   dimnames(host_env_weightings) <- list(individual_names, allele_names, loci_names)
      
      colnames(HostPopulation) <- paste("Host", seq_len(ncol(HostPopulation)), sep = "_")
      rownames(coalescence_checker_matrix) <- paste("Host", seq_len(ncol(HostPopulation)), sep = "_")
    
    
    new_gen <- Create_OffSpring_Pop_cpp(
      host_pop = HostPopulation,
      n_micro = N_Microbes,
      env_pool = env_used,
      envpoolsize = envpoolsize,
      X = XY[1],
      fixed_envpool = fixed_envpool,
      microbe_trait_list = traitpool_microbes,
      selection_parameter_microbes = selection_parameter_microbes,
      host_microbe_optima = host_microbe_optima_used,
      N_Species = N_Species,
      env_condition = env_cond_val_used[per_host_bac_gens],
      weightings=host_env_weightings_used
    )
    
    HostPopulation <- new_gen$Child
    BrayDiv <- vegdist(t(HostPopulation), method = "bray")
    div_pergen <- calc_div(t(new_gen$Child)) / log(calc_rich(new_gen$Child))
    species_richess<-specnumber(t(HostPopulation))
    new_gen$HostFitness_Abs <- (hostfitness_abs)
    new_gen$nomicrobiomehost_fitnessvector <- nomicrobiomehost_fitnessvector
    new_gen$BrayDiv <- BrayDiv
    new_gen$HostMicrobeOptima <- host_microbe_optima
    
    gen_data[[generation]] <- list(
      div_pergen, new_gen$HostFitness_Abs, new_gen$microbefitness, env_used,
      env_fits, host_microbe_optima, new_gen$nomicrobiomehost_fitnessvector, new_gen$microbe_samplingprob,
      microbiome_importances, new_gen$BrayDiv, coalescence_checker_matrix, selection_parameter_hosts,
      host_microbe_optima_used, microbiome_importances_used,host_env_weightings_used,host_phenotype,num_unique_parents,host_phenotype_nomicrobiome,microbiome_importances,mean_microbial_trait_vals,species_richess,Ne_new,hwe_result,optima_Observed_heterozygosity,optima_Expected_heterozygosity,
      importances_Observed_heterozygosity,importances_Expected_heterozygosity,optima_alleles_per_locus,importance_alleles_per_locus,Ne_old
    )
    
    names(gen_data[[generation]]) <- c(
      "Diversity", "HostFitness", "MicrobeFitness", "env_used", "env_fits",
      "host_microbe_optima", "nomicrobiomehost_fitnessvector",
      "microbe_samplingprob", "microbiome_importances", "BrayDiv", "CoalCheck", "sel_strengths",
      "host_microbe_optima_used", "microbiome_importances_used","host_env_weightings_used","host_phenotype","num_unique_parents","host_phenotype_nomicrobiome", "microbiome_importances","mean_microbial_trait_vals","species_richess","Ne",
      "hwe_result","optima_Observed_heterozygosity","optima_Expected_heterozygosity",
      "importances_Observed_heterozygosity","importances_Expected_heterozygosity","optima_alleles_per_locus","importance_alleles_per_locus","Ne_old"
    ) # ,"HostPrefOptima","HostMicrobeOptima","host_env_weightings_used)
    if(lineage_track){
      gen_data[[generation]]$lineage_data <- lineage_tracker}
  }
  new_gen$GenData <- gen_data
  # print(dim(table(gen_data[[1]]$host_microbe_optima[,,10])))
  
  #  if(!is.na(generation_data_file)){
  #    readr::write_rds(gen_data,generation_data_file,compress = "gz")
  #  }
  temp_list[[paste0("X", XY[1], "_Y", XY[2], "_EnvCont", XY[3], "_HostSel", XY[4], "_MicrobeSel", XY[5])]] <- new_gen
  temp_list
  # })
}

simulate_mating_with_pairs <- function(genotypes, parent_pairs) {
  num_matings <- nrow(parent_pairs)
  num_loci <- dim(genotypes)[3]
  offspring <- array(dim = c(num_matings, 2, num_loci)) # Adjusted dimensions
  
  for (i in 1:num_matings) {
    # Extract the parent indices from the specified pairs
    parent1 <- parent_pairs[i, 1]
    parent2 <- parent_pairs[i, 2]
    
    # Create offspring
    for (locus in 1:num_loci) {
      # Randomly sample an allele from each parent
      offspring[i, 1, locus] <- sample(genotypes[parent1, , locus], 1)
      offspring[i, 2, locus] <- sample(genotypes[parent2, , locus], 1)
    }
  }
  return(offspring)
}


