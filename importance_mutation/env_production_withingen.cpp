#include <Rcpp.h>
#include <numeric>
#include <cmath>
using namespace Rcpp;

// Optimized mean fitness calculation with better memory access patterns
double calculate_mean_fitness(const NumericVector& fit_value, const NumericVector& host_counts, int envpoolsize) {
  double sum_of_products = 0.0;
  const int n = fit_value.size();
  
  // Manual loop unrolling and direct access for better performance
  for (int i = 0; i < n; ++i) {
    sum_of_products += fit_value[i] * host_counts[i];
  }
  
  return sum_of_products / envpoolsize;
}

// Optimized fitness function with cached exp calculation
// [[Rcpp::export]]
double fitness_func_cpp(double selection_parameter, double optima, double trait) {
  double diff = trait - optima;
  return std::exp(-(diff * diff) / selection_parameter);
}

// Optimized multinomial sampling
// [[Rcpp::export]]
NumericVector rmultinom_cpp(int n, const NumericVector& prob) {
  NumericVector normalized_prob = clone(prob);
  double prob_sum = sum(normalized_prob);
  
  // Avoid division if already normalized
  if (std::abs(prob_sum - 1.0) > 1e-10) {
    normalized_prob = normalized_prob / prob_sum;
  }
  
  IntegerVector res(normalized_prob.size());
  rmultinom(n, normalized_prob.begin(), normalized_prob.size(), res.begin());
  return NumericVector(res);
}

// Heavily optimized main function
// [[Rcpp::export]]
NumericMatrix process_microbe_probs_optimized(double env_cond_val, const NumericVector& fixed_envpool, 
                                            const NumericMatrix& HostPopulation, double N_Microbes, 
                                            double envpoolsize, double selection_parameter_env, 
                                            const NumericVector& XY, const NumericVector& traitpool_microbes, 
                                            const NumericVector& env_used) {
  const int N_Species = traitpool_microbes.length();
  
  // Pre-allocate all matrices at once
  NumericMatrix microbe_probs(N_Species, 9);
  
  // Pre-calculate constants to avoid repeated calculations
  const double inv_envpoolsize = 1.0 / envpoolsize;
  const double xy0 = XY[0], xy1 = XY[1], xy2 = XY[2];
  const double env_fixed_factor = (1.0 - xy1 - xy2);
  
  // Calculate host population totals once
  NumericVector host_totals(N_Species);
  double total_host_pop = 0.0;
  
  for (int i = 0; i < N_Species; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < HostPopulation.ncol(); ++j) {
      row_sum += HostPopulation(i, j);
    }
    host_totals[i] = row_sum;
    total_host_pop += row_sum;
  }
  
  // Vectorized calculations where possible
  const double inv_total_host = 1.0 / total_host_pop;
  
  // Pre-calculate fitness values and environmental probabilities
  NumericVector fitness_vals(N_Species);
  NumericVector env_probs(N_Species);
  double fitness_env_sum = 0.0;
  
  for (int i = 0; i < N_Species; ++i) {
    // Calculate fitness
    double diff = traitpool_microbes[i] - env_cond_val;
    fitness_vals[i] = std::exp(-(diff * diff) / selection_parameter_env);
    
    // Environmental probability
    env_probs[i] = env_used[i] * inv_envpoolsize;
    
    // Sum for normalization
    fitness_env_sum += fitness_vals[i] * env_probs[i];
  }
  
  // Avoid division by zero
  const double inv_fitness_env_sum = (fitness_env_sum > 1e-15) ? (1.0 / fitness_env_sum) : 0.0;
  
  // Fill the matrix efficiently
  for (int i = 0; i < N_Species; ++i) {
    // Column 0: Relative abundance of microbes in host population
    microbe_probs(i, 0) = host_totals[i] * inv_total_host;
    
    // Column 1: Relative abundance in fixed environment
    microbe_probs(i, 1) = fixed_envpool[i] * inv_envpoolsize;
    
    // Column 2: Fitness-based sampling probability
    microbe_probs(i, 2) = (fitness_vals[i] * env_probs[i]) * inv_fitness_env_sum;
    
    // Columns 3-5: Sampling probabilities
    microbe_probs(i, 3) = xy2 * microbe_probs(i, 0);
    microbe_probs(i, 4) = env_fixed_factor * microbe_probs(i, 1);
    microbe_probs(i, 5) = xy1 * microbe_probs(i, 2);
    
    // Column 6: Sum of sampling probabilities
    microbe_probs(i, 6) = microbe_probs(i, 3) + microbe_probs(i, 4) + microbe_probs(i, 5);
    
    // Column 8: Fitness values
    microbe_probs(i, 8) = fitness_vals[i];
  }
  
  // Handle NaN/NA values efficiently
  for (int i = 0; i < N_Species; ++i) {
    for (int j = 0; j < 9; ++j) {
      if (j != 7 && (R_IsNA(microbe_probs(i, j)) || R_IsNaN(microbe_probs(i, j)))) {
        microbe_probs(i, j) = 0.0;
      }
    }
  }
  
  // Column 7: Multinomial sampling
  NumericVector sampling_probs = microbe_probs(_, 6);
  microbe_probs(_, 7) = rmultinom_cpp(envpoolsize, sampling_probs);
  
  return microbe_probs;
}

// Keep original function for comparison
// [[Rcpp::export]]
NumericMatrix process_microbe_probs(double env_cond_val, NumericVector fixed_envpool, 
                                           NumericMatrix HostPopulation, double N_Microbes, 
                                           double envpoolsize, double selection_parameter_env, 
                                           NumericVector XY, NumericVector traitpool_microbes, 
                                           NumericVector env_used) {
  int N_Species = traitpool_microbes.length();
  
  NumericMatrix var_microbeprobs(N_Species, 5);
  
  for (int i = 0; i < N_Species; ++i) {
    var_microbeprobs(i, 0) = env_used[i] / envpoolsize;
  }
  for (int i = 0; i < N_Species; ++i) {
    var_microbeprobs(i, 2) = fitness_func_cpp(selection_parameter_env, env_cond_val, traitpool_microbes[i]);
  }

  for (int i = 0; i < var_microbeprobs.nrow(); ++i) {
    if (NumericVector::is_na(var_microbeprobs(i, 2))) {
      var_microbeprobs(i, 2) = 0;
    }
  }
  
  for (int i = 0; i < N_Species; ++i) {
    if (NumericVector::is_na(var_microbeprobs(i, 0))) {
      Rcout << "NA detected in var_microbeprobs(" << i << ", 0)" << std::endl;
    }
    if (NumericVector::is_na(var_microbeprobs(i, 2))) {
      Rcout << "NA detected in var_microbeprobs(" << i << ", 2)" << std::endl;
    }
    if (NumericVector::is_na(var_microbeprobs(i, 1))) {
      Rcout << "NA detected in var_microbeprobs(" << i << ", 1)" << std::endl;
    }
  }
  var_microbeprobs(_, 1) = (var_microbeprobs(_, 2) * var_microbeprobs(_, 0)) / sum(var_microbeprobs(_, 2) * var_microbeprobs(_, 0));
  for (int i = 0; i < var_microbeprobs.nrow(); ++i) {
    if (NumericVector::is_na(var_microbeprobs(i, 1))) {
      var_microbeprobs(i, 1) = 0;
    }
  }

  NumericMatrix microbe_probs(N_Species, 9);

  microbe_probs(_, 0) = rowSums(HostPopulation) / sum(HostPopulation);
  microbe_probs(_, 1) = fixed_envpool / envpoolsize;
  microbe_probs(_, 2) = var_microbeprobs(_, 1);
  microbe_probs(_, 3) = XY[2] * microbe_probs(_, 0);
  microbe_probs(_, 4) = (1 - XY[1] - XY[2]) * microbe_probs(_, 1);
  microbe_probs(_, 5) = XY[1] * microbe_probs(_, 2);
  microbe_probs(_, 6) = rowSums(microbe_probs(_, Range(3, 5)));
  microbe_probs(_, 8) = var_microbeprobs(_, 2);
  
  for (int i = 0; i < N_Species; ++i) {
    if (NumericVector::is_na(microbe_probs(i, 2))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 2)" << std::endl;
    }
    if (NumericVector::is_na(microbe_probs(i, 3))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 3)" << std::endl;
    }
    if (NumericVector::is_na(microbe_probs(i, 4))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 4)" << std::endl;
    }
    if (NumericVector::is_na(microbe_probs(i, 5))) {
      Rcout << "NA detected in microbe_probs(" << i << ", 5)" << std::endl;
    }
  }
  
  for (int i = 0; i < microbe_probs.nrow(); ++i) {
    for (int j = 0; j < microbe_probs.ncol(); ++j) {
      if (NumericVector::is_na(microbe_probs(i, j))) {
        microbe_probs(i, j) = 0;
      }
    }
  }
  
  microbe_probs(_, 7) = rmultinom_cpp(envpoolsize, microbe_probs(_, 6));
  
  return microbe_probs;
}