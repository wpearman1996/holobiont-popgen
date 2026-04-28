#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double fitness_func_bacgen_cpp(double selection_parameter, double optima1, double optima2, double trait, double weighting) {
  // Calculate the mean of optima1 and optima2
  double meanval = ((weighting * optima1) + ((1 - weighting) * optima2));
  
  // Calculate fitness function
  double q = exp(pow((trait - meanval), 2) / -selection_parameter);
  
  return q;
}


// Function to simulate multinomial distribution in C++
// This is a simple version, replace with your actual implementation if needed

// [[Rcpp::export]]
NumericVector rmultinom_cpp(int n, NumericVector prob) {
  SEXP res = Rf_allocVector(INTSXP, prob.size());
  prob = prob / sum(prob); // Adjust probabilities to sum to 1
  rmultinom(n, prob.begin(), prob.size(), INTEGER(res));
  return NumericVector(res);
}


// [[Rcpp::export]]
List offspring_loopfunc(NumericMatrix host_pop, NumericVector ENV_sampling_probability, NumericVector host_microbe_optima,
                             double X, double selection_parameter_microbes, double env_condition, NumericVector microbe_trait_list,
                             int n_micro, StringVector microbe_names, NumericVector host_weighting) {
  
  int N_Species = host_pop.nrow();
  int n_hosts = host_pop.ncol();
  
  // Pre-allocate output vectors
  NumericVector fitness_microbes(n_hosts);
  NumericVector weighted_samplingprob(n_hosts);
  
  // Pre-calculate constants
  double one_minus_X = 1.0 - X;
  double inv_n_micro = 1.0 / n_micro;
  
  // Pre-allocate working vectors (reuse across iterations)
  NumericVector microbial_fitness(N_Species);
  NumericVector sampling_probability(N_Species);
  NumericVector combined_prob(N_Species);
  
  for (int i = 0; i < n_hosts; i++) {
    double weighting_used = host_weighting[i];
    double host_optimum = host_microbe_optima[i];
    
    // Calculate microbial fitness for all species (this is expensive, so do it efficiently)
    for (int z = 0; z < N_Species; ++z) {
      microbial_fitness[z] = fitness_func_bacgen_cpp(selection_parameter_microbes, host_optimum, 
                                                    env_condition, microbe_trait_list[z], weighting_used);
    }
    
    // Calculate combined probabilities and sampling probabilities in one pass
    double sum_sampling_prob = 0.0;
    for (int z = 0; z < N_Species; ++z) {
      // Equivalent to: parprob = parents[z] / n_micro
      double parprob = host_pop(z, i) * inv_n_micro;
      
      // Equivalent to: ParProvided = X * parprob, EnvProvided = (1-X) * envprob
      // combined_prob = ParProvided + EnvProvided
      combined_prob[z] = X * parprob + one_minus_X * ENV_sampling_probability[z];
      
      // sampling_probability = microbial_fitness * combined_prob
      sampling_probability[z] = microbial_fitness[z] * combined_prob[z];
      sum_sampling_prob += sampling_probability[z];
    }
    
    // Normalize sampling probabilities
    double inv_sum = 1.0 / sum_sampling_prob;
    for (int z = 0; z < N_Species; ++z) {
      sampling_probability[z] *= inv_sum;
    }
    
    // Sample microbes
    NumericVector sampled_counts = rmultinom_cpp(n_micro, sampling_probability);
    
    // Update host population and calculate metrics in one pass
    double sum_weighted_fitness = 0.0;
    double sum_weighted_sampling = 0.0;
    double sum_counts = 0.0;
    
    for (int z = 0; z < N_Species; ++z) {
      double count = sampled_counts[z];
      host_pop(z, i) = count;  // Update host_pop
      
      sum_counts += count;
      sum_weighted_fitness += microbial_fitness[z] * count;
      sum_weighted_sampling += sampling_probability[z] * count;
    }
    
    // Calculate final metrics
    fitness_microbes[i] = sum_weighted_fitness / sum_counts;
    weighted_samplingprob[i] = sum_weighted_sampling / sum_counts;
  }
  
  // Optimized NA check - check by column (better cache locality)
  for (int j = 0; j < n_hosts; j++) {
    for (int i = 0; i < N_Species; i++) {
      if (Rcpp::NumericMatrix::is_na(host_pop(i, j))) {
        Rcpp::stop("Error: NA values detected in host_pop.");
      }
    }
  }
  
  return List::create(Named("host_pop") = host_pop,
                      Named("fitness_microbes") = fitness_microbes,
                      Named("weighted_samplingprob") = weighted_samplingprob);
}
