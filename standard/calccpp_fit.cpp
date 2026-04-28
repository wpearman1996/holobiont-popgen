#include <Rcpp.h>
using namespace Rcpp;
// 
// // Fitness function
// double fitness_func_cpp(double selection_parameter, double optima, double trait) {
//   double q = exp(pow((trait - optima), 2) / -selection_parameter);
//   return q;
// }
// 
// // [[Rcpp::export]]
// NumericVector calculate_host_fitness_cpp(NumericMatrix HostPopulation, NumericVector traitpool_microbes, NumericVector microbiome_importances, NumericVector host_preference_optima, double selection_parameter_hosts, double env_cond_val_used) {
// 
//   int n = HostPopulation.ncol();
//   NumericVector hostfitness(n);
// 
//   for(int host = 0; host < n; host++) {
//     NumericVector pop = HostPopulation(_, host);
//     double mean_microbial_trait_val = sum(traitpool_microbes * pop) / sum(pop);
//     double composite_host_trait = ((mean_microbial_trait_val * microbiome_importances[host]) + (host_preference_optima[host] * (1 - microbiome_importances[host])));
// 
//     hostfitness[host] = fitness_func_cpp(selection_parameter_hosts, env_cond_val_used, composite_host_trait);
//   }
// 
//   return hostfitness;
// }
// 
// 
// 
// // [[Rcpp::export]]
// NumericVector generate_phenotype(NumericMatrix HostPopulation, NumericVector traitpool_microbes, NumericVector microbiome_importances, NumericVector host_preference_optima) {
// 
//   int n = HostPopulation.ncol();
//   NumericVector composite_host_trait(n);
// 
//   for(int host = 0; host < n; host++) {
//     NumericVector pop = HostPopulation(_, host);
//     double mean_microbial_trait_val = sum(traitpool_microbes * pop) / sum(pop);
//     composite_host_trait[host] = ((mean_microbial_trait_val * microbiome_importances[host]) + (host_preference_optima[host] * (1 - microbiome_importances[host])));
// 
//   }
//   return composite_host_trait;
// }
 
 
 // Fitness function with direct importance penalty
 double fitness_func_cpp(double selection_parameter, double optima, double trait, double importance, double cost_factor) {
   double q = exp(pow((trait - optima), 2) / -selection_parameter);
   //double q = exp(-(pow((trait - optima), 2)) / (2*pow(selection_parameter,2)));
//   double penalty = (pow(cost_factor, importance));
//   Rcpp::Rcout << "Received penalty: " << q * penalty << std::endl;
   return q; // * cost_factor;
 }
 
 inline double myAbs(double x) {
   return sqrt(x * x);
 }
 
 // [[Rcpp::export]]
 List calculate_host_fitness_cpp(NumericMatrix HostPopulation, NumericVector traitpool_microbes, 
                                 NumericVector microbiome_importances, NumericVector host_preference_optima, 
                                 double selection_parameter_hosts, double env_cond_val_used) {
   
   int n = HostPopulation.ncol();
   NumericVector hostfitness(n);
   NumericVector mean_microbial_trait_vals(n);
   NumericVector differences(n);
   NumericVector cost_factors_diffs(n);
   
   // Precompute mean microbial trait values for each host
   for (int host = 0; host < n; host++) {
     NumericVector pop = HostPopulation(_, host);
     mean_microbial_trait_vals[host] = sum(traitpool_microbes * pop) / sum(pop);
   }
   
   // Identify maximum absolute difference between mean_microbial_trait_vals and host_preference_optima
   NumericVector diff_vec = mean_microbial_trait_vals - host_preference_optima;
   NumericVector diff_abs(diff_vec.size());
   for (int i = 0; i < diff_vec.size(); i++) {
     diff_abs[i] = myAbs(diff_vec[i]);
   }
   
   double max_abs_diff = max(diff_abs);
   
   // Compute fitness for each host and store cost factor
   for (int host = 0; host < n; host++) {
     double mean_microbial_trait_val = mean_microbial_trait_vals[host];
     double composite_host_trait = ((mean_microbial_trait_val * microbiome_importances[host]) + 
                                    (host_preference_optima[host] * (1 - microbiome_importances[host])));
     double diff = host_preference_optima[host] - mean_microbial_trait_val;
     double cost_factor_adjusted = myAbs(diff)/max_abs_diff; //10*tanh(myAbs(diff));
     cost_factor_adjusted = exp(-cost_factor_adjusted);
     cost_factor_adjusted = (pow(cost_factor_adjusted,microbiome_importances[host]));
     cost_factors_diffs[host] = cost_factor_adjusted;
  //  cost_factor_adjusted = (cost_factor_adjusted);
     // if difference is high, and microbiome importance is high - then the decline it fitness will be high
     // if difference 
     
     differences[host] = myAbs(diff);
//     Rcpp::Rcout << "Host " << host << " cost_factor_adjusted: " << cost_factor_adjusted << std::endl;
//     Rcpp::Rcout << "Host " << host << " max_abs_diff: " << max_abs_diff  << std::endl;
//     Rcpp::Rcout << "Host " << host << " Difference: " << ((host_preference_optima[host] - mean_microbial_trait_val))   << std::endl;
     hostfitness[host] = fitness_func_cpp(selection_parameter_hosts, env_cond_val_used, composite_host_trait, 
                                          microbiome_importances[host], cost_factor_adjusted);
   }
   
   return List::create(Named("hostfitness") = hostfitness, 
                       Named("differences") = differences,
                       Named("mean_microbial_trait_vals") = mean_microbial_trait_vals);
 }
 
 // [[Rcpp::export]]
 NumericVector generate_phenotype(NumericMatrix HostPopulation, NumericVector traitpool_microbes, 
                                          NumericVector microbiome_importances, NumericVector host_preference_optima) {
   
   int n = HostPopulation.ncol();
   NumericVector composite_host_trait(n);
   
   for(int host = 0; host < n; host++) {
     NumericVector pop = HostPopulation(_, host);
     double mean_microbial_trait_val = sum(traitpool_microbes * pop) / sum(pop);

     composite_host_trait[host] = ((mean_microbial_trait_val * microbiome_importances[host]) + 
                                    (host_preference_optima[host] * (1 - microbiome_importances[host])));
   }


   return composite_host_trait;
 }
 
 
 
 // Simplified fitness calculation without cost factor
 
#include <Rcpp.h>
 using namespace Rcpp;
 
 // Simple fitness function - just the basic Gaussian selection
 double fitness_func_simple(double selection_parameter, double optima, double trait) {
   return exp(pow((trait - optima), 2) / -selection_parameter);
 }
 
 // Fitness function with dependency penalty only
 // [[Rcpp::export]]
 double fitness_func_with_dependency_penalty(double selection_parameter, double optima, 
                                             double composite_trait, double importance, 
                                             double host_trait, double microbe_trait,
                                             double penalty_strength = 1.0) {
   
   // Base fitness from composite trait
   double fitness_composite = exp(pow((composite_trait - optima), 2) / -selection_parameter);
   
   // Fitness if host used only its own trait
   double fitness_host_only = exp(pow((host_trait - optima), 2) / -selection_parameter);
   
   // Microbiome benefit
   double microbiome_benefit = fitness_composite - fitness_host_only;
   double normalized_benefit = std::max(0.0, microbiome_benefit);
   
   // Dependency penalty: penalize high importance when benefit is low
   double dependency_penalty = importance * (1.0 - normalized_benefit) * penalty_strength;
   
   // Apply only the dependency penalty
   return fitness_composite * exp(-dependency_penalty);
 }
 
 // Simplified host fitness calculation
 // [[Rcpp::export]]
 List calculate_host_fitness_simplified(NumericMatrix HostPopulation,
                                        NumericVector traitpool_microbes,
                                        NumericVector microbiome_importances,
                                        NumericVector host_preference_optima,  // Can be renamed to host_traits
                                        double selection_parameter_hosts,
                                        double env_cond_val_used,
                                        double penalty_strength = 1.0,
                                        bool use_penalty = true) {
   
   int n = HostPopulation.ncol();
   NumericVector hostfitness(n);
   NumericVector mean_microbial_trait_vals(n);
   NumericVector dependency_penalties(n);
   
   // Calculate mean microbial trait values
   for (int host = 0; host < n; host++) {
     NumericVector pop = HostPopulation(_, host);
     mean_microbial_trait_vals[host] = sum(traitpool_microbes * pop) / sum(pop);
   }
   
   // Calculate fitness for each host
   for (int host = 0; host < n; host++) {
     double host_trait = host_preference_optima[host];  // Now just the host's own trait
     double microbe_trait = mean_microbial_trait_vals[host];
     double importance = microbiome_importances[host];
     
     // Composite trait (weighted average)
     double composite_trait = (microbe_trait * importance) + (host_trait * (1 - importance));
     
     if (use_penalty) {
       hostfitness[host] = fitness_func_with_dependency_penalty(
         selection_parameter_hosts, env_cond_val_used, composite_trait,
         importance, host_trait, microbe_trait, penalty_strength
       );
       
       // Calculate penalty for diagnostics
       double fitness_no_penalty = exp(pow((composite_trait - env_cond_val_used), 2) / -selection_parameter_hosts);
       dependency_penalties[host] = fitness_no_penalty - hostfitness[host];
       
     } else {
       // Just basic fitness from composite trait
       hostfitness[host] = exp(pow((composite_trait - env_cond_val_used), 2) / -selection_parameter_hosts);
       dependency_penalties[host] = 0.0;
     }
   }
   
   return List::create(Named("hostfitness") = hostfitness,
                       Named("mean_microbial_trait_vals") = mean_microbial_trait_vals,
                       Named("dependency_penalties") = dependency_penalties);
 }
 