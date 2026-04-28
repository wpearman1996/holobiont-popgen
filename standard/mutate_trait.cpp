#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector mutate_trait_cpp(NumericVector trait, 
                               double mutation_rate,
                               double mutation_range_low, 
                               double mutation_range_high,
                               bool TraitIsOptima = false,
                               bool BinaryMutate = false,
                               bool DrawFromRange = false) {
  
  int trait_length = trait.size();
  NumericVector result = clone(trait); // Make a copy
  
  // Generate random numbers for mutation decision
  NumericVector random_numbers(trait_length);
  for(int i = 0; i < trait_length; i++) {
    random_numbers[i] = R::runif(0, 1);
  }
  
  // Count how many elements will mutate
  int mutate_count = 0;
  for(int i = 0; i < trait_length; i++) {
    if(random_numbers[i] < mutation_rate) {
      mutate_count++;
    }
  }
  
  // If nothing to mutate, return the original trait
  if (mutate_count == 0) {
    return result;
  }
  
  // Create a vector to store indices of elements that should mutate
  IntegerVector mutate_indices(mutate_count);
  int counter = 0;
  
  // Find indices of elements that should mutate
  for(int i = 0; i < trait_length; i++) {
    if(random_numbers[i] < mutation_rate) {
      mutate_indices[counter++] = i;
    }
  }
  
  if (!BinaryMutate) {
    if (TraitIsOptima) {
      // Generate new alleles from a uniform distribution
      for (int i = 0; i < mutate_count; i++) {
        result[mutate_indices[i]] = R::runif(mutation_range_low, mutation_range_high);
      }
    } else {
      if (DrawFromRange) {
        // Draw new alleles randomly from the full range
        for (int i = 0; i < mutate_count; i++) {
          result[mutate_indices[i]] = R::runif(mutation_range_low, mutation_range_high);
        }
      } else {
        // Mutate each allele within a range of its current value
        for (int i = 0; i < mutate_count; i++) {
          int idx = mutate_indices[i];
          double existing_value = result[idx];
          result[idx] = R::runif(existing_value - 0.15, existing_value + 0.15);
        }
      }
    }
  } else {
    // For binary mutation, we just have two possible values
    double value1 = mutation_range_low;
    double value2 = mutation_range_high;
    
    // Apply the mutation by swapping values
    for (int i = 0; i < mutate_count; i++) {
      int idx = mutate_indices[i];
      double existing_value = result[idx];
      
      if (existing_value == value1) {
        result[idx] = value2;
      } else {
        result[idx] = value1;
      }
    }
  }
  
  return result;
}