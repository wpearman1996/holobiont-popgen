#include <Rcpp.h>
#include <algorithm>  // for std::max and std::min
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector mutate_trait_cpp(NumericVector trait, 
                               double mutation_rate,
                               double mutation_range_low, 
                               double mutation_range_high,
                               bool TraitIsOptima = false,
                               bool BinaryMutate = false,
                               bool DrawFromRange = false,
                               double step_size = 0.15) {  // fix: valid default value
  
  int trait_length = trait.size();
  NumericVector result = clone(trait);
  
  NumericVector random_numbers(trait_length);
  for(int i = 0; i < trait_length; i++) {
    random_numbers[i] = R::runif(0, 1);
  }
  
  int mutate_count = 0;
  for(int i = 0; i < trait_length; i++) {
    if(random_numbers[i] < mutation_rate) {
      mutate_count++;
    }
  }
  
  if (mutate_count == 0) {
    return result;
  }
  
  IntegerVector mutate_indices(mutate_count);
  int counter = 0;
  
  for(int i = 0; i < trait_length; i++) {
    if(random_numbers[i] < mutation_rate) {
      mutate_indices[counter++] = i;
    }
  }
  
  if (!BinaryMutate) {
    if (TraitIsOptima) {
      for (int i = 0; i < mutate_count; i++) {
        result[mutate_indices[i]] = R::runif(mutation_range_low, mutation_range_high);
      }
    } else {
      if (DrawFromRange) {
        for (int i = 0; i < mutate_count; i++) {
          result[mutate_indices[i]] = R::runif(mutation_range_low, mutation_range_high);
        }
      } else {
        // Step mutation with clamping to valid range
        for (int i = 0; i < mutate_count; i++) {
          int idx = mutate_indices[i];
          double existing_value = result[idx];
          double new_value = R::runif(existing_value - step_size, existing_value + step_size);
          result[idx] = std::max(mutation_range_low, std::min(mutation_range_high, new_value));
        }
      }
    }
  } else {
    double value1 = mutation_range_low;
    double value2 = mutation_range_high;
    
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