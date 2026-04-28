#include <Rcpp.h>
#include <unordered_map>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simulate_mating_with_pairs_cpp(NumericVector genotypes, 
                                             CharacterMatrix parent_pairs,
                                             CharacterVector individual_names,
                                             int num_individuals, 
                                             int num_loci) {
  int num_matings = parent_pairs.nrow();
  
  // Create a map from individual names to indices
  std::unordered_map<std::string, int> name_to_index;
  for (int i = 0; i < individual_names.size(); ++i) {
    name_to_index[as<std::string>(individual_names[i])] = i;
  }
  
  // Create an offspring array with the desired dimensions
  NumericVector offspring(num_matings * 2 * num_loci);
  offspring.attr("dim") = Dimension(num_matings, 2, num_loci);
  
  for (int i = 0; i < num_matings; i++) {
    // Get parent indices from their names
    int parent1 = name_to_index[as<std::string>(parent_pairs(i, 0))];
    int parent2 = name_to_index[as<std::string>(parent_pairs(i, 1))];
    
    for (int locus = 0; locus < num_loci; locus++) {
      // Randomly sample an allele from each parent
      int allele1 = R::runif(0, 1) < 0.5 ? 0 : 1; // Choose allele 1 or 2 from parent1
      int allele2 = R::runif(0, 1) < 0.5 ? 0 : 1; // Choose allele 1 or 2 from parent2
      
      offspring(i + num_matings * (0 + 2 * locus)) = 
        genotypes(parent1 + num_individuals * (allele1 + 2 * locus));
      
      offspring(i + num_matings * (1 + 2 * locus)) = 
        genotypes(parent2 + num_individuals * (allele2 + 2 * locus));
    }
  }
  
  return offspring;
}
