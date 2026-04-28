#!/usr/bin/env Rscript

# =====================================================================
# HOLOBIONT SIMULATION - SETUP SCRIPT FOR SLURM ARRAY
# This script generates the job mapping file needed for array execution
# =====================================================================

setwd("/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/Importance_Mutation")

cat("HOLOBIONT SIMULATION - SLURM ARRAY SETUP\n")
cat("=========================================\n\n")

# =====================================================================
# CONFIGURABLE PARAMETERS (MUST MATCH MAIN SCRIPT)
# =====================================================================

# Evolutionary Parameters
VERTICAL_INHERITANCE <- c(0, 0.1, 0.2, 0.5, 0.9, 1)     
ENV_SELF_SEEDING <- c(0.9)                
HOST_CONTRIBUTION <- c(0.01)              
HOST_SEL_STRENGTH <- c(0.1,Inf)
MICROBE_SEL_STRENGTH <- c(0.1,Inf)
ENV_SEL_STRENGTH <- c(0.1,Inf)

# Genetic Architecture Parameters
LOCI_COUNTS <- c(100)           
IMPORTANCE_VALUES <- c(0)

WEIGHT_VALUES <- c(0.25)                   
MUTATION_RATE <- c(1e-4)      

# Simulation Run Parameters
N_REPS_TOTAL <- 10
OUTPUT_DIR <- "./results_importance_mutation_standard_rerun"
COMBINATION_REGISTRY_FILE <- file.path(OUTPUT_DIR, "combination_registry.RDS")
JOB_MAPPING_FILE <- file.path(OUTPUT_DIR, "job_mapping.RDS")

# Environmental Parameters
ENV_FILE_PATH <- "/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/env_list_rerun.RDS"
ENV_OBJECT_NAME <- "env_list_rerun"

# Control parameters
OVERWRITE_EXISTING <- FALSE  # Set TRUE to re-run all simulations
REPLICATES_TO_RUN <- 1:N_REPS_TOTAL  # Can specify subset: c(1, 3, 5)

# =====================================================================
# LOAD LIBRARIES
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(digest)
})

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

create_combo_hash <- function(combo_row) {
  combo_string <- paste(
    combo_row["X"], combo_row["Z"], combo_row["Y"], 
    combo_row["S_host"], combo_row["S_microbe"], combo_row["S_env"],
    combo_row["EnvType"], combo_row["Weighting"], combo_row["Importance"],
    combo_row["MutationRate"], combo_row["Genetics"],
    sep = "_"
  )
  return(digest::digest(combo_string, algo = "md5"))
}

create_combo_hashes <- function(combo_df) {
  hash_list <- vector("character", nrow(combo_df))
  for (i in seq_len(nrow(combo_df))) {
    combo_string <- paste(
      combo_df$X[i], combo_df$Z[i], combo_df$Y[i], 
      combo_df$S_host[i], combo_df$S_microbe[i], combo_df$S_env[i],
      combo_df$EnvType[i], combo_df$Weighting[i], combo_df$Importance[i],
      combo_df$MutationRate[i], combo_df$Genetics[i],
      sep = "_"
    )
    hash_list[i] <- digest::digest(combo_string, algo = "md5")
  }
  return(hash_list)
}

# =====================================================================
# CREATE OUTPUT DIRECTORY
# =====================================================================

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", OUTPUT_DIR, "\n\n")

# =====================================================================
# LOAD ENVIRONMENTAL CONDITIONS
# =====================================================================

cat("Loading environmental conditions...\n")
sec_varyring_envs_multigen <- read_rds(ENV_FILE_PATH)
sec_varyring_envs_multigen<-sec_varyring_envs_multigen#[1:3]
if (!is.list(sec_varyring_envs_multigen)) {
  sec_varyring_envs_multigen <- list(as.matrix(sec_varyring_envs_multigen, ncol = 1))
  names(sec_varyring_envs_multigen) <- "env_1gen"
}
cat("Loaded", length(sec_varyring_envs_multigen), "environmental condition(s)\n\n")

# =====================================================================
# CREATE/LOAD COMBINATION REGISTRY
# =====================================================================

cat("Creating parameter combinations...\n")

# Generate all combinations
genetic_names <- paste0("genetics_nloci_", LOCI_COUNTS)
importance_names <- paste0("importances_init", IMPORTANCE_VALUES)
weighting_names <- paste0("weightings_", WEIGHT_VALUES)

all_combinations <- expand.grid(
  X = VERTICAL_INHERITANCE,
  Z = ENV_SELF_SEEDING,
  Y = HOST_CONTRIBUTION,
  S_host = HOST_SEL_STRENGTH,
  S_microbe = MICROBE_SEL_STRENGTH,
  S_env = ENV_SEL_STRENGTH,
  EnvType = names(sec_varyring_envs_multigen),
  Weighting = weighting_names,
  Importance = importance_names,
  MutationRate = MUTATION_RATE,
  Genetics = genetic_names,
  stringsAsFactors = FALSE
)

# Apply filtering
all_combinations <- all_combinations %>%
  filter(S_host == S_microbe, S_host == S_env)

cat("Total combinations after filtering:", nrow(all_combinations), "\n")

# Load or create registry
if (file.exists(COMBINATION_REGISTRY_FILE)) {
  cat("Loading existing combination registry...\n")
  existing_registry <- read_rds(COMBINATION_REGISTRY_FILE)
  cat("Existing combinations:", nrow(existing_registry), "\n")
  
  # Check for new combinations
  all_combinations$combo_hash <- create_combo_hashes(all_combinations)
  new_mask <- !all_combinations$combo_hash %in% existing_registry$combo_hash
  new_combinations <- all_combinations[new_mask, ]
  
  if (nrow(new_combinations) > 0) {
    cat("Found", nrow(new_combinations), "new combinations to add\n")
    max_id <- max(existing_registry$combo_id)
    new_combinations$combo_id <- seq(from = max_id + 1, length.out = nrow(new_combinations))
    new_combinations$date_added <- as.character(Sys.time())
    
    registry <- bind_rows(existing_registry, new_combinations)
    cat("Updated registry now has", nrow(registry), "combinations\n")
  } else {
    cat("No new combinations to add\n")
    registry <- existing_registry
  }
} else {
  cat("Creating new combination registry...\n")
  all_combinations$combo_id <- seq_len(nrow(all_combinations))
  all_combinations$combo_hash <- create_combo_hashes(all_combinations)
  all_combinations$date_added <- as.character(Sys.time())
  registry <- all_combinations
  cat("Created registry with", nrow(registry), "combinations\n")
}

# Save registry
write_rds(registry, COMBINATION_REGISTRY_FILE, compress = "gz")
cat("Saved combination registry to:", COMBINATION_REGISTRY_FILE, "\n\n")

# =====================================================================
# CREATE JOB MAPPING
# =====================================================================

cat("Creating job mapping for SLURM array...\n")

job_mapping <- expand.grid(
  combo_id = registry$combo_id,
  replicate = REPLICATES_TO_RUN,
  stringsAsFactors = FALSE
) %>%
  arrange(combo_id, replicate) %>%
  mutate(task_id = row_number())

# Filter out jobs that already exist (unless OVERWRITE_EXISTING is TRUE)
if (!OVERWRITE_EXISTING) {
  cat("Checking for existing output files...\n")
  
  job_mapping$file_exists <- sapply(seq_len(nrow(job_mapping)), function(i) {
    file_path <- file.path(OUTPUT_DIR, 
                           paste0("rep_", job_mapping$replicate[i], 
                                  "_simulation_iteration_", job_mapping$combo_id[i], ".RDS"))
    file.exists(file_path)
  })
  
  existing_count <- sum(job_mapping$file_exists)
  cat("Found", existing_count, "existing output files\n")
  
  # Keep only jobs that need to run
  job_mapping <- job_mapping %>%
    filter(!file_exists) %>%
    select(-file_exists) %>%
    mutate(task_id = row_number())  # Renumber task IDs
  
  cat("Jobs remaining to run:", nrow(job_mapping), "\n")
} else {
  cat("OVERWRITE mode: will re-run all simulations\n")
}

if (nrow(job_mapping) == 0) {
  cat("\n")
  cat("No jobs to run! All simulations are complete.\n")
  cat("Set OVERWRITE_EXISTING = TRUE in this script to re-run all simulations.\n")
  quit(save = "no", status = 0)
}

# Save job mapping
write_rds(job_mapping, JOB_MAPPING_FILE, compress = "gz")
cat("Saved job mapping to:", JOB_MAPPING_FILE, "\n\n")

# =====================================================================
# PRINT SUMMARY
# =====================================================================

cat(paste(rep("=", 70), collapse=""), "\n")
cat("SETUP COMPLETE - READY FOR SLURM SUBMISSION\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("Summary:\n")
cat("- Total parameter combinations:", nrow(registry), "\n")
cat("- Replicates per combination:", length(REPLICATES_TO_RUN), "\n")
cat("- Total jobs to run:", nrow(job_mapping), "\n")
cat("- Array indices: 1-", nrow(job_mapping), "\n\n")

cat("Job distribution:\n")
jobs_per_combo <- job_mapping %>%
  count(combo_id, name = "n_replicates") %>%
  count(n_replicates, name = "n_combinations")
for (i in seq_len(nrow(jobs_per_combo))) {
  cat(sprintf("- %d combinations with %d replicates\n", 
              jobs_per_combo$n_combinations[i], jobs_per_combo$n_replicates[i]))
}

cat("\n")
cat("Next steps:\n")
cat("1. Review the SLURM submission script (holobiont_array.sl)\n")
cat("2. Adjust walltime, memory, and other resources if needed\n")
cat("3. Submit with: sbatch holobiont_array.sl\n")
cat("4. Monitor with: squeue -u $USER\n")
cat("5. After completion, run the analysis script\n\n")

# Save summary stats for easy reference
summary_stats <- list(
  n_combinations = nrow(registry),
  n_replicates = length(REPLICATES_TO_RUN),
  n_jobs = nrow(job_mapping),
  array_range = c(1, nrow(job_mapping)),
  setup_date = Sys.time(),
  overwrite_mode = OVERWRITE_EXISTING
)
write_rds(summary_stats, file.path(OUTPUT_DIR, "setup_summary.RDS"))

cat("Setup complete!\n")
