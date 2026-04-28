#!/usr/bin/env Rscript

# =====================================================================
# HOLOBIONT SIMULATION - SLURM ARRAY VERSION
# Each array task runs ONE simulation (combo_id + replicate pair)
# =====================================================================

setwd("/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov")

# Get SLURM array task ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No task ID provided. This script must be run via SLURM array job.")
}

TASK_ID <- as.numeric(args[1])
cat("Starting SLURM array task:", TASK_ID, "\n")
cat("Hostname:", Sys.info()["nodename"], "\n")
cat("PID:", Sys.getpid(), "\n\n")

# =====================================================================
# CONFIGURABLE PARAMETERS
# =====================================================================

# Basic Population Parameters
Host_PopSize <- 200         
MicrobePopSize <- 10^6        
EnvPoolSize <- 10^8           
N_Species <- 400              
N_Traits <- 25                
nhost_gens <- 3000            

# Evolutionary Parameters
VERTICAL_INHERITANCE <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.9, 1)      
ENV_SELF_SEEDING <- c(0.9)                
HOST_CONTRIBUTION <- c(0.01)              
HOST_SEL_STRENGTH <- c(0.1)
MICROBE_SEL_STRENGTH <- c(0.1)
ENV_SEL_STRENGTH <- c(0.1)

# Genetic Architecture Parameters
LOCI_COUNTS <- c(100)           
IMPORTANCE_VALUES <- c(0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.1,
                       0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.30,0.35,0.4,0.45,0.5)

WEIGHT_VALUES <- c(0.25)                   
MUTATION_RATE <- c(1e-4)      

# Simulation Run Parameters
N_REPS_TOTAL <- 10
OUTPUT_DIR <- "./results_rerun_standard_11Feb"
SAVE_COMPRESSED <- TRUE
COMBINATION_REGISTRY_FILE <- file.path(OUTPUT_DIR, "combination_registry.RDS")
JOB_MAPPING_FILE <- file.path(OUTPUT_DIR, "job_mapping.RDS")

# Environmental Parameters
ENV_FILE_PATH <- "/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/env_list_rerun.RDS"
ENV_OBJECT_NAME <- "env_list_rerun"
USE_ENV_TYPES <- c("1gen")
ENV_PATTERNS <- c("nochange", "sudden", "incmean", "incvar")

# =====================================================================
# LOAD REQUIRED LIBRARIES
# =====================================================================

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(readr)
  library(abind)
  library(vegan)
  library(digest)
})

# Load C++ functions
Rcpp::sourceCpp("/nesi/project/uoa04039/Simulation_Source_Code/withingen_process.cpp")
Rcpp::sourceCpp("/nesi/project/uoa04039/Simulation_Source_Code/calccpp_fit.cpp")
Rcpp::sourceCpp("/nesi/project/uoa04039/Simulation_Source_Code/env_production_withingen.cpp")
Rcpp::sourceCpp("/nesi/project/uoa04039/Simulation_Source_Code/simulate_mating_cpp.cpp")
Rcpp::sourceCpp("/nesi/project/uoa04039/NewSourceCode_Faster/create_offspring_loop.cpp")
Rcpp::sourceCpp("/nesi/project/uoa04039/NewSourceCode_Faster/mutate_trait.cpp")
source("/nesi/project/uoa04039/NewSourceCode_Faster/corefunctions.R")

# =====================================================================
# LOAD JOB MAPPING
# =====================================================================

if (!file.exists(JOB_MAPPING_FILE)) {
  stop("Job mapping file not found: ", JOB_MAPPING_FILE, 
       "\nPlease run the setup script first to generate job mappings.")
}

job_mapping <- read_rds(JOB_MAPPING_FILE)
cat("Loaded job mapping with", nrow(job_mapping), "total jobs\n")

# Get this task's assignment
if (TASK_ID > nrow(job_mapping)) {
  stop("Task ID ", TASK_ID, " exceeds total number of jobs (", nrow(job_mapping), ")")
}

this_job <- job_mapping[TASK_ID, ]
COMBO_ID <- this_job$combo_id
REPLICATE <- this_job$replicate

cat("Task assignment:\n")
cat("- Combination ID:", COMBO_ID, "\n")
cat("- Replicate:", REPLICATE, "\n\n")

# =====================================================================
# LOAD COMBINATION PARAMETERS
# =====================================================================

if (!file.exists(COMBINATION_REGISTRY_FILE)) {
  stop("Combination registry not found: ", COMBINATION_REGISTRY_FILE)
}

registry <- read_rds(COMBINATION_REGISTRY_FILE)
combo_params <- registry[registry$combo_id == COMBO_ID, ]

if (nrow(combo_params) == 0) {
  stop("Combination ID ", COMBO_ID, " not found in registry")
}

cat("Combination parameters:\n")
print(combo_params)
cat("\n")

# =====================================================================
# CHECK IF OUTPUT ALREADY EXISTS
# =====================================================================

output_file <- file.path(OUTPUT_DIR, 
                         paste0("rep_", REPLICATE, "_simulation_iteration_", COMBO_ID, ".RDS"))

if (file.exists(output_file)) {
  cat("Output file already exists:", output_file, "\n")
  cat("Skipping simulation (set OVERWRITE=TRUE in setup script to re-run)\n")
  quit(save = "no", status = 0)
}

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

interpolate_vector <- function(vec, N) {
  new_length <- (length(vec) - 1) * (N + 1) + 1
  interpolated_data <- approx(x = 1:length(vec), y = vec, n = new_length)$y
  return(interpolated_data)
}

expand_vector <- function(nhostgen, nmicrogen, env_vec) {
  newvec <- interpolate_vector(env_vec, nmicrogen - 1)
  newvec <- newvec[1:(nmicrogen * nhostgen)]
  newvec <- matrix(newvec, nrow = nhostgen, ncol = nmicrogen, byrow = TRUE)
  return(newvec)
}

extract_field <- function(gen_obj, field, stat = "mean") {
  if (field == "env_fits") {
    vec <- gen_obj$env_fits[, 9]
  } else {
    vec <- gen_obj[[field]]
  }
  
  if (is.null(vec)) return(NA)
  vec[is.na(vec)] <- 1  
  
  if (stat == "mean") {
    return(mean(vec, na.rm = TRUE))
  } else if (stat == "var") {
    return(var(vec, na.rm = TRUE))
  } else {
    stop("Unknown stat type. Use 'mean' or 'var'.")
  }
}

# =====================================================================
# INITIALIZE POPULATIONS AND ENVIRONMENTS
# =====================================================================

cat("Initializing populations...\n")

# Initialize host-microbe population
InitPopulation <- matrix(nrow = N_Species, ncol = Host_PopSize)
InitPopulation <- apply(InitPopulation, MARGIN = 2, FUN = function(x) {
  rmultinom(n = 1, size = MicrobePopSize, prob = rep(1 / N_Species, N_Species))
})
rownames(InitPopulation) <- paste("Microbe", 1:N_Species, sep = "_")
colnames(InitPopulation) <- paste("Host", 1:Host_PopSize, sep = "_")

# Initialize environmental pool
EnvPool <- rep(EnvPoolSize / N_Species, N_Species)
names(EnvPool) <- paste("Microbe", 1:N_Species, sep = "_")
fixed_envpool <- EnvPool

# Initialize microbe traits
traitpool_microbes <- runif(N_Species, -2.5, 2.5)
names(traitpool_microbes) <- paste("Microbe", 1:N_Species, sep = "_")

# Load environmental conditions
cat("Loading environmental conditions...\n")
sec_varyring_envs_multigen <- read_rds(ENV_FILE_PATH)
sec_varyring_envs_multigen<-sec_varyring_envs_multigen[1:3]
if (!is.list(sec_varyring_envs_multigen)) {
  sec_varyring_envs_multigen <- list(as.matrix(sec_varyring_envs_multigen, ncol = 1))
  names(sec_varyring_envs_multigen) <- "env_1gen"
}

# Initialize genetic arrays
cat("Initializing genetic arrays...\n")
genetic_list <- list()
for (i in seq_along(LOCI_COUNTS)) {
  host_microbe_optima <- lapply(1:LOCI_COUNTS[i], function(x) {
    matrix(sample(c(-2.5, 2.5), 2 * Host_PopSize, replace = TRUE), 
           ncol = 2, nrow = Host_PopSize)
  })
  host_microbe_optima <- abind(host_microbe_optima, along = 3)
  individual_names <- paste("Host", 1:Host_PopSize, sep = "_")
  loci_names <- paste("Locus", 1:LOCI_COUNTS[i], sep = "_")
  allele_names <- c("Allele_1", "Allele_2")
  dimnames(host_microbe_optima) <- list(individual_names, allele_names, loci_names)
  genetic_list[[i]] <- host_microbe_optima
}
names(genetic_list) <- paste0("genetics_nloci_", LOCI_COUNTS)

# Initialize importance arrays
importances_list <- list()
for (i in seq_along(IMPORTANCE_VALUES)) {
  microbiome_importances <- lapply(1:1, function(x) {
    matrix(sample(c(IMPORTANCE_VALUES[i], IMPORTANCE_VALUES[i]), 
                  2 * Host_PopSize, replace = TRUE),
           ncol = 2, nrow = Host_PopSize)
  })
  microbiome_importances <- abind(microbiome_importances, along = 3)
  individual_names <- paste("Host", 1:Host_PopSize, sep = "_")
  loci_names <- paste("Locus", 1:1, sep = "_")
  allele_names <- c("Allele_1", "Allele_2")
  dimnames(microbiome_importances) <- list(individual_names, allele_names, loci_names)
  importances_list[[i]] <- microbiome_importances
}
names(importances_list) <- paste0("importances_init", IMPORTANCE_VALUES)

# Initialize weightings arrays
weightings_list <- list()
for (i in seq_along(WEIGHT_VALUES)) {
  microbiome_weightings <- lapply(1:50, function(x) {
    matrix(sample(c(WEIGHT_VALUES[i], WEIGHT_VALUES[i]), 
                  2 * Host_PopSize, replace = TRUE),
           ncol = 2, nrow = Host_PopSize)
  })
  microbiome_weightings <- abind(microbiome_weightings, along = 3)
  individual_names <- paste("Host", 1:Host_PopSize, sep = "_")
  loci_names <- paste("Locus", 1:50, sep = "_")
  allele_names <- c("Allele_1", "Allele_2")
  dimnames(microbiome_weightings) <- list(individual_names, allele_names, loci_names)
  weightings_list[[i]] <- microbiome_weightings
}
names(weightings_list) <- paste0("weightings_", WEIGHT_VALUES)

# =====================================================================
# RUN SIMULATION
# =====================================================================

cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("STARTING SIMULATION\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Combination:", COMBO_ID, "| Replicate:", REPLICATE, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

start_time <- Sys.time()

# Extract parameters
env_matrix <- sec_varyring_envs_multigen[[as.character(combo_params$EnvType)]]
host_microbe_optima <- genetic_list[[as.character(combo_params$Genetics)]]
microbiome_importances <- importances_list[[as.character(combo_params$Importance)]]
host_weightings_i <- weightings_list[[as.character(combo_params$Weighting)]]

loci_count_genetics <- as.numeric(gsub("genetics_nloci_", "", combo_params$Genetics))
mutation_value <- combo_params$MutationRate

# Convert to numeric vector for C++ function
params <- as.numeric(c(combo_params$X, combo_params$Z, combo_params$Y,
                       combo_params$S_host, combo_params$S_microbe, 
                       combo_params$S_env))

cat("Running simulation with parameters:\n")
cat("- Vertical inheritance (X):", params[1], "\n")
cat("- Env self-seeding (Z):", params[2], "\n")
cat("- Host contribution (Y):", params[3], "\n")
cat("- Selection strengths:", params[4], params[5], params[6], "\n")
cat("- Generations:", nrow(env_matrix), "\n")
cat("- Loci:", loci_count_genetics, "\n")
cat("- Mutation rate:", mutation_value, "\n\n")

# Run simulation
tryCatch({
  sim_result <- lapply_wrapper_CPP(
    XY = params,
    HostPopulation = InitPopulation,
    N_Microbes = MicrobePopSize,
    envpoolsize = EnvPoolSize,
    env_pool = EnvPool,
    generations = nrow(env_matrix),
    fixed_envpool = fixed_envpool,
    selection_parameter_hosts = params[4],
    selection_parameter_microbes = params[5],
    selection_parameter_env = params[6],
    traitpool_microbes = traitpool_microbes,
    host_microbe_optima = host_microbe_optima,
    env_cond_val = env_matrix,
    N_Species = N_Species,
    microbiome_importances = microbiome_importances,
    per_host_bac_gens = ncol(env_matrix),
    self_seed_prop = 0.98,
    generation_data_file = NA,
    mutation_rate_optima = mutation_value,
    mutation_rate_importance = mutation_value,
    print_currentgen = TRUE,
    importances_mutation = FALSE, 
    optima_mutate = TRUE,
    nloci = loci_count_genetics,
    nloci_importance = 1,
    host_env_weightings = host_weightings_i,
    lineage_track = FALSE,
    Test_HWE = TRUE
  )
  
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("SIMULATION COMPLETED SUCCESSFULLY\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Duration:", round(duration, 2), "minutes\n")
  cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  
  # Save result
  cat("\nSaving results to:", output_file, "\n")
  if (SAVE_COMPRESSED) {
    write_rds(sim_result, output_file, compress = "gz")
  } else {
    write_rds(sim_result, output_file)
  }
  
  file_size <- file.info(output_file)$size / 1024^2
  cat("File size:", round(file_size, 2), "MB\n")
  cat(paste(rep("=", 70), collapse=""), "\n\n")
  
  # Clean up memory
  rm(sim_result)
  gc()
  
  cat("Task completed successfully!\n")
  quit(save = "no", status = 0)
  
}, error = function(e) {
  cat("\n")
  cat(paste(rep("!", 70), collapse=""), "\n")
  cat("ERROR IN SIMULATION\n")
  cat(paste(rep("!", 70), collapse=""), "\n")
  cat("Error message:", e$message, "\n")
  cat("Traceback:\n")
  print(sys.calls())
  cat(paste(rep("!", 70), collapse=""), "\n\n")
  
  # Log error
  error_log <- file.path(OUTPUT_DIR, "error_log.txt")
  error_msg <- sprintf("[%s] Task %d (Combo %d, Rep %d): %s\n",
                       Sys.time(), TASK_ID, COMBO_ID, REPLICATE, e$message)
  cat(error_msg, file = error_log, append = TRUE)
  
  quit(save = "no", status = 1)
})
