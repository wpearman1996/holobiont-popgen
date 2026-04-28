#!/bin/bash
#SBATCH --job-name=extract_seff
#SBATCH --array=1-720
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=uoa04039
#SBATCH --time=01:00:00
#SBATCH --output=logs/seff_combo_%a.out
#SBATCH --error=logs/seff_combo_%a.err

mkdir -p logs
ml load R

Rscript - "${SLURM_ARRAY_TASK_ID}" << 'EOF'

library(dplyr)
library(readr)

args     <- commandArgs(trailingOnly = TRUE)
combo_id <- as.integer(args[1])
cat(sprintf("Processing combo_id: %d\n", combo_id))

OUTPUT_DIR  <- "./results_importance_mutation_standard_rerun/"
EXTRACT_DIR <- file.path(OUTPUT_DIR, "extracted_seff")
dir.create(EXTRACT_DIR, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(EXTRACT_DIR, paste0("combo_", combo_id, "_seff.RDS"))

REPS      <- 1:10
MAX_GEN   <- 3000
EARLY_MAX <- 1999
LATE_MIN  <- 2001

ENV_FILE <- "/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/env_list_rerun.RDS"

registry <- read_rds(file.path(OUTPUT_DIR, "combination_registry.RDS")) %>%
  mutate(
    selection_type = ifelse(S_host == Inf, "Neutral", "Selective"),
    importance_num = as.numeric(gsub("importances_init", "", Importance)),
    NGen           = as.numeric(gsub(".*gen(\\d+).*", "\\1", EnvType))
  ) %>%
  select(combo_id, selection_type, importance_num, NGen, EnvType, X, S_host)

combo_meta <- registry[registry$combo_id == combo_id, ]

env_list   <- read_rds(ENV_FILE)
env_matrix <- env_list[[as.character(combo_meta$EnvType)]]
n_micro_gens <- ncol(env_matrix)

summary_rows    <- list()
trajectory_rows <- list()


for (rep in REPS) {
  file_path <- file.path(OUTPUT_DIR,
                         paste0("rep_", rep, "_simulation_iteration_", combo_id, ".RDS"))

  tryCatch({
    dat    <- read_rds(file_path)[[1]]$GenData
    n_gens <- min(length(dat), MAX_GEN, nrow(env_matrix))

    gen_df <- do.call(rbind, lapply(1:n_gens, function(g) {
      hf  <- dat[[g]]$HostFitness
      nmf <- dat[[g]]$nomicrobiomehost_fitnessvector
      imp <- dat[[g]]$microbiome_importances_used
      G   <- dat[[g]]$host_phenotype_nomicrobiome  # pure genetic phenotype (I=0)
      E   <- env_matrix[g, n_micro_gens]

      s_eff_with <- s_eff_without <- NA_real_
      m_imp <- if (!is.null(imp) && length(imp) > 0) mean(imp, na.rm = TRUE) else NA_real_

      can_compute <- !is.null(G) && !is.null(hf) && !is.null(nmf) &&
                     length(G) > 2 && length(hf) == length(G)

      if (can_compute) {
        d2 <- (G - E)^2

        # Direct calculation: s_eff = -(G - E)² / ln(W)
        # Per host, then average across population
        # Exclude hosts where W ≈ 1 (d2 ≈ 0) or W ≈ 0 (ln → -Inf)
        valid_with <- hf > 1e-10 & hf < (1 - 1e-10) & d2 > 1e-10
        if (sum(valid_with) > 5) {
          s_eff_with <- mean(-d2[valid_with] / log(hf[valid_with]), na.rm = TRUE)
        }

        valid_without <- nmf > 1e-10 & nmf < (1 - 1e-10) & d2 > 1e-10
        if (sum(valid_without) > 5) {
          s_eff_without <- mean(-d2[valid_without] / log(nmf[valid_without]), na.rm = TRUE)
        }
      }

      c(generation    = g,
        mean_imp      = m_imp,
        s_eff_with    = s_eff_with,
        s_eff_without = s_eff_without)
    }))

    gen_df <- as.data.frame(gen_df)

    # Ratio: how much larger is effective s with microbiome vs without?
    gen_df$s_ratio <- gen_df$s_eff_with / gen_df$s_eff_without

    # Summary
    pm <- function(col, phase) {
      if (phase == "early") mean(gen_df[[col]][gen_df$generation <= EARLY_MAX], na.rm = TRUE)
      else                  mean(gen_df[[col]][gen_df$generation >= LATE_MIN],  na.rm = TRUE)
    }

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      combo_id            = combo_id,
      rep                 = rep,
      early_imp           = pm("mean_imp", "early"),
      evolved_imp         = pm("mean_imp", "late"),
      early_s_eff_with    = pm("s_eff_with", "early"),
      late_s_eff_with     = pm("s_eff_with", "late"),
      early_s_eff_without = pm("s_eff_without", "early"),
      late_s_eff_without  = pm("s_eff_without", "late"),
      early_s_ratio       = pm("s_ratio", "early"),
      late_s_ratio        = pm("s_ratio", "late")
    )

    trajectory_rows[[length(trajectory_rows) + 1]] <- gen_df %>%
      mutate(combo_id = combo_id, rep = rep)

    rm(dat); gc()
    cat(sprintf("  Combo %d rep %d: %d gens\n", combo_id, rep, n_gens))

  }, error = function(e) {
    cat(sprintf("  ERROR combo %d rep %d: %s\n", combo_id, rep, e$message))
  })
}

summary_df    <- bind_rows(summary_rows)    %>% left_join(combo_meta, by = "combo_id")
trajectory_df <- bind_rows(trajectory_rows) %>% left_join(combo_meta, by = "combo_id")

write_rds(list(summary = summary_df, trajectory = trajectory_df), out_file, compress = "gz")
cat(sprintf("Saved: %s\n", out_file))

EOF
