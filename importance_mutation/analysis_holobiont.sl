#!/bin/bash
#SBATCH --job-name=extract_varratio
#SBATCH --array=1-720
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=uoa04039
#SBATCH --time=01:00:00
#SBATCH --output=logs/varratio_combo_%a.out
#SBATCH --error=logs/varratio_combo_%a.err

mkdir -p logs
ml load R

Rscript - "${SLURM_ARRAY_TASK_ID}" << 'EOF'

library(dplyr)
library(readr)

args     <- commandArgs(trailingOnly = TRUE)
combo_id <- as.integer(args[1])
cat(sprintf("Processing combo_id: %d\n", combo_id))

OUTPUT_DIR  <- "./results_importance_mutation_standard_rerun/"
EXTRACT_DIR <- file.path(OUTPUT_DIR, "extracted_varratio")
dir.create(EXTRACT_DIR, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(EXTRACT_DIR, paste0("combo_", combo_id, "_varratio.RDS"))

REPS      <- 1:10
MAX_GEN   <- 3000
EARLY_MAX <- 1999
LATE_MIN  <- 2001

registry <- read_rds(file.path(OUTPUT_DIR, "combination_registry.RDS")) %>%
  mutate(
    selection_type = ifelse(S_host == Inf, "Neutral", "Selective"),
    importance_num = as.numeric(gsub("importances_init", "", Importance)),
    NGen           = as.numeric(gsub(".*gen(\\d+).*", "\\1", EnvType))
  ) %>%
  select(combo_id, selection_type, importance_num, NGen, EnvType, X)

combo_meta <- registry[registry$combo_id == combo_id, ]

summary_rows    <- list()
trajectory_rows <- list()

for (rep in REPS) {
  file_path <- file.path(OUTPUT_DIR,
                         paste0("rep_", rep, "_simulation_iteration_", combo_id, ".RDS"))

  tryCatch({
    dat    <- read_rds(file_path)[[1]]$GenData
    n_gens <- min(length(dat), MAX_GEN)

    # Full per-generation extraction
    gen_df <- do.call(rbind, lapply(1:n_gens, function(g) {
      hf  <- dat[[g]]$HostFitness
      nmf <- dat[[g]]$nomicrobiomehost_fitnessvector
      imp <- dat[[g]]$microbiome_importances_used

      c(generation   = g,
        sd_actual    = if (!is.null(hf)  && length(hf)  > 1) sd(hf,  na.rm = TRUE) else NA_real_,
        sd_nomicro   = if (!is.null(nmf) && length(nmf) > 1) sd(nmf, na.rm = TRUE) else NA_real_,
        mean_imp     = if (!is.null(imp) && length(imp) > 0) mean(imp, na.rm = TRUE) else NA_real_,
        mean_fitness = if (!is.null(hf)  && length(hf)  > 0) mean(hf,  na.rm = TRUE) else NA_real_,
        mean_nomicro = if (!is.null(nmf) && length(nmf) > 0) mean(nmf, na.rm = TRUE) else NA_real_)
    }))

    gen_df                <- as.data.frame(gen_df)
    gen_df$var_ratio      <- gen_df$sd_nomicro / (gen_df$sd_actual + 1e-10)
    gen_df$rescue_effect  <- gen_df$mean_fitness - gen_df$mean_nomicro

    # Summary row
    early_vr <- mean(gen_df$var_ratio[gen_df$generation <= EARLY_MAX], na.rm = TRUE)
    late_vr  <- mean(gen_df$var_ratio[gen_df$generation >= LATE_MIN],  na.rm = TRUE)

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      combo_id              = combo_id,
      rep                   = rep,
      early_var_ratio       = early_vr,
      late_var_ratio        = late_vr,
      delta_var_ratio       = late_vr - early_vr,
      evolved_imp           = mean(gen_df$mean_imp[gen_df$generation      >= LATE_MIN],  na.rm = TRUE),
      early_imp             = mean(gen_df$mean_imp[gen_df$generation      <= EARLY_MAX], na.rm = TRUE),
      late_fitness          = mean(gen_df$mean_fitness[gen_df$generation  >= LATE_MIN],  na.rm = TRUE),
      early_fitness         = mean(gen_df$mean_fitness[gen_df$generation  <= EARLY_MAX], na.rm = TRUE),
      late_nomicro_fitness  = mean(gen_df$mean_nomicro[gen_df$generation  >= LATE_MIN],  na.rm = TRUE),
      early_nomicro_fitness = mean(gen_df$mean_nomicro[gen_df$generation  <= EARLY_MAX], na.rm = TRUE),
      late_rescue           = mean(gen_df$rescue_effect[gen_df$generation >= LATE_MIN],  na.rm = TRUE),
      early_rescue          = mean(gen_df$rescue_effect[gen_df$generation <= EARLY_MAX], na.rm = TRUE)
    )

    # Full trajectory
    trajectory_rows[[length(trajectory_rows) + 1]] <- data.frame(
      combo_id      = combo_id,
      rep           = rep,
      generation    = gen_df$generation,
      var_ratio     = gen_df$var_ratio,
      sd_actual     = gen_df$sd_actual,
      sd_nomicro    = gen_df$sd_nomicro,
      mean_imp      = gen_df$mean_imp,
      mean_fitness  = gen_df$mean_fitness,
      mean_nomicro  = gen_df$mean_nomicro,
      rescue_effect = gen_df$rescue_effect
    )

    rm(dat); gc()
    cat(sprintf("  Combo %d rep %d: %d gens extracted\n", combo_id, rep, n_gens))

  }, error = function(e) {
    cat(sprintf("  ERROR combo %d rep %d: %s\n", combo_id, rep, e$message))
  })
}

# Attach metadata and save as a list with both summary and trajectory
summary_df    <- bind_rows(summary_rows)    %>% left_join(combo_meta, by = "combo_id")
trajectory_df <- bind_rows(trajectory_rows) %>% left_join(combo_meta, by = "combo_id")

out <- list(summary = summary_df, trajectory = trajectory_df)
write_rds(out, out_file, compress = "gz")
cat(sprintf("Saved: %s (summary: %d rows | trajectory: %d rows)\n",
            out_file, nrow(summary_df), nrow(trajectory_df)))

EOF