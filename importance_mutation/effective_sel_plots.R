#!/usr/bin/env Rscript
library(dplyr)
library(readr)

OUTPUT_DIR  <- "./Importance_Mutation//results_importance_mutation_standard_rerun/"
EXTRACT_DIR <- file.path(OUTPUT_DIR, "extracted_seff")

files <- list.files(EXTRACT_DIR, pattern = "combo_.*_seff\\.RDS$", full.names = TRUE)
cat(sprintf("Found %d files\n", length(files)))

all_summary    <- list()
all_trajectory <- list()

for (f in files) {
  tryCatch({
    d <- read_rds(f)
    all_summary[[length(all_summary) + 1]]       <- d$summary
    all_trajectory[[length(all_trajectory) + 1]] <- d$trajectory
  }, error = function(e) cat(sprintf("Error: %s\n", e$message)))
}

summary_df    <- bind_rows(all_summary)
trajectory_df <- bind_rows(all_trajectory)

write_rds(summary_df,    file.path(OUTPUT_DIR, "seff_summary.RDS"),    compress = "gz")
write_rds(trajectory_df, file.path(OUTPUT_DIR, "seff_trajectory.RDS"), compress = "gz")

cat(sprintf("Summary: %d rows | Trajectory: %d rows\n", nrow(summary_df), nrow(trajectory_df)))

# ---- Diagnostics ----
cat("\n=== EFFECTIVE s BY IMPORTANCE ===\n")
cat("Both use (G - E)² as predictor; difference is fitness with vs without microbiome\n")
cat("s_eff_without should ≈ nominal s (0.1) — sanity check\n")
cat("s_eff_with should increase with importance — microbiome buffers genetics\n")
cat("s_ratio > 1 = microbiome weakens selection on genetics\n\n")

summary_df %>%
  group_by(importance_num) %>%
  summarise(
    nominal_s         = mean(as.numeric(S_host), na.rm = TRUE),
    evol_s         = mean(as.numeric(S_host), na.rm = TRUE),
    early_s_with      = mean(early_s_eff_with, na.rm = TRUE),
    early_s_without   = mean(early_s_eff_without, na.rm = TRUE),
    late_s_with       = mean(late_s_eff_with, na.rm = TRUE),
    late_s_without    = mean(late_s_eff_without, na.rm = TRUE),
    early_s_ratio     = mean(early_s_ratio, na.rm = TRUE),
    late_s_ratio      = mean(late_s_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(importance_num) %>%
  print(n = 30)

cat("\n=== EFFECTIVE s BY IMPORTANCE x VERTICAL INHERITANCE ===\n")
summary_df %>%
  group_by(importance_num, X) %>%
  summarise(
    early_s_with    = mean(early_s_eff_with, na.rm = TRUE),
    early_s_without = mean(early_s_eff_without, na.rm = TRUE),
    early_s_ratio   = mean(early_s_ratio, na.rm = TRUE),
    late_s_ratio    = mean(late_s_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(importance_num, X) %>%
  print(n = 80)



traj_sum<-trajectory_df %>%
  filter(generation %in% 2500:3000) %>%
  group_by(importance_num,selection_type,EnvType,X,rep) %>%
  summarise(
    nominal_s         = mean(as.numeric(S_host), na.rm = TRUE),
    evol_imp         = mean(as.numeric(mean_imp), na.rm = TRUE),
    s_with      = mean(s_eff_with, na.rm = TRUE),
    s_without   = mean(s_eff_without, na.rm = TRUE),
    s_ratio      = mean(s_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(importance_num) %>%
  mutate(X = as.character(X))



traj_sum %>%
  filter(selection_type %in% "Selective") %>%
  filter(EnvType %in% c("step_gen1", "step_gen5", 
                        "step_gen50")) %>%
  filter(!(X == "1" & EnvType == "step_gen1")) %>%
  ggplot() +
  geom_smooth(method="loess",se = F)+
  aes(x = evol_imp, y = s_with) +
  geom_point() +
  scale_colour_brewer(palette = "Spectral") +
  theme_minimal() +
  facet_wrap(vars(EnvType))



traj_sum %>%
  filter(selection_type == "Selective") %>%
  filter(EnvType %in% c("step_gen1", "step_gen5", "step_gen50")) %>%
  filter(!(X == "1" & EnvType == "step_gen1")) %>%
  ggplot(aes(x = evol_imp, y = s_with)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_colour_brewer(palette = "Spectral", name = "Vertical\nInheritance (X)") +
  labs(x = "Evolved Importance", y = "Effective Selection (with microbiome)") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  ) +
  facet_wrap(vars(EnvType))

make_eff_plot <- function(df) {
  ggplot(df, aes(x = evol_imp, y = s_with)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, colour = "black") +
    labs(x = "Mean evolved importance (last 500 gens)", y = "Effective Genotypic Selection") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

eff_data <- traj_sum %>%
  filter(selection_type == "Selective") %>%
  filter(EnvType %in% c("step_gen1", "step_gen5", "step_gen50")) %>%
  filter(!(X == "1" & EnvType == "step_gen1"))

eff_1 <- make_eff_plot(subset(eff_data, EnvType == "step_gen1"))
eff_2 <- make_eff_plot(subset(eff_data, EnvType == "step_gen5"))
eff_3 <- make_eff_plot(subset(eff_data, EnvType == "step_gen50"))

eff_sel_comb<-(eff_1 + eff_2 + eff_3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("eff_sel_plot.png",
       eff_sel_comb, width = 18, height = 7, dpi = 300)
cat("Saved: 04e_late_varratio_and_rescue.png\n")