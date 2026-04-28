library(dplyr)
library(readr)
library(ggplot2)

OUTPUT_DIR  <- "/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/Importance_Mutation//results_importance_mutation_standard_rerun/"
EXTRACT_DIR <- file.path(OUTPUT_DIR, "extracted_varratio")
PLOT_DIR    <- file.path(OUTPUT_DIR, "addiction_variance_plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

# =====================================================================
# LOAD AND COMBINE ALL COMBO FILES
# =====================================================================

combo_files <- list.files(EXTRACT_DIR, pattern = "combo_.*_varratio\\.RDS",
                          full.names = TRUE)
cat(sprintf("Found %d combo files\n", length(combo_files)))

all_data <- lapply(combo_files, read_rds)

var_summary    <- bind_rows(lapply(all_data, `[[`, "summary"))    %>%
  filter(is.finite(evolved_imp), is.finite(delta_var_ratio)) %>%
  mutate(NGen_label = paste0("NGen = ", NGen)) %>%
  filter(EnvType %in% c("step_gen1","step_gen5","step_gen50"))

var_trajectory <- bind_rows(lapply(all_data, `[[`, "trajectory")) %>%
  filter(is.finite(var_ratio)) %>%
  mutate(NGen_label = paste0("NGen = ", NGen)) %>%
  filter(EnvType %in% c("step_gen1","step_gen5","step_gen50"))

cat(sprintf("Summary rows: %d | Trajectory rows: %d\n",
            nrow(var_summary), nrow(var_trajectory)))




library(patchwork)

# Compute late-generation summaries from trajectory
late_ratio <- var_trajectory %>%
  filter(importance_num == 0) %>%
  group_by(combo_id, rep, NGen, X, EnvType, selection_type) %>%
  filter(!(X == 1 & EnvType == "step_gen1")) %>%
  filter(selection_type =="Selective") %>%
  filter(generation > 2499) %>%
  summarise(
    late_var_ratio  = mean(var_ratio[generation >= max(generation) - 500],  na.rm = TRUE),
    evolved_imp     = mean(mean_imp[generation  >= max(generation) - 500],  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(late_var_ratio), is.finite(evolved_imp)) %>%
  mutate(NGen_label = paste0("NGen = ", NGen))

# Compute late rescue effect from summary (host fitness - nomicrobiome fitness)
late_rescue <- var_summary %>%
  filter(importance_num == 0) %>%
  filter(selection_type =="Selective") %>%
  filter(!(X == 1 & EnvType == "step_gen1")) %>%
  mutate(
    rescue_effect = late_fitness - late_nomicro_fitness,
    NGen_label    = paste0("NGen = ", NGen)
  ) %>%
  filter(is.finite(rescue_effect), is.finite(evolved_imp))

imp_traj <- var_trajectory %>%
  filter(importance_num == 0) %>%
  filter(selection_type == "Selective") %>%
  group_by(generation, X, EnvType, selection_type, NGen_label) %>%
  summarise(
    mean_importance = mean(mean_imp, na.rm = TRUE),
    mean_rescue_effect = mean(rescue_effect, na.rm = TRUE),
    
    .groups = "drop"
  )

# Shared aesthetics
sel_colours <- c("Selective" = "#E63946", "Neutral" = "#888888")
base_theme  <- theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        panel.grid.minor = element_blank())

p_late <- ggplot(late_ratio,
                 aes(x = evolved_imp, y = late_var_ratio,group=NGen_label,shape=NGen_label,col=NGen_label)) +
#  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  labs(
  #  title = "Variance ratio vs evolved importance",
    x     = "Mean evolved importance (last 500 gens)",
    y     = "Standard deviation ratio"
  ) +
  base_theme

p_rescue <- ggplot(late_rescue[late_rescue$EnvType=="step_gen50",],
                   aes(x = evolved_imp, y = rescue_effect)) +
 # geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9,col="black") +
#  facet_wrap(~ EnvType,scales="free") +
  labs(
  #  title = "Rescue effect vs evolved importance",
    x     = "Mean evolved importance (last 500 gens)",
    y     = "Microbiome Rescue Effect"
  ) +
  base_theme


imp_traj_summary <- var_trajectory %>%
  filter(importance_num == 0) %>%
  filter(selection_type == "Selective") %>%
  group_by(generation, X, EnvType, selection_type, NGen_label) %>%
  summarise(
    mean_importance = mean(mean_imp, na.rm = TRUE),
    sd_importance = sd(mean_imp, na.rm = TRUE),
    se_importance = sd(mean_imp, na.rm = TRUE) / sqrt(n()),
    quant_importance = quantile(mean_imp,0.8, na.rm = TRUE), 
    n_reps = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_sd = mean_importance - sd_importance,
    upper_sd = mean_importance + sd_importance,
    lower_se = mean_importance - se_importance,
    upper_se = mean_importance + se_importance,
    lower_quant = pmax(0, mean_importance - quant_importance),
    upper_quant = mean_importance + quant_importance
    
  )


imp_traj_summary

p_imp_traj <-imp_traj_summary %>%
  filter (EnvType == "step_gen50") %>%ggplot() +
                     aes(x = generation, y = mean_importance,
                         colour = factor(X), fill = factor(X), group = factor(X)) +
  geom_ribbon(aes(ymin = lower_se, ymax = upper_se), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_brewer(palette = "Spectral", name = expression(Vertical ~ Inheritance ~ (italic(X)))) +
  scale_fill_brewer(palette = "Spectral", name = expression(Vertical ~ Inheritance ~ (italic(X)))) +
  #facet_grid( ~ EnvType) +
  labs(
    x = "Generation",
    y = "Mean Microbial Importance") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )



# 1. Update the legends and labels in the individual plots
library(ggplot2)
library(patchwork)
library(tagger)
# 1. Update Labels and Legend for top plots
# Using & to ensure changes apply across the combined object later
p_top <- (p_late + p_rescue) + 
  plot_layout(guides = "collect") & 
  scale_color_discrete(name = "Microbial Generation Count", labels = c("1", "5", "50")) &
  scale_shape_discrete(name = "Microbial Generation Count", labels = c("1", "5", "50"))
# Assuming your EnvType levels are "Type1", "Type2", "Type3"
# Replace these strings with your actual factor levels
traj_1 <- subset(imp_traj, EnvType == levels(factor(EnvType))[1])
traj_2 <- subset(imp_traj, EnvType == levels(factor(EnvType))[2])
traj_3 <- subset(imp_traj, EnvType == levels(factor(EnvType))[3])
make_traj_plot <- function(df) {
  ggplot(df, aes(x = generation, y = mean_importance, 
                 colour = factor(X), group = factor(X))) +
    geom_line(linewidth = 0.7, alpha = 1) +
    scale_colour_brewer(palette = "Spectral", name = "Vertical\nInheritance (X)") +
    labs(x = "Generation", y = "Mean Microbial Importance") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    ) +
    geom_vline(xintercept=2000,linetype = "dashed",colour = "grey")
}

  
final_plot <- (p_imp_traj + geom_vline(xintercept = 2000, linetype = "dashed", colour = "grey") |
                                  p_rescue +
                                 scale_color_discrete(name = "Microbial Generation Count", labels = c("1", "5", "50")) +
                                   scale_shape_discrete(name = "Microbial Generation Count", labels = c("1", "5", "50")) |
                                 eff_3) +
      plot_layout(guides = "collect") +
       plot_annotation(tag_levels = 'a') &
      theme(plot.tag = element_text(face = "bold"),
                  legend.position = "left")
# final_plot <- (p_imp_traj +  geom_vline(xintercept=2000,linetype = "dashed",colour = "grey") | p_rescue + 
#                  scale_color_discrete(name = "Microbial Generation Count", labels = c("1", "5", "50")) +
#                  scale_shape_discrete(name = "Microbial Generation Count", labels = c("1", "5", "50"))| eff_3) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = 'a') & 
#   theme(plot.tag = element_text(face = "bold"))

final_plot
ggsave(file.path("full_combined_plot.png"),
       final_plot, width = 15, height = 6, dpi = 1000)


m_var <- lm((late_var_ratio) ~ evolved_imp * X + EnvType + rep,
            data = late_ratio[late_ratio$selection_type =="Selective",])
summary(m_var)
anova(m_var)

m_rescue <- lm(rescue_effect ~ evolved_imp * X  + EnvType + rep,
               data = late_rescue[late_rescue$selection_type =="Selective",])
summary(m_rescue)


anova(m_rescue)



ggsave(file.path(PLOT_DIR, "04e_late_varratio_and_rescue.png"),
       p_combined, width = 18, height = 10, dpi = 300)
cat("Saved: 04e_late_varratio_and_rescue.png\n")





library(broom)
library(officer)

# Tidy your model results
model_tidy <- tidy(m_rescue)

# Create and format the table
library(knitr)
library(broom)

# This creates a simple table you can copy-paste into Word
kable(tidy(your_model_name), format = "pipe", digits = 3)






importance_summary_table <- late_ratio %>%
  group_by(X, NGen) %>%
  summarise(
    mean_evolved_importance = mean(evolved_imp, na.rm = TRUE),
    sd_evolved_importance   = sd(evolved_imp, na.rm = TRUE),
    n_replicates            = n(),
    .groups = "drop"
  ) %>%
  # Optional: Rename NGen for clarity in the final table
  mutate(Gen_Count = paste0(NGen, " generations")) %>%
  select(X, Gen_Count, mean_evolved_importance, sd_evolved_importance, n_replicates)








library(ggh4x)  # You may need to install this: install.packages("ggh4x")

supp_p_imp_traj <- imp_traj_summary %>%
  ggplot() +
  aes(x = generation, y = mean_importance,
      colour = factor(X), fill = factor(X), group = factor(X)) +
  geom_ribbon(aes(ymin = lower_se, ymax = upper_se), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_brewer(palette = "Spectral", name = "Vertical\nInheritance (X)") +
  scale_fill_brewer(palette = "Spectral", name = "Vertical\nInheritance (X)") +
  facet_grid(~ EnvType) +
  labs(
    x = "Generation",
    y = "Mean Microbial Importance"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Add alphabetical labels to facets
supp_p_imp_traj <- supp_p_imp_traj + 
  ggh4x::facetted_pos_scales(x = list(a = scale_x_continuous())) +
  ggh4x::force_panelsizes(rows = NULL, cols = NULL)

supp_p_imp_traj <- egg::tag_facet(supp_p_imp_traj, tag_pool = letters)
ggsave("supp_supp_p_imp_traj.png", supp_p_imp_traj, width = 10, height = 4,dpi=800)
