library(dplyr)
setwd("/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/")
RERUN_DATA_EXTRACTION = FALSE

library(data.table)
REGISTRY_FILE <- readr::read_rds("./results_rerun_standard_11Feb/combination_registry.RDS")
sim_files <- list.files(path = "./results_rerun_standard_11Feb//",pattern="*RDS")
sim_files <- sim_files[grepl("simulation_iteration",sim_files)]
N_REPS <- unique(stringr::word(sim_files,2,sep="_"))

clone_reg_file<-function(REGISTRY,N_REPS){
  reg_clone_list<-list()
  for(i in seq_along(N_REPS)){
    reg_clone_list[[i]] <- REGISTRY
    reg_clone_list[[i]]$Replicate_Number <- i
  }
  rbindlist(reg_clone_list)
}
cloned_file<-clone_reg_file(REGISTRY_FILE,N_REPS)
cloned_file$FileName<-paste0("rep_",cloned_file$Replicate_Number,"_simulation_iteration_",cloned_file$combo_id,".RDS")
#cloned_file<-cloned_file[cloned_file$EnvType %in% c("step_gen1","step_gen5","step_gen50"),]
calc_div <- function(population) {
  relabund <- population / rowSums(population,na.rm = T)
  -rowSums(relabund * log(relabund),na.rm = T)
}
calc_rich <- function(population) {
  population[population > 0] <- TRUE
  colSums(population)
}
gens_of_int <- 2000:3000
#filename<-cloned_file$FileName[5]
output_dir="./results_rerun_standard_11Feb//"
process_raw_sim_file<-function(filename,gens_of_int,output_dir){
  filename_output<-paste0(output_dir,gsub(".RDS","_individ_data.RDS",filename))
  if(file.exists(filename_output)){
    return(NULL)
  }
  x_sub<-readr::read_rds(paste0(output_dir,filename))
  x_sub <- x_sub[[1]]$GenData[gens_of_int]
  gendata<-lapply(x_sub,function(x){
    df <- data.frame(
      HostAlphaDiversity = x$Diversity,
      HostGenotype = x$host_microbe_optima_used,
      HostFitness = x$HostFitness,
      AvgMicrobialFitness = x$MicrobeFitness,
      HostPhenotype = x$host_phenotype,
      SpeciesRichness = x$species_richess,
      EnvDiv = calc_div(t(as.matrix(x$env_used))) / log(calc_rich((as.matrix(x$env_used)))),
      NumUniqueParents = length(unique(x$num_unique_parents)),
      BrayDiv = mean(x$BrayDiv),
      ObservedHeterozygosity = x$optima_Observed_heterozygosity,
      mean_microbial_trait_vals = x$mean_microbial_trait_vals,
      Ne = x$Ne,
      AverageNumAllelesPerLocus = mean(x$optima_alleles_per_locus),
      NumPolymorphicLoci = length(x$optima_alleles_per_locus[x$optima_alleles_per_locus > 1]),
      stringsAsFactors = FALSE
    )
    
  })
  readr::write_rds(gendata,filename_output)
  rm(gendata)
  gc()
  return(NULL)
}

if(RERUN_DATA_EXTRACTION){
  N_CORES<-128
  
  library(parallel)
  mclapply(cloned_file$FileName, 
           function(f) process_raw_sim_file(f, 1800:2600, "./results_indepth_reverse_updown//"),
           mc.cores = N_CORES)
}
# 
cloned_file$OutputFileName <-paste0(output_dir,gsub(".RDS","_individ_data.RDS",cloned_file$FileName))
files_genned<-list.files(path="./results_rerun_standard_11Feb/",pattern="*_individ_data.RDS")
raw_files_genned<-list.files(path="./results_rerun_standard_11Feb/",pattern="*RDS")
raw_files_genned<-raw_files_genned[grepl("rep_",raw_files_genned)]
raw_files_genned<-raw_files_genned[!grepl("_individ_data.RDS",raw_files_genned)]

files_genned<-paste0("./results_rerun_standard_11Feb//",files_genned)
cloned_file$FileComplete<-cloned_file$OutputFileName %in% files_genned
missing_indices <- which(!cloned_file$OutputFileName %in% files_genned)
clonedfile_nonrun<-cloned_file[!cloned_file$OutputFileName %in% files_genned,]

cloned_file$OutputFileName[!cloned_file$OutputFileName %in% files_genned]
job_mapping<-read_rds("./results_rerun_standard_11Feb/job_mapping.RDS")
REGISTRY_FILE[REGISTRY_FILE$combo_hash %in% clonedfile_nonrun$combo_hash,]

# 
# indiv_data<-lapply(files_genned[1:5],read_rds)
# write_rds(indiv_data,"indiv_data.RDS")
# indiv_data <- lapply(indiv_data,function(x){
#  rbindlist(x[99:501])
# })
# write_rds(indiv_data,"indiv_data_gensub.RDS")


#indiv_data<-readr::read_rds("indiv_data.RDS")
# Keep generations 1900-2100
#gen_range <- 1900:2100
#indiv_data<-readr::read_rds("indiv_data_gensub.RDS")
# indiv_data<-lapply(indiv_data,function(x){
#   rbindlist(x)
# })
selected_gen_indices <- 99:701
gens_of_int<-1800:2600
selected_generations <- gens_of_int[selected_gen_indices]  # 

# for(i in seq_along(indiv_data_gensub)){
#   registry_row <- cloned_file[i, ]
#   indiv_data_gensub[[i]]<-cbind(indiv_data_gensub[[i]],registry_row)
#   indiv_data_gensub[[i]]$Generation<-sort(rep(selected_generations,200))
# }

#indiv_data[[1]]$Generation
#write_rds(indiv_data,"indiv_data_withgen_subs.RDS")
# 
# 
# #indiv_data<-read_rds("indiv_data_withgen_subs.RDS")
# 
# # Now combine with registry (also fixed the gen_df variable name error)
# combined_data <- rbindlist(indiv_data)
# colnames(combined_data_all)
# write_rds(combined_data,"combined_data.RDS")
#combined_data_all <- readRDS("/nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/combined_data.RDS")
# combined_data_all$ImportanceNumeric = gsub("importances_init","",combined_data_all$Importance)
# combined_data_all$ImportanceNumeric<-(as.numeric(combined_data_all$ImportanceNumeric))*2
# combined_data_all$ImportanceLabel<-paste0("Importance=",combined_data_all$ImportanceNumeric)
# #####

#combined_data_all<- combined_data_all %>% filter(Generation %in% c(1900,2100)) %>%
#  filter(EnvType == "step_gen5")

aggregate_replicates <- function(data) {
  # Convert to tibble if it's a data.table
  data <- as_tibble(data)
  
  # Define grouping variables
  grouping_vars <- c("Generation", "EnvType", "Importance", 
                     "X")
  
  # Add registry columns if registry exists
  if (exists("registry")) {
    grouping_vars <- unique(c(colnames(registry), grouping_vars))
  }
  
  # Get grouping columns that actually exist in the data
  grouping_cols <- intersect(grouping_vars, names(data))
  
  # Get numeric columns (excluding grouping columns)
  remaining_cols <- setdiff(names(data), grouping_cols)
  numeric_cols <- remaining_cols[sapply(data[remaining_cols], is.numeric)]
  
  # Aggregate the data
  result <- data %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      # Calculate means for all numeric columns
      across(all_of(numeric_cols), 
             list(mean = ~ mean(.x, na.rm = TRUE),
                  sd = ~ sd(.x, na.rm = TRUE)), 
             .names = "{.col}_{.fn}"),
      
      # Count replicates
      n_replicates = n(),
      .groups = "drop"
    )
  
  return(result)
}
combined_data_all_aggregated <- aggregate_replicates(combined_data_all)
#write_rds(combined_data_all_aggregated,"combined_data_all_aggregated.RDS")
combined_data_all_aggregated<-readr::read_rds("combined_data_all_aggregated.RDS")
table(combined_data_all_aggregated$EnvType)
combined_data_all_aggregated <- combined_data_all_aggregated %>%
  filter(EnvType %in% c("step_gen1","step_gen5","step_gen50"))
combined_data_all_aggregated$X_Fac<-as.factor(combined_data_all_aggregated$X)
combined_data_all_aggregated$GenFac<-as.factor(combined_data_all_aggregated$Generation)
combined_data_all_aggregated$ImportanceNumeric = gsub("importances_init","",combined_data_all_aggregated$Importance)
combined_data_all_aggregated$ImportanceNumeric<-(as.numeric(combined_data_all_aggregated$ImportanceNumeric))*2
combined_data_all_aggregated$ImportanceLabel<-paste0("Importance=",combined_data_all_aggregated$ImportanceNumeric)
combined_data_all_aggregated$X_label = ifelse(combined_data_all_aggregated$X == 0, "No Inheritance (X=0)",
                                              ifelse(combined_data_all_aggregated$X == 0.1, "Very Low Inheritance (X=0.1)",
                                                     ifelse(combined_data_all_aggregated$X == 0.2, "Low Inheritance (X=0.2)",
                                                            ifelse(combined_data_all_aggregated$X == 0.3, "Low-Moderate Inheritance (X=0.3)",
                                                                   ifelse(combined_data_all_aggregated$X == 0.4, "Moderate Inheritance (X=0.4)",
                                                                          ifelse(combined_data_all_aggregated$X == 0.5, "Moderate Inheritance (X=0.5)",
                                                                                 ifelse(combined_data_all_aggregated$X == 0.9, "High Inheritance (X=0.9)",
                                                                                        ifelse(combined_data_all_aggregated$X == 1, "Complete Inheritance (X=1)",
                                                                                               paste("X =", combined_data_all_aggregated$X)))))))))






















library(ggplot2)
library(dplyr)
library(patchwork)


# ============================================================================
# THEME (keeping only what's needed)
# ============================================================================

presentation_theme <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0, margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, color = "gray30", hjust = 0, margin = margin(b = 15)),
      axis.title = element_text(size = base_size, face = "bold", color = "black"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 15),
      legend.key.size = unit(1.2, "lines"),
      strip.background = element_rect(fill = "gray95", color = "gray80", size = 0.5),
      strip.text = element_text(size = base_size, face = "bold", color = "black", margin = margin(t = 5, b = 5)),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_line(color = "gray95", size = 0.2),
      panel.spacing = unit(1, "lines"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    )
}
get_env_value <- function(generation, env_type) {
  env_matrix <- env_list[[env_type]]
  if (is.matrix(env_matrix)) {
    env_values <- env_matrix[, 1]  # First microbial generation of each host generation
  } else {
    env_values <- env_matrix
  }
  if (generation <= length(env_values)) {
    return(env_values[generation])
  } else {
    return(NA)
  }
}

# ============================================================================
# PREPARE DATA
# ============================================================================

# Load aggregated data

# Filter and prepare heatmap data
gens_of_int <- 1900:3000
label_gens <- gens_of_int[seq(1, length(gens_of_int), by = 50)]

heatmap_data <- combined_data_all_aggregated %>%
  #  filter(Generation %in% gens_of_int,
  #         EnvType %in% c("step_gen1", "step_gen5", "step_gen50")) %>%
  mutate(
    # Ensure proper factor ordering
    ImportanceLabel = factor(ImportanceLabel, 
                             levels = unique(ImportanceLabel)[order(as.numeric(gsub(".*=([0-9.]+).*", "\\1", unique(ImportanceLabel))))]),
    ImportanceLabel = factor(stringr::str_wrap(ImportanceLabel, width = 20),
                             levels = stringr::str_wrap(levels(ImportanceLabel), width = 20)),
    X_label = factor(X_label, 
                     levels = unique(X_label)[order(as.numeric(gsub(".*=([0-9.]+).*", "\\1", unique(X_label))))])
  )

# Set faceting variable
facet_by_x <- FALSE
if (facet_by_x) {
  facet_formula <- facet_grid(vars(X_label), vars(EnvType))
  y_var <- sym("ImportanceLabel")
  y_label <- "Microbiome Importance"
} else {
  facet_formula <- facet_grid(vars(ImportanceLabel), vars(EnvType))
  y_var <- sym("X_label")
  y_label <- "X Value"
}
heatmap_data <- heatmap_data %>% filter(Generation %in% 1900:3000)
# Create dashed lines data for environmental changes
dashed_lines_data <- heatmap_data %>%
  filter(EnvType == "step_gen50") %>%
  dplyr::select(EnvType) %>%
  distinct() %>%
  mutate(lines = list(2000)) %>%
  tidyr::unnest(lines) %>%
  dplyr::rename(generation_line = lines)

# ============================================================================
# NOTE: ADD ENVIRONMENTAL VALUES
# You'll need to add env_value to your data, e.g.:
env_list <- readRDS("env_list_rerun.RDS")
env_list <- env_list[1:3]

heatmap_data <- heatmap_data %>%
  rowwise() %>%
  mutate(env_value = get_env_value(Generation, EnvType)) %>%
  ungroup() %>%
  filter(ImportanceNumeric %in% c(0,0.02,0.1,0.2,0.52,0.8,1))
# ============================================================================

# Create environmental plot (assuming you've added env_value)
# If env_value is not available, you can skip the envplot or use EnvDiv_mean instead

envplot <- heatmap_data %>% 
  filter(ImportanceNumeric == 0, X == 0) %>%
  ggplot() +
  aes(x = factor(Generation), y = env_value) +  # Using EnvDiv as proxy
  geom_line(aes(group = 1), colour = "#112446", linewidth = 0.8) +
  geom_vline(data = dashed_lines_data,
             aes(xintercept = match(generation_line, sort(unique(heatmap_data$Generation)))),
             linetype = "dashed", color = "black", alpha = 0.7) +
  scale_x_discrete(breaks = label_gens, labels = label_gens, expand = c(0, 0)) +
  # facet_wrap(vars(EnvType), ncol = 6) +
  labs(x = NULL, y = "Environment") +
  presentation_theme() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_blank(),  # Changed from element_text(angle = 0)
    strip.background.x = element_blank(),  # Added to remove strip background
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    legend.position = "right"
  )

# ============================================================================
# CREATE HEATMAP FUNCTION
# ============================================================================

create_heatmap <- function(data, fill_var, fill_name, title = NULL, label_facets = TRUE) {
  
  # Determine faceting variables based on facet_by_x
  if (facet_by_x) {
    row_var <- "X_label"
    col_var <- "EnvType"
  } else {
    row_var <- "ImportanceLabel"
    col_var <- "EnvType"
  }
  
  # Create facet labels if requested
  if (label_facets) {
    facet_labels <- data %>%
      distinct(across(all_of(c(row_var, col_var)))) %>%
      arrange(across(all_of(c(row_var, col_var)))) %>%
      mutate(
        label = letters[row_number()],
        Generation = min(data$Generation),
        # Use the first level of the y-axis factor for positioning
        !!as_label(y_var) := levels(data[[as_label(y_var)]])[1]
      )
  }
  
  p <- ggplot(data, aes(x = factor(Generation), y = !!y_var, fill = !!sym(fill_var))) +
    geom_tile(width = 2, linewidth = 0) +
    geom_vline(data = dashed_lines_data,
               aes(xintercept = match(generation_line, sort(unique(data$Generation)))),
               linetype = "dashed", color = "black", alpha = 0.7)
  
  # Add facet labels if requested
  if (label_facets) {
    p <- p + 
      geom_text(data = facet_labels, 
                aes(label = label), 
                x = -Inf, y = Inf, 
                hjust = 1, vjust = 0,
                size = 5, fontface = "bold",
                inherit.aes = FALSE)
  }
  
  p <- p +
    scale_fill_gradient2(
      low = "#4575b4",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = mean(range(data[[fill_var]], na.rm = TRUE)),
      limits = range(data[[fill_var]], na.rm = TRUE),
      name = fill_name
    ) +
    scale_x_discrete(breaks = label_gens, labels = label_gens, expand = c(0, 0)) +
    labs(title = title, x = "Generation", y = NULL) +  # Changed y_label to NULL
    facet_formula +
    coord_cartesian(clip = "off") +
    presentation_theme() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 9),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x = unit(0.25, "cm"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
      strip.text.y = element_text(angle = 0),
      strip.text.x = element_blank(),  # Top facet labels already removed
      strip.background.x = element_blank(),  # Added: removes top strip background
      plot.margin = margin(t = 10, r = 5, b = 5, l = 15)
    )
  
  return(p)
}
create_heatmap <- function(data, fill_var, fill_name, title = NULL, label_facets = TRUE, 
                           show_envtype_labels = FALSE) {
  
  # Determine faceting variables based on facet_by_x
  if (facet_by_x) {
    row_var <- "X_label"
    col_var <- "EnvType"
  } else {
    row_var <- "ImportanceLabel"
    col_var <- "EnvType"
  }
  
  # Create facet labels if requested
  if (label_facets) {
    facet_labels <- data %>%
      distinct(across(all_of(c(row_var, col_var)))) %>%
      arrange(across(all_of(c(row_var, col_var)))) %>%
      mutate(
        label = letters[row_number()],
        Generation = min(data$Generation),
        # Use the first level of the y-axis factor for positioning
        !!as_label(y_var) := levels(data[[as_label(y_var)]])[1]
      )
  }
  
  # Create EnvType labels if requested - ONLY for top row
  if (show_envtype_labels) {
    # Calculate midpoint of x-axis for centering
    n_generations <- length(unique(data$Generation))
    mid_x <- ceiling(n_generations / 2)
    
    envtype_labels <- data %>%
      distinct(EnvType) %>%
      arrange(EnvType) %>%
      mutate(
        Generation = mid_x,
        # CRITICAL: Set to FIRST level of row variable (top row only)
        !!sym(row_var) := levels(data[[row_var]])[1],
        # Use the first level of the y-axis factor for positioning
        !!as_label(y_var) := levels(data[[as_label(y_var)]])[1]
      )
  }
  
  p <- ggplot(data, aes(x = factor(Generation), y = !!y_var, fill = !!sym(fill_var))) +
    geom_tile(width = 2, linewidth = 0) +
    geom_vline(data = dashed_lines_data,
               aes(xintercept = match(generation_line, sort(unique(data$Generation)))),
               linetype = "dashed", color = "black", alpha = 0.7)
  
  # Add facet labels if requested
  if (label_facets) {
    p <- p + 
      geom_text(data = facet_labels, 
                aes(label = label), 
                x = -Inf, y = Inf, 
                hjust = 1, vjust = 0,
                size = 5, fontface = "bold",
                inherit.aes = FALSE)
  }
  
  # Add EnvType labels if requested
  if (show_envtype_labels) {
    p <- p + 
      geom_text(data = envtype_labels,
                aes(x = Generation, y = Inf, label = EnvType),
                vjust = -1,      # Position above plot
                hjust = 0.5,       # Center horizontally
                size = 5,          # Text size
                fontface = "italic",
                inherit.aes = FALSE)
  }
  
  p <- p +
    scale_fill_gradient2(
      low = "#4575b4",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = mean(range(data[[fill_var]], na.rm = TRUE)),
      limits = range(data[[fill_var]], na.rm = TRUE),
      name = fill_name
    ) +
    scale_x_discrete(breaks = label_gens, labels = label_gens, expand = c(0, 0)) +
    labs(title = title, x = "Generation", y = NULL) +
    facet_formula +
    coord_cartesian(clip = "off") +
    presentation_theme() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 9),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x = unit(0.25, "cm"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
      strip.text.y = element_text(angle = 0),
      strip.text.x = element_blank(),
      strip.background.x = element_blank(),
      plot.margin = margin(t = if(show_envtype_labels) 25 else 10, r = 5, b = 5, l = 15)
    )
  
  return(p)
}
# ============================================================================
# CREATE ALL HEATMAPS
# ============================================================================

richness_heatmap <- create_heatmap(heatmap_data, "SpeciesRichness_mean", "Species Richness      ",show_envtype_labels = T)
diversity_heatmap <- create_heatmap(heatmap_data, "HostAlphaDiversity_mean", "Alpha Diversity",show_envtype_labels = T)
microbebray_heatmap <- create_heatmap(heatmap_data, "BrayDiv_mean", "Bray-Curtis Dissimilarity",show_envtype_labels = T)
hostfitness_heatmap <- create_heatmap(heatmap_data, "HostFitness_mean", "Host Fitness",show_envtype_labels = T)
ne_heatmap <- create_heatmap(heatmap_data, "Ne_mean", "Effective Population Size",show_envtype_labels=T)


ggsave("SupFig1_richness_heatmap.png", richness_heatmap, width = 16, height = 10)
ggsave("SupFig2_diversity_heatmap.png", diversity_heatmap, width = 16, height = 10)
ggsave("SupFig3_microbebray_heatmap.png", microbebray_heatmap, width = 16, height = 10)
ggsave("SupFig4_hostfitness_heatmap.png", hostfitness_heatmap, width = 16, height = 10)
ggsave("SupFig5_ne_heatmap.png", ne_heatmap, width = 16, height = 10)


envplot_single <- heatmap_data %>% 
  filter(ImportanceNumeric == 0, X == 0) %>%
  ggplot() +
  aes(x = factor(Generation), y = env_value, color = EnvType) +  # Color by EnvType
  geom_line(aes(group = EnvType), linewidth = 0.8) +
  geom_vline(data = dashed_lines_data,
             aes(xintercept = match(generation_line, sort(unique(heatmap_data$Generation)))),
             linetype = "dashed", color = "black", alpha = 0.7) +
  scale_x_discrete(breaks = label_gens, labels = label_gens, expand = c(0, 0)) +
  labs(x = "Generation", y = "Environment") +
  presentation_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9),  # Changed from element_blank()
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # Changed from element_blank()
    axis.ticks.length.x = unit(0.25, "cm"),  # Added for consistency with heatmap
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    legend.position = "none"
  )


richness_withenv <- (richness_heatmap | ((plot_spacer() / envplot_single / plot_spacer()) +   plot_layout(heights = c(5,4,5), guides = "collect"))) + 
  plot_layout(widths = c(10, 2), guides = "collect")

ggsave("SupFig1_richnesswithenv_heatmap.png", richness_withenv, width = 16, height = 10)

ne_heatmap <- create_heatmap(heatmap_data, "Ne_mean", "Effective Population Size",label_facets=T)
polymorphicsites_heatmap <- create_heatmap(heatmap_data, "NumPolymorphicLoci_mean", "Number of Polymorphic Loci")
microbefitness_heatmap <- create_heatmap(heatmap_data, "AvgMicrobialFitness_mean", "Microbe Fitness")
HostGenotype_sd <- create_heatmap(heatmap_data, "HostGenotype_sd", "HostGenotype_sd")
diversity_sd <- create_heatmap(heatmap_data, "BrayDiv_sd", "BrayDiv_sd")
HostGenotype_sd
envdiv <- create_heatmap(heatmap_data, "EnvDiv_mean", "EnvDiv_mean")

# ============================================================================
# COMBINE PLOTS
# ============================================================================
envplot_single <- heatmap_data %>% 
  filter(ImportanceNumeric == 0, X == 0) %>%
  ggplot() +
  aes(x = factor(Generation), y = env_value, color = EnvType) +  # Color by EnvType
  geom_line(aes(group = EnvType), linewidth = 0.8) +
  geom_vline(data = dashed_lines_data,
             aes(xintercept = match(generation_line, sort(unique(heatmap_data$Generation)))),
             linetype = "dashed", color = "black", alpha = 0.7) +
  scale_x_discrete(breaks = label_gens, labels = label_gens, expand = c(0, 0)) +
  labs(x = NULL, y = "Environment") +
  presentation_theme() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9),  # Changed from element_blank()
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # Changed from element_blank()
    axis.ticks.length.x = unit(0.25, "cm"),  # Added for consistency with heatmap
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
    legend.position = "none"
  )


hostfitness_combined_plot <- hostfitness_heatmap +guides(fill = guide_colorbar(title = "Host Fitness                      ")) +
  inset_element(envplot_single, 
                left = 0.84,   # Adjust these values to position
                bottom = 0.65, 
                right = 1, 
                top = 0.8,
                align_to = "full")


ggsave("hostfitness_combined_plot.png", hostfitness_combined_plot, width = 16, height = 10)
microbefitness_combined_plot <- microbefitness_heatmap +guides(fill = guide_colorbar(title = "Microbe Fitness              ")) +
  inset_element(envplot_single, 
                left = 0.84,   # Adjust these values to position
                bottom = 0.65, 
                right = 1, 
                top = 0.8,
                align_to = "full")
ggsave("microbefitness_combined_plot.png", microbefitness_combined_plot, width = 16, height = 10)

braydiv_combined_plot <- microbebray_heatmap + 
  inset_element(envplot_single, 
                left = 0.84,   # Adjust these values to position
                bottom = 0.65, 
                right = 1, 
                top = 0.8,
                align_to = "full")
ggsave("braydiv_combined_plot.png", braydiv_combined_plot, width = 16, height = 10)

richness_combined_plot <- richness_heatmap +guides(fill = guide_colorbar(title = "Species Richness              ")) +
  inset_element(envplot_single, 
                left = 0.84,   # Adjust these values to position
                bottom = 0.65, 
                right = 1, 
                top = 0.8,
                align_to = "full")

ggsave("richness_combined_plot.png", richness_combined_plot, width = 16, height = 10)

diversity_combined_plot <- diversity_heatmap  +guides(fill = guide_colorbar(title = "Alpha Diversity               ")) +
  inset_element(envplot_single, 
                left = 0.84,   # Adjust these values to position
                bottom = 0.65, 
                right = 1, 
                top = 0.8,
                align_to = "full")
ggsave("diversity_combined_plot.png", diversity_combined_plot, width = 16, height = 10)


# ============================================================================
# SAVE PLOTS
# ============================================================================


cat("All plots saved successfully!\n")

