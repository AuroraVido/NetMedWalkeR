# Clear environment
rm(list = ls())
graphics.off()

# Load packages
library(igraph)
library(dplyr)
library(magrittr)
library(data.table)
library(pROC)
library(future)
library(future.apply)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Load configuration and helper functions
source("config.R")
source("functions.R")

# Set seed from config
set.seed(SEED_VALUE)

# ==============================================================================
# 1. DATA PREPARATION & TOPOLOGY CHECK
# ==============================================================================
cat("--- 1. Loading and preparing data ---\n")

# Read input files defined in config
disease_genes_df <- read.delim(FILE_DISEASE_GENES)
interactome_df   <- fread(FILE_INTERACTOME)
drug_targets_df  <- read.delim(FILE_DRUG_TARGETS)

# Ensure ID formats
interactome_df$gene1 <- as.character(interactome_df$gene1)
interactome_df$gene2 <- as.character(interactome_df$gene2)

# Build graph and filter for Largest Connected Component (LCC)
full_graph <- igraph::graph_from_data_frame(interactome_df, directed = FALSE) %>%
  igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)

interactome_graph <- get_largest_connected_component(full_graph)
all_network_nodes <- V(interactome_graph)$name

cat("Interactome LCC loaded.\n")
cat(" - Nodes:", vcount(interactome_graph), "\n")
cat(" - Edges:", ecount(interactome_graph), "\n")

# Check Scale-Free Topology (Exploratory)
cat("Generating Scale-Free Topology Plot...\n")
check_scale_free_topology(interactome_graph)

# ==============================================================================
# 2. LOOCV OPTIMIZATION
# ==============================================================================
cat("\n--- 2. Optimizing RWR Damping Parameter ---\n")

loo_genes <- intersect(disease_genes_df$GeneID, all_network_nodes)

# Initialize parallel processing
plan(multisession, workers = N_CORES)

# Run LOOCV across the damping grid defined in config
mean_ranks <- future_sapply(DAMPING_GRID, function(d) {
  run_loocv_for_damping(d, interactome_graph, loo_genes)
}, future.seed = TRUE)

optimal_idx <- which.min(mean_ranks)
optimal_damping <- DAMPING_GRID[optimal_idx]

# ==============================================================================
# 3. RWR EXECUTION & DRUG PREPARATION
# ==============================================================================
cat("\n--- 3. Running Final RWR & Preparing Targets ---\n")

# Setup seed vector with disease genes
p0 <- rep(0, vcount(interactome_graph))
names(p0) <- V(interactome_graph)$name
p0[loo_genes] <- 1 / length(loo_genes)

# Execute RWR
rwr_scores <- igraph::page_rank(
  interactome_graph,
  algo = "prpack",
  directed = FALSE,
  damping = optimal_damping,
  personalized = p0,
  weights = E(interactome_graph)$weight 
)$vector

rwr_results_df <- data.frame(
  gene = names(rwr_scores),
  rwr_score = rwr_scores,
  stringsAsFactors = FALSE
)

# Log transform scores for statistical testing
epsilon <- min(rwr_results_df$rwr_score[rwr_results_df$rwr_score > 0]) / 10
rwr_results_df$log_score <- log10(rwr_results_df$rwr_score + epsilon)

# Prepare Drug Targets list
drug_targets_clean <- drug_targets_df %>%
  dplyr::rename(drug_id = DrugID, target_id = TargetID) %>%
  mutate(target_id = as.character(target_id)) %>%
  filter(target_id %in% all_network_nodes)

drug_list <- split(drug_targets_clean$target_id, drug_targets_clean$drug_id)

# ==============================================================================
# 3. NETWORK VISUALIZATION (Top 1000 nodes based on RWR)
# ==============================================================================
cat("\n--- 3. Generating Subgraph Visualization ---\n")

# Se non li hai, installali con install.packages(c("ggraph", "tidygraph"))
library(ggraph)
library(tidygraph)

# 3.1 Vettore di restart (Uniforme sui seed genes)
seeds <- intersect(disease_genes_df$GeneID, all_network_nodes)
restart_vec <- numeric(vcount(interactome_graph))
names(restart_vec) <- V(interactome_graph)$name
restart_vec[seeds] <- 1 / length(seeds)

# 3.2 Esegui RWR per trovare i nodi topologici più vicini
# Usiamo page_rank di igraph, che è il motore sotto al cofano per il RWR
rwr_res <- page_rank(interactome_graph, 
                     algo = "prpack", 
                     v = V(interactome_graph),
                     damping = optimal_damping, 
                     personalized = restart_vec)

# 3.3 Estrai i top 1000 nodi (inclusi sempre i seed genes per sicurezza)
top_nodes <- names(sort(rwr_res$vector, decreasing = TRUE)[1:2000])
nodes_to_keep <- unique(c(seeds, top_nodes))

# 3.4 Crea il sottografo
sub_graph <- induced_subgraph(interactome_graph, vids = nodes_to_keep)

# 3.5 Aggiungi attributi per il plot (Seed vs Top Nodes)
V(sub_graph)$node_type <- ifelse(V(sub_graph)$name %in% seeds, "Disease Seed", "Other Nodes")

cat("Subgraph created.\n")
cat(" - Nodes:", vcount(sub_graph), "\n")
cat(" - Edges:", ecount(sub_graph), "\n")

# ==============================================================================
# 4. STATISTICAL VALIDATION (PERMUTATIONS)
# ==============================================================================
cat("\n--- 4. Running Degree-Preserving Permutations ---\n")

# Bin nodes by degree
node_info <- data.frame(
  node = all_network_nodes,
  degree = degree(interactome_graph),
  log_rwr_score = rwr_results_df$log_score[match(all_network_nodes, rwr_results_df$gene)],
  stringsAsFactors = FALSE
)

breaks <- unique(quantile(node_info$degree, probs = seq(0, 1, by = 0.05)))
if(length(breaks) < 5) breaks <- 5
node_info$degree_bin <- cut(node_info$degree, breaks = breaks, include.lowest = TRUE, labels = FALSE)

# Create lookup bins
bins_list <- split(node_info$log_rwr_score, node_info$degree_bin)

# Run permutations in parallel
cat("Processing", length(drug_list), "drugs...\n")
results_list <- future_lapply(names(drug_list), function(d_id) {
  run_log_permutation(d_id, drug_list[[d_id]], bins_list, node_info, n_perms = N_PERMS)
}, future.seed = TRUE)

proximity_results <- do.call(rbind, results_list)

# ==============================================================================
# 5. NETWORK SEPARATION (s_AB)
# ==============================================================================
cat("\n--- 5. Calculating Network Separation ---\n")

# Pre-calculate internal disease distance
d_BB_matrix <- igraph::distances(interactome_graph, v = loo_genes, to = loo_genes, weights = NA)
mean_d_BB <- mean(d_BB_matrix[upper.tri(d_BB_matrix)])

sep_results_list <- future_lapply(names(drug_list), function(d_id) {
  res <- calc_network_separation(drug_list[[d_id]], loo_genes, interactome_graph, mean_d_BB)
  res$drug_id <- d_id
  return(res)
}, future.seed = TRUE)

separation_df <- do.call(rbind, sep_results_list)

# ==============================================================================
# 6. DISEASE CONTEXT SCORING
# ==============================================================================
cat("\n--- 6. Integrating Disease Context ---\n")

disease_db <- read.delim(FILE_DISEASE_DB, stringsAsFactors = FALSE)
# Ensure columns exist, adjust if file format differs
if(!"GeneID" %in% colnames(disease_db)) colnames(disease_db)[2] <- "GeneID" 
disease_db$GeneID <- as.character(disease_db$GeneID)

acm_gene_set <- unique(disease_genes_df$GeneID)

# Calculate similarity with other diseases
disease_similarity <- disease_db %>%
  group_by(disease) %>%
  summarise(
    jaccard_score = calc_jaccard(GeneID, acm_gene_set)
  ) %>%
  filter(jaccard_score > 0) %>%
  # Exclude the disease itself using config variable
  filter(!grepl(DISEASE_TO_EXCLUDE, disease, ignore.case = TRUE)) 

# Weight genes based on context
weighted_gene_data <- disease_db %>%
  inner_join(disease_similarity %>% select(disease, jaccard_score), by = "disease")

gene_relevance_df <- weighted_gene_data %>%
  group_by(GeneID) %>%
  summarise(raw_score = sum(jaccard_score)) %>%
  ungroup() %>%
  mutate(context_weight = raw_score / max(raw_score))

# Score drugs
context_scores_list <- future_lapply(names(drug_list), function(d_id) {
  s <- calc_context_score(drug_list[[d_id]], gene_relevance_df)
  data.frame(drug_id = d_id, context_score = s)
}, future.seed = TRUE)

context_scores_df <- do.call(rbind, context_scores_list)

# ==============================================================================
# 7. MERGING & RANKING
# ==============================================================================
cat("\n--- 7. Merging Results and Ranking ---\n")

# Merge proximity, separation, and context scores
full_analysis <- proximity_results %>%
  left_join(separation_df, by = "drug_id") %>%
  left_join(context_scores_df, by = "drug_id") %>%
  mutate(padj = p.adjust(emp_p_value, method = "BH"))

# Define Topology Class and order by Context Score (removing Final Score calculation)
full_analysis <- full_analysis %>%
  mutate(
    Topology_Class = case_when(
      is.na(s_ab) ~ "Undefined",
      s_ab < 0 ~ "Overlap",
      s_ab >= 0 & s_ab < 1 ~ "Proximal",
      s_ab > 1 ~ "Distant",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(desc(context_score)) # Ordered by context score

# Save results
saveRDS(full_analysis, "full_drug_analysis_results.rds")
write.csv(full_analysis, "full_drug_analysis_results.csv", row.names = FALSE)

# ==============================================================================
# 7b. SELECTION OF SIGNIFICANT CANDIDATES
# ==============================================================================
cat("\n--- 7b. Defining Significant Candidates ---\n")

# Filter criteria:
# 1. Adjusted P-value < 0.05 (Significant proximity)
# 2. Context Score > 0 (Relevance to disease context)
significant_candidates <- full_analysis %>%
  filter(padj < 0.05 & context_score > 0) %>%
  arrange(desc(context_score))

cat(paste("Total drugs analyzed:", nrow(full_analysis), "\n"))
cat(paste("Significant candidates selected:", nrow(significant_candidates), "\n"))

# Preview the top hits
# Removed Final_Score from select() as it no longer exists
print(head(significant_candidates %>% select(drug_id, z_score, padj, context_score, s_ab), 15))

# Save significant candidates list
saveRDS(significant_candidates, "significant_candidates.rds")
write.csv(significant_candidates, "significant_candidates.csv", row.names = FALSE)

# ==============================================================================
# 8. PLOTTING SEPARATION vs Z-SCORE
# ==============================================================================
cat("\n--- 9. Plotting Network Map ---\n")

plot_data <- full_analysis %>%
  filter(z_score > -5 & s_ab < 2.0) %>%
  mutate(
    is_hit = ifelse(padj < 0.05 & z_score > 1.65 & context_score > 0, "Yes", "No"),
    label_text = ifelse(is_hit == "Yes" & rank(-z_score) <= 20, drug_id, NA)
  )

p_sep <- ggplot(plot_data, aes(x = s_ab, y = z_score)) +
  geom_point(data = subset(plot_data, is_hit == "No"), color = "grey85", alpha = 0.5) +
  #geom_point(data = subset(plot_data, is_hit == "Yes"), aes(color = Topology_Class), size = 1.5) +
  geom_point(data = subset(plot_data, is_hit == "Yes"), aes(color = Topology_Class, size = context_score), alpha = 0.8) +
  scale_size_continuous(range = c(1, 4), name = "Context Score") +
  geom_hline(yintercept = 1.65, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_text_repel(aes(label = label_text), max.overlaps = 14, size = 3) +
  annotate("text", x = -0.1, y = max(plot_data$z_score), label = "Overlap", color = "firebrick", fontface = "bold", size=3) +
  annotate("text", x = 0.15, y = max(plot_data$z_score), label = "Proximal", color = "dodgerblue", fontface = "bold", size=3) +
  annotate("text", x = 1.2, y = max(plot_data$z_score), label = "Distant", color = "orange", fontface = "bold", size=3) +
  scale_color_manual(values = c("Overlap" = "firebrick", "Proximal" = "dodgerblue", "Distant" = "orange"), name = "Topology Class") +
  theme_minimal() +
  theme(legend.position = "bottom")+
  labs(title = "Drug-Disease Network Separation", x = "Network Separation (s_ab)", y = "Z-Score")

ggsave("Network_Separation_Plot.png", p_sep, width = 10, height = 7)

# ==============================================================================
# 9. MODULE DEFINITION & EXPORT
# ==============================================================================
cat("\n--- 10. Extracting Disease Module ---\n")

# Determine threshold using Elbow method
sorted_scores <- sort(rwr_scores, decreasing = TRUE)
# Limit calculation to top 2000 to avoid long tail noise
elbow_limit <- min(2000, length(sorted_scores))
elbow_idx <- find_elbow_point(sorted_scores[1:elbow_limit])

cat(paste("Module Size:", elbow_idx, "genes\n"))

# Extract module
module_genes <- names(sorted_scores)[1:elbow_idx]
module_subgraph <- igraph::induced_subgraph(interactome_graph, module_genes)

# Detect Communities (Louvain)
communities <- igraph::cluster_louvain(module_subgraph)
module_data <- data.frame(
  node = V(module_subgraph)$name,
  community = communities$membership,
  stringsAsFactors = FALSE
)

# Export for Metascape analysis (or other online tools for pathways analysis)
out_dir <- "Metascape_Input_Files"
if(!dir.exists(out_dir)) dir.create(out_dir)

for (mod in unique(module_data$community)) {
  genes <- module_data$node[module_data$community == mod]
  write.table(genes, file.path(out_dir, paste0("Module_", mod, ".txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

cat("Analysis Complete. Check output folder for results.\n")