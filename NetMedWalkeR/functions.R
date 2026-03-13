# ==============================================================================
# FUNCTIONS LIBRARY
# ==============================================================================
# This file contains the core logic and helper functions
# Do not modify unless changing the underlying algorithms

# Load necessary libraries inside functions or ensure they are loaded in main
library(igraph)
library(dplyr)
library(ggplot2)
library(pROC)

# --- Graph Preparation ---

get_largest_connected_component <- function(graph) {
  # Decompose graph into connected components
  components <- igraph::components(graph)
  # Identify the largest component
  lcc_id <- which.max(components$csize)
  # Extract node indices
  lcc_nodes <- which(components$membership == lcc_id)
  # Return subgraph
  return(igraph::induced_subgraph(graph, lcc_nodes))
}

check_scale_free_topology <- function(graph, filename = "Scale_Free_Plot.png") {
  # Calculate degrees
  d <- degree(graph)
  # Create frequency table
  dd <- as.data.frame(table(d))
  dd$d <- as.numeric(as.character(dd$d))
  dd$Freq <- as.numeric(dd$Freq)
  dd$Prob <- dd$Freq / sum(dd$Freq)
  
  # Plot
  p <- ggplot(dd, aes(x = d, y = Prob)) +
    geom_point(color = "darkgrey", alpha = 0.7) +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method = "lm", color = "firebrick", se = FALSE, linetype = "dashed") +
    labs(title = "Scale-Free Topology Check",
         x = "Degree (k)", y = "Probability P(k)") +
    theme_minimal()
  
  ggsave(filename, p, width = 8, height = 6)
  return(p)
}

# --- RWR & Optimization ---

run_loocv_for_damping <- function(d_val, graph, genes) {
  ranks <- numeric(length(genes))
  
  for (i in seq_along(genes)) {
    gene_out <- genes[i]
    seeds <- genes[-i]
    
    # Skip if no seeds available
    if (length(seeds) == 0) {
      ranks[i] <- NA
      next
    }
    
    # Setup restart vector
    p0 <- rep(0, vcount(graph))
    names(p0) <- V(graph)$name
    p0[seeds] <- 1 / length(seeds)
    
    # Run RWR
    rwr_res <- igraph::page_rank(
      graph,
      algo = "prpack",
      directed = FALSE,
      damping = d_val,
      personalized = p0,
      weights = E(graph)$weight 
    )$vector
    
    # Calculate rank
    ranks[i] <- rank(-rwr_res, ties.method = "min")[gene_out]
  }
  return(mean(ranks, na.rm = TRUE))
}

# --- Statistical Validation (Permutations) ---

run_log_permutation <- function(drug_id, real_targets, bins_list, node_info_df, n_perms) {
  
  valid_targets <- real_targets[real_targets %in% node_info_df$node]
  n_t <- length(valid_targets)
  
  # Calculate observed score
  target_indices <- match(valid_targets, node_info_df$node)
  obs_score <- mean(node_info_df$log_rwr_score[target_indices])
  target_bins <- node_info_df$degree_bin[target_indices]
  
  # Perform random sampling
  perm_means <- numeric(n_perms)
  for(i in 1:n_perms) {
    random_scores <- sapply(target_bins, function(b) {
      pool <- bins_list[[as.character(b)]]
      if(length(pool) == 1) return(pool)
      sample(pool, 1)
    })
    perm_means[i] <- mean(random_scores)
  }
  
  # Calculate statistics
  mean_rand <- mean(perm_means)
  sd_rand <- sd(perm_means)
  if (sd_rand == 0) sd_rand <- 1e-10
  
  z_score <- (obs_score - mean_rand) / sd_rand
  emp_p_value <- (sum(perm_means >= obs_score) + 1) / (n_perms + 1)
  
  return(data.frame(
    drug_id = drug_id,
    log_proximity = obs_score,
    z_score = z_score,
    emp_p_value = emp_p_value,
    n_targets = n_t
  ))
}

# --- Network Separation ---

calc_network_separation <- function(drug_targets, disease_genes, graph, d_bb_val) {
  
  targets_A <- intersect(drug_targets, V(graph)$name)
  
  # Return NA if no targets in graph
  if (length(targets_A) == 0) return(data.frame(s_ab = NA, d_aa = NA, d_ab = NA))
  
  # Internal Drug Distance (d_AA)
  if (length(targets_A) > 1) {
    d_AA_matrix <- igraph::distances(graph, v = targets_A, to = targets_A, weights = NA)
    mean_d_AA <- mean(d_AA_matrix[upper.tri(d_AA_matrix)])
  } else {
    mean_d_AA <- 0 
  }
  
  # Drug-Disease Distance (d_AB)
  d_AB_matrix <- igraph::distances(graph, v = targets_A, to = disease_genes, weights = NA)
  
  # Handle infinite distances in disconnected components
  if(any(is.infinite(d_AB_matrix))) {
    max_finite <- max(d_AB_matrix[is.finite(d_AB_matrix)], na.rm = TRUE)
    d_AB_matrix[is.infinite(d_AB_matrix)] <- max_finite + 1
  }
  mean_d_AB <- mean(d_AB_matrix)
  
  # Calculate s_AB
  s_ab <- mean_d_AB - (mean_d_AA + d_bb_val) / 2
  
  return(data.frame(s_ab = s_ab, d_aa = mean_d_AA, d_ab = mean_d_AB))
}

# --- Context Scoring (Jaccard) ---

calc_jaccard <- function(genes_x, genes_y) {
  inter <- length(intersect(genes_x, genes_y))
  uni <- length(union(genes_x, genes_y))
  if (uni == 0) return(0)
  return(inter / uni)
}

calc_context_score <- function(drug_targets, relevance_df) {
  valid_targets <- intersect(drug_targets, relevance_df$GeneID)
  if(length(valid_targets) == 0) return(0)
  
  # Sum weights of hit targets
  score <- sum(relevance_df$context_weight[relevance_df$GeneID %in% valid_targets])
  
  # Normalize by number of targets
  return(score / length(drug_targets))
}

# --- Module Detection ---

find_elbow_point <- function(scores) {
  n <- length(scores)
  if(n < 2) return(1)
  coords <- cbind(seq_len(n), scores)
  line_vec <- coords[n, ] - coords[1, ]
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))
  vec_from_first <- t(t(coords) - coords[1, ])
  scalar_prod <- vec_from_first %*% line_vec_norm
  vec_proj <- t(line_vec_norm %*% t(scalar_prod))
  vec_dist <- vec_from_first - vec_proj
  dist_to_line <- sqrt(rowSums(vec_dist^2))
  return(which.max(dist_to_line))
}