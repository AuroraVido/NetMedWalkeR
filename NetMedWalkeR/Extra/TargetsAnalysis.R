# --- LIBRARIES ---
# Ensure all necessary libraries are loaded
library(ggplot2)
library(ggraph)
library(tidygraph)
library(dplyr)
library(tidyr)
library(boot)
library(scales)
library(corrplot)
library(future)
library(future.apply)
library(rlang) # Required for !!sym()


# --- SECTION 10: TOP CANDIDATE TARGET ANALYSIS ---
# ---------------------------------------------------------------------------
cat("\n\n--- 10. TOP CANDIDATE TARGET ANALYSIS ---\n")

if (nrow(significant_candidates) == 0) {
  cat("No significant drugs found. Skipping network analysis.\n")
} else {
  
  cat(paste("Starting network analysis on the targets of significant vs non-significant drugs...\n"))
  
  # --- 10.1 Setup: Identify genes of interest ---
  
  # 1. Targets of "Good" Drugs (padj < 0.05)
  sig_targets <- drug_targets_prepared %>%
    filter(drug_id %in% significant_candidates$drug_id) %>% 
    pull(target_id) %>%
    unique()
  
  # 2. Targets of "Bad" Drugs (padj > 0.8)
  non_sig_drugs <- full_analysis %>% filter(padj > 0.8) %>% pull(drug_id)
  nonsig_targets <- drug_targets_prepared %>%
    filter(drug_id %in% non_sig_drugs) %>%
    pull(target_id) %>%
    unique()
  
  # 3. Clean up overlap (Remove any "Bad" targets that are also "Good" targets)
  nonsig_targets <- setdiff(nonsig_targets, sig_targets)
  
  # Define ALL nodes needed for analysis
  nodes_of_interest <- unique(c(sig_targets, nonsig_targets))
  
  cat(paste(" - Target Genes to Analyze:", length(nodes_of_interest), "\n"))
  cat(paste("   (Significant:", length(sig_targets), "vs Non-Significant:", length(nonsig_targets), ")\n"))
  
  # --- 10.2 Parallelized Network Metrics Calculation ---
  cat("Calculating centrality and distance metrics (Parallelized)...\n")
  
  # FIX: Ensure at least 1 worker is used
  n_workers <- max(1, availableCores() - 1)
  plan(multisession, workers = n_workers) 
  
  # Split nodes into chunks for parallel processing
  node_chunks <- split(nodes_of_interest, ceiling(seq_along(nodes_of_interest) / (length(nodes_of_interest) / n_workers)))
  
  calc_metrics_chunk <- function(nodes_subset, graph, disease_genes) {
    # Check if nodes exist in the graph to avoid errors
    valid_nodes <- nodes_subset[nodes_subset %in% igraph::V(graph)$name]
    
    if(length(valid_nodes) == 0) return(NULL)
    
    d_mat <- igraph::distances(graph, v = valid_nodes, to = disease_genes)
    min_dist <- apply(d_mat, 1, min)
    deg <- igraph::degree(graph, v = valid_nodes)
    bet <- igraph::betweenness(graph, v = valid_nodes, directed = FALSE, normalized = TRUE)
    clust <- igraph::transitivity(graph, type = "local", vids = valid_nodes, isolates = "zero")
    
    return(data.frame(
      node = valid_nodes,
      degree = deg,
      betweenness = bet,
      clustering_coeff = clust,
      dist_to_disease = min_dist,
      stringsAsFactors = FALSE
    ))
  }
  
  # Run parallel calculation
  metrics_list <- future_lapply(node_chunks, function(chunk) {
    calc_metrics_chunk(chunk, interactome_graph, acm_genes_in_network)
  }, future.seed = TRUE)
  
  target_metrics_df <- do.call(rbind, metrics_list)
  
  # Add Global Metrics (Eigenvector is fast enough to run on single core usually, but good to add here)
  eigen_global <- igraph::eigen_centrality(interactome_graph, directed = FALSE)$vector
  target_metrics_df$eigenvector <- eigen_global[target_metrics_df$node]
  
  # Assign Groups
  target_metrics_df <- target_metrics_df %>%
    mutate(group = case_when(
      node %in% sig_targets ~ "Significant (padj < 0.05)",
      node %in% nonsig_targets ~ "Non-Significant (padj > 0.8)"
    ))
  
  cat("Metrics calculation complete.\n")
  
  # --- 10.3 UNIVARIATE ANALYSIS ---
  cat("\n--- A. Univariate Analysis Results (Wilcoxon) ---\n")
  
  metrics_to_test <- c("degree", "betweenness", "dist_to_disease", "eigenvector", "clustering_coeff")
  wilcox_summary <- data.frame()
  
  for(m in metrics_to_test) {
    sig_v <- target_metrics_df %>% filter(group == "Significant (padj < 0.05)") %>% pull(!!sym(m))
    nonsig_v <- target_metrics_df %>% filter(group == "Non-Significant (padj > 0.8)") %>% pull(!!sym(m))
    
    if(length(sig_v) > 0 && length(nonsig_v) > 0) {
      wt <- wilcox.test(sig_v, nonsig_v)
      
      row <- data.frame(
        Metric = m,
        Median_Sig = median(sig_v, na.rm=TRUE),
        Median_NonSig = median(nonsig_v, na.rm=TRUE),
        P_Value = wt$p.value,
        Signif = case_when(
          wt$p.value < 0.001 ~ "***",
          wt$p.value < 0.01 ~ "**",
          wt$p.value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
      wilcox_summary <- rbind(wilcox_summary, row)
    }
  }
  
  print(wilcox_summary)
  
  # Plots
  cat("\nGenerating Boxplots...\n")
  for(metric_name in metrics_to_test) {
    scale_type <- "linear"
    if(metric_name == "degree") scale_type <- "log1p"
    if(metric_name %in% c("betweenness", "eigenvector")) scale_type <- "pseudolog"
    
    plot_data <- target_metrics_df
    
    # Temporary scaling for visualization
    if(scale_type == "log1p") plot_data[[metric_name]] <- log1p(plot_data[[metric_name]])
    if(scale_type == "pseudolog") {
      eps <- min(plot_data[[metric_name]][plot_data[[metric_name]] > 0], na.rm=T)/10
      plot_data[[metric_name]] <- log10(plot_data[[metric_name]] + eps)
    }
    
    p <- ggplot(plot_data, aes(x = group, y = .data[[metric_name]], fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
      labs(title = paste(metric_name), 
           y = ifelse(scale_type!="linear", paste(scale_type, metric_name), metric_name),
           x = "") +
      theme_minimal() +
      theme(legend.position = "none")
    
    ggsave(paste0("Univariate_", metric_name, ".png"), p, width = 5, height = 5)
  }
  cat("Univariate plots saved.\n")
  
  # --- 10.4 CORRELATION ANALYSIS ---
  cat("\n--- B. Multicollinearity Check ---\n")
  cor_data <- target_metrics_df %>% select(degree, betweenness, clustering_coeff, dist_to_disease, eigenvector)
  M <- cor(cor_data, method = "spearman", use = "complete.obs")
  
  # Save plot
  png("Target_Metrics_Correlation.png", width = 1800, height = 1800, res = 300)
  corrplot(M, method = "color", type = "upper", addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE, mar=c(0,0,2,0))
  dev.off()
  cat("Correlation matrix saved.\n")
  
  # --- 10.5 DATA PREPARATION FOR REGRESSION ---
  cat("\n--- Preparing Data for Logistic Regression (Scaling) ---\n")
  
  # Create the scaled dataset ONCE to ensure consistency between univariate and multivariate models
  model_data <- target_metrics_df %>%
    # Define the binary outcome variable (Y): 1 if Significant, 0 if Non-Significant
    mutate(is_significant = ifelse(group == "Significant (padj < 0.05)", 1, 0)) %>%
    mutate(
      # --- PREDICTOR TRANSFORMATION & SCALING ---
      # corresponding to the equation terms: deg_log, betw, eigen, clust, dist
      
      # 1. Degree: Apply log1p first to correct the power-law skewness, then scale (Z-score).
      #    This creates the "deg_log" term in your equation.
      degree_scaled = as.numeric(scale(log1p(degree))), 
      
      # 2. Other Metrics: Apply Z-score scaling directly (Mean = 0, SD = 1).
      #    This allows comparison of beta coefficients (effect size per 1 SD increase).
      betweenness_scaled = as.numeric(scale(betweenness)),
      dist_scaled = as.numeric(scale(dist_to_disease)),
      eigen_scaled = as.numeric(scale(eigenvector)),
      clust_scaled = as.numeric(scale(clustering_coeff))
    )
  
  # Define the list of predictors to iterate over
  predictors <- c("degree_scaled", "betweenness_scaled", "dist_scaled", "eigen_scaled", "clust_scaled")
  
  # --- 10.5a UNIVARIATE LOGISTIC REGRESSION ---
  cat("\n--- C1. Univariate Logistic Regression (Unadjusted Odds Ratios) ---\n")
  
  uni_results <- data.frame()
  
  # Loop through each predictor to calculate Unadjusted Odds Ratios (Single variable models)
  for (pred in predictors) {
    
    # Dynamically create the formula: "is_significant ~ predictor"
    f <- as.formula(paste("is_significant ~", pred))
    
    # Fit the GLM (Generalized Linear Model)
    # family = binomial(link = "logit") specifies we are modeling the log-odds (logit)
    uni_model <- glm(f, data = model_data, family = binomial(link = "logit"))
    
    # Extract coefficients
    coef_summ <- summary(uni_model)$coefficients
    or <- exp(coef_summ[2, "Estimate"])       # Exponentiate coefficient to get Odds Ratio
    ci <- exp(confint.default(uni_model)[2,]) # Calculate 95% Confidence Interval
    p_val <- coef_summ[2, "Pr(>|z|)"]         # Extract P-value
    
    # Store results
    row <- data.frame(
      Predictor = pred,
      Type = "Univariate",
      Odds_Ratio = round(or, 3),
      Lower_CI = round(ci[1], 3),
      Upper_CI = round(ci[2], 3),
      P_Value = p_val
    )
    uni_results <- rbind(uni_results, row)
  }
  
  # Add significance stars for the table
  uni_results <- uni_results %>%
    mutate(Signif = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )) %>%
    arrange(P_Value)
  
  print(uni_results)
  
  # --- 10.5b MULTIVARIATE LOGISTIC REGRESSION ---
  cat("\n--- C2. Multivariate Logistic Regression (Adjusted Odds Ratios) ---\n")
  
  # This is the implementation of your LaTeX Equation:
  # logit(p) = beta0 + beta1*deg + beta2*betw + beta3*dist + beta4*eigen + beta5*clust
  logit_model <- glm(
    is_significant ~ degree_scaled + betweenness_scaled + dist_scaled + eigen_scaled + clust_scaled,
    data = model_data,
    family = binomial(link = "logit") # Uses the logit link function
  )
  
  cat("Logistic Regression Summary (Multivariate):\n")
  print(summary(logit_model))
  
  cat("\n--- Multivariate Results Table ---\n")
  
  # Extract Adjusted Odds Ratios (controlling for all other variables)
  coefficients <- summary(logit_model)$coefficients
  or_estimates <- exp(coefficients[, "Estimate"])
  conf_int <- exp(confint.default(logit_model))
  
  multi_table <- data.frame(
    Predictor = rownames(coefficients),
    Type = "Multivariate",
    Odds_Ratio = or_estimates,
    Lower_CI = conf_int[, 1],
    Upper_CI = conf_int[, 2],
    P_Value = coefficients[, "Pr(>|z|)"]
  ) %>%
    filter(Predictor != "(Intercept)") %>% # Remove intercept from final table
    mutate(Signif = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )) %>%
    mutate(across(c(Odds_Ratio, Lower_CI, Upper_CI), \(x) round(x, 3)))
  
  print(multi_table)
  
  # --- 10.5c VIF CHECK ---
  # Check for Multicollinearity to ensure coefficients are stable
  cat("\n--- Checking Multicollinearity (VIF) ---\n")
  vif_values <- car::vif(logit_model)
  print(vif_values)
  
  high_vif <- names(vif_values)[vif_values > 5]
  if (length(high_vif) > 0) {
    cat("\n⚠️ WARNING: High Multicollinearity detected (VIF > 5).\n")
  } else {
    cat("\n✅ VIF check passed.\n")
  }
}


# --- DIAGNOSTIC PLOTS: DISTRIBUTION CHECK ---
cat("\n--- Generating Diagnostic Plots for Transformed Variables ---\n")

# Reshape data to long format for ggplot (visualization purposes)
plot_data_long <- model_data %>%
  select(degree_scaled, betweenness_scaled, dist_scaled, eigen_scaled, clust_scaled) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

# 1. Histogram + Density Plot
# Visual check to see if log1p(degree) looks Gaussian and if others are well-behaved
p_hist <- ggplot(plot_data_long, aes(x = Value)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red", linetype = "dashed", size = 0.8) +
  facet_wrap(~Metric, scales = "free") +
  labs(title = "Distribution of Scaled Predictors",
       subtitle = "Blue: Actual Data Density | Red Dashed: Perfect Normal Distribution",
       x = "Z-Score (Scaled Value)", y = "Density") +
  theme_minimal()

ggsave("Distribution_Check_Histograms.png", p_hist, width = 10, height = 6)
print(p_hist)

# 2. Q-Q Plots (Quantile-Quantile)
# Check for normality deviations (points should lie on the red line)
p_qq <- ggplot(plot_data_long, aes(sample = Value)) +
  stat_qq(color = "darkgrey", size = 0.5) +
  stat_qq_line(color = "red", size = 1) +
  facet_wrap(~Metric, scales = "free") +
  labs(title = "Q-Q Plots for Normality Check",
       subtitle = "Points on the red line = Normal Distribution",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

ggsave("Distribution_Check_QQ.png", p_qq, width = 10, height = 6)
print(p_qq)

cat("Diagnostic plots saved: 'Distribution_Check_Histograms.png' and 'Distribution_Check_QQ.png'.\n")

