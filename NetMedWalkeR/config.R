# ==============================================================================
# CONFIGURATION FILE
# ==============================================================================
# This file contains all input parameters and file paths
# Modify this section to adapt the analysis to a different disease

# --- 1. Input File Paths ---
# Path to the list of disease-associated genes (e.g. ACM genes)
FILE_DISEASE_GENES <- "ACM genes.txt"

# Path to the protein-protein interaction network
FILE_INTERACTOME   <- "interactome_mapped.txt"

# Path to the drug-target database
FILE_DRUG_TARGETS  <- "farmaci_filtrati_clean.txt"

# Path to the disease ontology/phenotype database (e.g. Phenopedia)
FILE_DISEASE_DB    <- "Phenopedia.txt"

# --- 2. Analysis Parameters ---
# Global seed for reproducibility
SEED_VALUE <- 42

# Number of permutations for statistical validation
N_PERMS <- 10000

# Number of CPU cores to use (leave 1 core free for OS)
# Set to specific number (e.g., 4) if preferred
N_CORES <- parallel::detectCores() - 1

# --- 3. Context & Disease Specifics ---
# Name of the specific disease to exclude from Context Scoring 
# (To prevent self-loops in Jaccard similarity)
# Example: "Arrhythmogenic Right Ventricular Dysplasia" for ACM
DISEASE_TO_EXCLUDE <- "Arrhythmogenic Right Ventricular Dysplasia"

# --- 4. RWR Optimization ---
# Grid of damping factors to test during Leave-One-Out Cross-Validation
DAMPING_GRID <- seq(0.15, 0.85, by = 0.05)
