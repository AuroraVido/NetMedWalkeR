# NetMedWalkeR: Network-Based Drug Repurposing Pipeline

## Overview
**This repository contains an R-based computational pipeline for network-based drug repurposing developed for my thesis. While originally applied to Arrhythmogenic Cardiomyopathy (ACM), this pipeline is highly adaptable to any other disease**. By simply replacing the input list of disease-associated genes, the algorithm recalculates topological and context-specific relevance scores to identify novel therapeutic candidates for the new target condition. It integrates protein-protein interaction (PPI) networks, disease-gene associations, and drug-target databases to compute context-specific relevance scores and extract disease modules. 

## Repository Structure

### Scripts
* **`config.R`**: Configuration file containing global parameters, file paths, seed values for reproducibility, and specific disease context settings.
* **`functions.R`**: Core algorithmic library. Contains functions for graph topology analysis, Largest Connected Component (LCC) extraction, Jaccard-based context scoring, and significance testing.
* **`main.R`**: The main execution script. It loads the data, prepares the topological graph, executes the Random Walk with Restart (RWR) algorithm, evaluates statistical separation, and extracts the final disease module for downstream pathway analysis (e.g., via Metascape).
**CRITICAL: Do not clear your R workspace or restart your R session after execution of `main.R` if you want to run `Extra/TargetsAnalysis.R`.**
* **`Extra/TargetsAnalysis.R`**: Supplementary script for downstream network and statistical analysis of the top candidate drug targets (*long execution time depending on your CPU capabilities*)

### Data Files (Inputs)
* **`ACM genes.txt` (or custom disease genes)**: List of genes associated with the disease of interest. Replace this file with the gene list of any other pathology to adapt the pipeline.
* **`interactome_mapped.txt`**: The baseline protein-protein interaction (PPI) network downloaded from STRING (nodes and weighted edges). You can replace this with any other custom interactome if needed.
* **`farmaci_filtrati_clean.txt`**: Cleaned database mapping drugs to their specific target genes, sourced from DrugBank. You can substitute this with your own curated list of drugs and corresponding targets.
* **`Phenopedia.txt`**: Disease ontology and phenotype database used for context scoring.

## Prerequisites
To run this pipeline, you need **R** installed along with the following primary libraries:
* `igraph` (Network manipulation and topological calculations)
* `dplyr`, `tidyr`, `data.table` (Data wrangling)
* `pROC` (Statistical evaluations)
* `ggplot2`, `ggrepel` (Visualization)
* `future`, `future.apply` (Parallel processing)

Install missing packages using `install.packages("package_name")`.

## Usage
1. Clone this repository to your local machine.
2. Open `config.R` to update the input file names (e.g., if analyzing a new disease) and ensure the file paths and CPU core allocations (`N_CORES`) match your system environment.
3. Execute `main.R`.

*Note: The script utilizes parallel processing for statistical permutation tests (`N_PERMS = 10000`). Execution time will vary heavily depending on your CPU capabilities.*

## Output
Running `main.R` will generate:
1. **Topological and Validation Plots**: Plot of the network scale-free distribution and plot of the distribution of drugs based on the Network Separation Metric.
2. **Candidate Drugs List (.csv)**: A final dataset containing the statistically significant drug candidates identified for repurposing.
