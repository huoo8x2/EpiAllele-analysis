#!/usr/bin/env Rscript
#
# =============================================================================
#
# Description:
#   This script performs comprehensive quality control (QC) and preprocessing
#   for single-nucleus RNA-seq (snRNA-seq) data. It takes Cell Ranger
#   output and generates a cleaned, normalized, and clustered Seurat object.
#
# Analysis Workflow:
#   1.  Load raw and filtered count matrices from Cell Ranger output.
#   2.  Create an initial Seurat object.
#   3.  Filter low-quality cells based on gene count, UMI count, and
#       mitochondrial content.
#   4.  Perform standard Seurat processing workflow: normalization, scaling, PCA,
#       and clustering (UMAP).
#   5.  Save the final, clean Seurat object to an .rds file for downstream
#       analysis.
#   6.  Run SoupX to estimate and remove ambient RNA contamination.
#   7.  Run DoubletFinder to identify and remove potential doublets.
#
# Required Arguments:
#   --sample / -s           The name/identifier for the sample (e.g., "WT").
#   --project_path / -p     The root path for the project directory. Output will be saved under
#                           'project_path/results/snRNAseq'.
#   --matrix_path / -m      Parent directory containing Cell Ranger output for ALL samples.
#                           The script expects a structure like:
#                           'matrix_path/<sample_name>/outs/raw_feature_bc_matrix'.
#
# Example:
#   Rscript 03_runQC.R \
#     --sample WT \
#     --project_path ~/project/EpiAllele \
#     --matrix_path ~/project/EpiAllele/result/snRNAseq/matrix
#
# Output:
#   Three subdirectories (soupx, qc, doubletfinder) are created, each containing
#   results from a specific stage of the QC process.
#
#   1. qc/
#      - An .rds file of the Seurat object after standard filtering 
#        (low-quality cells were filtered based on the following criteria: 
#        fewer than 200 or more than the 99th percentile of detected gene counts, 
#        more than 20,000 UMI counts, and more than 20% of mitochondrial content).
#      - Violin plots showing nFeature_RNA, nCount_RNA, and mitochondrial
#        percentage for QC assessment.
#      - UMAPs showing Immunoglobulin (IG) and Hemoglobin (HB) gene expression
#   2. soupx/
#      - An .rds file of the Seurat object after SoupX ambient RNA correction.
#      - Change plots of IG and HB gene expression after performing SoupX
#   3. doubletfinder/
#      - A UMAP plot visualizing the detected doublets.
#      - The final cleaned Seurat object saved as an .rds file after
#        doublet removal. The main output for downstream analysis.
#
# =============================================================================


# Define and parse arguments ----
suppressWarnings(library("getopt", quietly = T, character.only = T))

## Define the argument specifications (spec matrix) ----
spec <- matrix(c(
  'help',         'h', 0, 'logical',   'Show this help message and exit',
  'sample',       's', 1, "character", 'Sample identifier to process (e.g., "WT")',
  'project_path', 'p', 1, "character", 'Root path for the project',
  'matrix_path',  'm', 1, "character", 'Parent directory for CellRanger outputs'
), byrow=TRUE, ncol=5)


## Create a funtion to generate help information ----
show_custom_help <- function(spec) {
  cat("\n=============================================================================\n")
  cat("\nDescription:\n")
  cat("  This script performs comprehensive QC and preprocessing for snRNA-seq data\n")
  
  
  cat("\n-----------------------------------------------------------------------------\n")
  # This line automatically prints the usage and arguments from the 'spec' matrix
  cat(getopt(spec, usage = TRUE))
  cat("-----------------------------------------------------------------------------\n")
  
  cat("\nExample:\n")
  cat("  Rscript your_script_name.R \\\n")
  cat("    --sample WT \\\n")
  cat("    --project_path ~/project/EpiAllele \\\n")
  cat("    --matrix_path ~/project/EpiAllele/result/snRNAseq/matrix\n")
  
  cat("\nOutput:\n")
  cat("  Three subdirectories (qc, soupx, doubletfinder) are created, each containing\n")
  cat("  results from each QC step.\n\n")
  
  cat("  1. qc/\n")
  cat("     - An .rds file of the Seurat object after standard filtering.\n")
  cat("       (Criteria: gene counts < 200 or > 99th percentile, UMI > 20,000,\n")
  cat("        and mitochondrial content > 20%).\n")
  cat("     - Violin plots showing nFeature_RNA, nCount_RNA, and mitochondrial\n")
  cat("       percentage for QC assessment.\n")
  cat("     - UMAPs showing Immunoglobulin (IG) and Hemoglobin (HB) gene expression.\n\n")
  
  cat("  2. soupx/\n")
  cat("     - An .rds file of the Seurat object after SoupX ambient RNA correction.\n")
  cat("     - Change plots of IG and HB gene expression after performing SoupX.\n\n")
  
  cat("  3. doubletfinder/\n")
  cat("     - A UMAP plot visualizing the detected doublets.\n")
  cat("     - The final, cleaned Seurat object saved as an .rds file after\n")
  cat("       doublet removal. This is the main output for downstream analysis.\n")
  
  cat("\n=============================================================================\n\n")
}


## Parse the arguments and call the custom help function when needed ----
opt <- getopt(spec=spec)

if (!is.null(opt$help) || is.null(opt$result_path) || is.null(opt$species)) {
  # Call our custom help function
  show_custom_help(spec)
  q(status = 1) # Exit the script
}

sample <- opt$sample
proj_path <- opt$project_path
matrix_path <- opt$matrix_path


# library packages and load QC.funtions ----
#remotes::install_github('lzmcboy/DoubletFinder_204_fix')
#remotes::install_github('MaxMeieran/DoubletFinder')
for (p in c("tidyverse", "Seurat", "SoupX", "ggsci", "ggpubr", 
            "DoubletFinder", "ggrepel", "clustree", "homologene")) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}

source(paste0(proj_path, "/script/snRNAseq/03_QC_functions.R"))

# quality check ----
result.path <- paste0(proj_path, "/result/snRNAseq/qc/")
if (!dir.exists(result.path)) {
  dir.create(result.path)
}
runQC(sample, matrix_path, result.path)

# Soupx to remove ambient RNA ----
result.path <- paste0(proj_path, "/result/snRNAseq/soupx/")
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

RunSoupX(sample, matrix_path, result.path)

# DoubletFinder to detect and remove doublets ----
soupx.path <- paste0(proj_path, "/result/snRNAseq/soupx/")
result.path <- paste0(proj_path, "/result/snRNAseq/doubletfinder/")
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

RunDoubletFinder(sample, data.path, result.path)


# remove doublets and perform standard cluster steps ----
doubletfinder.path <- paste0(proj_path, "/result/snRNAseq/doubletfinder/")
result.path <- paste0(proj_path, "/result/snRNAseq/clustering/")
if (!dir.exists(result.path)) {
  dir.create(result.path)
}

RunCluster(sample, data.path, result.path, 
           checkAAV = T, Findmarker = T)


