#!/bin/Rscript
library(getopt)
library(stringr)
source("/data02/hukaijie/EpiAllele/final_script/chromatin_profiling/05_modified_trackplot.R")
sgRNA_peaks = "/data02/hukaijie/EpiAllele/data/sgRNA_info/sgRNA.bed"
source("/data02/hukaijie/EpiAllele/final_script/RcolorSet.R")


spec <- matrix(c(
  'merge_path', 'm', 1, "character",
  'rep_path', 'r', 1, "character",
  'results_path', 'o', 1, "character",
  'trackplot_name', 'n', 1, "character",
  'group', 'g', 1, "character",
  'patterns', 'p', 1, "character",
  'genome', 'e', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec =spec)
merge_path <- opt$merge_path
rep_path <- opt$rep_path
result.path <- opt$results_path
trackplot_name <- opt$trackplot_name
group <- opt$group
patterns <- strsplit(opt$patterns, ",")[[1]]
genome <- opt$genome

bw_files <- c()
for (pattern in patterns) {
  rep_bw_files <- list.files(rep_path, pattern = paste0(group, ".*", pattern), full.names = TRUE)
  merge_bw_files <- list.files(merge_path, pattern = paste0(group, ".*", pattern), full.names = TRUE)
  bw_files <- c(bw_files, rep_bw_files, merge_bw_files)
}

bigWigs <- read_coldata(bws = bw_files, build = genome, 
                        sample_names = sub("(.*).bw", "\\1", basename(bw_files)))
bigWigs$modification <- gsub("_rep[1-3]", "", bigWigs$bw_sample_names)
bigWigs$group <- sub(".*\\.(.*)", "\\1", bigWigs$bw_sample_names)

pos1 <- "chr14:55177378-55208581"
pos2 <- "chr14:55194524-55208581"
pos3 <- "chr14:55199278-55208581"
region_list <- c(pos1, pos2, pos3)
bigWigs$color <- colorset[["rep_track_col"]][bigWigs$group]

retry_track_extract <- function(colData, loci, max_attempts = 5, delay = 5) {
  attempt <- 1
  while (attempt <= max_attempts) {
    try({
      t <- track_extract(colData = colData, loci = loci)
      return(t)
    }, silent = TRUE)
    message(paste("Attempt", attempt, "failed. Retrying in", delay, "seconds..."))
    Sys.sleep(delay)
    attempt <- attempt + 1
  }
  stop("All attempts to extract track data failed.")
}

pdf(paste0(result.path, "/", trackplot_name, ".pdf"), width = 10, height = length(bw_files)*1)
for (region in region_list) {
  t <- retry_track_extract(colData = bigWigs, loci = region)
  track_plot(summary_list = t, 
             show_axis = TRUE,
             y_min = 0,
             condition = "modification",
             groupScaleByCondition = TRUE,
             col = bigWigs$color, 
             show_ideogram = TRUE, peaks = sgRNA_peaks)
}
dev.off()