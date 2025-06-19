#!/bin/Rscript
library(getopt)
library(Cairo)
source("/data02/hukaijie/EpiAllele/final_script/chromatin_profiling/05_modified_trackplot.R")
sgRNA_peaks = "/data02/hukaijie/EpiAllele/data/sgRNA_info/sgRNA.bed"
source("/data02/hukaijie/EpiAllele/final_script/RcolorSet.R")

spec <- matrix(c(
  'data_path', 'd', 1, "character",
  'results_path', 'r', 1, "character",
  'trackplot_name', 'n', 1, "character",
  'groups', 'g', 1, "character",
  'pattern', 'p', 1, "character",
  'genome', 'o', 1, "character",
  "type", "t", 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec =spec)
data.path <- opt$data_path
result.path <- opt$results_path
trackplot_name <- opt$trackplot_name
groups <- strsplit(opt$groups, ",")[[1]]
pattern <- opt$pattern
genome <- opt$genome
type <- opt$type



bw_files <- c()
for (group in groups) {
  bw_file <- list.files(data.path, pattern = paste0(group, ".*", pattern), full.names = TRUE)
  bw_files <- c(bw_files, bw_file)
}

bigWigs <- read_coldata(bws = bw_files, build = genome, sample_names = groups)
bigWigs$modification <- sub("(.*?)_.*.", "\\1", bigWigs$bw_sample_names)
bigWigs$group <- sub(".*_(.*)", "\\1", groups)
if (type %in% c("track_invitro", "track_invivo")) {
   bigWigs$color <- colorset[[type]][bigWigs$group]
} else {
   stop(paste0("type should be track_invitro or track_invivo"))
}


pos1 <- "chr14:55177378-55208581"
pos2 <- "chr14:55194524-55208581"
pos3 <- "chr14:55199278-55208581"
region_list <- c(pos1, pos2, pos3)
height <- ifelse(length(bw_files) < 4, 4, length(bw_files))

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

pdf(paste0(result.path, "/", trackplot_name, ".pdf"), width = 10, height = height)
for (region in region_list) {
  t <- retry_track_extract(colData = bigWigs, loci = region)
  t$data <- lapply(t$data, function(df) {
    df$max <- pmax(0, df$max)  # 把小于 0 的值截断为 0
    df
  })
  attr(t$data, "meta") <- attr(retry_track_extract(colData = bigWigs, loci = region)$data, "meta")
  track_plot(summary_list = t, 
             show_axis=T,
             y_min=0,
             condition="modification",
             groupScaleByCondition=T,
             col = bigWigs$color, 
             show_ideogram = T, peaks =sgRNA_peaks)         
}
dev.off()
