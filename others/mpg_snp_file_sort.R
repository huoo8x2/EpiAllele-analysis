#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c("tidyverse", "stringr", "getopt")) {
  if (!require(p, character.only = T)) {
    print(paste0("Please install ", p))
  }else{
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}


spec <- matrix(c(
  'snp_file', 's', 2, "character",
  'out_file', 's', 2, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec =spec)
snp_file <- opt$snp_file
out_file <- opt$out_file


base_list <- list(R=c("A", "G"),
                  Y=c("A", "C"),
                  K=c("G", "T"),
                  S=c("C", "G"),
                  H=c("A", "C", "T"),
                  B=c("C", "G", "T"),
                  V=c("A", "C", "G"),
                  N=c("A", "C", "G", "T"))


SNP <- read.csv(snp_file, header = T, check.names = F)
SNP <- SNP %>% select(!`Group A`&!`Group B`) %>% 
  mutate(`DBA/2J`=str_replace(`DBA/2J`, "\\|null", ""),
         `C57BL/6J`=str_replace(`C57BL/6J`, "\\|null", "")) %>% 
  mutate(across(everything(), ~ str_replace_all(., "^\\s*$", NA_character_))) %>% 
  drop_na()

for (i in 1:nrow(SNP)) {
  dba_base <- SNP$`DBA/2J`[i]
  c57_base <- SNP$`C57BL/6J`[i]
  if (dba_base %in% names(base_list)) {
    SNP$`DBA/2J`[i] <- paste(base_list[[dba_base]], collapse = "/")
  }
  if (c57_base %in% names(base_list)) {
    SNP$`C57BL/6J`[i] <- paste(base_list[[c57_base]], collapse = "/")
  }
  
  dba_base <- strsplit(SNP$`DBA/2J`[i], "/")[[1]]
  c57_base <- strsplit(SNP$`C57BL/6J`[i], "/")[[1]]
  
  if (length(intersect(dba_base, c57_base)) > 0) {
    SNP$`DBA/2J`[i] <- NA
    SNP$`C57BL/6J`[i] <- NA
  }
}
SNP <- drop_na(SNP)
write.csv(SNP, out_file)



