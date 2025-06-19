#!/usr/bin/env Rscript
# Define and parse arguments ----
suppressWarnings(library("getopt", quietly = T, character.only = T))

## Define the argument specifications (spec matrix) ----
spec <- matrix(c(
  'help', 'h', 0, 'logical', 'Show this help message and exit',
  'path_RawCGdata', 'c', 1, "character", 'Path for call CG results',
  'path_Sample_info_tab', 'i', 1, "character", 'Path for sampleinfo file',
  'path_Result',  'r', 1, "character", 'Path for output',
  'species', 's', 1, "character", 'species, human or mouse',
  'w_dml', 'w', 1, "integer", 'Width for dml plot',
  'h_dml', 'v', 1, "integer", 'Height for dml plot',
), byrow=TRUE, ncol=5)



## Create a funtion to generate help information ----
## =========================================================================
show_custom_help <- function(spec) {
  cat("\n=============================================================================\n")
  
  cat("\nDescription:\n")
  cat("  This script performs a comprehensive Differentially Methylated Region (DMR)\n")
  cat("  analysis from raw methylation call data. It identifies significant DMRs,\n")
  cat("  annotates them to genomic features, and generates visualizations including\n")
  cat("  a Manhattan plot for genome-wide overview.\n")
  
  cat("\n-----------------------------------------------------------------------------\n")
  # This line automatically prints the usage and arguments from the 'spec' matrix
  cat(getopt(spec, usage = TRUE))
  cat("-----------------------------------------------------------------------------\n")
  
  cat("\nExample:\n")
  cat("  Rscript your_DMR_script.R \\\n")
  cat("    --path_RawCGdata /data02/hukaijie/EpiAllele/result/20250212_wgbs/CG \\\n")
  cat("    --path_Sample_info_tab /data02/hukaijie/EpiAllele/result/20250212_wgbs/D5_wgbs_sampleinfo.txt \\\n")
  cat("    --path_Result /data02/hukaijie/EpiAllele/result/20250212_wgbs/DMR \\\n")
  cat("    --species mouse \\\n")
  cat("    --w_dml 7 --h_dml 3\n")
  
  cat("\nOutput:\n")
  cat("  All results are saved in the directory specified by --path_Result.\n")
  cat("  The primary outputs include:\n\n")
  
  cat("  1. DMR Result Tables:\n")
  cat("     - A plain text or .csv file containing statistically significant DMRs\n")
  cat("       with columns for location (chr, start, end), strand, p-value,\n")
  cat("       q-value (FDR), and methylation difference between groups.\n\n")
  
  cat("  2. Annotated DMRs:\n")
  cat("     - An extended table where each DMR is annotated with its relationship\n")
  cat("       to genomic features (e.g., promoter, exon, intron, intergenic) and\n")
  cat("       the nearest gene(s), based on the specified species.\n\n")
  
  cat("  3. Visualizations (in .png or .pdf format):\n")
  cat("     - Manhattan Plot: Displays the genome-wide significance (-log10(p-value))\n")
  cat("       of differentially methylated sites across all chromosomes.\n")

  cat("\n=============================================================================\n\n")
}

## Parse the arguments and call the custom help function when needed ----
opt <- getopt(spec=spec)

if (!is.null(opt$help) || is.null(opt$result_path) || is.null(opt$species)) {
  # Call our custom help function
  show_custom_help(spec)
  q(status = 1) # Exit the script
}

path_RawCGdata = opt$path_RawCGdata
path_Sample_info_tab = opt$path_Sample_info_tab
path_Result = opt$path_Result
species = opt$species
w_dml = opt$w_dml
h_dml = opt$h_dml


# library packages and load QC.funtions ----
for (p in c("data.table", "tidyverse", "magrittr", "TxDb.Hsapiens.UCSC.hg38.refGene", 
            "TxDb.Mmusculus.UCSC.mm39.refGene", "ChIPseeker", "DSS", "bsseq", "patchwork","ggpubr")) {
  if (!require(p, character.only = T)) {
    stop(paste0("Please install ", p))
  } else {
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}




# QC plot ----
# files_use = list.files(path_RawCGdata)
# list_QC = lapply(files_use,function(f){
 
#  dat_raw = fread(file.path(path_RawCGdata,f),data.table = F) %>% 
#    set_colnames(c('chr','coordinate','strand','context','unconverted','converted')) %>% 
#    mutate(chr = as.character(chr),coordinate = as.integer(coordinate)) %>% 
#    mutate(all = unconverted+converted) %>% 
#    mutate(ratio = unconverted/all) %>% 
#    mutate(context = factor(context,levels = c('CH','CG')))
 
#  p_ratio = ggplot(dat_raw ,aes(context,ratio,color = context))+
#    geom_boxplot()+
#    # geom_jitter()+
#    theme_classic(base_size = 14)+
#    labs(x = NULL)+
#    theme(legend.position = 'none')+
#    ggtitle(f)
 
#  return(p_ratio)
 
# })

# png(paste0(path_Result,'/','QC.png'),width = 8,height = 6,res = 300,units = 'in')
# print(wrap_plots(list_QC))
# dev.off()

# data processing ----
tab_SampleInfo = fread(path_Sample_info_tab,colClasses = list(character = "Sample_name",character = "Sample_group"))

file_path = paste0(tab_SampleInfo$Sample_name,'.cg')

allDat <- lapply(file_path,function(f){
  # f="C-1.cg";
  print(f);

  tmp=fread(file.path(path_RawCGdata,f))
  colnames(tmp)<-c('chr','coordinate','strand','context','unconverted','converted')
  tmp$chr=as.character(tmp$chr)
  tmp$coordinate=as.integer(tmp$coordinate)
  tmp$all=as.numeric(tmp$unconverted)+as.numeric(tmp$converted)
  tmp$methy_ratio=as.numeric(tmp$unconverted)/as.numeric(tmp$all)

  newtmp<-tmp[,c('chr','coordinate','all','unconverted')]

  colnames(newtmp)=c('chr', 'pos' ,'N' ,'X')

  return(newtmp)
})


BSobj <- makeBSseqData(allDat,tab_SampleInfo$Sample_name)

if (species=='human') {
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.refGene
  
  anndb ="org.Hs.eg.db"
  uniq_chr = paste0('chr',c(1:22,'X','Y'))
  
}else if(species=='mouse'){
  txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
  
  anndb ="org.Mm.eg.db"
  uniq_chr = paste0('chr',c(1:19,'X','Y'))
}else{
  stop(paste0(species, " must be human or mouse"))
}

groups <- relevel(factor(tab_SampleInfo$Sample_group), ref = "sg")
comparisons <- combn(unique(levels(groups)), 2, simplify = F)
for (comp in comparisons) {
  dmlTest <- DMLtest(BSobj,
                    group1=tab_SampleInfo %>% dplyr::filter(.,Sample_group==comp[1]) %>% pull(Sample_name),
                    group2=tab_SampleInfo %>% dplyr::filter(.,Sample_group==comp[2]) %>% pull(Sample_name),
                    smoothing=T,smoothing.span = 500)
  dmrs <- callDMR(dmlTest, delta = 0.2, p.threshold = 0.01)

  saveRDS(dmlTest, paste0(path_Result, '/', comp[1], ".vs.", comp[2], '.dmlTest.rds'))
  saveRDS(dmrs, paste0(path_Result, '/', comp[1], ".vs.", comp[2], '.dmrs.rds'))

  se_dmrs = data.frame(dmrs)

  dmlTest_df <- data.frame(dmlTest) %>% 
    filter(.,chr%in%uniq_chr) %>% 
    mutate(chr = factor(chr,levels= uniq_chr)) %>% 
    mutate(logP = ifelse(diff>0, -log10(fdr), log10(fdr))) %>% 
    mutate(group = ifelse((fdr<0.05 & diff>0.2), "UP",
                          ifelse((fdr < 0.05 & diff < -0.2), "DOWN", "NOT-SIG"))) %>% 
    mutate(group = ifelse((chr=="chr14" &
                            pos>=55169378 &
                            pos<=55214384), "Myh6-related", group))
  
  dml_main <- dmlTest_df %>% filter(group != "Myh6-related")
  dml_myh6 <- dmlTest_df %>% filter(group == "Myh6-related")

  ## whole genome plot (pvalue) ----
  p1 <- ggplot() +
    geom_point(data = dml_main,
              aes(x = pos, y = logP, color = group, alpha=abs(logP), size=abs(logP))) +
    geom_point(data = dml_myh6,
              aes(x = pos, y = logP, color = group,
              alpha=abs(logP), size=abs(logP))) + 
    scale_size(range = c(0.5, 1), guide = "none") +
    scale_alpha(range = c(0.5, 1), guide = "none") +
    scale_color_manual(values = c("DOWN" = 'steelblue', 
                                  "UP" = 'firebrick',
                                  "NOT-SIG" = "gray",
                                  "Myh6-related"="green4")) +
    facet_wrap(~chr, nrow = 1, scales = 'free_x') +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 3), nrow = 1)) +
    labs(x = NULL, y = '-log10(pvalue)')
  
  ## whole genome dotplot (diff) ----
  p2 <- ggplot() +
    geom_point(data = dml_main,
              aes(x = pos, y = diff, color = group, alpha=abs(logP), size=abs(logP))) +
    geom_point(data = dml_myh6,
              aes(x = pos, y = diff, color = group,
              alpha=abs(logP), size=abs(logP))) + 
    scale_size(range = c(0.5,1),guide = "none")+
    scale_alpha(range = c(0.5,1),guide = "none")+
    scale_color_manual(values = c("DOWN"='steelblue', 
                                  "UP"='firebrick',
                                  "NOT-SIG"="gray",
                                  "Myh6-related"="green4"))+
    facet_wrap(~chr,nrow=1,scales = 'free_x')+
    theme_classic(base_size = 8)+
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    guides(color=guide_legend(override.aes = list(size=3),nrow = 1))+
    theme(legend.position = "bottom")+
    labs(x=NULL,y='DNA methylation difference')+
    geom_hline(yintercept = 0.1, linetype = 'dashed')+
    geom_hline(yintercept = -0.1, linetype = 'dashed')

  
  png(paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.dmlTest.p_val.png'),width = w_dml, height=h_dml, res = 300, units = 'in')
  print(p1)
  dev.off()

  png(paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.dmlTest.diff.png'),width = w_dml, height=h_dml, res = 300, units = 'in')
  print(p2)
  dev.off()

  # pdf(paste0(path_Result,'/', comp[1], ".vs.", comp[2], '.dmlTest.pdf'),width = w_dml,height=h_dml)
  # print(p1)
  # print(p2)
  # dev.off()

  # ## volcano plot ----
  # dmls_in_dmrs <- purrr::map_dfr(1:nrow(se_dmrs), function(i) {
  #   dmr <- se_dmrs[i, ]
  #   dmlTest_df %>%
  #     filter(
  #       chr == dmr$chr,
  #       pos >= dmr$start,
  #       pos <= dmr$end
  #     )
  # })
  # dmls_in_dmrs <- dmls_in_dmrs %>% 
  #   mutate(CG_pos=paste0(chr, ": ", pos),
  #         logP_real=abs(logP))
  # Myh6_related_CGs <- dmls_in_dmrs %>% 
  #   filter(group=="Myh6-related") %>% 
  #   arrange(desc(abs(diff))) %>%
  #   head(5) %>% 
  #   pull(CG_pos)
  # p <- ggscatter(dmls_in_dmrs, x = "diff", y = "logP_real",
  #                 color = "group", 
  #                 size = 1.6,
  #                 shape = 19,
  #                 label = "gene",
  #                 label.select = Myh6_related_CGs, 
  #                 font.label = c(14,'black'), 
  #                 repel = T,
  #                 xlab = "DNA methylation difference", 
  #                 ylab = "-log10(P.value)",
  #                 title = paste0(comp[1], ".vs.", comp[2])) +
  #   geom_vline(xintercept = c(-0.2, 0.2), lty = 4, col = "black", lwd = 0.8) +
  #   geom_hline(yintercept = -log10(0.01), lty = 4, col = "black", lwd = 0.8) +
  #   scale_color_manual(values = c("UP" = 'firebrick', 
  #                               "DOWN" = 'steelblue',
  #                               "NOT-SIG" = "gray",
  #                               "Myh6-SIG" = "limegreen")) +
  #   theme_classic() +
  #   theme(axis.title = element_text(size = 15),
  #         axis.text = element_text(size = 14),
  #         legend.title = element_text(size = 14),
  #         legend.text = element_text(size = 13),
  #         plot.title = element_text(size = 20, face = "bold"))
  
  # png(paste0(path_Result,'/',comp[1], ".vs.", comp[2], ".dmls.in.DMRs.volcano.png"),
  #     width = 6, height = 5, res = 300, units = 'in')
  # print(p)
  # dev.off()

  # pdf(paste0(path_Result,'/',comp[1], ".vs.", comp[2], ".dmlTest.volcano.pdf"),
  #     width = 6, height = 5)
  # print(p)
  # dev.off()


  ## DMR annotation ----
  if(nrow(se_dmrs)!=0){
    
    write.table(se_dmrs,paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.dmrs_info.txt'),
                quote = F,row.names = F,col.names = T,sep = '\t')
    write.csv(se_dmrs, paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.dmrs_info.csv'), 
              quote = F,row.names = F)
    gr_DMR <- GRanges(seqnames=se_dmrs$chr,
                      ranges=IRanges(start=as.integer(se_dmrs$start), 
                                      end=as.integer(se_dmrs$end)))
    
    anno_df <- annotatePeak(gr_DMR, tssRegion=c(-500, 500), TxDb=txdb, annoDb=anndb,
                            columns = c('ENSEMBL',"SYMBOL",'GENENAME'),
                            level = "transcript") %>% 
      as.data.frame() %>% 
      mutate(anno_type = str_split(annotation,' \\(',simplify = T)[,1]) %>% 
      na.omit()

    write.table(anno_df,paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.anno_dmrs.bed'), 
                quote = F,row.names = F,col.names = T,sep = '\t')
    write.csv(anno_df, paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.anno_dmrs.csv'))
    
  }else if (nrow(se_dmrs)==0){
    file.create(paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.dmrs_info.txt'))
    # file.create(paste0(path_Result,'/','dmrs_CON_EXP.bed'))
    file.create(paste0(path_Result,'/',comp[1], ".vs.", comp[2], '.anno_dmrs.bed'))
  }
}



