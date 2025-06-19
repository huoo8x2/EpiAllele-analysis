#!/usr/bin/env Rscript

## {{{ Install required libraries if required
for (p in c("data.table", "tidyverse", "DSS", "bsseq", "patchwork",
            "ggpubr", "getopt")) {
  if (!require(p, character.only = T)) {
    print(paste0("Please install ", p))
  }else{
    suppressMessages(library(p, quietly = T, character.only = T))
  }
}
## }}}

##############arguments from command line#######################
spec <- matrix(c(
  'CG_path', 'c', 2, "character",
  'sample_info', 's', 2, "character",
  'result_path', 'o', 2, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec =spec)
path_RawCGdata <- opt$CG_path
path_Sample_info_tab <- opt$sample_info
path_Result <- opt$result_path


############# data processing ######################
tab_SampleInfo = fread(path_Sample_info_tab,colClasses = list(character = "Sample_name",character = "Sample_group"))

for (type in c("c57", "dba", "total")) {
  sampleInfo <- tab_SampleInfo %>% mutate(Sample_name=paste0(Sample_name, ".", type))
  file_path = paste0(sampleInfo$Sample_name,'.cg')
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

  BSobj <- makeBSseqData(allDat, sampleInfo$Sample_name)
  BSobj <- BSobj[rowMeans(getCoverage(BSobj)) >= 3, ]

  groups <- relevel(factor(sampleInfo$Sample_group), ref = "sg")
  comparisons <- combn(unique(levels(groups)), 2, simplify = F)

  for (comp in comparisons) {
    dmlTest <- DMLtest(BSobj,
                      group1=sampleInfo %>% dplyr::filter(.,Sample_group==comp[1]) %>% pull(Sample_name),
                      group2=sampleInfo %>% dplyr::filter(.,Sample_group==comp[2]) %>% pull(Sample_name),
                      smoothing=F)

    dmlTest_df <- data.frame(dmlTest) %>% 
      mutate(logP = ifelse(diff>0, -log10(fdr), log10(fdr))) %>% 
      mutate(group = ifelse((fdr<0.05 & diff>0.2), "UP",
                            ifelse((fdr< 0.05 & diff < -0.2), "DOWN", "NOT-SIG")))
    
    ## whole genome dot plot ----
    pdf(paste0(path_Result, "/", comp[1], ".vs.", comp[2], ".", type, ".dotplot.pdf"), width = 6, height = 4)
    p <- ggplot(dmlTest_df, aes(pos, diff))+
      geom_point(aes(col=group))+
      theme_bw()+
      scale_color_manual(values=c("UP"="red3", "DOWN"="blue4", "NOT-SIG"="gray"))+
      geom_hline(yintercept = 0)+
      labs(x=NULL,y='DNA methylation difference')+
      geom_hline(yintercept = 0.2, linetype = 'dashed')+
      geom_hline(yintercept = -0.2, linetype = 'dashed')
    print(p)

    p <- ggplot(dmlTest_df, aes(pos, logP))+
      geom_point(aes(col=group))+
      theme_bw()+
      scale_color_manual(values=c("UP"="red3", "DOWN"="blue4", "NOT-SIG"="gray"))+
      geom_hline(yintercept = 0)+
      labs(x=NULL,y='-log10(pvalue)')+
      geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
      geom_hline(yintercept = log10(0.05), linetype = 'dashed')
    print(p)
    dev.off()

    write.csv(dmlTest_df, file = paste0(path_Result, "/", comp[1], ".vs.", comp[2], ".", type, ".dmlTest.csv"), row.names = F)
  }
}
