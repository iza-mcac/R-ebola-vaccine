rm(list = ls())
pkgs <- c('tidyverse', 'SummarizedExperiment','ggplot2', 'CEMiTool','DESeq2')
sapply(pkgs, require,character.only = T)

basedir <- '~/Documents/Projects/CSBL/Ebola/Ebola_lncRNA/'
setwd(basedir)

inputdir <- 'data_preprocessing/SExperiments/'
SEs <- list.files(inputdir, pattern = 'RData')

for(SE_name in SEs){ # E_name <- SEs[2]
  message('\n\n\t',SE_name,'\n')
  load(file = file.path(inputdir, SE_name))
  dataset <- metadata(SE)$Cohort
  outdir = file.path("./CoExpression/CEMiTool/",dataset)

  SE <- SummarizedExperiment(assays = as.matrix(assay(SE)),
                             colData = colData(SE),
                             rowData = rowData(SE),
                             metadata = metadata(SE))

  dds <- DESeqDataSet(SE, design =~1)
  norm.data <- vst(dds) %>%
    assay %>% as.data.frame()
  message('Normalized Exprs (VST):')
  print(norm.data[1:5,1:4])

  message('Pheno Data:')
  pheno <- colData(SE) %>% as.data.frame %>%
    rownames_to_column('SampleName') %>%
    select(SampleName, sampling_day, treatment) %>%
    mutate(Class = paste0('D',sampling_day),
           Class = paste0(Class,'_',treatment)) %>%
    select(SampleName, Class)
  print(head(pheno))

  # gmt
  # int_df =

  # run cemitool
  message('\nRunning CEMiTool...')
  cem <- cemitool(norm.data, pheno,
                  # gmt_in = gmt,
                  # interactions=int_df,
                  filter=TRUE, plot=TRUE, verbose=TRUE)

  message('\nSaving Report...')
  generate_report(cem, directory=paste0(outdir,'/Report'),force=T)

  message('\nSaving Tables...')
  write_files(cem, directory=paste0(outdir,'/Tables'),force=T)

  message('\nSaving Plots...')
  save_plots(cem, "all", directory=paste0(outdir,'/Plots'),force=T)
  message('\nDONE!\n\n\n\n')
}
