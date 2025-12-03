library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggplot2)
options(stringsAsFactors = F)

#Load Geneva antibodies - lncRNA
geneva_antibody_lncRNA <- readRDS("Corr_analysis/geneva_antibody_lncRNA.RDS")
USA_antibody_lncRNA <- readRDS("Corr_analysis/USA_antibody_lncRNA.RDS")

geneva <- geneva_antibody_lncRNA %>%
  dplyr::select(1,2,3,7,5,8)

USA <- USA_antibody_lncRNA %>%
  dplyr::filter(antibody_type=="Ebola Zaire Virus GP IgG Antibody Titer") %>%
  dplyr::mutate(batch=3) %>%
  dplyr::select(1,2,3,5,6,7)

colnames(USA) <- colnames(geneva)

USA_and_geneva <- rbind(USA,geneva)

p <- USA_and_geneva %>%
  dplyr::filter(lncRNA %in% c("FAM225A"),day==14) %>%
  ggplot(aes(x=D1_DE,y=log2(IgGTitter+1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()
pdf(file = "Corr_analysis/figures/FAM225A_IgG_D14.pdf",
    width = 2.7,height = 2.7)
print(p)
dev.off()

p <- USA_and_geneva %>%
  dplyr::filter(lncRNA %in% c("DOCK8-AS1"),day==168) %>%
  ggplot(aes(x=D1_DE,y=log2(IgGTitter+1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()
pdf(file = "Corr_analysis/figures/DOCK8-AS1_IgG_D168.pdf",
    width = 2.7,height = 2.7)
print(p)
dev.off()

p <- USA_and_geneva %>%
  dplyr::filter(lncRNA %in% c("LEF1-AS1"),day==730) %>%
  ggplot(aes(x=D1_DE,y=log2(IgGTitter+1)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()
pdf(file = "Corr_analysis/figures/LEF1-AS1_IgG_D730.pdf",
    width = 2.7,height = 2.7)
print(p)
dev.off()
