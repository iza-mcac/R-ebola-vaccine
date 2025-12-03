library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(metafor)
library(igraph)

#load lncRNA x IgG meta analysis data
lnc_ab_meta <- fread("data/Metafor_results/Metafor_IgG_perDay/RandomEffects_Model_Corr_lncRNAs_IgG_perDay_AllCohorts_v3.csv")

lnc_ab_meta <- lnc_ab_meta %>%
  dplyr::mutate(sig=ifelse(pval<0.05,yes = "sig",no = "not.sig"),
         dir=ifelse(beta>0,yes = "positive",no = "negative"),
         sig_dir=ifelse(sig=="sig" & dir=="positive",
                        yes = "positive",
                        no = ifelse(sig=="sig",
                                    yes = "negative",
                                    no = "not.sig"))) %>%
  dplyr::filter(!duplicated(paste0(Gene,day)))

write.csv(lnc_ab_meta,file = "data/lncRNA_antibody_metacorrelation_results.csv",
          row.names = F)

p <- lnc_ab_meta %>%
  dplyr::filter(day!=0) %>%
  ggplot(aes(x = beta,y = -log(pval),color=sig_dir))+
  geom_point(size=3)+
  geom_text_repel(data = lnc_ab_meta %>%
                    dplyr::filter(day!=0,
                           sig=="sig"),
                  aes(label=Gene),nudge_y = 1.2,
                  size=3,force = 1,color="black")+
  theme_bw()+
  geom_hline(yintercept = -log(0.05))+
  scale_color_manual(values = c("blue","grey","red3"))+
  facet_wrap(facets = ~day,nrow = 1)
pdf(file = "Corr_analysis/figures/lncRNA_IgG_metacor_volcano_all_days.pdf",
    width = 8,height = 2.25)
print(p)
dev.off()

#Plot results for matching days: 28, 180 and 360
p <- lnc_ab_meta %>%
  dplyr::filter(day %in% c(28,168,365)) %>%
  ggplot(aes(x = beta,y = -log(pval),color=sig_dir))+
  geom_point(size=3)+
  geom_text_repel(data = lnc_ab_meta %>%
                    dplyr::filter(day %in% c(28,168,365),
                           sig=="sig"),
                  aes(label=Gene),nudge_y = 1.2,
                  size=3,force = 1,color="black")+
  theme_bw()+
  geom_hline(yintercept = -log(0.05))+
  scale_color_manual(values = c("blue","grey","red3"))+
  facet_wrap(facets = ~day,nrow = 1)
p

usa <- USA_antibody_RNAseq_volunteers %>%
  filter(LBTEST=="Ebola Zaire Virus GP IgG Antibody Titer")
table(usa$day)
table(geneva_antibody$day)

cor_individual <- fread("data/files_Durão/data/lncRNA_antibody_correlation_individual_studies_results.csv")
