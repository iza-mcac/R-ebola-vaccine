library(biomaRt) 
library(tidyverse)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", 
                host = 'https://www.ensembl.org')

txtogenesv37 <-read_csv("txtogenesv37.csv")

######

DEA_b1_day01_res_joined <- DEA_b1_day01_res %>% left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b1_day01_res_joined_na <-DEA_b1_day01_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b1_day01_res_joined_na$X,
                        mart = mart)

reannot_b1_d1 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b1_day01_res_reannot <- DEA_b1_day01_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b1_d1$external_gene_name)

write_csv2(DEA_b1_day01_res_reannot, "DEGs_day1_batch1_reannot.csv")

######################

DEA_b1_day02_res_joined <- DEA_b1_day02_res %>% left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b1_day02_res_joined_na <-DEA_b1_day02_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b1_day02_res_joined_na$X,
                        mart = mart)

reannot_b1_d2 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b1_day02_res_reannot <- DEA_b1_day02_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b1_d2$external_gene_name)

write_csv2(DEA_b1_day02_res_reannot, "DEGs_day2_batch1_reannot.csv")

######################

DEA_b1_day03_res_joined <- DEA_b1_day03_res %>% left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b1_day03_res_joined_na <-DEA_b1_day03_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b1_day03_res_joined_na$X,
                        mart = mart)

reannot_b1_d3 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b1_day03_res_reannot <- DEA_b1_day03_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b1_d3$external_gene_name)

write_csv2(DEA_b1_day03_res_reannot, "DEGs_day3_batch1_reannot.csv")

######################

DEA_b1_day07_res_joined <- DEA_b1_day07_res %>%
  left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b1_day07_res_joined_na <-DEA_b1_day07_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b1_day07_res_joined_na$X,
                        mart = mart)

reannot_b1_d7 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b1_day07_res_reannot <- DEA_b1_day07_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b1_d7$external_gene_name)

write_csv2(DEA_b1_day07_res_reannot, "DEGs_day7_batch1_reannot.csv")


######################

DEA_b2_day01_res_joined <- DEA_b2_day01_res %>%
  left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b2_day01_res_joined_na <-DEA_b2_day01_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b2_day01_res_joined_na$X,
                        mart = mart)

reannot_b2_d1 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b2_day01_res_reannot <- DEA_b2_day01_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b2_d1$external_gene_name)

write_csv2(DEA_b2_day01_res_reannot, "DEGs_day1_batch2_reannot.csv")

######################

DEA_b2_day02_res_joined <- DEA_b2_day02_res %>%
  left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b2_day02_res_joined_na <-DEA_b2_day02_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b2_day02_res_joined_na$X,
                        mart = mart)

reannot_b2_d2 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b2_day02_res_reannot <- DEA_b2_day02_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b2_d2$external_gene_name)

write_csv2(DEA_b2_day02_res_reannot, "DEGs_day2_batch2_reannot.csv")

######################

DEA_b2_day03_res_joined <- DEA_b2_day03_res %>%
  left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b2_day03_res_joined_na <-DEA_b2_day03_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b2_day03_res_joined_na$X,
                        mart = mart)

reannot_b2_d3 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b2_day03_res_reannot <- DEA_b2_day03_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b2_d3$external_gene_name)

write_csv2(DEA_b2_day03_res_reannot, "DEGs_day3_batch2_reannot.csv")

######################

DEA_b2_day07_res_joined <- DEA_b2_day07_res %>%
  left_join(txtogenesv37, by=c("X"="gene_name"))

DEA_b2_day07_res_joined_na <-DEA_b2_day07_res_joined %>%
  filter(is.na(transcript_type)) %>% select(X)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_b2_day07_res_joined_na$X,
                        mart = mart)

reannot_b2_d7 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_b2_day07_res_reannot <- DEA_b2_day07_res %>% 
  filter(X %in% human_lncRNAs_v38$gene_name |
           X %in% reannot_b2_d7$external_gene_name)

write_csv2(DEA_b2_day07_res_reannot, "DEGs_day7_batch2_reannot.csv",
           row.names = FALSE)


######################

Eboplus_DEGsD1 <- read_delim("Eboplus_DEGsD1.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


DEA_day1_eboplus_res_joined <- Eboplus_DEGsD1 %>%
  left_join(txtogenesv37, by=c("genes"="gene_name"))

DEA_day1_eboplus_res_joined_na <-DEA_day1_eboplus_res_joined %>%
  filter(is.na(transcript_type)) %>% select(genes)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_day1_eboplus_res_joined_na$genes,
                        mart = mart)

reannot_eboplus_d1 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_day1_eboplus_res_reannot <- Eboplus_DEGsD1 %>% 
  filter(genes %in% human_lncRNAs_v38$gene_name |
           genes %in% reannot_eboplus_d1$external_gene_name)

write.csv(DEA_day1_eboplus_res_reannot, "DEGs_day1_eboplus_reannot_v2.csv",
           row.names = FALSE)

######################

Eboplus_DEGsD2 <- read_delim("Eboplus_DEGsD2.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

DEA_day2_eboplus_res_joined <- Eboplus_DEGsD2 %>%
  left_join(txtogenesv37, by=c("genes"="gene_name"))

DEA_day2_eboplus_res_joined_na <-DEA_day2_eboplus_res_joined %>%
  filter(is.na(transcript_type)) %>% select(genes)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEA_day2_eboplus_res_joined_na$genes,
                        mart = mart)

reannot_eboplus_d2 <- genes %>% filter(gene_biotype == "lncRNA")

DEA_day2_eboplus_res_reannot <- Eboplus_DEGsD2 %>% 
  filter(genes %in% human_lncRNAs_v38$gene_name |
           genes %in% reannot_eboplus_d2$external_gene_name)

write.csv(DEA_day2_eboplus_res_reannot, "DEGs_day2_eboplus_reannot_v2.csv",
           row.names = FALSE)

######################

Eboplus_DEGsD3 <- read_delim("Eboplus_DEGsD3.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


DEGs_day3eboplus_res_joined <- Eboplus_DEGsD3 %>%
  left_join(txtogenesv37, by=c("genes"="gene_name"))

DEGs_day3eboplus_res_joined_na <-DEGs_day3eboplus_res_joined %>%
  filter(is.na(transcript_type)) %>% select(genes)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEGs_day3eboplus_res_joined_na$genes,
                        mart = mart)

reannot_eboplus_d3 <- genes %>% filter(gene_biotype == "lncRNA")

DEGs_day3eboplus_res_reannot <- Eboplus_DEGsD3 %>% 
  filter(genes %in% human_lncRNAs_v38$gene_name |
           genes %in% reannot_eboplus_d3$external_gene_name)

write.csv(DEGs_day3eboplus_res_reannot, "DEGs_day3_eboplus_reannot_v2.csv",
           row.names = FALSE)

######################

Eboplus_DEGsD7 <- read_delim("Eboplus_DEGsD7.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

DEGs_day7eboplus_joined <- Eboplus_DEGsD7 %>%
  left_join(txtogenesv37, by=c("genes"="gene_name"))

DEGs_day7eboplus_joined_na <-DEGs_day7eboplus_joined %>%
  filter(is.na(transcript_type)) %>% select(genes)


genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype"),
                        filters = c("external_synonym"),
                        values = DEGs_day7eboplus_joined_na$genes,
                        mart = mart)

reannot_eboplus_d7 <- genes %>% filter(gene_biotype == "lncRNA")

DEGs_day7eboplus_res_reannot <- Eboplus_DEGsD7 %>% 
  filter(genes %in% human_lncRNAs_v38$gene_name |
           genes %in% reannot_eboplus_d7$external_gene_name)

write.csv(DEGs_day7eboplus_res_reannot, "DEGs_day7_eboplus_reannot_v2.csv", 
           row.names = FALSE)

#####################################
################################################
#VEBCON

vebcon_d1<- rownames_to_column(vebcon_d1)
vebcon_d3<- rownames_to_column(vebcon_d3)
vebcon_d7<- rownames_to_column(vebcon_d7)
#######################

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", 
                host = 'https://www.ensembl.org')



#################################

att_df <-listAttributes(mart)

genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype", "ensembl_gene_id"),
                        filters = c("ensembl_gene_id"),
                        values = vebcon_d1$rowname,
                        mart = mart)

vebcon_d1_all <-vebcon_d1 %>% left_join(genes, by = c("rowname"="ensembl_gene_id"))

write_csv2(vebcon_d1_all, "Vebcon_d1_all.csv")

vebcon_d1_lncRNA <-vebcon_d1_all %>% filter(gene_biotype == "lncRNA")

#

txtogenesnodot <-txtogenesv37

txtogenesnodot$gene_id <- gsub("\\..*","",txtogenesnodot$gene_id)

#

vebcon_d1_dic <-vebcon_d1_lncRNA %>% left_join(txtogenesnodot, by = c("rowname"="gene_id"))

write.csv(DEGs_day7eboplus_res_reannot, "DEGs_day7_eboplus_reannot_v2.csv", 
          row.names = FALSE)

####

genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype", "ensembl_gene_id"),
                        filters = c("ensembl_gene_id"),
                        values = vebcon_d3$rowname,
                        mart = mart)

vebcon_d3_all <-vebcon_d3 %>% left_join(genes, by = c("rowname"="ensembl_gene_id"))

vebcon_d3_lncRNA <-vebcon_d3_all %>% filter(gene_biotype == "lncRNA")

vebcon_d3_dic <-vebcon_d3_lncRNA %>% left_join(txtogenesnodot, by = c("rowname"="gene_id"))

write_csv2(vebcon_d3_all, "Vebcon_d3_all.csv")

####
#########################

genes <- biomaRt::getBM(attributes = c("external_synonym",
                                       "external_gene_name",
                                       "gene_biotype", "ensembl_gene_id"),
                        filters = c("ensembl_gene_id"),
                        values = vebcon_d7$rowname,
                        mart = mart)

vebcon_d7_all <-vebcon_d7 %>% left_join(genes, by = c("rowname"="ensembl_gene_id"))

vebcon_d7_lncRNA <-vebcon_d7_all %>% filter(gene_biotype == "lncRNA")

vebcon_d7_dic <-vebcon_d7_lncRNA %>% left_join(txtogenesnodot, by = c("rowname"="gene_id"))


write_csv2(vebcon_d7_all, "Vebcon_d7_all.csv")






