library("groHMM")
library("GenomicFeatures")
library("org.Hs.eg.db")
library("edgeR")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(Rsamtools)
library(dplyr)
library(data.table)
library(ggplot2)
library(scales)

.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

options(stringsAsFactors = F)

#Load PRJEB22484 reads metadata
#Paper source for data:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6280751/
#Methods paper:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522910/
files <- list.files(path = "GROseq_lymphoblasts/data",
                    pattern = "DOCK8-AS1_DOCK8.bam",
                    full.names = T)

tps <- gsub(list.files(path = "GROseq_lymphoblasts/data",
                  pattern = "DOCK8-AS1_DOCK8.bam",
                  full.names = F),pattern = "_.*",replacement = "")

tps_order <- tps[c(1,7,4,6,9,10,12,2,3,5,8,11)]

depths <- fread("GROseq_lymphoblasts/data/counts_per_lib.txt")
lib_depths <- data.frame(lib=tps,
                         depth=depths)

DOCK8_positions <- c(214854,465259)
DOCK8_AS1_positions <- c(212824,215893)
ENSG00000235880_positions <- c(267966,273002)
gene_lenghts <- c(DOCK8_positions[2]-DOCK8_positions[1],
                  DOCK8_AS1_positions[2]-DOCK8_AS1_positions[1],
                  ENSG00000235880_positions[2]-ENSG00000235880_positions[1])

gene_lenghts <- gene_lenghts/1000

for(i in 1:length(files)){
  fil <- files[i]
  bam <- Rsamtools::scanBam(fil)
  bam_field <- names(bam[[1]])
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  
  pos_freq <- as.data.frame(table(bam_df$pos))
  pos_freq$Var1 <- as.numeric(as.character(pos_freq$Var1))
  
  pos_all <- pos_strand_df <- data.frame(pos=bam_df$pos,
                              strand=bam_df$flag)
  
  pos_all <- pos_all %>%
    dplyr::mutate(timepoint=tps[i],
           gene=ifelse(pos>=DOCK8_AS1_positions[1] & pos<=DOCK8_AS1_positions[2] & strand==16,
                       yes = "DOCK8-AS1",
                       no = ifelse(pos>=DOCK8_positions[1] & pos<=DOCK8_positions[2] & strand==0,
                                   yes = "DOCK8",
                                   no = ifelse(pos>=ENSG00000235880_positions[1] & pos<=ENSG00000235880_positions[2] & strand==16,
                                               yes="ENSG00000235880",
                                               no = ""))))
  
  depth <- depths$V1[i]
  pm_factor <- depth/1000000
  
  pos_freq_plot <- pos_freq %>%
    dplyr::left_join(pos_strand_df,by = c("Var1"="pos")) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::mutate(Freq=ifelse(strand==16,-1*Freq,Freq),
           timepoint=tps[i],
           gene=ifelse(Var1>=DOCK8_AS1_positions[1] & Var1<=DOCK8_AS1_positions[2] & strand==16,
                       yes = "DOCK8-AS1",
                       no = ifelse(Var1>=DOCK8_positions[1] & Var1<=DOCK8_positions[2] & strand==0,
                                   yes = "DOCK8",
                                   no = ifelse(Var1>=ENSG00000235880_positions[1] & Var1<=ENSG00000235880_positions[2] & strand==16,
                                               yes="ENSG00000235880",
                                               no = ""))))
  
  counts <- pos_freq_plot %>%
    dplyr::group_by(gene,timepoint) %>%
    dplyr::summarise(counts=sum(abs(Freq))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(gene!="") %>%
    dplyr::mutate(length=gene_lenghts,lib_factor=pm_factor,
           RPKM=(counts/lib_factor)/length)
  
  
  test <- subset(bam_df, rname == 'chr9')
  chr9_neg <- bam_df[bam_df$rname == 'chr9' &
                      apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos']
  
  chr9_pos <- bam_df[bam_df$rname == 'chr9' &
                      apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos']
  
  chr9_neg_density <- density(chr9_neg)
  chr9_pos_density <- density(chr9_pos)
  chr9_neg_density$y <- chr9_neg_density$y * -1
  
  tp <- files[i]
  df_plot <- data.frame(density_pos=chr9_pos_density$y,
                      position_pos=chr9_pos_density$x,
                      density_neg=chr9_neg_density$y,
                      position_neg=chr9_neg_density$x,
                      timepoint=tps[i],
                      genes="DOCK8-AS1_DOCK8")
  
  if(i==1){
    df_plot_all <- df_plot
    counts_all <- counts
    pos_freq_plot_all <- pos_freq_plot
    pos_df_all <- pos_all
  }else{
    df_plot_all <- rbind(df_plot_all,df_plot)
    counts_all <- rbind(counts_all,counts)
    pos_freq_plot_all <- rbind(pos_freq_plot_all,pos_freq_plot)
    pos_df_all <- rbind(pos_df_all,pos_all)
  }
}

#Alternative 01 - density plot with 512 positions interval (bins)
p <- df_plot_all %>%
  #filter(timepoint=="GM0h") %>%
  mutate(timepoint=factor(timepoint,levels = tps_order)) %>%
  ggplot()+
  geom_area(aes(x = position_pos,y=density_pos),fill="red")+
  geom_area(aes(x = position_neg,y=density_neg),fill="blue")+
  theme_bw()+
  scale_x_continuous(labels = comma)+
  facet_wrap(facets = ~timepoint,ncol = 2)
pdf(file = "GROseq_lymphoblasts/figures/DOCK8-AS1_DOCK8_read_density_timepoints.pdf",
    width = 18.92,height = 2.03)
print(p)
dev.off()

#Alternative 02 -histogram of reads with 100 positions interval (bins)
p <- pos_df_all %>%
  mutate(timepoint=factor(timepoint,levels = tps_order)) %>%
  ggplot(aes(group=strand,fill=as.character(strand)))+
  geom_histogram(aes(x = pos),bins = 100)+
  theme_bw()+
  scale_x_continuous(labels = comma)+
  facet_wrap(facets = ~timepoint)
p

#Alternative 03 - dotplot with one point in each position with the number of
#counts in each position
library(ggrastr)
p <- pos_freq_plot_all %>%
  dplyr::mutate(timepoint=factor(timepoint,tps_order),
         gene=ifelse(gene=="","not_annotated",gene)) %>%
  ggplot(aes(x=Var1,y=Freq,color=gene))+
  geom_point(size=1)+
  theme_bw()+
  scale_x_continuous(labels = comma,limits = c(200000,500000),
                     breaks = c(200000,300000,400000,500000))+
  facet_wrap(facets = ~timepoint,ncol = 2)
pdf(file = "GROseq_lymphoblasts/figures/DOCK8-AS1_DOCK8_read_dotplot_timepoints.pdf",
    width = 11.4,height = 9.5)
print(p)
dev.off()

p <- pos_freq_plot_all %>%
  dplyr:::filter(timepoint=="GM24h") %>%
  dplyr:::mutate(gene=ifelse(gene=="","not_annotated",gene)) %>%
  ggplot(aes(x=Var1,y=Freq,color=gene))+
  #geom_point_rast(size=1,raster.dpi = 1200)+
  geom_point(size=1)+
  theme_bw()+
  scale_x_continuous(labels = comma,limits = c(200000,500000),
                     breaks = c(200000,300000,400000,500000))
pdf(file = "GROseq_lymphoblasts/figures/DOCK8-AS1_DOCK8_read_dotplot_24h.pdf",
    width = 9.33,height = 2.46)
print(p)
dev.off()

#Plot genomic region with gene structures
library(ggbio)
library(Homo.sapiens)
class(Homo.sapiens)
data(genesymbol, package = "biovizBase")
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75

p.ideo <- Ideogram(genome = "hg19")+
  xlim(GRanges("chr9", IRanges(1,138394717)))
p.ideo

gr <- GRanges(seqnames = 9, IRanges(212824, 465259), strand = "*")

p <- autoplot(ensdb, GRangesFilter(gr), names.expr = "gene_name")+
  theme_bw()+
  scale_x_continuous(labels = comma)
pdf(file = "GROseq_lymphoblasts/figures/DOCK8-AS1_DOCK8_genomic_range.pdf",
    width = 4.25*2,height = 3.52*2)
print(p)
dev.off()


#Number of counts of each gene per timepoint
p <- counts_all %>%
  mutate(timepoint=factor(timepoint,levels = tps_order)) %>%
  ggplot(aes(x=timepoint,y = RPKM,color=gene))+
  geom_point()+
  geom_line(aes(x=as.numeric(timepoint)))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
pdf(file = "GROseq_lymphoblasts/figures/DOCK8-AS1_DOCK8_RPKM_timepoints.pdf",
    width = 6.88,height = 2.55)
print(p)
dev.off()

#Correlation between DOCK8-AS1 and DOCK8 using all timepoints
df_plot <- reshape2::dcast(data = counts_all,formula = timepoint~gene,fun.aggregate = sum,value.var = "RPKM")
cor.test(x = df_plot$DOCK8,y = df_plot$`DOCK8-AS1`)
p <- df_plot %>%
  ggplot(aes(x = `DOCK8-AS1`,y = DOCK8))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()
pdf(file = "GROseq_lymphoblasts/figures/DOCK8-AS1_DOCK8_correlation.pdf",
    width = 5.18,height = 2.55)
print(p)
dev.off()


x <- fread("http://ds.biogps.org/dataset/csv/BDS_00001/gene/157983/")
