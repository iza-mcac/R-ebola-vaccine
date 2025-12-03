library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(metafor)
library(igraph)

#Load lncRNA - genes meta analysis
lnc_gene_meta <- fread("data/Metafor_results/Metafor_genes/RandomEffects_Model_Corr_lncRNAs_Genes_AllCohorts.csv")

lnc_gene_meta <- lnc_gene_meta %>%
  dplyr::select(-1) %>%
  dplyr::filter(!duplicated(.)) %>%
  mutate(sig=ifelse(p.adj<0.05,yes = "sig",no = "not.sig"),
         dir=ifelse(beta>0,yes = "positive",no = "negative"),
         sig_dir=ifelse(sig=="sig" & dir=="positive",
                        yes = "positive",
                        no = ifelse(sig=="sig",
                                    yes = "negative",
                                    no = "not.sig")))

write.csv(lnc_gene_meta,file = "Corr_analysis/lncRNA_RNA_metacorrelation_results.csv",
          row.names = F)

#Volcano plots for the meta-correlation results for each lncRNA
p <- lnc_gene_meta %>%
  ggplot(aes(x = beta,y = -log(p.adj),color=sig_dir))+
  geom_point(size=1)+
  # geom_text_repel(data = lnc_ab_meta %>%
  #                   filter(sig=="sig"),
  #                 aes(label=Gene),
  #                 size=3,force = 1,color="black")+
  theme_bw()+
  geom_hline(yintercept = -log(0.05))+
  scale_color_manual(values = c("blue","grey","red3"))+
  facet_wrap(facets = ~LncRNA)
p

#Plot bar plots for the number of sig correlation for each lncRNA
df_plot <- lnc_gene_meta %>%
  filter(sig=="sig") %>%
  mutate(value=1) %>%
  group_by(LncRNA,dir) %>%
  summarise(correlated_genes=sum(value)) %>%
  ungroup() %>%
  mutate(correlated_genes=ifelse(dir=="negative",
                                 yes = correlated_genes*-1,
                                 no = correlated_genes)) %>%
  arrange(desc(correlated_genes)) %>%
  mutate(LncRNA=factor(LncRNA,levels = rev(unique(LncRNA))))

p <- df_plot %>%
  ggplot(aes(x=LncRNA,y = correlated_genes,fill=dir))+
  geom_col()+
  geom_text(data = df_plot %>%
              filter(dir=="positive"),
            aes(label = abs(correlated_genes)), vjust = -0.5)+
  geom_text(data = df_plot %>%
              filter(dir=="negative"),
            aes(label = abs(correlated_genes)), vjust = 1.5)+
  scale_fill_manual(values = c("blue","red3"))+
  scale_y_continuous(breaks = c(-1000,-500,0,500,1000),
                     limits = c(-1100,1100),
                     labels = c(1000,500,0,500,1000))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p


lncRNAs <- unique(lnc_gene_meta$LncRNA)

reactome_df <- readRDS("data/support_files/all_level_reactome.RDS")
colnames(reactome_df)[1] <- "term"
BTM_df <- clusterProfiler::read.gmt("data/support_files/BTM_for_GSEA_20131008.gmt")

reactome_BTM_df <- rbind(reactome_df,BTM_df)

reactome_BTM_fgsea <- split(reactome_BTM_df$gene,
                            reactome_BTM_df$term)

universe <- unique(c(reactome_df$gene,
                     BTM_df$gene,
                     lnc_gene_meta$gene))

for(i in 1:length(lncRNAs)){
  lnc <- lncRNAs[i]
  df <- lnc_gene_meta %>%
    filter(LncRNA==lnc)
  
  
  #fgsea (Reactome and BTMs)
  ranks <- df$beta
  names(ranks) <- df$gene
  
  fgsea_res <- fgsea::fgseaSimple(pathways = reactome_BTM_fgsea,
                            stats = ranks,nperm = 1000)
  fgsea_res$lncRNA <- lnc
  
  #Overrepresentation analysis (Reactome and BTMs)
  cor_pos <- df %>%
    filter(sig_dir=="positive") %>%
    pull(gene)
  
  ora_pos_res <- fgsea::fora(pathways = reactome_BTM_fgsea,
                             genes = cor_pos,
                             universe = universe)
  
  ora_pos_res$dir <- "positive"
  
  cor_neg <- df %>%
    filter(sig_dir=="negative") %>%
    pull(gene)
  
  ora_neg_res <- fgsea::fora(pathways = reactome_BTM_fgsea,
                             genes = cor_neg,
                             universe = universe)
  
  ora_neg_res$dir <- "negative"
  
  ora_res <- rbind(ora_pos_res,
                   ora_neg_res)
  
  ora_res$lncRNA <- lnc
  
  if(i==1){
    fgsea_all <- fgsea_res
    ora_all <- ora_res
  }else{
    fgsea_all <- rbind(fgsea_all,
                       fgsea_res)
    ora_all <- rbind(ora_all,
                     ora_res)
  }
}

# fgsea_all_2 <- fgsea_all %>%
#   mutate(db=ifelse(pathway %in% unique(reactome_df$term),
#                    yes = "Reactome",no = "BTM"))
# 
# top_sig <- fgsea_all_2 %>%
#   filter(padj<0.05,abs(NES)>3) %>%
#   pull(pathway) %>%
#   unique()
# 
# p <- fgsea_all_2 %>%
#   filter(pathway %in% top_sig) %>%
#   arrange(desc(NES)) %>%
#   mutate(pathway=factor(pathway,levels = rev(unique(pathway))),
#          is.sig=ifelse(padj < 0.05,yes = "yes",no = "no")) %>%
#   ggplot(aes(x=lncRNA,y = pathway,color=NES,size=-log(padj),alpha=is.sig))+
#   geom_point()+
#   theme_bw()+
#   scale_color_gradient2(low = "blue",mid = "white",high = "red3")+
#   scale_alpha_manual(values = c(0.1,1))+
#   theme(axis.text.x = element_text(angle = 90,hjust = 1))+
#   facet_wrap(facets = ~db,scales = "free_y",ncol = 1)
# p

ora_all_2 <- ora_all %>%
  mutate(db=ifelse(pathway %in% unique(reactome_df$term),
                   yes = "Reactome",no = "BTM"))

top_sig <- ora_all_2 %>%
  group_by(dir,db) %>%
  top_n(n = 10,wt = -log(padj)) %>%
  ungroup() %>%
  pull(pathway) %>%
  unique()

all_sig_BTM <- ora_all_2 %>%
  filter(padj<0.05,db=="BTM") %>%
  pull(pathway) %>%
  unique()

all_sig_reactome <- ora_all_2 %>%
  filter(padj<0.05,db=="Reactome") %>%
  pull(pathway) %>%
  unique()

all_sig <- ora_all_2 %>%
  filter(padj<0.05) %>%
  pull(pathway) %>%
  unique()

p <- ora_all_2 %>%
  filter(pathway %in% all_sig_BTM) %>%
  arrange(desc(-log(padj))) %>%
  mutate(pathway=factor(pathway,levels = rev(unique(pathway))),
         is.sig=ifelse(padj < 0.05,yes = "yes",no = "no")) %>%
  filter(is.sig=="yes") %>%
  ggplot(aes(x=lncRNA,y = pathway,
             color=dir,size=-log(padj)))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("blue","red"))+
  scale_alpha_manual(values = c(0.1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  facet_wrap(facets = ~db,scales = "free_y",ncol = 2)
p

#Path-path network
genes_sig <- ora_all_2 %>%
  filter(pathway %in% all_sig_BTM) %>%
  pull(overlapGenes) %>%
  unlist() %>% unique()

reactome_BTM_df$term <- as.character(reactome_BTM_df$term)
path_sig_genes <- reactome_BTM_df %>%
  filter(term %in% all_sig_BTM,
         gene %in% genes_sig)

path_sig_fgsea <- split(path_sig_genes$gene,
                        path_sig_genes$term)


path_sig_enr <- lapply(X = path_sig_fgsea,
                       FUN = fgsea::fora,
                       pathways=path_sig_fgsea,
                       universe=universe)

for(i in 1:length(path_sig_enr)){
  df <- path_sig_enr[[i]]
  df$term <- names(path_sig_enr)[i]
  df_out <- df %>%
    dplyr::mutate(score=-log(padj)) %>%
    dplyr::select(term,pathway,score) %>%
    dplyr::filter(term!=pathway) %>%
    dplyr::rename(from=1,to=2,weight=3)
  
  if(i==1){
    term_term_network <- df_out
  }else{
    term_term_network <- rbind(term_term_network,df_out)
  }
}

term_term_network <- term_term_network %>%
  group_by(grp = paste(pmax(from, to), pmin(from, to), sep = "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp) %>%
  as.data.frame()

#term_term_network$weight[is.infinite(term_term_network$weight)] <- 654.7855

#create igraph object and run louvain algorithm to detect clusters
term_term_ggraph <- term_term_network %>%
  filter(weight > -log(0.05))

g <- graph_from_data_frame(term_term_ggraph,directed = F,
                           vertices = unique(c(term_term_ggraph$from,
                                               term_term_ggraph$to)))
louv <- cluster_louvain(graph = g,weights = term_term_ggraph$weight)
mods <- data.frame(node=louv$names,mod=louv$membership)
degree <- strength(graph = g,vids = V(g),mode = "all",weights = E(g)$weight)
nodes <- mods %>%
  mutate(degree = degree) %>%
  mutate(node=ifelse(str_sub(node,
                             start = nchar(node),
                             end = nchar(node))==" ",
                     no=node,
                     yes=str_sub(node,
                                 start = 1,
                                 end = nchar(node)-1)
  ))

nodes_annot <- data.frame(mod=c(1,2,3,4,5),
                          class=c("Antiviral IFN signaling",
                                  "Monocytes and inflammatory response",
                                  "T cells",
                                  "Antigen presentation",
                                  "Platellets"))

edges <- term_term_network %>%
  rename(Source=1,Target=2)

nodes_out <- nodes %>%
  dplyr::mutate(Label=node) %>%
  dplyr::rename(Id=node) %>%
  dplyr::select(Id,Label,everything())

# write.csv(edges,file = "data/Metafor_results/term_term_edges_BTM.csv",
#           row.names = F)
# 
# write.csv(nodes_out,file = "data/Metafor_results/term_term_nodes_BTM.csv",
#           row.names = F)


#update igraph object with degree, annot and labels
V(g)$degree <- nodes$degree
V(g)$mod <- as.character(nodes$mod)

#remove smaller component (contains only 3 terms)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g2 <- igraph::induced_subgraph(g, vert_ids)

####Figure 3A network####
#plot term-term network with mod annotations
library(ggraph)
p <- ggraph(graph = g2,layout = "auto")+
  geom_edge_link(aes(color=log(weight)))+
  geom_node_point(aes(size=degree,
                      fill=as.character(mod)),
                  color="black",pch=21)+
  #geom_node_text(aes(label = V(g)$name),repel=TRUE)+
  #scale_fill_manual(values = pal)+
  scale_edge_color_continuous(low = "white",high = "grey40")+
  theme_void()
print(p)

#Plot only selected pathways
selected_pathways <- c("type I interferon response (M127)",
                       "enriched in T cells (I) (M7.0)",
                       "enriched in monocytes (II) (M11.0)")

p <- ora_all_2 %>%
  filter(pathway %in% selected_pathways) %>%
  arrange(desc(-log(padj))) %>%
  mutate(pathway=factor(pathway,levels = rev(unique(pathway))),
         lncRNA=factor(lncRNA,levels = rev(c("LEF1-AS1","FAM225A",
                                             "DOCK8-AS1","FOXN3-AS1",
                                             "CYTOR","DLEU2","LINC00189",
                                             "SERPINB9P1","C8orf31"))),
         is.sig=ifelse(padj < 0.05,yes = "yes",no = "no")) %>%
  filter(is.sig=="yes") %>%
  ggplot(aes(x=pathway,y = lncRNA,
             color=dir,size=-log(padj)))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("blue","red3"))+
  scale_alpha_manual(values = c(0.1,1))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  facet_wrap(facets = ~db,scales = "free_y",ncol = 2)
p

#Plot lncRNA-gene networks from top BTM modules
selected_BTM_mod_genes <- BTM_df %>%
  filter(term %in% selected_pathways)

selected_lncRNA_genes_cors <- lnc_gene_meta %>%
  filter(sig=="sig",
         gene %in% unique(selected_BTM_mod_genes$gene))

edges <- selected_lncRNA_genes_cors %>%
  mutate(score=abs(beta*-log(p.adj))) %>%
  select(LncRNA,gene,dir,beta,p.adj,score) %>%
  rename(Source=1,Target=2,weight=6)

nodes <- data.frame(Id=unique(c(edges$Source,
                                edges$Target)),
                    Label="",
                    Class=c(rep("lncRNA",11),
                            rep("gene",171))) %>%
  left_join(selected_BTM_mod_genes,by=c("Id"="gene")) %>%
  mutate(term=ifelse(is.na(term),"lncRNA",as.character(term)))

write.csv(edges,file = "data/lnc_selected_BTM_genes_edges.csv",
          row.names = F)

write.csv(nodes,file = "data/lnc_selected_BTM_genes_nodes.csv",
          row.names = F)

library(ggalluvial)

alluvial_df <- edges %>%
  group_by(Source,dir) %>%
  top_n(n = 5,wt = abs(beta)) %>%
  arrange(desc(abs(beta)),.by_group = T) %>%
  ungroup() %>%
  select(Source,Target,dir) %>%
  left_join(selected_BTM_mod_genes,by=c("Target"="gene")) %>%
  mutate(term=as.character(term),
         term=factor(term,levels = unique(term)),
         Source=factor(Source,levels = unique(Source)),
         Target=factor(Target,levels = unique(Target)),
         dir=factor(dir,levels = c("negative","positive")))

p <- alluvial_df %>%
  filter(Source %in% c("FAM225A","DOCK8-AS1","LEF1-AS1")) %>%
  ggplot(aes(axis1 = Source, axis2 = Target, axis3 = term)) +
  # geom_alluvium(aes(fill = dir),
  #               width = 0, knot.pos = 0, reverse = FALSE) +

  geom_flow(aes(fill=dir),width = 0,reverse=F)+
  geom_stratum(width = 1/8, reverse = F) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            reverse = F,size=2.5) +
  scale_x_continuous(breaks = 1:3)+
  scale_fill_manual(values = c("blue","red3"))+
  theme_void()+
  theme(legend.position = "top")
pdf(file = "Corr_analysis/figures/alluvial_lef_fam_dock.pdf",
    width = 9.41,height = 7.41)
print(p)
dev.off()

