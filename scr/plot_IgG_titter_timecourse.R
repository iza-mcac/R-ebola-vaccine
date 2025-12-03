geneva_plot <- geneva_antibody %>%
  dplyr::select(volunteer_id,day,IgGTitter) %>%
  dplyr::mutate(dataset="Ebovac")

usa_plot <- USA_antibody_RNAseq_volunteers %>%
  dplyr::filter(LBTEST=="Ebola Zaire Virus GP IgG Antibody Titer") %>%
  dplyr::select(Volunteer.ID,day,LBSTRESN_2) %>%
  dplyr::mutate(dataset="Eboplus")

colnames(usa_plot) <- colnames(geneva_plot)

df_plot <- rbind(geneva_plot,usa_plot)

p <- df_plot %>%
  mutate(day=ifelse(day==168,yes = 180,
                    no=ifelse(day==360,yes = 365,
                              no=day)),
         dataset=factor(dataset,levels = c("Ebovac","Eboplus"))) %>%
  ggplot(aes(x=day,y = log2(IgGTitter+1)))+
  geom_jitter()+
  geom_smooth()+
  theme_bw()+
  scale_x_continuous(breaks = c(0,7,14,28,56,84,180,365,730))+
  facet_wrap(facets = ~dataset,ncol = 1)
pdf(file = "Corr_analysis/figures/IgG_titter_timecourse.pdf",
    width = 5.08,height = 5.08)
print(p)
dev.off()


