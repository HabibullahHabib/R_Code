args<-commandArgs(TRUE);root<-args[1];setwd(file.path(root,'full'))
suppressPackageStartupMessages(library(ggplot2))
d<-read.csv('forest_plot_data.csv');d$label<-factor(d$label,levels=c(as.character(1:40),'Pooled (NS)','Pooled (Sh)','Case deletion'))
d$kind<-ifelse(d$label=='Pooled (NS)','Proposed without shrinkage',ifelse(d$label=='Pooled (Sh)','Proposed with shrinkage',ifelse(d$label=='Case deletion','Case deletion','Individual imputation')))
kinds<-c('Individual imputation','Proposed without shrinkage','Proposed with shrinkage','Case deletion');d$kind<-factor(d$kind,levels=kinds)
stopifnot(nrow(d)==172L,sum(d$kind=='Individual imputation')==160L,sum(d$kind=='Proposed without shrinkage')==4L,sum(d$kind=='Proposed with shrinkage')==4L)
input<-read.csv('analysis_input_incomplete.csv');actual<-100*mean(is.na(input$y))
p<-ggplot(d,aes(x=estimate,y=label,colour=kind,shape=kind))+
 geom_vline(xintercept=0,linetype=2,linewidth=.4,colour='grey35')+
 geom_segment(aes(x=lower,xend=upper,yend=label,linewidth=kind),show.legend=FALSE)+
 geom_point(aes(size=kind))+
 scale_colour_manual(name='Result',values=c('Individual imputation'='#A0A0A0','Proposed without shrinkage'='#0072B2','Proposed with shrinkage'='#D55E00','Case deletion'='#333333'))+
 scale_shape_manual(name='Result',values=c('Individual imputation'=16,'Proposed without shrinkage'=15,'Proposed with shrinkage'=18,'Case deletion'=17))+
 scale_linewidth_manual(values=c('Individual imputation'=.25,'Proposed without shrinkage'=.85,'Proposed with shrinkage'=.85,'Case deletion'=.45),guide='none')+
 scale_size_manual(values=c('Individual imputation'=.8,'Proposed without shrinkage'=2.6,'Proposed with shrinkage'=3,'Case deletion'=2),guide='none')+
 scale_y_discrete(labels=function(x){x[x=='Pooled (NS)']<-'Proposed without shrinkage';x[x=='Pooled (Sh)']<-'Proposed with shrinkage';x})+
 facet_wrap(~contrast,scales='free_x',ncol=4)+
 labs(title='IBScovars: regression imputation - individual imputations and proposed procedures',subtitle=sprintf('40 regression imputations; 30%% target MAR (%.2f%% realized missingness)',actual),x='Treatment minus control: point estimate and simultaneous 95% confidence interval',y='Imputed dataset / proposed procedure',caption='Highlighted rows show the final proposed procedures using all 40 regression imputations. Grey numbered rows are individual imputations.\nRubin combination is used for both procedures; this is not an equal-variance assumption. Highlighting does not indicate superiority.')+
 theme_minimal(base_size=10)+theme(panel.grid.minor=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',legend.text=element_text(size=9),axis.text.y=element_text(size=8),plot.title=element_text(face='bold',size=14),plot.caption=element_text(size=10,hjust=0))+
 guides(colour=guide_legend(override.aes=list(size=c(1.5,3,3.5,2.5))))
ggsave('comparison_graph.png',p,width=16,height=11,dpi=180,bg='white');ggsave('comparison_graph.pdf',p,width=16,height=11,bg='white')
write.csv(d,'highlighted_forest_plot_data.csv',row.names=FALSE)
writeLines(c('Reference: Jiang, Zhou and Liang (2023), Figure 1, p.13; adjacent application paragraph describes pooled point/interval estimates across imputed datasets.','The reference row Pooled is the Rubin-combined MI result, not a chosen individual imputation and not an equal-variance covariance assumption.','This figure labels the final procedures directly: blue square = proposed without shrinkage, orange diamond = proposed with shrinkage.','Grey numbered rows = individual regression imputations; black triangle = naive case deletion; dashed vertical line = zero difference.','Our y-axis and legend use proposed with/without shrinkage instead of Pooled. Colour indicates procedure, not significance or superiority.'),'FOREST_PLOT_EXPLANATION.txt')
