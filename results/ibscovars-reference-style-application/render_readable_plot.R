args<-commandArgs(TRUE);root<-args[1];setwd(file.path(root,'full'))
suppressPackageStartupMessages(library(ggplot2))
plot_df<-read.csv('forest_plot_data.csv');plot_df$label<-factor(plot_df$label,levels=c(as.character(1:40),'Pooled (NS)','Pooled (Sh)','Case deletion'));plot_title_method<-'propensity-score imputation'
expressions<-parse(file.path(root,'original_script.R'))
for(e in expressions)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('p'))){eval(e,.GlobalEnv);break}
stopifnot(nrow(plot_df)==172L)
ggsave('comparison_graph.png',p,width=10,height=14,dpi=200,bg='white')
ggsave('comparison_graph.pdf',p,width=10,height=14,bg='white')
