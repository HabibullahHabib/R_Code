library(grid)
table_files <- list.files(pattern='^table[23]_.*rds$',recursive=TRUE)
pdf('main_tables.pdf',width=16,height=10,onefile=TRUE)
for(f in table_files) {
 d <- readRDS(f)
 # Split adjustment blocks so ten-treatment tables remain legible.
 groups <- split(seq.int(2,ncol(d)),ceiling(seq_along(seq.int(2,ncol(d)))/11))
 for(cols in groups) for(start in seq(1,nrow(d),by=20)) {
  rows <- start:min(nrow(d),start+19); x <- d[rows,c(1,cols),drop=FALSE]
  x[] <- lapply(x,function(v) {if(is.numeric(v)) v <- format(round(v,2),trim=TRUE); v[is.na(v)] <- ''; as.character(v)})
  grid.newpage(); grid.text(sub('.rds$','',f),.5,.96,gp=gpar(fontsize=10))
  widths <- c(.24,rep(.72/(ncol(x)-1),ncol(x)-1)); xpos <- cumsum(widths)-widths/2
  for(j in seq_len(ncol(x))) {
   grid.text(names(x)[j],xpos[j],.88,gp=gpar(fontsize=7,fontface='bold'))
   for(i in seq_len(nrow(x))) grid.text(x[[j]][i],xpos[j],.83-(i-1)*.038,gp=gpar(fontsize=if(j==1)5.5 else 7))
  }
 }
}
dev.off()
comparison <- sample_size_comparison_df
pdf('simulation_plots.pdf',width=12,height=8,onefile=TRUE)
for(v in unique(comparison$variance_method)) for(s in unique(comparison$scenario)) for(mr_value in missing_rates) for(metric in c('Joint_many_coverage','Mean_many_CI_length')) {
 z <- comparison[comparison$variance_method==v & comparison$scenario==s & comparison$mr==mr_value,]
 if(!nrow(z) || !any(is.finite(z[[metric]]))) next
 plot(range(z$n_per_arm),range(c(z[[metric]],if(metric=='Joint_many_coverage')95),finite=TRUE),type='n',xlab='Sample size per arm',ylab=metric,main=paste(v,s,paste0('Missingness ',100*mr_value,'%'),sep='\n'),cex.main=.8)
 if(metric=='Joint_many_coverage') abline(h=95,lty=2,col='grey40')
 methods <- unique(z$method)
 for(i in seq_along(methods)) {a <- z[z$method==methods[i],];a<-a[order(a$n_per_arm),];lines(a$n_per_arm,a[[metric]],type='b',col=i,lwd=2)}
 legend('bottomright',legend=methods,col=seq_along(methods),lty=1,cex=.6,bty='n')
}
dev.off()
stopifnot(nzchar(Sys.which('pdfunite')))
stopifnot(system2(Sys.which('pdfunite'),c('main_tables.pdf','simulation_plots.pdf','simulation_report.pdf'))==0L)
files <- list.files('checkpoints',pattern='completion.csv',recursive=TRUE,full.names=TRUE)
completion <- dplyr::bind_rows(lapply(files,read.csv))
stopifnot(nrow(completion)==length(sample_sizes)*length(scenarios)*length(missing_rates),sum(completion$attempted)==if(.h$mode=='full')2000L else length(sample_sizes)*10L,length(table_files)==4*length(sample_sizes))
write.csv(completion,'completion_summary.csv',row.names=FALSE)
for(name in c('warnings.csv','errors.csv')) if(!file.exists(name)) writeLines(if(name=='warnings.csv')'condition,job,warning' else 'condition,job,error',name)
writeLines(c(paste('Treatments:',n_treatments),paste('Replicates per setting:',if(.h$mode=='full')'28 or 29; see allocation_plan.csv' else '1'),paste('Total jobs:',sum(completion$attempted)),paste('Successful:',sum(completion$successful)),paste('Failed:',sum(completion$failed)),'Original uploaded source retained unchanged. Runner overrides n_sim only; all seeds and 40 imputations retained.','Raw CSVs are losslessly gzip-compressed. Supplementary plots show cross-sample-size coverage and mean CI length.'),'run_summary.txt')

if(.h$mode=='full')file.copy(file.path(.h$root,'allocation_plan.csv'),'allocation_plan.csv',overwrite=TRUE)
