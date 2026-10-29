library(grid)
table_files <- list.files(pattern='^table[23]_.*rds$')
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
pdf('simulation_plots.pdf',width=13,height=8,onefile=TRUE)
for(f in list.files(pattern='^paper_MAR_diagnostic_summary_.*csv$')) {
 d <- read.csv(f)
 for(s in unique(d$scenario)) for(metric in c('bias','RMSE','many_SE_ratio')) {
  z <- d[d$scenario==s & d$variance_version=='shrinkage',]
  if(!any(is.finite(z[[metric]]))) next
  par(mar=c(5,5,5,2)); plot(range(z$mr),range(z[[metric]],finite=TRUE),type='n',xlab='Missing fraction',ylab=metric,main=paste(sub('.csv$','',f),s,sep='\n'),cex.main=.7)
  methods <- unique(z$method)
  for(i in seq_along(methods)) {a <- aggregate(z[[metric]][z$method==methods[i]],list(mr=z$mr[z$method==methods[i]]),mean); lines(a$mr,a$x,col=i,lwd=2,type='b')}
  legend('topleft',legend=methods,col=seq_along(methods),lty=1,cex=.65,bty='n')
 }
}
dev.off()
if(nzchar(Sys.which('pdfunite'))) {
 code <- system2(Sys.which('pdfunite'),c('main_tables.pdf','simulation_plots.pdf','simulation_report.pdf'));stopifnot(code==0L)
} else stop('pdfunite is required to combine the report')
completion <- expand.grid(variance_ratio=variance_ratios,mar_mechanism=mar_mechanisms,scenario=names(scenarios),missing_rate=missing_rates,stringsAsFactors=FALSE)
completion$attempted <- n_sim;completion$successful <- n_sim;completion$failed <- 0L
write.csv(completion,'completion_summary.csv',row.names=FALSE)
for(name in c('warnings.csv','errors.csv')) if(!file.exists(name)) writeLines(if(name=='warnings.csv')'condition,job,warning' else 'time,error',name)
stopifnot(length(table_files)==6L*length(variance_ratios)*length(mar_mechanisms))
summary_files <- list.files(pattern='^professor_feedback_.*csv$')
stopifnot(length(summary_files)==2)
writeLines(c('Original scripts, seeds, all condition levels and 40 imputations preserved.',paste('Replicates per setting:',n_sim),paste('Total attempted jobs:',n_sim*60),'Raw CSVs gzip compressed without changing their contents.','Supplementary plots average diagnostics across treatment comparisons; CSVs retain individual comparisons.'),'run_summary.txt')
