library(grid)
files<-list.files(pattern='\\.csv$',recursive=TRUE,full.names=TRUE)
tables<-files[!grepl('validation_raw|progress|warnings|errors',files)]
pdf('main_tables.pdf',width=16,height=10,onefile=TRUE)
for(f in tables){
 d<-read.csv(f,check.names=FALSE);if(!nrow(d))next
 for(cs in split(seq_len(ncol(d)),ceiling(seq_len(ncol(d))/9)))for(start in seq(1,nrow(d),20)){
  x<-d[start:min(start+19,nrow(d)),cs,drop=FALSE]
  x[]<-lapply(x,function(a){if(is.numeric(a))a<-format(round(a,3),trim=TRUE);a[is.na(a)]<-'';as.character(a)})
  grid.newpage();grid.text(basename(f),.5,.95,gp=gpar(fontsize=11));pos<-seq(.07,.93,length.out=ncol(x))
  for(j in seq_len(ncol(x))){grid.text(names(x)[j],pos[j],.88,gp=gpar(fontsize=6,fontface='bold'));for(i in seq_len(nrow(x)))grid.text(x[[j]][i],pos[j],.82-(i-1)*.037,gp=gpar(fontsize=6))}
 }
}
dev.off()
plots<-list.files(pattern='\\.pdf$',recursive=TRUE,full.names=TRUE);plots<-plots[!grepl('main_tables|simulation_report|simulation_plots',plots)]
stopifnot(length(plots)>0,system2(Sys.which('pdfunite'),c(shQuote(plots),'simulation_plots.pdf'))==0L)
stopifnot(system2(Sys.which('pdfunite'),c('main_tables.pdf','simulation_plots.pdf','simulation_report.pdf'))==0L)
expected<-if(.h$mode=='full')1000 else if(grepl('ibscovars',root))4 else 1
stopifnot(length(unique(validation_raw$rep))==validation_reps)
counts<-validation_raw|>dplyr::distinct(missing_rate,rep)
stopifnot(nrow(counts)==expected)
completion<-counts|>dplyr::count(missing_rate,name='successful_experiments');completion$attempted<-completion$successful_experiments;completion$failed<-0L
write.csv(completion,'completion_summary.csv',row.names=FALSE)
for(name in c('warnings.csv','errors.csv'))if(!file.exists(name))writeLines(if(name=='warnings.csv')'stage,rep,warning' else 'time,error',name)
writeLines(c(paste('Total repeated-deletion experiments:',expected),paste('Imputations:',m_imputations),'All three imputation methods evaluated in each experiment; both variance versions share imputations.','Original primary analysis retained. Original source and seeds unchanged; only validation_reps overridden.','Recovery rates describe coverage of fixed complete-sample benchmark contrasts, not population-parameter coverage.'),'run_summary.txt')

write.csv(if(exists('actg'))actg else ibs,'analysis_input_data.csv',row.names=FALSE)
writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
