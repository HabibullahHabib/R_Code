ok<-vapply(results,function(z)isTRUE(z$ok),logical(1));success<-results[ok];stopifnot(length(success)>0)
paired<-bind_rows(lapply(success,`[[`,'paired'));intervals<-bind_rows(lapply(success,`[[`,'intervals'))
errors<-bind_rows(lapply(results[!ok],function(z)data.frame(sim=z$sim,error=z$error,missing_seed=z$missing_seed,imputation_seed=z$imputation_seed,qseed=z$qseed)))
if(!ncol(errors))errors<-data.frame(sim=integer(),error=character(),missing_seed=integer(),imputation_seed=integer(),qseed=integer())
warnings<-bind_rows(lapply(results,function(z)if(length(z$warnings))data.frame(sim=z$sim,warning=z$warnings)else NULL));if(!ncol(warnings))warnings<-data.frame(sim=integer(),warning=character())
write.csv(errors,'errors.csv',row.names=FALSE);write.csv(warnings,'warnings.csv',row.names=FALSE)
write.csv(paired,'paired_repetitions.csv',row.names=FALSE);write.csv(intervals,'raw_intervals.csv',row.names=FALSE)
widths<-c('Non-shrinkage'=mean(paired$mean_width_nonshrinkage),'Shrinkage'=mean(paired$mean_width_shrinkage))
summary<-data.frame(procedure=names(widths),successful_repetitions=length(success),mean_interval_width=as.numeric(widths),row.names=NULL)
if(app=='IBScovars'){
 summary$joint_benchmark_recovery_pct<-100*c(mean(paired$joint_recovery_nonshrinkage),mean(paired$joint_recovery_shrinkage))
 summary$mean_achieved_missing_pct<-100*mean(paired$achieved_missing_rate)
}
write.csv(summary,'comparison_summary.csv',row.names=FALSE)
mc<-function(x,label){n<-length(x);avg<-mean(x);se<-if(n>1)sd(x)/sqrt(n)else NA_real_;data.frame(metric=label,estimate=avg,monte_carlo_se=se,mc95_lower=avg-1.96*se,mc95_upper=avg+1.96*se)}
paired_comparison<-bind_rows(mc(paired$mean_width_nonshrinkage-paired$mean_width_shrinkage,'Width reduction: non-shrinkage minus shrinkage'),mc(paired$width_saving_pct,'Mean percentage width reduction'),mc(100*as.numeric(paired$shrinkage_narrower),'Repetitions in which shrinkage is narrower (%)'))
if(app=='IBScovars')paired_comparison<-bind_rows(paired_comparison,mc(100*(as.numeric(paired$joint_recovery_shrinkage)-as.numeric(paired$joint_recovery_nonshrinkage)),'Joint benchmark recovery difference: shrinkage minus non-shrinkage (percentage points)'))
write.csv(paired_comparison,'paired_comparison.csv',row.names=FALSE)
per_contrast<-intervals|>group_by(variance_method,contrast)|>summarise(repetitions=n(),mean_estimate=mean(estimate),mean_width=mean(width),.groups='drop')
if(app=='IBScovars')per_contrast<-left_join(per_contrast,intervals|>group_by(variance_method,contrast)|>summarise(benchmark_recovery_pct=100*mean(benchmark_recovered),.groups='drop'),by=c('variance_method','contrast'))
write.csv(per_contrast,'per_contrast_comparison.csv',row.names=FALSE)
completion<-data.frame(dataset=app,attempted=n_experiments,successful=sum(ok),failed=sum(!ok),imputations_per_repetition=40,imputation_method='Regression',missing_design=missing_design)
write.csv(completion,'completion_summary.csv',row.names=FALSE)
colors<-c('#777777','#2878B5')
draw_graph<-function(){
 par(mfrow=c(1,if(app=='IBScovars')2 else 1),mar=c(6,5,5,2),oma=c(5,0,1,0))
 if(app=='IBScovars'){
  recovery<-summary$joint_benchmark_recovery_pct
  x<-barplot(recovery,names.arg=summary$procedure,col=colors,ylim=c(0,110),ylab='Repetitions recovering every benchmark contrast (%)',main='Did the intervals recover the known benchmark?',cex.main=.85,border=NA)
  text(x,recovery+3,labels=paste0(format(round(recovery,2),nsmall=2),'%'))
 }
 x<-barplot(widths,names.arg=names(widths),col=colors,ylim=c(0,max(widths)*1.22),ylab=if(app=='ACTG175')'Average interval width (CD4 count units)'else'Average interval width (response units)',main='How wide were the intervals?',cex.main=.95,border=NA)
 text(x,widths+max(widths)*.04,labels=format(round(widths,if(app=='ACTG175')2 else 3),nsmall=if(app=='ACTG175')2 else 3))
 gain<-mean(paired$width_saving_pct)
 mtext(paste0('Shrinkage was ',format(round(abs(gain),2),nsmall=2),'% ',if(gain>=0)'narrower'else'wider',' on average across repetitions.'),side=1,line=3,cex=.8)
 mtext(paste0(app,': ',n_experiments,' repetitions; 40 regression imputations each; identical imputations for both procedures.'),outer=TRUE,side=3,line=-.4,cex=.85,font=2)
 if(app=='ACTG175')note<-'Same patients and natural missing values in every repetition; only imputation draws change.\nWidth comparisons assess imputation stability and precision; recovery of unknown outcomes cannot be measured.' else note<-'New 30% expected MAR missingness in each repetition. Benchmark = original complete-data contrasts.\nBenchmark recovery is not population-parameter coverage; higher recovery alone does not establish superiority.'
 mtext(note,outer=TRUE,side=1,line=1,cex=.75)
}
pdf('comparison_graph.pdf',width=if(app=='IBScovars')14 else 10,height=7);draw_graph();dev.off()
png('comparison_graph.png',width=if(app=='IBScovars')1960 else 1400,height=980,res=140);draw_graph();dev.off()
library(grid)
pdf('main_tables.pdf',width=14,height=9)
for(item in list(list(title='Shrinkage versus non-shrinkage',data=summary),list(title='Paired comparison (Monte Carlo uncertainty)',data=paired_comparison),list(title='Results by treatment versus control',data=per_contrast))){
 d<-item$data
 for(start in seq(1,nrow(d),20)){
  z<-d[start:min(nrow(d),start+19),,drop=FALSE];z[]<-lapply(z,function(x){if(is.numeric(x))x<-format(round(x,4),trim=TRUE);x[is.na(x)]<-'';as.character(x)})
  grid.newpage();grid.text(paste(app,item$title,sep=': '),.5,.95,gp=gpar(fontsize=13));pos<-seq(.08,.92,length.out=ncol(z))
  for(j in seq_len(ncol(z))){grid.text(names(z)[j],pos[j],.84,gp=gpar(fontsize=7,fontface='bold'));for(i in seq_len(nrow(z)))grid.text(paste(strwrap(z[[j]][i],width=if(j==1)32 else 23),collapse='\n'),pos[j],.75-(i-1)*.06,gp=gpar(fontsize=7))}
 }
}
dev.off();stopifnot(system2(Sys.which('pdfunite'),c('main_tables.pdf','comparison_graph.pdf','simulation_report.pdf'))==0)
writeLines(c(paste('Dataset:',app),paste('Attempted repetitions:',n_experiments),paste('Successful:',sum(ok)),paste('Failed:',sum(!ok)),'Regression imputation only; 40 imputations per repetition. Only shrinkage and non-shrinkage are compared.',missing_design,'ACTG175 retains its original natural missingness and patients: its repetitions measure variation across repeated MI analyses, not new patient samples or new missingness patterns.','IBScovars repeats empirical covariate-driven MAR: gamma=log(2), same in all arms, calibrated to 30% expected missingness, original group assignments/covariate/outcomes retained before masking.','Both procedures share each repetition\'s data, 40 imputations and point estimates.','Known complete-data benchmark recovery is only available for IBScovars; it is not population-parameter coverage.','No equal-variance, PMM, propensity, or cross-dataset comparisons are performed.'),'run_summary.txt')
writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
