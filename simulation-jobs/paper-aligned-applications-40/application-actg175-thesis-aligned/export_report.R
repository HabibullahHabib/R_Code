grDevices::graphics.off()
comparison<-function(method_label){
 z<-primary_long|>dplyr::filter(method==method_label,variance_method%in%c('Equal-variance','Shrinkage'))|>dplyr::select(method,contrast,estimate,variance_method,lower,upper,width,significant)
 z|>tidyr::pivot_wider(names_from=variance_method,values_from=c(lower,upper,width,significant))
}
thesis_primary<-comparison('Regression MI');reference_method<-comparison('Propensity MI')
write.csv(thesis_primary,'THESIS_PRIMARY_regression_equal_vs_shrinkage.csv',row.names=FALSE)
write.csv(reference_method,'REFERENCE_METHOD_propensity_equal_vs_shrinkage.csv',row.names=FALSE)
write.csv(primary_long,'ALL_METHODS_full_precision_intervals.csv',row.names=FALSE)
for(info in list(list(data=thesis_primary,prefix='THESIS_PRIMARY',method='Regression MI'),list(data=reference_method,prefix='REFERENCE_METHOD',method='Propensity MI'))){
 gt::gtsave(gt::tab_header(gt::gt(info$data),title=paste(info$method,'equal-variance vs shrinkage'),subtitle='Adaptation to the supplied dataset; 40 imputations requested by user (documents use 20).'),paste0(info$prefix,'_comparison.html'))
 z<-primary_long|>dplyr::filter(method==info$method,variance_method%in%c('Equal-variance','Shrinkage'))
 p<-ggplot2::ggplot(z,ggplot2::aes(x=contrast,y=estimate,ymin=lower,ymax=upper,shape=variance_method))+ggplot2::geom_hline(yintercept=0,linetype=2)+ggplot2::geom_pointrange(position=ggplot2::position_dodge(width=.5))+ggplot2::coord_flip()+ggplot2::theme_minimal()+ggplot2::labs(title=paste(info$method,'equal-variance vs shrinkage'),subtitle='40 imputations; real-data adaptation, not replication of the original paper dataset',y='Treatment minus control (simultaneous 95% CI)',x=NULL)
 .h$saveplot(paste0(info$prefix,'_forest.png'),p,width=11,height=7,dpi=300)
}
writeLines(c('Primary procedure follows thesis section 4.1: regression MI; common-variance versus unequal-variance shrinkage intervals.','User requested 40 imputations; both uploaded documents use 20 in their real-data applications.','Reference article section 6 uses propensity MI; that method is reported separately as a reference-method adaptation.','Original corrected-script non-shrinkage means unshrunk unequal variance and is retained only as a supplementary comparator.','Equal-variance intervals use the per-imputation degrees-of-freedom-weighted pooled arm variance, Rubin pooling and the same multivariate-t trace-based df construction.','ACTG175 and IBScovars are different datasets from the original 100-patient dose-response application. Covariates, natural missingness/MAR setup and unequal arm sizes follow the corrected user scripts; these are adaptations, not exact paper replications.','No repeated-deletion simulations were added. The 10000-imputation runs are preserved separately and superseded by this user-requested 40-imputation analysis.'),'METHOD_ALIGNMENT.txt')
source(file.path(.h$root,'render_report.R'))
writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
completion<-data.frame(analysis=list.dirs('checkpoints',recursive=FALSE),imputations=40L)
stopifnot(nrow(completion)==if(grepl('actg175',.h$root))3L else 12L)
write.csv(completion,'completion_summary.csv',row.names=FALSE)
writeLines('Completed real-data adaptation with 40 imputations per method/analysis.','run_summary.txt')
