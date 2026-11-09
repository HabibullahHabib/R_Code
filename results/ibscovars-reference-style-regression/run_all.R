args<-commandArgs(TRUE);.h<-new.env();.h$root<-normalizePath(args[1]);.h$out<-file.path(.h$root,'full');dir.create(.h$out,recursive=TRUE,showWarnings=FALSE);setwd(.h$out)
.h$log<-function(file,rows){p<-file.path(.h$out,file);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
tryCatch(withCallingHandlers({
 for(e in parse(file.path(.h$root,'executed_script.R')))eval(e,.GlobalEnv)
 grDevices::graphics.off()
 write.csv(pooled_res,'computed_pooled_results_full_precision.csv',row.names=FALSE)
 write.csv(cc_res,'computed_case_deletion_full_precision.csv',row.names=FALSE)
 write.csv(per_imp_df,'individual_imputation_intervals.csv',row.names=FALSE)
 write.csv(plot_df,'forest_plot_data.csv',row.names=FALSE)
 write.csv(dat_incomplete,'analysis_input_incomplete.csv',row.names=FALSE)
 write.csv(dat_full,'analysis_input_complete.csv',row.names=FALSE)
 writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
 saveRDS(list(imputations=comps,pooled=pooled_res,case_deletion=cc_res,per_imputation=per_imp_df,data=dat_incomplete),'analysis_snapshot.rds')
 stopifnot(m_imputations==40L,imputation_method=='regression',length(comps)==40L,nrow(dat_full)==369L,nrow(pooled_res)==8L,nrow(cc_res)==4L,nrow(per_imp_df)==160L,nrow(main_table_df)==12L)
 sh<-pooled_res[pooled_res$method=='With shrinkage',];ns<-pooled_res[pooled_res$method=='Without shrinkage',]
 stopifnot(identical(sh$contrast,ns$contrast),identical(sh$estimate,ns$estimate))
 writeLines(paste('Completed as supplied: 40',imputation_method,'imputations; target 30% MAR; achieved',100*achieved_missing,'percent.'),'COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
