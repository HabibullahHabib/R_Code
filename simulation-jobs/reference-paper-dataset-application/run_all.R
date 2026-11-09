args<-commandArgs(TRUE);.h<-new.env();.h$root<-normalizePath(args[1]);.h$out<-file.path(.h$root,'full');dir.create(.h$out,recursive=TRUE,showWarnings=FALSE);setwd(.h$out)
.h$log<-function(file,rows){p<-file.path(.h$out,file);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
.h$saveplot<-function(filename,plot,...){ggplot2::ggsave(filename,plot,...);ggplot2::ggsave(paste0(tools::file_path_sans_ext(filename),'.pdf'),plot,...)}
.h$transform<-function(e){if(!is.call(e))return(e);if(identical(e[[1]],as.name('ggsave')))e[[1]]<-quote(.h$saveplot);for(i in seq_along(e)[-1])if(is.call(e[[i]]))e[[i]]<-.h$transform(e[[i]]);e}
tryCatch(withCallingHandlers({
 for(e in parse(file.path(.h$root,'original_script.R')))eval(.h$transform(e),.GlobalEnv)
 grDevices::graphics.off()
 write.csv(results_long,'computed_results_full_precision.csv',row.names=FALSE)
 write.csv(dat,'analysis_input_data.csv',row.names=FALSE)
 writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
 saveRDS(list(results=results_long,imputations=comps,mean_alpha=pooled$mean_alpha,alpha_by_imputation=pooled$alpha_by_imputation),'analysis_snapshot.rds')
 stopifnot(nrow(dat)==100L,sum(is.na(dat$y))==32L,length(comps)==40L,nrow(results_long)==8L)
 sh<-results_long[results_long$Procedure=='With shrinkage',];ns<-results_long[results_long$Procedure=='Without shrinkage',]
 stopifnot(identical(sh$Comparison,ns$Comparison),identical(sh$Estimate,ns$Estimate))
 writeLines('Completed application as supplied: 40 treatment-restricted ABB imputations. No repeated-simulation loop added.','COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
