args<-commandArgs(TRUE);.h<-new.env();.h$root<-normalizePath(args[1]);.h$out<-file.path(.h$root,'full');dir.create(.h$out,recursive=TRUE,showWarnings=FALSE);setwd(.h$out)
.h$log<-function(file,rows){p<-file.path(.h$out,file);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
.h$validation<-function(){
 dir.create('checkpoints',showWarnings=FALSE)
 spec<-list(source=unname(tools::md5sum(file.path(.h$root,'executed_script.R'))),n=n_validation,m=m_imputations,methods=all_methods,seed=master_seed,input=digest::digest(dat_full))
 if(file.exists('protocol.rds'))stopifnot(identical(spec,readRDS('protocol.rds')))else saveRDS(spec,'protocol.rds')
 cl<-parallel::makeCluster(4,port=as.integer(Sys.getenv('SIM_CLUSTER_PORT')),outfile='worker_console.log');on.exit(parallel::stopCluster(cl))
 parallel::clusterEvalQ(cl,{library(multcomp);library(mvtnorm);library(dplyr);library(ggplot2)})
 parallel::clusterExport(cl,setdiff(ls(.GlobalEnv,all.names=TRUE),'.h'),envir=.GlobalEnv)
 fun<-function(i){w<-character();v<-withCallingHandlers(validation_one(i),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')});list(value=v,warnings=w,rep=i)};environment(fun)<-.GlobalEnv
 all<-vector('list',n_validation)
 for(start in seq(1L,n_validation,10L)){
  ids<-start:min(start+9L,n_validation);p<-file.path('checkpoints',sprintf('batch_%05d.rds',ceiling(start/10)))
  if(file.exists(p))z<-readRDS(p)else{z<-parallel::parLapplyLB(cl,ids,fun,chunk.size=1L);saveRDS(z,paste0(p,'.tmp'));stopifnot(file.rename(paste0(p,'.tmp'),p))}
  all[ids]<-z;validation_list[ids]<<-lapply(z,`[[`,'value');writeLines(as.character(length(z)),sub('rds$','count',p))
  for(a in z)for(w in a$warnings).h$log('warnings.csv',data.frame(stage='validation',rep=a$rep,warning=w))
  saved<-sum(vapply(list.files('checkpoints',pattern='count$',full.names=TRUE),function(p)as.integer(readLines(p)),integer(1)))
  current<-all[seq_len(max(ids))];skipped<-sum(vapply(current,function(z)is.null(z$value),logical(1)))
  write.csv(data.frame(saved_repetitions=saved,total_repetitions=n_validation,processed_repetitions=max(ids),skipped_repetitions=skipped,updated=as.character(Sys.time())),'progress.csv.tmp',row.names=FALSE);file.rename('progress.csv.tmp','progress.csv');cat('Validation checkpoints:',max(ids),'/',n_validation,'\n')
 }
 direct<-validation_one(1L);stopifnot(identical(direct,all[[1]]$value));writeLines('Original serial validation repetition 1 matches its checkpointed parallel result exactly.','FIDELITY_PASSED.txt')
 omitted<-data.frame(rep=which(vapply(all,function(z)is.null(z$value),logical(1))));write.csv(omitted,'omitted_validation_repetitions.csv',row.names=FALSE)
}
.h$transform<-function(e){
 if(!is.call(e))return(e)
 if(identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('bb')))return(quote(.h$validation()))
 for(i in seq_along(e)[-1])if(is.call(e[[i]]))e[[i]]<-.h$transform(e[[i]])
 e
}
tryCatch(withCallingHandlers({
 for(e in parse(file.path(.h$root,'executed_script.R')))eval(.h$transform(e),.GlobalEnv)
 grDevices::graphics.off();stopifnot(identical(all_methods,'REG'),m_imputations==40L)
 write.csv(forest_df,'forest_plot_data.csv',row.names=FALSE);write.csv(reg_pooled,'computed_pooled_results_full_precision.csv',row.names=FALSE);write.csv(dat_primary,'analysis_input_incomplete.csv',row.names=FALSE);write.csv(dat_full,'analysis_input_complete.csv',row.names=FALSE)
 write.csv(validation_raw,'validation_raw.csv',row.names=FALSE)
 saveRDS(list(primary_pools=primary_pools,primary=primary_reg_table,complete=complete_benchmark,case_deletion=case_deletion,forest=forest_df),'analysis_snapshot.rds')
 writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
 writeLines('Completed Litter application: 40 REG imputations; regression only; 500 validation repetitions requested.','COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(stage='main',rep=NA,warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
