args<-commandArgs(TRUE);root<-normalizePath(args[1]);mode<-if(length(args)>1)args[2] else 'full'
out<-file.path(root,mode);dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(out)
.h<-new.env();.h$root<-root;.h$out<-out;.h$mode<-mode
.h$log<-function(name,rows){p<-file.path(.h$out,name);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
.h$saveplot<-function(filename,plot,...){ggplot2::ggsave(filename,plot,...);ggplot2::ggsave(sub('\\.png$','.pdf',filename),plot,...);invisible(filename)}
.h$validation_loop<-function(body){
 tag<-if(exists('r_idx',.GlobalEnv))paste0('MR',round(100*mr)) else 'MR30'
 folder<-file.path(.h$out,'checkpoints',tag);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 spec<-list(source_md5=unname(tools::md5sum(file.path(.h$root,'original_script.R'))),reps=validation_reps,m=m_imputations,seed=master_seed,tag=tag)
 meta<-file.path(folder,'setting.rds');if(file.exists(meta))stopifnot(identical(readRDS(meta),spec))else saveRDS(spec,meta)
 cl<-parallel::makeCluster(4,port=as.integer(Sys.getenv('SIM_CLUSTER_PORT')),outfile=file.path(.h$out,'worker_console.log'));on.exit(parallel::stopCluster(cl))
 parallel::clusterEvalQ(cl,{library(mvtnorm);library(dplyr);library(tidyr);library(tibble);library(ggplot2)})
 vars<-setdiff(ls(.GlobalEnv,all.names=TRUE),c('.h'))
 parallel::clusterExport(cl,vars,envir=.GlobalEnv)
 fun<-function(b,expr){
  e<-new.env(parent=.GlobalEnv);e$b<-b;e$raw_validation<-list();e$cc<-1L;warnings<-character()
  result<-withCallingHandlers({eval(expr,e);e$raw_validation},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  list(results=result,warnings=warnings,rep=b)
 };environment(fun)<-.GlobalEnv
 batches<-split(seq_len(validation_reps),ceiling(seq_len(validation_reps)/5))
 for(i in seq_along(batches)){
  ids<-batches[[i]];p<-file.path(folder,sprintf('batch_%05d.rds',i))
  if(file.exists(p))z<-readRDS(p) else {
   z<-parallel::parLapplyLB(cl,ids,fun,expr=body,chunk.size=1L)
   saveRDS(z,paste0(p,'.tmp'));stopifnot(file.rename(paste0(p,'.tmp'),p))
  }
  writeLines(as.character(length(z)),sub('rds$','count',p))
  for(a in z){
   for(res in a$results){raw_validation[[cc]]<<-res;cc<<-cc+1L}
   for(w in a$warnings).h$log('warnings.csv',data.frame(stage=tag,rep=a$rep,warning=w))
  }
  counts<-list.files(file.path(.h$out,'checkpoints'),pattern='batch_.*count$',recursive=TRUE,full.names=TRUE)
  saved<-sum(vapply(counts,function(p)as.integer(readLines(p)),integer(1)))
  progress<-data.frame(saved_experiments=saved,total_experiments=if(.h$mode=='full')1000 else if(grepl('ibscovars',.h$root))4 else 1,condition=tag,updated=as.character(Sys.time()))
  write.csv(progress,file.path(.h$out,'progress.csv.tmp'),row.names=FALSE);file.rename(file.path(.h$out,'progress.csv.tmp'),file.path(.h$out,'progress.csv'))
  message('Saved validation experiments: ',saved,' | ',tag)
 }
 invisible(NULL)
}
.h$transform<-function(e){
 if(!is.call(e))return(e)
 if(identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('b'))){return(substitute(.h$validation_loop(quote(BODY)),list(BODY=e[[4]])))}
 if(identical(e[[1]],as.name('ggsave'))||identical(e[[1]],quote(ggplot2::ggsave)))e[[1]]<-quote(.h$saveplot)
 for(i in seq_along(e)[-1])if(is.call(e[[i]]))e[[i]]<-.h$transform(e[[i]])
 e
}
tryCatch(withCallingHandlers({
 for(e in parse(file.path(root,'original_script.R'))){
  if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('validation_reps'))) e[[3]]<-if(.h$mode=='smoke')1L else if(grepl('ibscovars',root))250L else 1000L
  eval(.h$transform(e),.GlobalEnv)
 }
 source(file.path(root,'export_report.R'))
 writeLines(paste('Completed',Sys.time()),'COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(stage='main',rep=NA,warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(time=as.character(Sys.time()),error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
