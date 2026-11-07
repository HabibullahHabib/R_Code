args<-commandArgs(TRUE);root<-normalizePath(args[1]);out<-file.path(root,'full');dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(out)
.h<-new.env();.h$root<-root;.h$out<-out
.h$log<-function(file,rows){p<-file.path(.h$out,file);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
.h$saveplot<-function(filename,plot,...){ggplot2::ggsave(filename,plot,...);ggplot2::ggsave(paste0(tools::file_path_sans_ext(filename),'.pdf'),plot,...)}
.h$make<-function(dat,method,m=10000L,base_seed,pmm_donors=5L,propensity_groups=3L){
 stopifnot(m==10000L)
 tag<-paste0(method,'_',base_seed);folder<-file.path(.h$out,'checkpoints',tag);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 spec<-list(source=unname(tools::md5sum(file.path(.h$root,'original_script.R'))),data=digest::digest(dat),method=method,m=m,base_seed=base_seed,pmm_donors=pmm_donors,propensity_groups=propensity_groups)
 meta<-file.path(folder,'setting.rds');if(file.exists(meta))stopifnot(identical(spec,readRDS(meta)))else saveRDS(spec,meta)
 cl<-parallel::makeCluster(2,port=as.integer(Sys.getenv('SIM_CLUSTER_PORT')),outfile=file.path(.h$out,'worker_console.log'));on.exit(parallel::stopCluster(cl))
 parallel::clusterEvalQ(cl,{library(dplyr);library(tibble)})
 vars<-ls(.GlobalEnv);vars<-vars[vapply(mget(vars,.GlobalEnv),is.function,logical(1))];parallel::clusterExport(cl,vars,envir=.GlobalEnv)
 fun<-function(mm,dat,method,base_seed,pmm_donors,propensity_groups){
  ss<-base_seed+mm-1L;warnings<-character()
  value<-withCallingHandlers(switch(method,regression=impute_regression_once(dat,seed=ss),pmm=impute_pmm_once(dat,seed=ss,donors=pmm_donors),propensity=impute_propensity_once(dat,seed=ss,n_groups=propensity_groups),stop('Unknown imputation method.')),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  list(value=value,warnings=warnings,index=mm)
 };environment(fun)<-.GlobalEnv
 comps<-vector('list',m)
 for(start in seq(1L,m,100L)){
  ids<-start:min(start+99L,m);p<-file.path(folder,sprintf('batch_%05d.rds',ceiling(start/100)))
  if(file.exists(p))z<-readRDS(p)else{
   z<-parallel::parLapplyLB(cl,ids,fun,dat=dat,method=method,base_seed=base_seed,pmm_donors=pmm_donors,propensity_groups=propensity_groups,chunk.size=1L)
   saveRDS(z,paste0(p,'.tmp'));stopifnot(file.rename(paste0(p,'.tmp'),p))
  }
  writeLines(as.character(length(z)),sub('rds$','count',p))
  comps[ids]<-lapply(z,`[[`,'value')
  for(a in z)for(w in a$warnings).h$log('warnings.csv',data.frame(analysis=tag,imputation=a$index,warning=w))
  if(start==1L){
   original<-.h$original_make(dat,method,m=3L,base_seed=base_seed,pmm_donors=pmm_donors,propensity_groups=propensity_groups)
   stopifnot(identical(original,comps[1:3]));writeLines('Exact original/parallel match for the first three imputations.',file.path(folder,'FIDELITY_PASSED.txt'))
  }
  paths<-list.files(file.path(.h$out,'checkpoints'),pattern='batch_.*count$',recursive=TRUE,full.names=TRUE)
  saved<-sum(vapply(paths,function(p)as.integer(readLines(p)),integer(1)))
  total<-if(grepl('actg175',.h$root))30000L else 120000L
  write.csv(data.frame(saved_imputations=saved,total_imputations=total,analysis=tag,analysis_imputations=max(ids),updated=as.character(Sys.time())),'progress.csv.tmp',row.names=FALSE);file.rename('progress.csv.tmp','progress.csv')
  message(tag,': ',max(ids),' / ',m,' imputations; ',saved,' / ',total,' overall')
 }
 comps
}
.h$transform<-function(e){
 if(is.character(e))return(gsub('m=40','m=10000',e,fixed=TRUE))
 if(!is.call(e))return(e)
 if(identical(e[[1]],as.name('ggsave'))||identical(e[[1]],quote(ggplot2::ggsave)))e[[1]]<-quote(.h$saveplot)
 for(i in seq_along(e)[-1])if(is.call(e[[i]])||is.character(e[[i]]))e[[i]]<-.h$transform(e[[i]])
 e
}
tryCatch(withCallingHandlers({
 for(e in parse(file.path(root,'original_script.R'))){
  if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('m_imputations')))e[[3]]<-10000L
  .h$is_make<-is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('make_imputations'))
  eval(.h$transform(e),.GlobalEnv)
  if(.h$is_make){.h$original_make<-make_imputations;make_imputations<-.h$make}
 }
 source(file.path(.h$root,'export_10000.R'))
 writeLines(paste('Completed all application analyses with 10000 imputations:',Sys.time()),'COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(analysis='main',imputation=NA,warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(time=as.character(Sys.time()),error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
