args<-commandArgs(TRUE);root<-normalizePath(args[1]);mode<-if(length(args)>1)args[2]else'full';app<-if(grepl('actg175',root))'ACTG175'else'IBScovars'
out<-file.path(root,mode);dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(out)
suppressPackageStartupMessages({library(mvtnorm);library(dplyr);library(tibble)})
# Load unchanged scientific functions only; do not execute the application's main analysis.
for(e in parse(file.path(root,'original_application.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e,.GlobalEnv)
data_env<-new.env()
if(app=='ACTG175'){
 data('ACTG175',package='BART',envir=data_env);d<-data_env$ACTG175
 study_data<-tibble(group=factor(d$arms,levels=sort(unique(d$arms))),x=as.numeric(scale(as.numeric(d$cd40))),y=as.numeric(d$cd496))
 stopifnot(!anyNA(study_data$x));benchmark<-NULL;missing_design<-'Natural missingness retained exactly; repetitions change imputation draws only'
}else{
 data('IBScovars',package='DoseFinding',envir=data_env);d<-data_env$IBScovars
 study_data<-tibble(group=factor(d$dose,levels=sort(unique(d$dose))),x=as.numeric(scale(as.numeric(d$gender))),y=as.numeric(d$resp))
 stopifnot(!anyNA(study_data));means<-tapply(study_data$y,study_data$group,mean);benchmark<-means[-1]-means[1];names(benchmark)<-paste0(names(means)[-1],' vs ',names(means)[1])
 missing_design<-'New Bernoulli 30% expected covariate-driven MAR pattern per repetition; gamma=log(2), empirical intercept calibration'
}
control_level<-levels(study_data$group)[1];m_imputations<-40L;n_experiments<-if(mode=='smoke')2L else 10000L;master_seed<-20261108L;dataset_offset<-if(app=='ACTG175')10000000L else 20000000L
mar_gamma<-log(2);mar_intercept<-if(app=='IBScovars')uniroot(function(a)mean(plogis(a+mar_gamma*study_data$x))-.30,c(-30,30),tol=1e-12)$root else NA_real_
spec<-list(dataset=app,input_hash=digest::digest(study_data),scientific_source=unname(tools::md5sum(file.path(root,'original_application.R'))),experiments=n_experiments,imputations=40L,method='regression',master_seed=master_seed,dataset_offset=dataset_offset,missing_design=missing_design,mar_gamma=if(app=='IBScovars')mar_gamma else NA_real_,mar_intercept=mar_intercept,comparison=c('Shrinkage','Non-shrinkage'))
dir.create('checkpoints',showWarnings=FALSE);if(file.exists('protocol.rds'))stopifnot(identical(spec,readRDS('protocol.rds')))else saveRDS(spec,'protocol.rds')
jsonlite::write_json(spec,'protocol.json',pretty=TRUE,auto_unbox=TRUE,na='null');write.csv(study_data,'analysis_input_data.csv',row.names=FALSE)
summary_data<-study_data|>group_by(group)|>summarise(n=n(),observed=sum(!is.na(y)),missing=sum(is.na(y)),missing_pct=100*mean(is.na(y)),.groups='drop');write.csv(summary_data,'dataset_summary.csv',row.names=FALSE)
if(!is.null(benchmark))write.csv(data.frame(contrast=names(benchmark),benchmark=as.numeric(benchmark)),'complete_data_benchmark.csv',row.names=FALSE)
run_experiment<-function(sim){
 warnings<-character();missing_seed<-master_seed+dataset_offset+sim;imputation_seed<-master_seed+dataset_offset+20000000L+sim*100L;qseed<-master_seed+dataset_offset+50000000L+sim*10L
 tryCatch(withCallingHandlers({
  dat<-study_data
  if(app=='IBScovars'){set.seed(missing_seed);mask<-rbinom(nrow(dat),1L,plogis(mar_intercept+mar_gamma*dat$x));dat$y[mask==1L]<-NA_real_}
  comps<-make_imputations(dat,method='regression',m=40L,base_seed=imputation_seed)
  intervals<-pool_many_to_one_both(comps,control_level=control_level,alpha=.05,qseed_shrink=qseed+1L,qseed_nonshrink=qseed+2L)
  sh<-intervals[intervals$variance_method=='Shrinkage',];ns<-intervals[intervals$variance_method=='Non-shrinkage',]
  stopifnot(identical(sh$contrast,ns$contrast),identical(sh$estimate,ns$estimate),all(is.finite(intervals$lower)),all(is.finite(intervals$upper)),all(intervals$width>0))
  intervals$sim<-sim
  if(!is.null(benchmark)){intervals$benchmark<-as.numeric(benchmark[intervals$contrast]);intervals$benchmark_recovered<-intervals$lower<=intervals$benchmark&intervals$upper>=intervals$benchmark}
  avgsh<-mean(sh$width);avgns<-mean(ns$width)
  paired<-data.frame(sim=sim,achieved_missing_rate=mean(is.na(dat$y)),mean_width_shrinkage=avgsh,mean_width_nonshrinkage=avgns,width_saving_pct=100*(avgns-avgsh)/avgns,shrinkage_narrower=avgsh<avgns,joint_recovery_shrinkage=if(is.null(benchmark))NA else all(sh$lower<=benchmark[sh$contrast]&sh$upper>=benchmark[sh$contrast]),joint_recovery_nonshrinkage=if(is.null(benchmark))NA else all(ns$lower<=benchmark[ns$contrast]&ns$upper>=benchmark[ns$contrast]))
  list(ok=TRUE,sim=sim,paired=paired,intervals=intervals,warnings=warnings,error=NA_character_,missing_seed=missing_seed,imputation_seed=imputation_seed,qseed=qseed)
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(ok=FALSE,sim=sim,warnings=warnings,error=conditionMessage(e),missing_seed=missing_seed,imputation_seed=imputation_seed,qseed=qseed))
}
tryCatch({
 cl<-parallel::makeCluster(2,port=as.integer(Sys.getenv('SIM_CLUSTER_PORT')),outfile='worker_console.log')
 parallel::clusterEvalQ(cl,{library(mvtnorm);library(dplyr);library(tibble)})
 vars<-setdiff(ls(.GlobalEnv),c('cl'));parallel::clusterExport(cl,vars,envir=.GlobalEnv)
 results<-vector('list',n_experiments)
 for(start in seq(1L,n_experiments,10L)){
  ids<-start:min(start+9L,n_experiments);p<-file.path('checkpoints',sprintf('batch_%05d.rds',ceiling(start/10)))
  if(file.exists(p))z<-readRDS(p)else{z<-parallel::parLapplyLB(cl,ids,run_experiment,chunk.size=1L);saveRDS(z,paste0(p,'.tmp'));stopifnot(file.rename(paste0(p,'.tmp'),p))}
  writeLines(as.character(length(z)),sub('rds$','count',p));results[ids]<-z
  saved<-sum(vapply(list.files('checkpoints',pattern='count$',full.names=TRUE),function(p)as.integer(readLines(p)),integer(1)))
  current<-results[seq_len(max(ids))];failed<-sum(vapply(current,function(z)!isTRUE(z$ok),logical(1)));nw<-sum(vapply(current,function(z)length(z$warnings),integer(1)))
  write.csv(data.frame(saved_experiments=saved,total_experiments=n_experiments,processed_experiments=max(ids),errors_recorded=failed,warnings_recorded=nw,updated=as.character(Sys.time())),'progress.csv.tmp',row.names=FALSE);file.rename('progress.csv.tmp','progress.csv')
  message(app,': ',max(ids),' / ',n_experiments,' repetitions saved')
 }
 parallel::stopCluster(cl)
 # Verify a serial experiment against the checkpointed parallel result.
 serial<-run_experiment(1L);stopifnot(identical(serial,results[[1]]));writeLines('Original scientific functions: serial and parallel repetition 1 match exactly. Shrinkage/non-shrinkage share all 40 imputations and point estimates.','FIDELITY_PASSED.txt')
 source(file.path(root,'export_report.R'))
 writeLines(paste('Completed',n_experiments,'repetitions at',Sys.time()),'COMPLETED.txt')
},error=function(e){writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
