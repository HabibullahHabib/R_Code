args <- commandArgs(TRUE)
root <- normalizePath(args[1]); mode <- if(length(args)>1) args[2] else 'full'
out <- file.path(root,mode); dir.create(out,recursive=TRUE,showWarnings=FALSE); setwd(out)
.h <- new.env()
.h$root <- root; .h$out <- out; .h$mode <- mode; .h$started <- Sys.time()
.h$log <- function(file, rows) {
 p <- file.path(.h$out,file); write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')
}
.h$checkpoint <- function(cl, X, fun, ...) {
 tag <- condition_tag
 folder <- file.path(.h$out,'checkpoints',tag); dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 spec <- list(md5=unname(tools::md5sum(file.path(.h$root,'original_script.R'))),n_sim=n_sim,seed=master_seed,tag=tag)
 meta <- file.path(folder,'setting.rds')
 if(file.exists(meta)) stopifnot(identical(readRDS(meta),spec)) else saveRDS(spec,meta)
 batches <- split(X,ceiling(seq_along(X)/20L)); result <- vector('list',length(X))
 for(b in seq_along(batches)) {
  ids <- batches[[b]]; p <- file.path(folder,sprintf('batch_%05d.rds',b))
  if(file.exists(p)) z <- readRDS(p) else {
   z <- parallel::parLapplyLB(cl,ids,fun,...,chunk.size=1L)
   saveRDS(z,paste0(p,'.tmp'),compress=TRUE); stopifnot(file.rename(paste0(p,'.tmp'),p))
  }
  result[ids] <- z
  for(j in seq_along(z)) if(length(attr(z[[j]],'captured_warnings'))) for(w in attr(z[[j]],'captured_warnings'))
   .h$log('warnings.csv',data.frame(condition=tag,job=ids[j],warning=w))
  saved <- length(list.files(file.path(.h$out,'checkpoints'),pattern='^batch_.*rds$',recursive=TRUE)) * if(.h$mode=='full')20L else 10L
  progress <- data.frame(saved_jobs=saved,total_jobs=if(.h$mode=='full')600000 else 60,condition=tag,condition_jobs=max(ids),updated=as.character(Sys.time()))
  write.csv(progress,file.path(.h$out,'progress.csv.tmp'),row.names=FALSE); file.rename(file.path(.h$out,'progress.csv.tmp'),file.path(.h$out,'progress.csv'))
  message('Checkpoint ',tag,': ',max(ids),' / ',length(X))
 }
 result
}
.h$save_table <- function(data,filename,...) {
 gt::gtsave(data,filename,...)
 write.csv(data$`_data`,sub('html$','csv',filename),row.names=FALSE)
 saveRDS(data$`_data`,sub('html$','rds',filename))
 invisible(data)
}
.h$csv <- function(x,file,...) {
 if(grepl('^raw_',file)) {
  con <- gzfile(paste0(file,'.gz'),'wt'); on.exit(close(con)); utils::write.csv(x,con,...)
 } else utils::write.csv(x,file,...)
}
.h$transform <- function(e) {
 if(!is.call(e)) return(e)
 if(identical(e[[1]],as.name('makeCluster'))) return(quote(parallel::makeCluster(n_cores,outfile=file.path(.h$out,'worker_console.log'))))
 if(identical(e[[1]],as.name('parLapplyLB'))) e[[1]] <- quote(.h$checkpoint)
 if(identical(e[[1]],as.name('gtsave'))) e[[1]] <- quote(.h$save_table)
 if(identical(e[[1]],as.name('write.csv'))) e[[1]] <- quote(.h$csv)
 if(identical(e[[1]],as.name('setdiff')) && length(e)==3 && identical(e[[3]],'cl')) e[[3]] <- c('cl','.h')
 for(i in seq_along(e)[-1]) if(is.call(e[[i]])) e[[i]] <- .h$transform(e[[i]])
 e
}
tryCatch(withCallingHandlers({
 expr <- parse(file.path(root,'original_script.R'))
 for(e in expr) {
  if(is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('n_cores'))) e[[3]] <- 4L
  if(.h$mode=='smoke' && is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('n_sim'))) e[[3]] <- 1L
  if(.h$mode=='smoke' && is.call(e) && identical(e[[1]],as.name('<-')) && as.character(e[[2]]) %in% c('variance_ratios','mar_mechanisms')) e[[3]] <- e[[3]][[2]]
  eval(.h$transform(e),.GlobalEnv)
  if(is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name('run_one_job'))) {
   original_run_one_job <- run_one_job
   run_one_job <- function(i) {
    warnings <- character()
    z <- withCallingHandlers(original_run_one_job(i),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    attr(z,'captured_warnings') <- warnings; z
   }
  }
 }
 writeLines(capture.output(sessionInfo()),'sessionInfo.txt')
 source(file.path(.h$root,'export_report.R'))
 writeLines(paste('Completed',Sys.time()),'COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(condition='main',job=NA,warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){
 .h$log('errors.csv',data.frame(time=as.character(Sys.time()),error=conditionMessage(e)))
 writeLines(conditionMessage(e),'FAILED.txt'); message('FAILED: ',conditionMessage(e));quit(status=1)
})
