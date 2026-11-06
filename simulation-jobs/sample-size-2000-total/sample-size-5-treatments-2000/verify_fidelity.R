args <- commandArgs(TRUE);root<-args[1]
for(e in parse(file.path(root,'original_script.R'))) {if(is.call(e)&&identical(e[[1]],as.name('for')))break;eval(e,.GlobalEnv)}
n_sim<-1L;jobs<-expand.grid(scen_idx=1:2,mr_idx=1:5,sim=1L,KEEP.OUT.ATTRS=FALSE)
z<-readRDS(file.path(root,'smoke/checkpoints/n20/batch_00001.rds'))
for(i in c(1L,3L)) {a<-run_one_job_safe(i);b<-z[[i]];attr(b,'captured_warnings')<-NULL;stopifnot(identical(a,b));cat('Original/checkpoint exact match:',i,'\n')}
writeLines('Exact original/checkpoint match for complete-data and MAR-imputation jobs.',file.path(root,'FIDELITY_PASSED.txt'))
