args <- commandArgs(TRUE); root <- args[1]
e <- parse(file.path(root,'original_script.R'))
for(x in e) {
 if(is.call(x) && identical(x[[1]],as.name('for'))) break
 eval(x,.GlobalEnv)
}
n_sim <- 1L
jobs <- expand.grid(scen_idx=seq_along(scenarios),mr_idx=seq_along(missing_rates),sim=1L,KEEP.OUT.ATTRS=FALSE)
z <- readRDS(file.path(root,'smoke/checkpoints/R2_control_MCAR_treatment_MAR/batch_00001.rds'))
for(i in c(1L,3L)) {
 a <- run_one_job(i);b <- z[[i]];attr(b,'captured_warnings') <- NULL
 stopifnot(identical(a,b));cat('Exact original/checkpoint match, job',i,'\n')
}
writeLines('PASS: exact original/checkpoint equality for complete data and MAR imputation jobs.',file.path(root,'FIDELITY_PASSED.txt'))
