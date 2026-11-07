args<-commandArgs(TRUE);root<-args[1];dir.create(file.path(root,'fidelity'),showWarnings=FALSE);setwd(file.path(root,'fidelity'))
for(e in parse(file.path(root,'original_script.R'))){if(is.call(e)&&identical(e[[1]],as.name('for')))break;eval(e,.GlobalEnv)}
if(grepl('actg175',root)){
 validation_base<-actg|>filter(!is.na(y))|>select(-missing_y)|>droplevels();truth<-benchmark_contrasts(validation_base,control_level)
 dat_miss<-impose_mcar_within_arm(validation_base,rate=.30,seed=master_seed+5000000L+1L)
 seed<-master_seed+10000000L+1000L+100L;scenario<-'Validation MCAR 30%';rate<-.30;tag<-'MR30'
}else{
 truth<-benchmark_contrasts(ibs,control_level);dat_miss<-impose_mcar_within_arm(ibs,rate=.10,seed=master_seed+50000000L+100000L+1L)
 seed<-master_seed+100000000L+1000000L+1000L+100L;scenario<-'Repeated MCAR 10%';rate<-.10;tag<-'MR10'
}
a<-make_mi_results(dat_miss,method='regression',method_label='Regression MI',m=m_imputations,control_level=control_level,outcome_covariates=outcome_covariates,propensity_covariates=propensity_covariates,pmm_donors=pmm_donors,propensity_ngroups=propensity_ngroups,scenario_label=scenario,missing_rate=rate,seed_base=seed)|>add_benchmark_metrics(truth)|>mutate(rep=1L)
b<-readRDS(file.path(root,'smoke/checkpoints',tag,'batch_00001.rds'))[[1]]$results[[1]]
stopifnot(identical(a,b));writeLines('PASS: original serial validation calculation exactly matches checkpointed worker output.',file.path(root,'FIDELITY_PASSED.txt'));cat('PASS: exact original/worker equality\n')
