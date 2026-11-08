.h$equal_pool<-function(comp_list,control_level,alpha=.05,seed=NULL){
 st<-get_completed_stats(comp_list);groups<-colnames(st$means);tr<-setdiff(groups,control_level);labels<-paste0(tr,' vs ',control_level);k<-length(tr);m<-st$m
 Q<-sweep(st$means[,tr,drop=FALSE],1,st$means[,control_level],'-');Qbar<-colMeans(Q);B<-stats::cov(Q)
 weights<-st$n-1;df_resid<-sum(weights)
 pooled<-as.numeric(st$vars%*%weights/df_resid)
 A<-matrix(1/st$n[control_level],k,k);diag(A)<-1/st$n[tr]+1/st$n[control_level]
 Ubar<-mean(pooled)*A;Tm<-Ubar+(1+1/m)*B;Tm<-(Tm+t(Tm))/2
 rm<-(1+1/m)*sum(diag(B%*%safe_solve(Ubar)))/k
 stopifnot(is.finite(rm),rm>=-1e-12,all(eigen(Tm,symmetric=TRUE,only.values=TRUE)$values>0))
 nu<-if(max(abs(B))<=1e-14||rm<=1e-14)1e6 else (m-1)*(1+1/rm)^2
 crit<-safe_qmvt_crit(make_corr_safe(Tm),df=nu,alpha=alpha,seed=seed);se<-sqrt(diag(Tm));lower<-Qbar-crit*se;upper<-Qbar+crit*se
 tibble::tibble(contrast=labels,estimate=as.numeric(Qbar),lower=as.numeric(lower),upper=as.numeric(upper),width=as.numeric(upper-lower),significant=as.numeric(lower)>0|as.numeric(upper)<0,variance_method='Equal-variance',shrink_alpha=NA_real_)
}
.h$pool<-function(comp_list,control_level,alpha=.05,qseed_shrink=NULL,qseed_nonshrink=NULL){
 original<-.h$original_pool(comp_list,control_level,alpha,qseed_shrink,qseed_nonshrink)
 equal<-.h$equal_pool(comp_list,control_level,alpha,qseed_nonshrink)
 stopifnot(isTRUE(all.equal(equal$estimate,original$estimate[original$variance_method=='Shrinkage'],check.attributes=FALSE)))
 dplyr::bind_rows(original,equal)
}
.h$make<-function(dat,method,m=40L,base_seed,pmm_donors=5L,propensity_groups=3L){
 stopifnot(m==40L);tag<-paste0(method,'_',base_seed);folder<-file.path(.h$out,'checkpoints',tag);dir.create(folder,recursive=TRUE,showWarnings=FALSE)
 spec<-list(source=unname(tools::md5sum(file.path(.h$root,'original_script.R'))),data=digest::digest(dat),method=method,m=m,base_seed=base_seed,pmm_donors=pmm_donors,propensity_groups=propensity_groups)
 meta<-file.path(folder,'setting.rds');if(file.exists(meta))stopifnot(identical(spec,readRDS(meta)))else saveRDS(spec,meta)
 cache<-file.path(folder,'imputations.rds');provenance<-'New computation from original scientific function'
 if(file.exists(cache))comps<-readRDS(cache) else {
  previous<-sub('thesis-aligned','corrected',.h$root,fixed=TRUE);old<-file.path(previous,'full/checkpoints',tag)
  if(file.exists(file.path(old,'batch_00001.rds'))){
   prior<-readRDS(file.path(old,'setting.rds'));compare<-prior;compare$m<-40L;stopifnot(identical(compare,spec))
   z<-readRDS(file.path(old,'batch_00001.rds'));comps<-lapply(z[1:40],`[[`,'value');provenance<-'First 40 seeded imputations from verified 10000-imputation checkpoint'
  }else comps<-.h$original_make(dat,method,m,base_seed,pmm_donors,propensity_groups)
  stopifnot(length(comps)==40L)
  check<-.h$original_make(dat,method,m=3L,base_seed=base_seed,pmm_donors=pmm_donors,propensity_groups=propensity_groups);stopifnot(identical(check,comps[1:3]))
  saveRDS(comps,paste0(cache,'.tmp'));stopifnot(file.rename(paste0(cache,'.tmp'),cache))
 }
 writeLines(provenance,file.path(folder,'PROVENANCE.txt'))
 writeLines('First three imputations exactly match the original serial scientific function.',file.path(folder,'FIDELITY_PASSED.txt'))
 comps
}
