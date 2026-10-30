old <- '/workspace/paper-feedback-5-treatments'; new <- '/workspace/paper-feedback-5-treatments-1000'
a <- expand.grid(scen_idx=1:2,mr_idx=1:5,sim=1:10000,KEEP.OUT.ATTRS=FALSE)
b <- expand.grid(scen_idx=1:2,mr_idx=1:5,sim=1:1000,KEEP.OUT.ATTRS=FALSE)
stopifnot(identical(a[1:nrow(b),],b))
for(folder in list.dirs(file.path(new,'full/checkpoints'),recursive=FALSE)) {
 spec <- readRDS(file.path(old,'full/checkpoints',basename(folder),'setting.rds'))
 stopifnot(spec$n_sim==10000,identical(spec$md5,unname(tools::md5sum(file.path(new,'original_script.R')))))
 spec$n_sim <- 1000
 saveRDS(spec,file.path(folder,'setting.rds'))
}
writeLines('Verified exact job-grid prefix and identical original source hash. First 1000 replicates retain their original seeds and results.',file.path(new,'CHECKPOINT_REUSE_VERIFIED.txt'))
