args<-commandArgs(TRUE);root<-normalizePath(args[1]);data('IBScovars',package='DoseFinding');d<-IBScovars
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20261115L)
ids<-sort(sample.int(nrow(d),100L,replace=FALSE));chosen<-data.frame(original_row_id=ids,d[ids,,drop=FALSE]);row.names(chosen)<-NULL
p<-file.path(root,'selected_patients.csv')
if(file.exists(p))stopifnot(identical(read.csv(p)$original_row_id,chosen$original_row_id),all(vapply(names(d),function(nm){a<-read.csv(p)[[nm]];b<-chosen[[nm]];if(is.factor(b)){a<-as.character(a);b<-as.character(b)};isTRUE(all.equal(a,b,check.attributes=FALSE))},logical(1))))else write.csv(chosen,p,row.names=FALSE)
stopifnot(nrow(chosen)==100L,!anyDuplicated(ids),length(unique(chosen$dose))==5L,!anyNA(chosen))
manifest<-list(dataset='DoseFinding::IBScovars',source_n=nrow(d),subset_n=100L,selection='One simple random sample without replacement; not stratified; frozen before analysis',selection_seed=20261115L,rng=RNGkind(),selected_original_rows=ids,patients_md5=unname(tools::md5sum(p)))
jsonlite::write_json(manifest,file.path(root,'selection_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
print(table(chosen$dose));cat('Selection frozen and verified.\n')
