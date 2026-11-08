args<-commandArgs(TRUE);root<-normalizePath(args[1]);out<-file.path(root,'full');dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(out)
.h<-new.env();.h$root<-root;.h$out<-out
.h$log<-function(file,rows){p<-file.path(.h$out,file);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
.h$saveplot<-function(filename,plot,...){ggplot2::ggsave(filename,plot,...);ggplot2::ggsave(paste0(tools::file_path_sans_ext(filename),'.pdf'),plot,...)}
.h$transform<-function(e){
 if(is.character(e))return(gsub('m=40','m=40',e,fixed=TRUE))
 if(!is.call(e))return(e)
 if(identical(e[[1]],as.name('ggsave'))||identical(e[[1]],quote(ggplot2::ggsave)))e[[1]]<-quote(.h$saveplot)
 for(i in seq_along(e)[-1])if(is.call(e[[i]])||is.character(e[[i]]))e[[i]]<-.h$transform(e[[i]])
 e
}
source(file.path(root,'alignment_helpers.R'))
tryCatch(withCallingHandlers({
 for(e in parse(file.path(root,'original_script.R'))){
  if(is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('m_imputations')))e[[3]]<-40L
  .h$assigned<-if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]]))as.character(e[[2]]) else ''
  eval(.h$transform(e),.GlobalEnv)
  if(.h$assigned=='make_imputations'){.h$original_make<-make_imputations;make_imputations<-.h$make}
  if(.h$assigned=='pool_many_to_one_both'){.h$original_pool<-pool_many_to_one_both;pool_many_to_one_both<-.h$pool}
  if(.h$assigned=='save_primary_outputs'){
   .h$save_original<-save_primary_outputs
   save_primary_outputs<-function(long_results,primary_table,output_dir,prefix,title,subtitle).h$save_original(long_results,primary_table,output_dir,paste0(prefix,'_SUPPLEMENTARY'),paste0(title,' (supplementary methods)'),subtitle)
  }
 }
 source(file.path(.h$root,'export_report.R'))
 writeLines(paste('Thesis-aligned adaptation completed, m=40:',Sys.time()),'COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(time=as.character(Sys.time()),warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(time=as.character(Sys.time()),error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
