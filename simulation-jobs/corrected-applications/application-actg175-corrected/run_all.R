args<-commandArgs(TRUE);root<-normalizePath(args[1]);mode<-if(length(args)>1)args[2] else 'as-written'
out<-file.path(root,mode);dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(out)
.h<-new.env();.h$root<-root;.h$out<-out
.h$log<-function(file,rows){p<-file.path(.h$out,file);write.table(rows,p,sep=',',row.names=FALSE,col.names=!file.exists(p),append=file.exists(p),qmethod='double')}
.h$saveplot<-function(filename,plot,...){ggplot2::ggsave(filename,plot,...);ggplot2::ggsave(sub('\\.png$','.pdf',filename),plot,...)}
.h$transform<-function(e){
 if(!is.call(e))return(e)
 if(identical(e[[1]],as.name('ggsave'))||identical(e[[1]],quote(ggplot2::ggsave)))e[[1]]<-quote(.h$saveplot)
 for(i in seq_along(e)[-1])if(is.call(e[[i]]))e[[i]]<-.h$transform(e[[i]])
 e
}
tryCatch(withCallingHandlers({
 for(e in parse(file.path(root,'original_script.R')))eval(.h$transform(e),.GlobalEnv)
 source(file.path(root,'export_report.R'))
 writeLines(paste('As-written analysis completed with',m_imputations,'imputations at',Sys.time()),'COMPLETED.txt')
},warning=function(w){.h$log('warnings.csv',data.frame(time=as.character(Sys.time()),warning=conditionMessage(w)));invokeRestart('muffleWarning')}),error=function(e){.h$log('errors.csv',data.frame(time=as.character(Sys.time()),error=conditionMessage(e)));writeLines(conditionMessage(e),'FAILED.txt');message('FAILED: ',conditionMessage(e));quit(status=1)})
