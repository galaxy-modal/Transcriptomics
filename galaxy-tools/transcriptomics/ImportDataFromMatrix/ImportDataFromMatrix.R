library(Biobase)
library(GEOquery)
library(GEOmetadb)
library(limma)
library(jsonlite)
library(affy)
library(dplyr)
library(affyPLM)

cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]
nbargs=length(cargs)

dataFile=cargs[[nbargs-9]]
normalization=cargs[[nbargs-8]]
conditionsFile=cargs[[nbargs-7]]
annotation=cargs[[nbargs-6]]
result_export_eset=cargs[[nbargs-5]]
result=cargs[[nbargs-4]]
result.path=cargs[[nbargs-3]]
result.template=cargs[[nbargs-2]]

dir.create(result.path, showWarnings = TRUE, recursive = FALSE)

data=as.matrix(read.table(file = dataFile))
conditions=read.table(file=conditionsFile,sep = "\t",row.names=1)
htmlfile=readChar(result.template, file.info(result.template)$size)

colnames(conditions)=c("source_name_ch1","description")
phenodata<-new("AnnotatedDataFrame",data=conditions)

eset=ExpressionSet(assayData=data,phenoData=phenodata,annotation=annotation)

if (normalization == "quantile") {
	eset <- normalize.ExpressionSet.quantiles(eset, transfn="log")
} else if (normalization == "log2") {
	exprs(eset) = log2(exprs(eset)) 
} 

boxplotnorm="boxplotnorm.png"
png(boxplotnorm,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(data.frame(exprs(eset)),las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOTNORM###",replacement = boxplotnorm, fixed = TRUE)
file.copy(boxplotnorm,result.path)

plotMAnorm="plotMAnorm.png"
nblines=length(colnames(data))%/%3 + as.numeric((length(colnames(data))%%3)!=0) 
png(plotMAnorm,width=800,height =300*nblines )
par(mfrow=c(nblines,3))
#for (i in 1:length(colnames(data))){
	MAplot(eset)
#}

dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###PLOTMANORM###",replacement = plotMAnorm, fixed = TRUE)
file.copy(plotMAnorm,result.path)
#write.table(tolower(c(condition1Name,condition2Name)),quote = FALSE,col.names = FALSE, row.names=FALSE,file=result_export_conditions)
#saveConditions=c(condition1Name,condition2Name)
save(eset,file=result_export_eset)
write(htmlfile,result)

#l=list()
#for(i in 1:length(esets))
#{
#	l[[paste("study",i,sep="")]]<-res[[i]]
#}
#l[["Meta"]]=res[[length(res)-1]]
#showVenn(res,file.path(temp.files.path,"venn.png"))
#writeLines(c("<h2>Venn diagram</h2>"),file.conn)
#writeLines(c("<img src='venn.png'><br/><br/>"),file.conn)
#writeLines(c("</body></html>"),file.conn)
#close(file.conn)