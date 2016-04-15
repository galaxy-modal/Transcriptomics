library(Biobase)
library(GEOquery)
library(GEOmetadb)
library(limma)
library(jsonlite)
library(affy)
library(dplyr)

cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]
nbargs=length(cargs)
celList=vector()
celFileNameList=vector()
for (i in seq(1,nbargs-7,2))
{
	celList=c(celList,cargs[[i]])
	celFileNameList=c(celFileNameList,cargs[[i+1]])
}


normalization=cargs[[nbargs-6]]
result_export_eset=cargs[[nbargs-5]]
result=cargs[[nbargs-4]]
result.path=cargs[[nbargs-3]]
result.template=cargs[[nbargs-2]]

dir.create(result.path, showWarnings = TRUE, recursive = FALSE)
for(i in 1:length(celList))
{
	file.copy(celList[i],paste0("./",celFileNameList[i]))
}

data <- ReadAffy(filenames=celFileNameList, celfile.path=".")
htmlfile=readChar(result.template, file.info(result.template)$size)

boxplot="boxplot.png"
png(boxplot,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(data,las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOT###",replacement = boxplot, fixed = TRUE)
file.copy(boxplot,result.path)

images="images.png"
nblines=length(celList)%/%4 + as.numeric((length(celList)%%4)!=0) 
png(images,width=800,height = 200*nblines)
par(mfrow=c(nblines,4))
image(data)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###IMAGES###",replacement = images, fixed = TRUE)
file.copy(images,result.path)


plotMA="plotMA.png"
nblines=length(celList)%/%3 + as.numeric((length(celList)%%3)!=0) 
png(plotMA,width=800,height =300*nblines )
par(mfrow=c(nblines,3))
MAplot(data)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###PLOTMA###",replacement = plotMA, fixed = TRUE)
file.copy(plotMA,result.path)


if (normalization == "rma") {
	eset <- rma(data)
} else if (normalization == "quantile") {
	eset = rma(data,background = FALSE,normalize = TRUE)
} else if (normalization == "background"){
	eset = rma(data,background = TRUE ,normalize = FALSE)
} else if (normaization == "log2") {
	eset = rma(data,background = FALSE ,normalize = FALSE)
}
	

boxplotnorm="boxplotnorm.png"
png(boxplotnorm,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(data.frame(exprs(eset)),las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOTNORM###",replacement = boxplotnorm, fixed = TRUE)
file.copy(boxplotnorm,result.path)

plotMAnorm="plotMAnorm.png"
nblines=length(celList)%/%3 + as.numeric((length(celList)%%3)!=0) 
png(plotMAnorm,width=800,height =300*nblines )
par(mfrow=c(nblines,3))
for (i in 1:length(celList)){
plotMA(eset,i)
}

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