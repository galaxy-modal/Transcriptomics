library(Biobase)
library(GEOquery)
library(GEOmetadb)
library(limma)
library(jsonlite)
library(affy)
library(dplyr)
#source("http://bioconductor.org/biocLite.R")
#library(DeSousa2013)


cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]
nbargs=length(cargs)
celList=vector()
celFileNameList=vector()
for (i in seq(1,nbargs-10,2))
{
	celList=c(celList,cargs[[i]])
	celFileNameList=c(celFileNameList,cargs[[i+1]])
}


targetFile=cargs[[nbargs-9]]
condition1=cargs[[nbargs-8]]
condition2=cargs[[nbargs-7]]
result_export_eset=cargs[[nbargs-6]]
result_export_conditions=cargs[[nbargs-5]]
result=cargs[[nbargs-4]]
result.path=cargs[[nbargs-3]]
result.template=cargs[[nbargs-2]]

print(targetFile)
file.copy(targetFile,"./targetFile.txt")

condition1Name="condition1"
condition2Name="condition2"

condition1_tmp <- strsplit(condition1,",")
condition1 <-unlist(condition1_tmp)
#condition1_vecstring = sub("^([^.]*).*", "\\1", condition1_tmp_vecstring) 

condition2_tmp <- strsplit(condition2,",")
condition2<-unlist(condition2_tmp)
#condition2_vecstring = sub("^([^.]*).*", "\\1", condition2_tmp_vecstring) 

conditions=c(condition1,condition2)

# We keep only the files matching the conditions
#filesToKeep=which(celFileNameList %in% conditions)
#celList=celList[filesToKeep]
#celFileNameList=celFileNameList[filesToKeep]

nbresult=1000
dir.create(result.path, showWarnings = TRUE, recursive = FALSE)
for(i in 1:length(celList))
{
	print(celList[i])
	print(celFileNameList[i])
	file.copy(celList[i],paste0("./",celFileNameList[i]))
}


#normalization<-function(data){
#	ex <- exprs(data)
#	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#	LogC <- (qx[5] > 100) ||
#			(qx[6]-qx[1] > 50 && qx[2] > 0) ||
#			(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
#	if (LogC) { ex[which(ex <= 0)] <- NaN
#		return (log2(ex)) } else {
#		return (ex)
#	}
#}
targets <- readTargets("targetFile.txt",path=".",sep="\t")
ab <- ReadAffy(filenames=targets[,1], celfile.path=".")
eset <- rma(ab)
#if (normalization=="rma") {
#	eset <- rma(ab)
#} else if (normalization=="quantile") {
#	eset=expresso()
#} else if (normalization=="log2") {
#		
#} else if (normalization=="none") {
#	
#}


#exprset=exprs(eset)
#sink("/dev/null")
#dir.create(result.file.path,recursive=TRUE)
eset=eset[,which(rownames(eset@phenoData@data) %in% conditions)]
eset@phenoData@data$source_name_ch1=""
#print(which(rownames(eset@phenoData@data) %in% condition1))
eset@phenoData@data$source_name_ch1[which(rownames(eset@phenoData@data) %in% condition1)]=condition1Name
eset@phenoData@data$source_name_ch1[which(rownames(eset@phenoData@data) %in% condition2)]=condition2Name
sml=paste0("G",as.numeric(tolower(as.character(pData(eset)["source_name_ch1"][,1]))!=condition1Name))
print(pData(eset))
fl <- as.factor(sml)
eset$description <- fl

design <- model.matrix(~ description + 0, eset)

colnames(design) <- levels(fl)
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nbresult)
annot <- annotation(eset)
mapping=read.csv("/galaxy-tools/transcriptomics/db/gplToBioc.csv",stringsAsFactors=FALSE)
gpl=mapping[which(mapping$bioc_package==annotation(eset)),]$gpl
gpl=gpl[1]
#getSQLiteFile()
#con <- dbConnect('SQLite','GEOmetadb.sqlite') 
#gpl=getBiocPlatformMap(con, bioc=annot)$gpl
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- add_rownames(tT, "ID")
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("Platform_SPOTID","ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Chromosome.annotation","GO.Function.ID"))
tT<-format(tT, digits=2, nsmall=2)
colnames(tT)=gsub(pattern = "\\.",replacement = "_",colnames(tT))
matrixtT=as.matrix(tT)
datajson=toJSON(matrixtT,pretty = TRUE)
#json=paste0("{\"data\":",datajson,"}")
#jsonfile=paste0(result.file.path,"/data.json")
#write(file = jsonfile, json)
#file.copy("./GEOAnalyse_tpl.html",result)

htmlfile=readChar(result.template, file.info(result.template)$size)
htmlfile=gsub(x=htmlfile,pattern = "###DATAJSON###",replacement = datajson, fixed = TRUE)
dir.create(result.path, showWarnings = TRUE, recursive = FALSE)

boxplot="boxplot.png"
png(boxplot,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(exprs(eset),las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOT###",replacement = boxplot, fixed = TRUE)
file.copy(boxplot,result.path)

histopvalue="histopvalue.png"

png(histopvalue,width=800,height = 400)
par(mfrow=c(1,2))
hist(fit2$F.p.value,nclass=100)
volcanoplot(fit2,coef=1,highlight=10)
htmlfile=gsub(x=htmlfile,pattern = "###HIST###",replacement = histopvalue, fixed = TRUE)
dev.off()
file.copy(histopvalue,result.path)
#file.conn=file(diag.html,open="w")

#writeLines(c("<html><body bgcolor='lightgray'>"),file.conn)


write.table(tolower(c(condition1Name,condition2Name)),quote = FALSE,col.names = FALSE, row.names=FALSE,file=result_export_conditions)
saveConditions=c(condition1Name,condition2Name)
save(eset,saveConditions,file=result_export_eset)
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