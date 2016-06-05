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

load(cargs[[nbargs-12]])
targetFile=cargs[[nbargs-11]]
condition1Name=cargs[[nbargs-10]]
condition1=cargs[[nbargs-9]]
condition2Name=cargs[[nbargs-8]]
condition2=cargs[[nbargs-7]]
nbresult=cargs[[nbargs-6]]
result_export_eset=cargs[[nbargs-5]]
result=cargs[[nbargs-4]]
result.path=cargs[[nbargs-3]]
result.template=cargs[[nbargs-2]]

file.copy(targetFile,"./targetFile.txt")

condition1_tmp <- strsplit(condition1,",")
condition1 <-unlist(condition1_tmp)


condition2_tmp <- strsplit(condition2,",")
condition2<-unlist(condition2_tmp)


conditions=c(condition1,condition2)

#nbresult=1000
dir.create(result.path, showWarnings = TRUE, recursive = FALSE)

targets <- readTargets("targetFile.txt",path=".",sep="\t")
print("passe")
eset=eset[,which(rownames(eset@phenoData@data) %in% conditions)]
eset@phenoData@data$source_name_ch1=""
eset@phenoData@data$source_name_ch1[which(rownames(eset@phenoData@data) %in% condition1)]=condition1Name
eset@phenoData@data$source_name_ch1[which(rownames(eset@phenoData@data) %in% condition2)]=condition2Name

sml=paste0("G",as.numeric(as.character(pData(eset)["source_name_ch1"][,1])!=condition1Name))
fl <- as.factor(sml)
eset$description <- fl
design <- model.matrix(~ description + 0, eset)

colnames(design) <- levels(fl)
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nbresult)

annot <- annotation(eset)
mapping=read.csv("/galaxy-tools/transcriptomics/db/gplToBioc.csv",stringsAsFactors=FALSE)
gpl=mapping[which(mapping$bioc_package==annotation(eset)),]$gpl
gpl=gpl[1]
annotation(eset)=gpl

platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

fData(eset)["ID"]=row.names(fData(eset))
fData(eset)=merge(x=fData(eset),y=ncbifd,all.x = TRUE, by = "ID")
colnames(fData(eset))[4]="ENTREZ_GENE_ID"
row.names(fData(eset))=fData(eset)[,"ID"]

# replace original platform annotation
tT <- add_rownames(tT, "ID")
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("Platform_SPOTID","ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title","Gene.ID","Chromosome.annotation","GO.Function.ID"))
tT<-format(tT, digits=2, nsmall=2)
colnames(tT)=gsub(pattern = "\\.",replacement = "_",colnames(tT))
matrixtT=as.matrix(tT)
datajson=toJSON(matrixtT,pretty = TRUE)

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



write.table(tolower(c(condition1Name,condition2Name)),quote = FALSE,col.names = FALSE, row.names=FALSE,file=result_export_conditions)
saveConditions=c(condition1Name,condition2Name)
save(eset,saveConditions,file=result_export_eset)
write(htmlfile,result)

