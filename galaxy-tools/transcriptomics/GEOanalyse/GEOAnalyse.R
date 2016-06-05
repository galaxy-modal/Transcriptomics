library(Biobase)
library(GEOquery)
library(limma)
library(jsonlite)

cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]

nbargs=length(cargs)
Rdata=cargs[[1]]
condition1=cargs[[2]]
condition2=cargs[[3]]
nbresult=cargs[[4]]
transformation=cargs[[5]]
result<-cargs[[6]]
result.path<-cargs[[7]]
export_rdata=cargs[[8]]
result.template=cargs[[9]]


load(Rdata)
gset=eset

normalization<-function(data){
	ex <- exprs(data)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	LogC <- (qx[5] > 100) ||
			(qx[6]-qx[1] > 50 && qx[2] > 0) ||
			(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	if (LogC) { ex[which(ex <= 0)] <- NaN
		return (log2(ex)) } else {
		return (ex)
	}
}

if (transformation=="auto"){
	exprs(eset)=normalization(eset)
} else if (transformation=="yes"){
	exprs(eset)=log2(exprs(eset))			
}

selected=c(which((tolower(as.character(pData(eset)["source_name_ch1"][,1]))==condition1)), 
		which(tolower(as.character(pData(eset)["source_name_ch1"][,1]))==condition2))

eset=eset[,selected]
sml=paste0("G",as.numeric(tolower(as.character(pData(eset)["source_name_ch1"][,1]))!=condition1))

fl <- as.factor(sml)
eset$description <- fl

design <- model.matrix(~ description + 0, eset)

colnames(design) <- levels(fl)
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nbresult)
gpl <- annotation(eset)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  

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
boxplot(exprs(gset),las=2,outline=FALSE,main="Raw data")
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###RAWBOXPLOT###",replacement = boxplot, fixed = TRUE)
file.copy(boxplot,result.path)

if (!identical(exprs(eset),exprs(gset))) {
	boxplotlog2="boxplotlog2.png"
	png(boxplotlog2,width=800,height = 400)
	par(mar=c(7,5,1,1))
	boxplot(exprs(eset),las=2,outline=FALSE,main="Log2 transform data")
	dev.off()
	htmlfile=gsub(x=htmlfile,pattern = "###LOG2BOXPLOT###",replacement = paste0("<img src='",boxplotlog2,"'>"), fixed = TRUE)
	file.copy(boxplotlog2,result.path) 
} else{
	htmlfile=gsub(x=htmlfile,pattern = "###LOG2BOXPLOT###",replacement = "", fixed = TRUE)
}

histopvalue="histopvalue.png"

png(histopvalue,width=800,height = 400)
par(mfrow=c(1,2))
hist(fit2$F.p.value,nclass=100)
volcanoplot(fit2,coef=1,highlight=10)
htmlfile=gsub(x=htmlfile,pattern = "###HIST###",replacement = histopvalue, fixed = TRUE)
dev.off()
file.copy(histopvalue,result.path)

save(eset,file=export_rdata)

write(htmlfile,result)

