library(Biobase)
library(GEOquery)
library(limma)
library(jsonlite)

#source("http://bioconductor.org/biocLite.R")
#library(DeSousa2013)


cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]

nbargs=length(cargs)
Rdata=cargs[[1]]
condition=cargs[[2]]

#tables<-cargs[[1]]
#tech<-cargs[[2]]
result<-cargs[[3]]
result.path<-cargs[[4]]
result.template=cargs[[5]]

nbresult=1000
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
exprs(gset)=normalization(gset)
#sink("/dev/null")
#dir.create(result.file.path,recursive=TRUE)
sml=paste0("G",as.numeric(tolower(as.character(pData(gset)["source_name_ch1"][,1]))!=condition))
fl <- as.factor(sml)
gset$description <- fl

design <- model.matrix(~ description + 0, gset)

colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nbresult)
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
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
boxplot(exprs(gset),las=2,outline=FALSE)
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