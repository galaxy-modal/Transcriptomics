
#---- lib ---------------
library(GEOquery)


cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]

GEOQueryID<-cargs[[1]]
GEOQueryData<-cargs[[2]]
GEOQueryRData<-cargs[[3]]
conditionFile<-cargs[[4]]
transformation<-cargs[[5]]

data1 = getGEO(GEOQueryID)
eset=data1[[1]]

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

matrixData=exprs(eset)
write.table(matrixData,col.names=NA,row.names=TRUE,sep="\t",file=GEOQueryData)
if (length(unique(tolower(pData(data1[[1]])["source_name_ch1"][,1])))>1)
{
	conditions=pData(data1[[1]])["source_name_ch1"]
	description=paste0(as.vector(pData(data1[[1]])["geo_accession"][,1]), " ",as.vector(pData(data1[[1]])["title"][,1]), " ", as.vector(conditions[,1]))
	#file=rownames(pData(data1[[1]]))
}	else
{
	conditions=pData(data1[[1]])["description"]
	description=paste0(as.vector(pData(data1[[1]])["geo_accession"][,1]), " ",as.vector(pData(data1[[1]])["title"][,1]), " ", as.vector(conditions[,1]))
	#file=rownames(pData(data1[[1]]))
}
conditions[,1]=tolower(conditions[,1])
pData(eset)["source_name_ch1"]=conditions

#write.table(unique(tolower(conditions[,1])),quote = FALSE,col.names = FALSE, row.names=FALSE,file=conditionFile)
write.table(cbind(conditions,description),quote = FALSE,col.names = FALSE, row.names=TRUE,file=conditionFile,sep="\t")
save(eset,conditions,file=GEOQueryRData)