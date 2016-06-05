
#---- lib ---------------
library(GEOquery)


cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]

GEOQueryID<-cargs[[1]]
GEOQueryData<-cargs[[2]]
GEOQueryRData<-cargs[[3]]
conditionFile<-cargs[[4]]

data1 = getGEO(GEOQueryID)
eset=data1[[1]]

write.table(exprs(data1[[1]]),file=GEOQueryData)
if (length(unique(tolower(pData(data1[[1]])["source_name_ch1"][,1])))>1)
{
	conditions=pData(data1[[1]])["source_name_ch1"]
	description=paste0(as.vector(pData(data1[[1]])["title"][,1]), " ", as.vector(conditions[,1]))
	#file=rownames(pData(data1[[1]]))
}	else
{
	conditions=pData(data1[[1]])["description"]
	description=paste0(as.vector(pData(data1[[1]])["title"][,1]), " ", as.vector(conditions[,1]))
	#file=rownames(pData(data1[[1]]))
}
conditions[,1]=tolower(conditions[,1])
pData(eset)["source_name_ch1"]=conditions

#write.table(unique(tolower(conditions[,1])),quote = FALSE,col.names = FALSE, row.names=FALSE,file=conditionFile)
write.table(cbind(conditions,description),quote = FALSE,col.names = FALSE, row.names=TRUE,file=conditionFile,sep="\t")
save(eset,conditions,file=GEOQueryRData)