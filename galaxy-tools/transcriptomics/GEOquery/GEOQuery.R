
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
}	else
{
	conditions=pData(data1[[1]])["description"]
}
pData(eset)["source_name_ch1"]=conditions
write.table(unique(tolower(conditions[,1])),quote = FALSE,col.names = FALSE, row.names=FALSE,file=conditionFile)
save(eset,conditions,file=GEOQueryRData)