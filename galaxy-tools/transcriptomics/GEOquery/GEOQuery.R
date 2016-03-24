
#---- lib ---------------
library(GEOquery)
#source("http://bioconductor.org/biocLite.R")
#library(DeSousa2013)

cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]

GEOQueryID<-cargs[[1]]
GEOQueryData<-cargs[[2]]
GEOQueryRData<-cargs[[3]]
conditionFile<-cargs[[4]]

data1 = getGEO(GEOQueryID)
eset=data1[[1]]
#names(GPLList(data1))
#gds <- getGEO(GEOQueryID)
#gpl_name=annotation(data1[[1]])
#gpl=getGEO(gpl_name)

write.table(exprs(data1[[1]]),file=GEOQueryData)
if (length(unique(tolower(pData(data1[[1]])["source_name_ch1"][,1])))>1)
{
	conditions=pData(data1[[1]])["source_name_ch1"]
	#conditionValue= gsub(" ", "_", lapply(pData(data1[[1]])["source_name_ch1"],as.character)[[1]])
}	else
{
	conditions=pData(data1[[1]])["description"]
	#conditionValue= gsub(" ", "_", lapply(pData(data1[[1]])["description"],as.character)[[1]])
}

write.table(unique(tolower(conditions[,1])),quote = FALSE,col.names = FALSE, row.names=FALSE,file=conditionFile)
save(eset,conditions,file=GEOQueryRData)