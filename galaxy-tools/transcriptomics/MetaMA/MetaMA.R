library(metaMA)
library(affy)
library(annaffy)
library(VennDiagram)
library(GEOquery)

#source("http://bioconductor.org/biocLite.R")
#library(DeSousa2013)


cargs<-commandArgs()
cargs<-cargs[(which(cargs=="--args")+1):length(cargs)]

nbargs=length(cargs)
rdataList=list()
conditionList=list()

for (i in seq(1,nbargs-5,2))
{
	Rdata=cargs[[i]]	
	condition=cargs[[i+1]]
	load(Rdata)
	rdataList=c(rdataList,(eset))
	conditionList=c(conditionList,condition)
}

#tables<-cargs[[1]]
#tech<-cargs[[2]]
result.html<-cargs[[nbargs-4]]
result.path<-cargs[[nbargs-3]]
#result.venn<-cargs[[nbargs-3]]
result.template<-cargs[[nbargs-2]]

#sink("/dev/null")
#dir.create(temp.files.path,recursive=TRUE)
#file.conn=file(diag.html,open="w")

#writeLines(c("<html><body bgcolor='lightgray'>"),file.conn)

print(result.path)

showVenn<-function(res,file)
{
	venn.plot<-venn.diagram(x = c(res[c(1:(length(res)-3))],meta=list(res$Meta)),
			filename = NULL, col = "black", 
			fill = c(1:(length(res)-2)),
			margin=0.05, alpha = 0.6)
	jpeg(file)
	grid.draw(venn.plot)
	dev.off()
}


library("org.Hs.eg.db")
x <- org.Hs.egUNIGENE
mapped_genes <- mappedkeys(x)
link <- as.list(x[mapped_genes])

probe2unigene<-function(expset){
	#construction of the map probe->unigene
	probes=rownames(exprs(expset))
	gene_id=fData(expset)[probes,"ENTREZ_GENE_ID"]
	unigene=link[gene_id]
	names(unigene)<-probes
	probe_unigene=unigene
}

unigene2probe<-function(map)
{
	suppressWarnings(x <- cbind(unlist(map), names(map)))
	unigene_probe=split(x[,2], x[,1])
}

convert2metaMA<-function(listStudies,mergemeth=mean)
{
	if (!(class(listStudies) %in% c("list"))) {
		stop("listStudies must be a list")
	}
	conv_unigene=lapply(listStudies,
			FUN=function(x) unigene2probe(probe2unigene(x)))
	id=lapply(conv_unigene,names)
	inter=Reduce(intersect,id)
	if(length(inter)<=0){stop("no common genes")}
	print(paste(length(inter),"genes in common"))
	esets=lapply(1:length(listStudies),FUN=function(i){
				l=lapply(conv_unigene[[i]][inter],
						FUN=function(x) exprs(listStudies[[i]])[x,,drop=TRUE])
				esetsgr=t(sapply(l,FUN=function(ll) if(is.null(dim(ll))){ll}
									else{apply(ll,2,mergemeth)}))
				esetsgr
			})
	return(list(esets=esets,conv.unigene=conv_unigene))
}

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


classes=list()
for (i in 1:length(conditionList))
{
	currentData=rdataList[[i]]
	currentCondition=conditionList[[i]]
		
	if (length(unique(tolower(pData(currentData)["source_name_ch1"][,1])))>1)
	{
		currentClasses=as.numeric(tolower(as.character(pData(currentData)["source_name_ch1"][,1]))==currentCondition)
		
	}	else
	{
		currentClasses=as.numeric(tolower(as.character(pData(currentData)["description"][,1]))==currentCondition)
		
	}
	classes=c(classes,list(currentClasses))
	#write(file="~/galaxy-modal/classes.txt",classes)
}
rdatalist <- lapply(rdataList, FUN=function(datalist) normalization(datalist))

conv=convert2metaMA(rdataList)
esets=conv$esets
conv_unigene=conv$conv.unigene

#write(file="~/galaxy-modal/esets.txt",length(esets))
#write(file="~/galaxy-modal/classes.txt",length(classes))
res=pvalcombination(esets=esets,classes=classes)
resIDDIRR=IDDIRR(res$Meta,res$AllIndStudies)
length(res$Meta)
Hs.Meta=rownames(esets[[1]])[res$Meta]
origId.Meta=lapply(conv_unigene,FUN=function(vec) as.vector(unlist(vec[Hs.Meta])))

gpllist <- lapply(rdataList, FUN=function(ann) annotation(ann))
platflist <- lapply(gpllist, FUN=function(gpl) getGEO(gpl, AnnotGPL=TRUE))
ncbifdlist <- lapply(platflist, FUN=function(data) data.frame(attr(dataTable(data), "table")))
ncbifdresult=lapply(1:length(origId.Meta), FUN=function(i) ncbifdlist[[i]][which(ncbifdlist[[i]]$ID %in% origId.Meta[[i]]),])
ncbidfannot=do.call(rbind,ncbifdresult)
ncbidfannot <- subset(ncbidfannot, select=c("Platform_SPOTID","ID","Gene.title","Gene.symbol","GenBank.Accession","Chromosome.location","Chromosome.annotation","GO.Function.ID"))

library(jsonlite)
matrixncbidfannot=as.matrix(ncbidfannot)
datajson=toJSON(matrixncbidfannot,pretty = TRUE)
summaryjson=toJSON(as.matrix(t(resIDDIRR)),pretty = TRUE)


#vennsplit=strsplit(result.venn,split="/")[[1]]
#venn=paste0("./",vennsplit[length(vennsplit)])


vennFilename="venn.png"
vennFile=file.path(result.path,vennFilename)
htmlfile=readChar(result.template, file.info(result.template)$size)
htmlfile=gsub(x=htmlfile,pattern = "###DATAJSON###",replacement = datajson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###SUMMARYJSON###",replacement = summaryjson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###VENN###",replacement = vennFilename, fixed = TRUE)
write(htmlfile,result.html)

library(VennDiagram)
flog.threshold(ERROR)

#venn.plot<-venn.diagram(x = c(res[c(1:(length(res)-3))],meta=list(res$Meta)),filename = v, col = "black", fill = c(1:(length(res)-2)), margin=0.05, alpha = 0.6,imagetype = "png")
dir.create(result.path, showWarnings = TRUE, recursive = FALSE)

showVenn<-function(liste,file)
{
	venn.plot<-venn.diagram(x = liste,
			filename = vennFilename, col = "black",
			fill = 1:length(liste)+1,
			margin=0.05, alpha = 0.6,imagetype = "png")
#	png(file);
#	grid.draw(venn.plot);
#	dev.off();
	
}

l=list()
for(i in 1:length(esets))
{
	l[[paste("study",i,sep="")]]<-res[[i]]
}
l[["Meta"]]=res[[length(res)-1]]
showVenn(l,vennFile)
file.copy(vennFilename,result.path)
#l=list()
#for(i in 1:length(esets))
#{
#	l[[paste("study",i,sep="")]]<-res[[i]]
#}
#l[["Meta"]]=res[[length(res)-1]]
#showVenn(res,result.venn)
#writeLines(c("<h2>Venn diagram</h2>"),file.conn)
#writeLines(c("<img src='venn.png'><br/><br/>"),file.conn)
#writeLines(c("</body></html>"),file.conn)
#close(file.conn)