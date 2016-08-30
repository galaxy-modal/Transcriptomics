FROM ubuntu:14.04

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list


RUN apt-get update && apt-get install -y \
    r-base \
   r-base-dev 

RUN apt-get install -y r-base-core libxml2-dev libcurl4-openssl-dev curl

RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" \
			-e "biocLite('GEOquery')"  \
			-e "biocLite('limma')" \
			-e "biocLite('affy')" \
			-e "biocLite('annaffy')" \
			-e "biocLite('org.Hs.eg.db')" \
			-e "biocLite('DESeq')" \
			-e "biocLite('HTSFilter')" \
			-e "install.packages('metaMA', repos='http://cran.us.r-project.org')" \
			-e "install.packages('metaRNASeq', repos='http://cran.us.r-project.org')" \
			-e "install.packages('jsonlite', repos='http://cran.us.r-project.org')" \
			-e "install.packages('VennDiagram', repos='http://cran.us.r-project.org')"

RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" \
			-e "biocLite('GEOmetadb')"
			
RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" \
			-e "biocLite('affyPLM')"

RUN Rscript -e "install.packages('dplyr', repos='http://cran.us.r-project.org')" \

RUN mkdir -p /galaxy-tools/transcriptomics/GEOquery
RUN mkdir -p /galaxy-tools/transcriptomics/GEOanalyse
RUN mkdir -p /galaxy-tools/transcriptomics/Analyse
RUN mkdir -p /galaxy-tools/transcriptomics/CustomAnalyse
RUN mkdir -p /galaxy-tools/transcriptomics/MetaMA
RUN mkdir -p /galaxy-tools/transcriptomics/MetaRNASeq
RUN mkdir -p /galaxy-tools/transcriptomics/AffyQCnormalization
RUN mkdir -p /galaxy-tools/transcriptomics/db

COPY transcriptomics/db/gplToBioc.csv /galaxy-tools/transcriptomics/db/
COPY transcriptomics/GEOquery/GEOQuery.R /galaxy-tools/transcriptomics/GEOquery/
#COPY transcriptomics/GEOanalyse/GEOAnalyse.R /galaxy-tools/transcriptomics/GEOanalyse/
#COPY transcriptomics/GEOanalyse/GEOAnalyse_tpl.html /galaxy-tools/transcriptomics/GEOanalyse/
COPY transcriptomics/GEOanalyse/ /galaxy-tools/transcriptomics/GEOanalyse/
#COPY transcriptomics/GEOanalyse/ /galaxy-tools/transcriptomics/GEOanalyse/
COPY transcriptomics/Analyse/ /galaxy-tools/transcriptomics/Analyse/
COPY transcriptomics/MetaMA/MetaMA.R /galaxy-tools/transcriptomics/MetaMA/
COPY transcriptomics/MetaMA/MetaMa_tpl.html /galaxy-tools/transcriptomics/MetaMA/
#COPY transcriptomics/MetaRNASeq/MetaRNASeq.R /galaxy-tools/transcriptomics/MetaRNASeq/
COPY transcriptomics/MetaRNASeq/MetaRNASeq2.R /galaxy-tools/transcriptomics/MetaRNASeq/
COPY transcriptomics/MetaRNASeq/MetaRNASeq_tpl.html /galaxy-tools/transcriptomics/MetaRNASeq/
COPY transcriptomics/CustomAnalyse/CustomAnalyse_tpl.html /galaxy-tools/transcriptomics/CustomAnalyse/
COPY transcriptomics/CustomAnalyse/CustomAnalyse.R /galaxy-tools/transcriptomics/CustomAnalyse/
COPY transcriptomics/AffyQCnormalization/AffyQCnormalization.R /galaxy-tools/transcriptomics/AffyQCnormalization/
COPY transcriptomics/AffyQCnormalization/AffyQCnormalization_tpl.html /galaxy-tools/transcriptomics/AffyQCnormalization/
COPY transcriptomics/ImportDataFromMatrix/ImportDataFromMatrix.R /galaxy-tools/transcriptomics/ImportDataFromMatrix/
COPY transcriptomics/ImportDataFromMatrix/ImportDataFromMatrix_tpl.html /galaxy-tools/transcriptomics/ImportDataFromMatrix/
COPY stderr_wrapper.py /galaxy-tools/

RUN chmod a+rw /usr/local/lib/R/site-library
ENV PATH $PATH:/GEOquery:/GEOanalyse:/MetaMA


