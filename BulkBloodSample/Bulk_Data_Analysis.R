#############################################################
####                   Load libraries                    ####
#############################################################
library(dplyr)
library(DESeq2)
library(tibble)
library(AnnotationHub)
library(ensembldb)
library(purrr)

#############################################################
####                     Load data                       ####
#############################################################
setwd('./data')
count_data<-read.csv('GSE161731_counts.csv',header=T,row.names = 1)
col_data<-read.csv('GSE161731_counts_key.csv',header=T,row.names = 1)
colnames(count_data)<-gsub('\\.','-',colnames(count_data))
colnames(count_data)<-gsub('X','',colnames(count_data))
(row.names(col_data)%in%colnames(count_data))

#############################################################
####                Use only COVID data                  ####
#############################################################
## coldata
col_data_covid<-col_data %>% dplyr::filter(cohort=='COVID-19')
col_data_covid$gender<-factor(col_data_covid$gender)
col_data_covid$race<-factor(col_data_covid$race)
col_data_covid$time_since_onset<-factor(col_data_covid$time_since_onset,levels=c('early','middle','late'))
col_data_covid$batch<-factor(col_data_covid$batch)
col_data_covid<-col_data_covid[order(col_data_covid$time_since_onset),]
## count data
count_data_covid<-count_data[,row.names(col_data_covid)]
all(rownames(col_data_covid) %in% colnames(count_data))
## Taking the average of same id and time
mean_round<-function(data){
  data<-mean(data)
  data<-round(data,0)
  data
}
col_data_covid$id_time<-paste0(col_data_covid$subject_id,'_',col_data_covid$time_since_onset)
count_data_covid_unique<-NULL
for(i in unique(col_data_covid$id_time)){
  id_time_index<-which(col_data_covid$id_time==i)
  sample_name<-row.names(col_data_covid)[id_time_index]
  if(length(sample_name)>1){
    id_time_count<-count_data_covid[,sample_name]
    new_col<-apply(id_time_count,1,mean_round)
  }
  else{
    new_col<-count_data_covid[,sample_name]
  }
  count_data_covid_unique<-cbind(count_data_covid_unique,new_col)
}
colnames(count_data_covid_unique)<-unique(col_data_covid$id_time)

col_data_covid_unique<-col_data_covid[!duplicated(col_data_covid[,'id_time']),]
row.names(col_data_covid_unique)<-col_data_covid_unique$id_time
col_data_covid_unique<-as.data.frame(col_data_covid_unique)
count_data_covid_unique<-as.data.frame(count_data_covid_unique)
## filter the low-expressed gene
gene_expressed_index<-which(apply(count_data_covid_unique,1,sum)>=100)
## index for pairwise comparison
index_early_middle <- which(col_data_covid_unique$time_since_onset %in% c("early", "middle")) 
index_early_late <- which(col_data_covid_unique$time_since_onset %in% c("early", "late")) 
index_middle_late <- which(col_data_covid_unique$time_since_onset %in% c("middle", "late")) 

#############################################################
####                    DESeq matrix                     ####
#############################################################
dds_time<-DESeqDataSetFromMatrix(countData = count_data_covid_unique[gene_expressed_index, index_early_late],
                                 colData = col_data_covid_unique[index_early_late, ],
                                 design = ~ batch + time_since_onset)
dds_time_lrt<-DESeq(dds_time, test="LRT", reduced = ~ batch)
# Extract results
res_time_lrt <- results(dds_time_lrt)
# Create a tibble for LRT results
res_time_lrt_tb <- res_time_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
# Subset to return genes with padj < 0.05
sigLRT_genes <- res_time_lrt_tb %>% 
  dplyr::filter(padj < 0.05)

#############################################################
####                  gene annotation                    ####
#############################################################
# Connect to AnnotationHub
ah <- AnnotationHub()
# Query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
# Extract annotations of interest
human_ens <- human_ens[["AH89426"]]
# Create a gene-level dataframe
annotations <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, symbol,entrezid)
annotations$entrezid <- map(annotations$entrezid,1) %>%  unlist()
## Merge the AnnotationHub dataframe with the results 
res_ids <- inner_join(sigLRT_genes, annotations, by=c("gene"="gene_id")) 





