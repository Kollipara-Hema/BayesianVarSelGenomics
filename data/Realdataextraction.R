# Import data using the data hub
# Data Hub List
# All datasets are available at https://xenabrowser.net/datapages/.
set.seed(321)
#install.packages("UCSCXenaTools")
library(UCSCXenaTools)


data(XenaData)

head(XenaData)

XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "KIRC") %>%
  XenaFilter(filterDatasets = "HiSeqV2_PANCAN")  -> df_todo

df_todo

XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download

#Prepare data into R for analysis.
# Data downloaded on: 09/13/2023
cli = XenaPrepare(xe_download)
class(cli)

#write.csv(cli,"mRNA_KIRC_sep2023.csv")



#*****************************************************
#Download the miRNA data 
#*****************************************************


XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "KIRC") %>%
  XenaFilter(filterDatasets = "miRNA_HiSeq_gene")  -> df_mRNA
df_mRNA


XenaQuery(df_mRNA) %>%
  XenaDownload() -> xe_my_download
tcgaHub_mRNA = XenaPrepare(xe_my_download)
#write.csv(tcgaHub_mRNA,"miRNA_sep2023.csv")



#*****************************************************
#Download the protein data 
#*****************************************************


XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "KIRC") %>%
  XenaFilter(filterDatasets = "RPPA")  -> df_pRNA
df_pRNA


XenaQuery(df_pRNA) %>%
  XenaDownload() -> xe_my_download
tcgaHub_pRNA = XenaPrepare(xe_my_download)
#write.csv(tcgaHub_mRNA,"miRNA_sep2023.csv")

###############
# cleaning data

common_cols<-intersect(colnames(cli),colnames(tcgaHub_mRNA))
mrna<-cli[,common_cols]
mirna<-tcgaHub_mRNA[,common_cols]
mirna<-na.omit(mirna)

Y<-t(mrna[,-1])
colnames(Y)<-mrna$sample
X<-t(mirna[,-1])
colnames(X)<-mirna$sample
X<-as.data.frame(X)
Y<-as.data.frame(Y)


clinical<-readr::read_tsv("~/Desktop/BayesianVarSelGenomics/data/clinical.tsv")


X$type<-substr(rownames(X), nchar(rownames(X))-1, nchar(rownames(X)))
X$case_submitter_id<-substr(rownames(X), 1, nchar(rownames(X))-3 )
X$rownames<-rownames(X)
X1<-merge(X,clinical[,c("case_submitter_id","age_at_index","gender")],by="case_submitter_id")
X<-X1[!duplicated(X1$rownames),]
X$sample_type<-ifelse(X$type=="11","Solid Tissue normalwrite",ifelse(X$type=="01","Primary Tumor","New primar"))
rownames(X)<-X$rownames

#write.csv(Y,"Y_data.csv")
#write.csv(X,"X_data.csv")
Y1<- Y[,c("VHL","MET","XIST")]


merged_data <- merge(Y1, X, by = 0, all = TRUE)
#Analysis for MET#

#### tumor data###
columns_to_exclude <- c("rownames", "age_at_index", "gender", "case_submitter_id","type","sample_type")
tumor_D<- subset(merged_data[c(merged_data$type == "01"& merged_data$gender == "male"), ], select = -which(names(merged_data) %in% columns_to_exclude))#&merged_data$gender == "female"

#tumor_Y<-merged_data[which(merged_data$type=="01"),c("VHL")]
X=as.matrix(tumor_D[,-c(1,2,3,4)])
Y=as.matrix(tumor_D[,3])
Y_scale<-Y-mean(Y)
X_scale<-scale(X)

write.csv(merged_data,"regressiondata.csv")



common_cols1<-intersect(colnames(cli),colnames(tcgaHub_pRNA$RPPA.gz))
length(common_cols1)
length(unique(common_cols1))
prd<-t(tcgaHub_pRNA$RPPA.gz)
prd<-t(tcgaHub_pRNA$RPPA.gz[,common_cols1])
dim(prd)
mrna<-t(cli[,common_cols1])
dim(mrna)


prd<-t(tcgaHub_pRNA$RPPA.gz[,common_cols1[-1]])
mrna<-t(cli[,common_cols1[-1]])
colnames(prd)<-tcgaHub_pRNA$RPPA.gz$sample
colnames(mrna)<-cli$sample
newdata<-merge(prd[,c("VEGFR2-R-V","SETD2-R-E")],mrna,by=row.names(mrna))

Y<-prd[,c("VEGFR2-R-V","SETD2-R-E")]
Y$case<-rownames(prd)
Y$case<-c(rownames(Y))
X<-as.data.frame(mrna)
X$case<-rownames(mrna)
newdata<-merge(Y,X,by="case")
Y<-as.matrix(pd$VEGFR2.R.V)





