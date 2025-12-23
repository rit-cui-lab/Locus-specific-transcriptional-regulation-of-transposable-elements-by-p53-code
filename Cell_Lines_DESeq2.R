#DEseq2 differential anyalisis different fold increases were also compared. We would eventually decide to use fold increase of 0.


library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90") #Adjustments to folders can be made.

filter<-read.table("IMR90_hr12_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90/DESeq2")


#control
#Comparison by Control

samples<-read.table("hr12_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","hr12_5_FU","control")) #comparison between treatment (hr12_5_FU) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]

#DMSO
#Comparison by DMSO

res_DMSO<-results(dds, contrast=c("condition","hr12_5_FU","hr12_DMSO")) #comparison between treatment (hr12_5_FU) and control (hr12_DMSO)
res_DMSO<-as.data.frame(res_DMSO)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_DMSO<-resSig_DMSO[,c("RE","con_1","con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","baseMean",
                            "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "con_1.x", "con_2.x", "hr12_DMSO_1.x", "hr12_DMSO_2.x","hr12_5_FU_1.x","hr12_5_FU_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "con_1", "con_2", "hr12_DMSO_1", "hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

write.table(resSig, "IMR90_hr12_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


F_vals<-c()
vals<-c()
i<-0
resSig_con <- subset(res_con, padj < 0.05)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig<-merge(resSig_con,resSig_DMSO, by="row.names")
resSig<-resSig[,c("Row.names", "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y")]

names(resSig)<-c("RE","con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj"
                 ,"DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")
while (i<4.5)
{
  temp<-resSig[abs(resSig$con_log2FoldChange) >= i,]
  temp<-temp[abs(temp$DMSO_log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#HCT116
setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

filter<-read.table("HCT116_hr12_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116/DESeq2")


#control
#Comparison by Control
samples<-read.table("hr12_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","hr12_5_FU","control")) #comparison between treatment (hr12_5_FU) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]

#DMSO

res_DMSO<-results(dds, contrast=c("condition","hr12_5_FU","hr12_DMSO")) #comparison between treatment (hr12_5_FU) and control (hr12_DMSO)
res_DMSO<-as.data.frame(res_DMSO)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_DMSO<-resSig_DMSO[,c("RE","con_1","con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","baseMean",
                            "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "con_1.x", "con_2.x", "hr12_DMSO_1.x", "hr12_DMSO_2.x","hr12_5_FU_1.x","hr12_5_FU_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "con_1", "con_2", "hr12_DMSO_1", "hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

write.table(resSig, "HCT116_hr12_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)



F_vals<-c()
vals<-c()
i<-0
resSig_con <- subset(res_con, padj < 0.05)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig<-merge(resSig_con,resSig_DMSO, by="row.names")
resSig<-resSig[,c("Row.names", "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y")]

names(resSig)<-c("RE","con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj"
                 ,"DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")
while (i<4.5)
{
  temp<-resSig[abs(resSig$con_log2FoldChange) >= i,]
  temp<-temp[abs(temp$DMSO_log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#hr6 DESeq2 RE
#DESeq2

library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

filter<-read.table("IMR90_hr6_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90/DESeq2")


#control
#Comparison by Control
samples<-read.table("hr6_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","hr6_5_FU","control")) #comparison between treatment (hr6_5_FU) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]

#DMSO
#Comparison by DMSO

res_DMSO<-results(dds, contrast=c("condition","hr6_5_FU","hr6_DMSO")) #comparison between treatment (hr6_5_FU) and control (hr6_DMSO)
res_DMSO<-as.data.frame(res_DMSO)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_DMSO<-resSig_DMSO[,c("RE","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2","baseMean",
                            "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "con_1.x", "con_2.x", "hr6_DMSO_1.x", "hr6_DMSO_2.x","hr6_5_FU_1.x","hr6_5_FU_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "con_1", "con_2", "hr6_DMSO_1", "hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

write.table(resSig, "IMR90_hr6_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


F_vals<-c()
vals<-c()
i<-0
resSig_con <- subset(res_con, padj < 0.05)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig<-merge(resSig_con,resSig_DMSO, by="row.names")
resSig<-resSig[,c("Row.names", "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y")]

names(resSig)<-c("RE","con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj"
                 ,"DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")
while (i<4.5)
{
  temp<-resSig[abs(resSig$con_log2FoldChange) >= i,]
  temp<-temp[abs(temp$DMSO_log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#HCT116
setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

filter<-read.table("HCT116_hr6_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116/DESeq2")


#control
#Comparison by Control
samples<-read.table("hr6_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","hr6_5_FU","control")) #comparison between treatment (hr6_5_FU) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]

#DMSO
#Comparison by DMSO

res_DMSO<-results(dds, contrast=c("condition","hr6_5_FU","hr6_DMSO")) #comparison between treatment (hr6_5_FU) and control (hr6_DMSO)
res_DMSO<-as.data.frame(res_DMSO)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_DMSO<-resSig_DMSO[,c("RE","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2","baseMean",
                            "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "con_1.x", "con_2.x", "hr6_DMSO_1.x", "hr6_DMSO_2.x","hr6_5_FU_1.x","hr6_5_FU_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "con_1", "con_2", "hr6_DMSO_1", "hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

write.table(resSig, "HCT116_hr6_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


F_vals<-c()
vals<-c()
i<-0
resSig_con <- subset(res_con, padj < 0.05)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig<-merge(resSig_con,resSig_DMSO, by="row.names")
resSig<-resSig[,c("Row.names", "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y")]

names(resSig)<-c("RE","con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj"
                 ,"DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")
while (i<4.5)
{
  temp<-resSig[abs(resSig$con_log2FoldChange) >= i,]
  temp<-temp[abs(temp$DMSO_log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#updated Nutlin IMr90 and HCT116
library("DESeq2")
library(VennDiagram)
#edit sample for DESeq2
#HCT116

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

filter<-read.table("HCT116_total_Nutlin_my_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116/DESeq2")

#control and DMSO Andrew

#control
#Comparison by Control
samples<-read.table("table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "HCT116_Nutlin_vs_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

F_vals<-c()
vals<-c()
i<-0
name<-"HCT116_Nutlin_vs_DMSO_padj_0.05_FC"

resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)

resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

library("DESeq2")
library(VennDiagram)
#edit sample for DESeq2
#IMR90

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

filter<-read.table("IMR90_total_Nutlin_my_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90/DESeq2")



#control
#Comparison by Control
samples<-read.table("total_control_sample_table.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","Nutlin","control")) #comparison between treatment (Nutlin) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]

#DMSO
#Comparison by DMSO

res_DMSO<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res_DMSO<-as.data.frame(res_DMSO)

resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "con_1.x", "con_2.x", "DMSO_1.x", "DMSO_2.x","Nutlin_1.x","Nutlin_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "con_1", "con_2", "DMSO_1", "DMSO_2","Nutlin_1","Nutlin_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

write.table(resSig, "IMR90_Nutlin_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"IMR90_nutlin_overlap_padj_0.05_FC"

resSig_con <- subset(res_con, padj < 0.05)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig<-merge(resSig_con,resSig_DMSO, by="row.names")
resSig<-resSig[,c("Row.names", "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y")]

names(resSig)<-c("RE","con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj"
                 ,"DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")


resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[c("RE","con_1", "con_2", "DMSO_1", "DMSO_2","Nutlin_1","Nutlin_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")]
while (i<4.5)
{
  temp<-resSig[abs(resSig$con_log2FoldChange) >= i,]
  temp<-temp[abs(temp$DMSO_log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines")

#Tables of combined telescope samples that have been filtered
Saos<-read.table("Saos-2_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
GM06170<-read.table("GM06170_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
GM00011<-read.table("GM00011_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF7_GSE47042<-read.table("MCF7_GSE47042_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
A549<-read.table("A549_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
Fibro<-read.table("Fibro_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF10A<-read.table("MCF10A_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF7_GSE86221<-read.table("MCF7_GSE86221_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
SJSA<-read.table("SJSA_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines/DESeq2")

#Saos-2

samples<-read.table("Saos-2_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =Saos , colData = samples, design = ~ condition) #Saves the data, sample table, and design of the DESeq2 comparison we used. condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","doxycyclin","control")) #comparison between treatment (doxycyclin) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,Saos, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","dox_1","dox_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "Saos-2_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"Saos-2_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,Saos, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","dox_1","dox_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#GM06170

samples<-read.table("GM06170_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =GM06170 , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Doxorubicin","control")) #comparison between treatment (Doxorubicin) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,GM06170, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","dox_1","dox_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "GM06170_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"GM06170_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,GM06170, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","dox_1","dox_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#GM00011

samples<-read.table("GM00011_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =GM00011 , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Doxorubicin","control")) #comparison between treatment (Doxorubicin) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,GM00011, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","dox_1","dox_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "GM00011_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"GM00011_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,GM00011, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","dox_1","dox_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#MCF7_GSE47042

samples<-read.table("MCF7_GSE47042_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )#Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =MCF7_GSE47042 , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","control")) #comparison between treatment (Nutlin) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,MCF7_GSE47042, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","Nutlin_1","Nutlin_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "MCF7_GSE47042_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"MCF7_GSE47042_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,MCF7_GSE47042, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","Nutlin_1","Nutlin_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#A549

samples<-read.table("A549_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =A549 , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","control")) #comparison between treatment (Nutlin) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,A549, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","con_3","Nutlin_1","Nutlin_2","Nutlin_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "A549_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"A549_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,A549, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","con_3","Nutlin_1","Nutlin_2","Nutlin_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#Fibro

samples<-read.table("Fibro_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =Fibro , colData = samples, design = ~ condition) #Saves the data, sample table, and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,Fibro, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nutlin_1","Nutlin_2","Nutlin_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "Fibro_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"Fibro_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,Fibro, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nutlin_1","Nutlin_2","Nutlin_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#MCF10A

samples<-read.table("MCF10A_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =MCF10A , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,MCF10A, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nutlin_1","Nutlin_2","Nutlin_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "MCF10A_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"MCF10A_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,MCF10A, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nutlin_1","Nutlin_2","Nutlin_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#MCF7_GSE86221

samples<-read.table("MCF7_GSE86221_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =MCF7_GSE86221 , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,MCF7_GSE86221, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "MCF7_GSE86221_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"MCF7_GSE86221_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,MCF7_GSE86221, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)


#SJSA

samples<-read.table("SJSA_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specificic version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =SJSA , colData = samples, design = ~ condition) #Saves the data, sample table, and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,SJSA, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "SJSA_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"SJSA_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,SJSA, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

filter<-read.table("HCT116_total_Nutlin_my_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116/DESeq2")

#control and DMSO Andrew

#control
#Comparison by Control
samples<-read.table("table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "HCT116_Nutlin_vs_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)



vals<-c()
F_vals<-c()
i<-0
name<-"HCT116_Nutlin_vs_DMSO_padj_0.05_FC"

resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)

resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

library("DESeq2")
library(VennDiagram)

#IMR90

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

filter<-read.table("IMR90_total_Nutlin_my_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90/DESeq2")

#control and DMSO Andrew

#control
#Comparison by Control
samples<-read.table("total_control_sample_table.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","Nutlin","control")) #comparison between treatment (Nutlin) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]

#DMSO
#Comparison by DMSO
res_DMSO<-results(dds, contrast=c("condition","Nutlin","DMSO"))
res_DMSO<-as.data.frame(res_DMSO)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","con_1","con_2","DMSO_1","DMSO_2","Nutlin_1","Nutlin_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "con_1.x", "con_2.x", "DMSO_1.x", "DMSO_2.x","Nutlin_1.x","Nutlin_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "con_1", "con_2", "DMSO_1", "DMSO_2","Nutlin_1","Nutlin_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

write.table(resSig, "IMR90_Nutlin_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"IMR90_nutlin_overlap_padj_0.05_FC"

resSig_con <- subset(res_con, padj < 0.05)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig<-merge(resSig_con,resSig_DMSO, by="row.names")
resSig<-resSig[,c("Row.names", "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y")]

names(resSig)<-c("RE","con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj"
                 ,"DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")


resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[c("RE","con_1", "con_2", "DMSO_1", "DMSO_2","Nutlin_1","Nutlin_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")]
while (i<4.5)
{
  temp<-resSig[abs(resSig$con_log2FoldChange) >= i,]
  temp<-temp[abs(temp$DMSO_log2FoldChange) >= i,]
  F_vals<-c(F_vals,nrow(temp))
  
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)


library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines")

#MV4_11

filter<-read.table("MV4_11_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/MV4_11_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/MV4_11_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/MV4_11_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#PAEC

filter<-read.table("PAEC_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/PAEC_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/PAEC_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/PAEC_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#PAEC

filter<-read.table("PAEC_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/PAEC_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specificic version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/PAEC_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/PAEC_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#IMR90_GSE139003

filter<-read.table("IMR90_GSE139003_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/IMR90_GSE139003_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
#FDR<-p.adjust(res_con$pvalue, method = "BH")
#res_con<-cbind(res_con,FDR)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1", "DMSO_2", "hr12_Nut_1", "hr12_Nut_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/IMR90_GSE139003_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/IMR90_GSE139003_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1", "DMSO_2", "hr12_Nut_1", "hr12_Nut_2","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#MCF10A_hr4

filter<-read.table("MCF10A_hr4_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/MCF10A_hr4_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/MCF10A_hr4_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/MCF10A_hr4_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)


#hNCCs

filter<-read.table("hNCCs_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/hNCCs_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/hNCCs_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/hNCCs_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#hSMCs
filter<-read.table("hSMCs_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

samples<-read.table("DESeq2/hSMCs_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Nutlin","DMSO")) #comparison between treatment (Nutlin) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/hSMCs_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)

vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/hSMCs_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","DMSO_3","Nut_1","Nut_2", "Nut_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#A431

filter<-read.table("A431_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered

samples<-read.table("DESeq2/A431_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res<-results(dds, contrast=c("condition","Dox","control")) #comparison between treatment (Dox) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","con_3","Dox_1","Dox_2", "Dox_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/A431_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/A431_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","con_3","Dox_1","Dox_2", "Dox_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

#M231

filter<-read.table("M231_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered
names(filter)<-c("con_1","con_2","con_3","5_FU_1","5_FU_2", "5_FU_3")

samples<-read.table("DESeq2/M231_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.




res<-results(dds, contrast=c("condition","5_FU","control")) #comparison between treatment (5_FU) and control (control)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","con_3","5_FU_1","5_FU_2", "5_FU_3","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/M231_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/M231_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","con_3","5_FU_1","5_FU_2", "5_FU_3","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)


#SJSA_GSE89807

filter<-read.table("SJSA_GSE89807_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered
names(filter)<-c("DMSO_1","DMSO_2","5_FU_1","5_FU_2" )

samples<-read.table("DESeq2/SJSA_GSE89807_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.




res<-results(dds, contrast=c("condition","5_FU","DMSO")) #comparison between treatment (5_FU) and control (DMSO)
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","DMSO_2","5_FU_1","5_FU_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/SJSA_GSE89807_DMSO_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/SJSA_GSE89807_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","DMSO_1","con_2","5_FU_1","5_FU_2","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)


#MRC5

filter<-read.table("MRC5_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered


samples<-read.table("DESeq2/MRC5_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData =filter , colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.




res<-results(dds, contrast=c("condition","UV_hr3","control"))
res<-as.data.frame(res)
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","UV_hr3_1","UV_hr3_2","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig <- resSig[order(resSig$padj),]

write.table(resSig, "DESeq2/MRC5_control_treatment_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


vals<-c()
F_vals<-c()
i<-0
name<-"DESeq2/MRC5_padj_0.05_FC"
resSig <- subset(res, padj < 0.05)
resSig<-cbind(RE=row.names(resSig),resSig)
resSig<-merge(resSig,filter, by.x="RE", by.y="row.names")
resSig<-resSig[,c("RE","con_1","con_2","UV_hr3_1","UV_hr3_2","baseMean",
                  "log2FoldChange","lfcSE","stat","pvalue","padj")]

while (i<4.5)
{
  temp<-resSig[abs(resSig$log2FoldChange) >= i,]
  
  F_vals<-c(F_vals,nrow(temp))
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  i<-i+0.5
}

temp<-data.frame("FC"=vals,"p_values"=p_vals,"padj"=F_vals)

library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

filter<-read.table("IMR90_hr24_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Table of combined telescope samples that have been filtered


#control and DMSO overlap hr24

#control
#Comparison by Control
samples<-read.table("DESeq2/hr24_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Sample table indicating the controls and treatments needs to follow the specific version of DESeq2 being used
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition) #Saves the data, sample table and design of the DESeq2 comparison we used condition DESeqDataSet object.
dds <- DESeq(dds) #Run DESeq2 differential anyalisis on prepped DESeqDataSet object.


res_con<-results(dds, contrast=c("condition","hr24_5_FU","control")) #comparison between treatment (hr24_5_FU) and control (control)
res_con<-as.data.frame(res_con)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","control_1","control_2","hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1","hr24_5_FU_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]



#DMSO
#Comparison by DMSO

res_DMSO<-results(dds, contrast=c("condition","hr24_5_FU","hr24_DMSO")) #comparison between treatment (hr24_5_FU) and control (hr24_DMSO)
res_DMSO<-as.data.frame(res_DMSO)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_DMSO<-resSig_DMSO[,c("RE","control_1","control_2","hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1","hr24_5_FU_2","baseMean",
                            "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#Merge control and DMSO comparison
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "control_1.x", "control_2.x", "hr24_DMSO_1.x", "hr24_DMSO_2.x","hr24_5_FU_1.x","hr24_5_FU_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "control_1", "control_2", "hr24_DMSO_1", "hr24_DMSO_2","hr24_5_FU_1","hr24_5_FU_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

nrow(resSig)

write.table(resSig, "DESeq2/IMR90_hr24_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)
