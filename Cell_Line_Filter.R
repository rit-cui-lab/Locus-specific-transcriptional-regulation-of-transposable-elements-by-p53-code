#merge the telescope results for each sample into one file via the final_count and contains at least five counts per million at least four times.

library("edgeR")
setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines")

#Saos-2
con_1<-read.table(file="Saos-2_con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="Saos-2_con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
dox_1<-read.table(file="Saos-2_p53_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
dox_1<-dox_1[,c("transcript","final_count")]
dox_2<-read.table(file="Saos-2_p53_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
dox_2<-dox_2[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,dox_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","dox_1")

final<-merge(final,dox_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","dox_1", "dox_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","dox_1", "dox_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "Saos-2_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


#GM06170

con_1<-read.table(file="GM06170_con_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="GM06170_con_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
dox_1<-read.table(file="GM06170_dox_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
dox_1<-dox_1[,c("transcript","final_count")]
dox_2<-read.table(file="GM06170_dox_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
dox_2<-dox_2[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,dox_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","dox_1")

final<-merge(final,dox_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","dox_1", "dox_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","dox_1", "dox_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,] 


write.table(final_temp, "GM06170_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


#GM00011

con_1<-read.table(file="GM00011_con_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="GM00011_con_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
dox_1<-read.table(file="GM00011_dox_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
dox_1<-dox_1[,c("transcript","final_count")]
dox_2<-read.table(file="GM00011_dox_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
dox_2<-dox_2[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,dox_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","dox_1")

final<-merge(final,dox_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","dox_1", "dox_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","dox_1", "dox_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "GM00011_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#MCF7_GSE47042
con_1<-read.table(file="MCF7_GSE47042_con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="MCF7_GSE47042_con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
nut_1<-read.table(file="MCF7_GSE47042_Nut_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_1<-nut_1[,c("transcript","final_count")]
nut_2<-read.table(file="MCF7_GSE47042_Nut_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_2<-nut_2[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","Nutlin_1")

final<-merge(final,nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","Nutlin_1", "Nutlin_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","Nutlin_1", "Nutlin_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "MCF7_GSE47042_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#A549

con_1<-read.table(file="A549_con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="A549_con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
con_3<-read.table(file="A549_con_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_3<-con_3[,c("transcript","final_count")]
nut_1<-read.table(file="A549_Nut_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_1<-nut_1[,c("transcript","final_count")]
nut_2<-read.table(file="A549_Nut_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_2<-nut_2[,c("transcript","final_count")]
nut_3<-read.table(file="A549_Nut_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_3<-nut_3[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final, con_3, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3")

final<-merge(final,nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","con_3","Nutlin_1")

final<-merge(final,nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","con_3","Nutlin_1", "Nutlin_2")

final<-merge(final,nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","con_3","Nutlin_1", "Nutlin_2", "Nutlin_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","con_3","Nutlin_1", "Nutlin_2", "Nutlin_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "A549_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#Fibro

DMSO_1<-read.table(file="Fibro_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="Fibro_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="Fibro_DMSO_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]
nut_1<-read.table(file="Fibro_Nut_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_1<-nut_1[,c("transcript","final_count")]
nut_2<-read.table(file="Fibro_Nut_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_2<-nut_2[,c("transcript","final_count")]
nut_3<-read.table(file="Fibro_Nut_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_3<-nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final,nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","DMSO_3","Nutlin_1")

final<-merge(final,nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","DMSO_3","Nutlin_1", "Nutlin_2")

final<-merge(final,nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","DMSO_3","Nutlin_1", "Nutlin_2", "Nutlin_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2","DMSO_3","Nutlin_1", "Nutlin_2", "Nutlin_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "Fibro_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#MCF10A
DMSO_1<-read.table(file="MCF10A_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="MCF10A_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="MCF10A_DMSO_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]
nut_1<-read.table(file="MCF10A_Nut_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_1<-nut_1[,c("transcript","final_count")]
nut_2<-read.table(file="MCF10A_Nut_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_2<-nut_2[,c("transcript","final_count")]
nut_3<-read.table(file="MCF10A_Nut_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_3<-nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final,nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","DMSO_3","Nutlin_1")

final<-merge(final,nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","DMSO_3","Nutlin_1", "Nutlin_2")

final<-merge(final,nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","DMSO_3","Nutlin_1", "Nutlin_2", "Nutlin_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2","DMSO_3","Nutlin_1", "Nutlin_2", "Nutlin_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "MCF10A_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#MCF7_GSE86221
DMSO_1<-read.table(file="MCF7_GSE86221_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="MCF7_GSE86221_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
nut_1<-read.table(file="MCF7_GSE86221_Nut_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_1<-nut_1[,c("transcript","final_count")]
nut_2<-read.table(file="MCF7_GSE86221_Nut_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_2<-nut_2[,c("transcript","final_count")]



#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final,nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","Nutlin_1")

final<-merge(final,nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","Nutlin_1", "Nutlin_2")



#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2","Nutlin_1", "Nutlin_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "MCF7_GSE86221_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#SJSA

DMSO_1<-read.table(file="SJSA_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="SJSA_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
nut_1<-read.table(file="SJSA_Nut_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_1<-nut_1[,c("transcript","final_count")]
nut_2<-read.table(file="SJSA_Nut_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
nut_2<-nut_2[,c("transcript","final_count")]



#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final,nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","Nutlin_1")

final<-merge(final,nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","Nutlin_1", "Nutlin_2")



#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2","Nutlin_1", "Nutlin_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "SJSA_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


library("edgeR")
setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines")


#IMR90_GSE139003
DMSO_1<-read.table(file="IMR90_GSE139003_DMSO_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="IMR90_GSE139003_DMSO_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]

hr12_Nut_1<-read.table(file="IMR90_GSE139003_hr12_Nut_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr12_Nut_1<-hr12_Nut_1[,c("transcript","final_count")]
hr12_Nut_2<-read.table(file="IMR90_GSE139003_hr12_Nut_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr12_Nut_2<-hr12_Nut_2[,c("transcript","final_count")]



#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final,hr12_Nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","hr12_Nut_1")

final<-merge(final,hr12_Nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2","hr12_Nut_1","hr12_Nut_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2","hr12_Nut_1","hr12_Nut_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "IMR90_GSE139003_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#MCF10A still work to do to get ready.
DMSO_1<-read.table(file="MCF10A_DMSO_hr4_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="MCF10A_DMSO_hr4_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="MCF10A_DMSO_hr4_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]

Nut_1<-read.table(file="MCF10A_Nut_hr4_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_1<-Nut_1[,c("transcript","final_count")]
Nut_2<-read.table(file="MCF10A_Nut_hr4_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_2<-Nut_2[,c("transcript","final_count")]
Nut_3<-read.table(file="MCF10A_Nut_hr4_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_3<-Nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final, Nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1")

final<-merge(final, Nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2")

final<-merge(final, Nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "MCF10A_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


#PAEC

DMSO_1<-read.table(file="PAEC_DMSO_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="PAEC_DMSO_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="PAEC_DMSO_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]

Nut_1<-read.table(file="PAEC_Nut_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_1<-Nut_1[,c("transcript","final_count")]
Nut_2<-read.table(file="PAEC_Nut_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_2<-Nut_2[,c("transcript","final_count")]
Nut_3<-read.table(file="PAEC_Nut_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_3<-Nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final, Nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1")

final<-merge(final, Nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2")

final<-merge(final, Nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "PAEC_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#hNCCs

DMSO_1<-read.table(file="hNCCs_DMSO_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="hNCCs_DMSO_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="hNCCs_DMSO_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]

Nut_1<-read.table(file="hNCCs_Nut_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_1<-Nut_1[,c("transcript","final_count")]
Nut_2<-read.table(file="hNCCs_Nut_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_2<-Nut_2[,c("transcript","final_count")]
Nut_3<-read.table(file="hNCCs_Nut_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_3<-Nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final, Nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1")

final<-merge(final, Nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2")

final<-merge(final, Nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "hNCCs_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


#hSMCs

DMSO_1<-read.table(file="hSMCs_DMSO_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="hSMCs_DMSO_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="hSMCs_DMSO_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]

Nut_1<-read.table(file="hSMCs_Nut_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_1<-Nut_1[,c("transcript","final_count")]
Nut_2<-read.table(file="hSMCs_Nut_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_2<-Nut_2[,c("transcript","final_count")]
Nut_3<-read.table(file="hSMCs_Nut_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_3<-Nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final, Nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1")

final<-merge(final, Nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2")

final<-merge(final, Nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "hSMCs_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#MV4_11

DMSO_1<-read.table(file="MV4_11_DMSO_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="MV4_11_DMSO_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]
DMSO_3<-read.table(file="MV4_11_DMSO_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_3<-DMSO_3[,c("transcript","final_count")]

Nut_1<-read.table(file="MV4_11_Nut_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_1<-Nut_1[,c("transcript","final_count")]
Nut_2<-read.table(file="MV4_11_Nut_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_2<-Nut_2[,c("transcript","final_count")]
Nut_3<-read.table(file="MV4_11_Nut_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Nut_3<-Nut_3[,c("transcript","final_count")]


#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, DMSO_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3")

final<-merge(final, Nut_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1")

final<-merge(final, Nut_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2")

final<-merge(final, Nut_3, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2", "DMSO_3", "Nut_1", "Nut_2", "Nut_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "MV4_11_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


library("edgeR")
setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines")


#HCT116_GSE158021
con_1<-read.table(file="HCT116_GSE158021_con_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="HCT116_GSE158021_con_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]

hr2_5_FU_1<-read.table(file="HCT116_GSE158021_5_FU_hr2_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr2_5_FU_1<-hr2_5_FU_1[,c("transcript","final_count")]
hr2_5_FU_2<-read.table(file="HCT116_GSE158021_5_FU_hr2_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr2_5_FU_2<-hr2_5_FU_2[,c("transcript","final_count")]



#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,hr2_5_FU_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr2_5_FU_1")

final<-merge(final,hr2_5_FU_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr2_5_FU_1","hr2_5_FU_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","hr2_5_FU_1","hr2_5_FU_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "HCT116_GSE158021_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#A431
con_1<-read.table(file="A431_untreated_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="A431_untreated_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
con_3<-read.table(file="A431_untreated_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_3<-con_3[,c("transcript","final_count")]

Dox_1<-read.table(file="A431_Dox_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Dox_1<-Dox_1[,c("transcript","final_count")]
Dox_2<-read.table(file="A431_Dox_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Dox_2<-Dox_2[,c("transcript","final_count")]
Dox_3<-read.table(file="A431_Dox_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
Dox_3<-Dox_3[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final, con_3, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3")

final<-merge(final, Dox_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3", "Dox_1")

final<-merge(final, Dox_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3", "Dox_1", "Dox_2")

final<-merge(final, Dox_3, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3", "Dox_1", "Dox_2", "Dox_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2", "con_3", "Dox_1", "Dox_2", "Dox_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "A431_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#SJSA_GSE89807
DMSO_1<-read.table(file="SJSA_GSE89807_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_1<-DMSO_1[,c("transcript","final_count")]
DMSO_2<-read.table(file="SJSA_GSE89807_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
DMSO_2<-DMSO_2[,c("transcript","final_count")]


treat_1<-read.table(file="SJSA_GSE89807_5_FU_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
treat_1<-treat_1[,c("transcript","final_count")]
treat_2<-read.table(file="SJSA_GSE89807_5_FU_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
treat_2<-treat_2[,c("transcript","final_count")]



#merge final
final<-merge(DMSO_1, DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2")

final<-merge(final, treat_1, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "5_FU_1")

final<-merge(final, treat_2, by="transcript",all = TRUE)
names(final)<-c("transcript","DMSO_1","DMSO_2", "5_FU_1", "5_FU_2")



#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("DMSO_1","DMSO_2", "5_FU_1", "5_FU_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "SJSA_GSE89807_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#MRC5
con_1<-read.table(file="MRC5_con_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="MRC5_con_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]

UV_hr3_1<-read.table(file="MRC5_UV_hr3_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
UV_hr3_1<-UV_hr3_1[,c("transcript","final_count")]
UV_hr3_2<-read.table(file="MRC5_UV_hr3_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
UV_hr3_2<-UV_hr3_2[,c("transcript","final_count")]



#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,UV_hr3_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","UV_hr3_1")

final<-merge(final,UV_hr3_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","UV_hr3_1","UV_hr3_2")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","UV_hr3_1","UV_hr3_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "MRC5_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#M231 script
 
library("edgeR")
setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/other_cell_lines")

#M231
con_1<-read.table(file="M231_untreated_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="M231_untreated_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
con_3<-read.table(file="M231_untreated_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_3<-con_3[,c("transcript","final_count")]

treat_1<-read.table(file="M231_5_FU_1_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
treat_1<-treat_1[,c("transcript","final_count")]
treat_2<-read.table(file="M231_5_FU_1_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
treat_2<-treat_2[,c("transcript","final_count")]
treat_3<-read.table(file="M231_5_FU_1_3.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
treat_3<-treat_3[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final, con_3, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3")

final<-merge(final, treat_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3", "5_FU_1")

final<-merge(final, treat_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3", "5_FU_1", "5_FU_2")

final<-merge(final, treat_3, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2", "con_3", "5_FU_1", "5_FU_2", "5_FU_3")


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2", "con_3", "5_FU_1", "5_FU_2", "5_FU_3")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "M231_filter.csv", sep=',', row.names = T, col.names = T, quote = F)
library("edgeR")

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/count_tables")

IMR90<-read.table(file="IMR90_DMSOvsNutlin_countsTable.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)

HcT116<-read.table(file="HCT116_p53WT_DMSOvsNutlin_countsTable.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

IMR90_DMSO_1<-read.table(file="IMR90_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
IMR90_DMSO_1<-IMR90_DMSO_1[,c("transcript","final_count")]
IMR90_DMSO_2<-read.table(file="IMR90_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
IMR90_DMSO_2<-IMR90_DMSO_2[,c("transcript","final_count")]
IMR90_Nutlin_1<-read.table(file="IMR90_Nutlin_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
IMR90_Nutlin_1<-IMR90_Nutlin_1[,c("transcript","final_count")]
IMR90_Nutlin_2<-read.table(file="IMR90_Nutlin_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
IMR90_Nutlin_2<-IMR90_Nutlin_2[,c("transcript","final_count")]
IMR90_con_1<-read.table(file="IMR90_con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
IMR90_con_1<-IMR90_con_1[,c("transcript","final_count")]
IMR90_con_2<-read.table(file="IMR90_con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
IMR90_con_2<-IMR90_con_1[,c("transcript","final_count")]

#merge final
final<-merge(IMR90_con_1, IMR90_con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,IMR90_DMSO_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","DMSO_1")

final<-merge(final,IMR90_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","DMSO_1", "DMSO_2")

final<-merge(final,IMR90_Nutlin_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","DMSO_1", "DMSO_2", "Nutlin_1")

final<-merge(final,IMR90_Nutlin_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","DMSO_1", "DMSO_2", "Nutlin_1", 
                "Nutlin_2")

#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","DMSO_1", "DMSO_2", "Nutlin_1", "Nutlin_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.


keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final_temp <- final[keep,]


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope/")
write.table(final_temp, "IMR90/IMR90_total_Nutlin_my_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#Nutlin and DMSO
final<-final[,c("DMSO_1", "DMSO_2", "Nutlin_1", "Nutlin_2")]

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final_temp <- final[keep,]


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope/")
write.table(final_temp, "IMR90/IMR90_total_Nutlin_my_filter_NutlinvsDMSO.csv", sep=',', row.names = T, col.names = T, quote = F)

HCT116_DMSO_1<-read.table(file="HTC116_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
HCT116_DMSO_1<-HCT116_DMSO_1[,c("transcript","final_count")]
HCT116_DMSO_2<-read.table(file="HTC116_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
HCT116_DMSO_2<-HCT116_DMSO_2[,c("transcript","final_count")]
HCT116_Nutlin_1<-read.table(file="HTC116_Nutlin_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
HCT116_Nutlin_1<-HCT116_Nutlin_1[,c("transcript","final_count")]
HCT116_Nutlin_2<-read.table(file="HTC116_Nutlin_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
HCT116_Nutlin_2<-HCT116_Nutlin_2[,c("transcript","final_count")]
#merge final
final<-merge(HCT116_DMSO_1, HCT116_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2")
final<-merge(final, HCT116_Nutlin_1, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1")
final<-merge(final, HCT116_Nutlin_2, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1","wt_Nutlin_2")
final<-merge(final, HCT116_Knockout_DMSO_1, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1","wt_Nutlin_2", "knockout_DMSO_1")
final<-merge(final, HCT116_Knockout_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1","wt_Nutlin_2", "knockout_DMSO_1",
                "knockout_DMSO_2")
final<-merge(final, HCT116_Knockout_Nutlin_1, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1","wt_Nutlin_2", "knockout_DMSO_1",
                "knockout_DMSO_2","knockout_Nutlin_1")
final<-merge(final, HCT116_Knockout_Nutlin_2, by="transcript",all = TRUE)
names(final)<-c("transcript","wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1","wt_Nutlin_2", "knockout_DMSO_1",
                "knockout_DMSO_2","knockout_Nutlin_1","knockout_Nutlin_2")


final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("wt_DMSO_1","wt_DMSO_2", "wt_Nutlin_1","wt_Nutlin_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)
final_temp <- final[keep,]

names(final)<-c("DMSO_1","DMSO_2","Nutlin_1","Nutlin_2")


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope/")
write.table(final_temp, "HCT116/HCT116_total_Nutlin_my_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#hr6 filter
#hr6 IMR90 filter
library("edgeR")

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

con_1<-read.table(file="con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]

hr6_DMSO_1<-read.table(file="hr6_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_DMSO_1<-hr6_DMSO_1[,c("transcript","final_count")]
hr6_DMSO_2<-read.table(file="hr6_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_DMSO_2<-hr6_DMSO_2[,c("transcript","final_count")]

hr6_5_FU_1<-read.table(file="hr6_5-FU_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_5_FU_1<-hr6_5_FU_1[,c("transcript","final_count")]
hr6_5_FU_2<-read.table(file="hr6_5-FU_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_5_FU_2<-hr6_5_FU_2[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,hr6_DMSO_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1")

final<-merge(final,hr6_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2")

final<-merge(final,hr6_5_FU_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1" )

final<-merge(final,hr6_5_FU_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2" )


#clean up IMR90

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final <- final[keep,]

write.table(final, "IMR90_hr6_5_FU_time_course_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

#HCT116 hr6

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

con_1<-read.table(file="con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]

hr6_DMSO_1<-read.table(file="hr6_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_DMSO_1<-hr6_DMSO_1[,c("transcript","final_count")]
hr6_DMSO_2<-read.table(file="hr6_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_DMSO_2<-hr6_DMSO_2[,c("transcript","final_count")]

hr6_5_FU_1<-read.table(file="hr6_5-FU_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_5_FU_1<-hr6_5_FU_1[,c("transcript","final_count")]
hr6_5_FU_2<-read.table(file="hr6_5-FU_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr6_5_FU_2<-hr6_5_FU_2[,c("transcript","final_count")]


#merge final
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2")

final<-merge(final,hr6_DMSO_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1")

final<-merge(final,hr6_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2")

final<-merge(final,hr6_5_FU_1, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1" )

final<-merge(final,hr6_5_FU_2, by="transcript",all = TRUE)
names(final)<-c("transcript","con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2" )


#clean up HCT116

final[is.na(final)] <- 0 #Turns any NAs into 0s


row.names(final)<-final$transcript
final<-final[,c("con_1","con_2","hr6_DMSO_1","hr6_DMSO_2","hr6_5_FU_1","hr6_5_FU_2")]
final<-subset(final, row.names(final)!="__no_feature")

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix witht he counts to keep list
final <- final[keep,]

write.table(final, "HCT116_hr6_5_FU_time_course_filter.csv", sep=',', row.names = T, col.names = T, quote = F)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

con_1<-read.table(file="con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
hr24_DMSO_1<-read.table(file="hr24_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_DMSO_1<-hr24_DMSO_1[,c("transcript","final_count")]
hr24_DMSO_2<-read.table(file="hr24_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_DMSO_2<-hr24_DMSO_2[,c("transcript","final_count")]
hr24_5_FU_1<-read.table(file="hr24_5-FU_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_5_FU_1<-hr24_5_FU_1[,c("transcript","final_count")]
hr24_5_FU_2<-read.table(file="hr24_5-FU_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_5_FU_2<-hr24_5_FU_2[,c("transcript","final_count")]

#merge unique
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2")

final<-merge(final, hr24_DMSO_1, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2", "hr24_DMSO_1")

final<-merge(final, hr24_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2", "hr24_DMSO_1","hr24_DMSO_2")

final<-merge(final, hr24_5_FU_1, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2","hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1")

final<-merge(final, hr24_5_FU_2, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2", "hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1", "hr24_5_FU_2")

#clean up

final[is.na(final)] <- 0 #Turns any NAs into 0s

row.names(final)<-final$transcript
final<-final[,c("control_1","control_2", "hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1", "hr24_5_FU_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "HCT116_hr24_5_FU_time_course_filter.csv", sep=',', row.names = T, col.names = T, quote = F)


#turn to hr24

library("DESeq2")
library(VennDiagram)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/HCT116")

filter<-read.table("HCT116_hr24_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

#setwd("C:/Users/dazer/Documents/Science/Bioinformatics/Bioinformatics_work/DE_REs/telescope/time_course/HCT116/DESeq2")

#control and DMSO overlap hr12

#control
samples<-read.table("DESeq2/hr24_table_samples.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
dds <- DESeqDataSetFromMatrix(countData = filter, colData = samples, design = ~ condition)
dds <- DESeq(dds)


res_con<-results(dds, contrast=c("condition","hr24_5_FU","control"))
res_con<-as.data.frame(res_con)
#FDR<-p.adjust(res_con$pvalue, method = "BH")
#res_con<-cbind(res_con,FDR)
resSig_con <- subset(res_con, padj < 0.05)
resSig_con<-cbind(RE=row.names(resSig_con),resSig_con)
resSig_con<-merge(resSig_con,filter, by.x="RE", by.y="row.names")
resSig_con<-resSig_con[,c("RE","control_1","control_2","hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1","hr24_5_FU_2","baseMean",
                          "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_con <- resSig_con[order(resSig_con$padj),]



#DMSO

res_DMSO<-results(dds, contrast=c("condition","hr24_5_FU","hr24_DMSO"))
res_DMSO<-as.data.frame(res_DMSO)
#FDR<-p.adjust(res_DMSO$pvalue, method = "BH")
#res_DMSO<-cbind(res_DMSO,FDR)
resSig_DMSO <- subset(res_DMSO, padj < 0.05)
resSig_DMSO<-cbind(RE=row.names(resSig_DMSO),resSig_DMSO)
resSig_DMSO<-merge(resSig_DMSO,filter, by.x="RE", by.y="row.names")
resSig_DMSO<-resSig_DMSO[,c("RE","control_1","control_2","hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1","hr24_5_FU_2","baseMean",
                            "log2FoldChange","lfcSE","stat","pvalue","padj")]
resSig_DMSO <- resSig_DMSO[order(resSig_DMSO$padj),]

#overlap
resSig<-merge(resSig_con,resSig_DMSO, by="RE")
resSig<-resSig[,c("RE", "control_1.x", "control_2.x", "hr24_DMSO_1.x", "hr24_DMSO_2.x","hr24_5_FU_1.x","hr24_5_FU_2.x",
                  "baseMean.x","log2FoldChange.x","lfcSE.x","stat.x","pvalue.x","padj.x","baseMean.y",
                  "log2FoldChange.y","lfcSE.y","stat.y","pvalue.y","padj.y" )]

names(resSig)<-c("RE", "control_1", "control_2", "hr24_DMSO_1", "hr24_DMSO_2","hr24_5_FU_1","hr24_5_FU_2",
                 "con_baseMean","con_log2FoldChange","con_lfcSE","con_stat","con_pvalue","con_padj",
                 "DMSO_baseMean","DMSO_log2FoldChange","DMSO_lfcSE","DMSO_stat","DMSO_pvalue","DMSO_padj")

nrow(resSig)

write.table(resSig, "DESeq2/HCT116_hr24_control_overlap_DMSO_padj_0.05.csv", sep=',', row.names = F,col.names = T, quote = F)


#IMR90

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/IMR90")

con_1<-read.table(file="con_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_1<-con_1[,c("transcript","final_count")]
con_2<-read.table(file="con_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
con_2<-con_2[,c("transcript","final_count")]
hr24_DMSO_1<-read.table(file="hr24_DMSO_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_DMSO_1<-hr24_DMSO_1[,c("transcript","final_count")]
hr24_DMSO_2<-read.table(file="hr24_DMSO_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_DMSO_2<-hr24_DMSO_2[,c("transcript","final_count")]
hr24_5_FU_1<-read.table(file="hr24_5-FU_1.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_5_FU_1<-hr24_5_FU_1[,c("transcript","final_count")]
hr24_5_FU_2<-read.table(file="hr24_5-FU_2.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)
hr24_5_FU_2<-hr24_5_FU_2[,c("transcript","final_count")]

#merge unique
final<-merge(con_1, con_2, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2")

final<-merge(final, hr24_DMSO_1, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2", "hr24_DMSO_1")

final<-merge(final, hr24_DMSO_2, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2", "hr24_DMSO_1","hr24_DMSO_2")

final<-merge(final, hr24_5_FU_1, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2","hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1")

final<-merge(final, hr24_5_FU_2, by="transcript",all = TRUE)
names(final)<-c("transcript","control_1","control_2", "hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1", "hr24_5_FU_2")

#clean up

final[is.na(final)] <- 0 #Turns any NAs into 0s

row.names(final)<-final$transcript
final<-final[,c("control_1","control_2", "hr24_DMSO_1","hr24_DMSO_2","hr24_5_FU_1", "hr24_5_FU_2")]
final<-subset(final, row.names(final)!="__no_feature") #remove the "__no_feature" row.

keep <- rowSums(cpm(final) > 5) >= 4 #Selection five counts per million at least four times
table(keep)

#subset the counts matrix with the counts to keep list
final_temp <- final[keep,]


write.table(final_temp, "IMR90_hr24_5_FU_time_course_filter.csv", sep=',', row.names = T, col.names = T, quote = F)
