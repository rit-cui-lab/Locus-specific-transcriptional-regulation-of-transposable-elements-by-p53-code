#all RE occurance

library(VennDiagram)


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope")


#padj
HCT116_Nut_padj<-read.table("HCT116/DESeq2/HCT116_Nutlin_vs_DMSO_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

HCT116_FU_padj<-read.table("time_course/HCT116/DESeq2/HCT116_hr12_control_overlap_DMSO_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

IMR90_Nut_padj<-read.table("IMR90/DESeq2/IMR90_Nutlin_control_overlap_DMSO_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

IMR90_FU_padj<-read.table("time_course/IMR90/DESeq2/IMR90_hr12_control_overlap_DMSO_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

Saos_2_padj<-read.table("other_cell_lines/DESeq2/Saos-2_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

GM06170_padj<-read.table("other_cell_lines/DESeq2/GM06170_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

GM00011_padj<-read.table("other_cell_lines/DESeq2/GM00011_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

MCF7_GSE47042_padj<-read.table("other_cell_lines/DESeq2/MCF7_GSE47042_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

A549_padj<-read.table("other_cell_lines/DESeq2/A549_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

Fibro_padj<-read.table("other_cell_lines/DESeq2/Fibro_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

MCF10A_padj<-read.table("other_cell_lines/DESeq2/MCF10A_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

MCF7_GSE86221_padj<-read.table("other_cell_lines/DESeq2/MCF7_GSE86221_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

SJSA_padj<-read.table("other_cell_lines/DESeq2/SJSA_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

MV4_11_padj<-read.table("other_cell_lines/DESeq2/MV4_11_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

PAEC_padj<-read.table("other_cell_lines/DESeq2/MV4_11_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

IMR90_GSE139003_padj<-read.table("other_cell_lines/DESeq2/IMR90_GSE139003_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

MCF10A_hr4_padj<-read.table("other_cell_lines/DESeq2/MCF10A_hr4_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

hNCCs_padj<-read.table("other_cell_lines/DESeq2/hNCCs_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

hSMCs_padj<-read.table("other_cell_lines/DESeq2/hSMCs_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

A431_padj<-read.table("other_cell_lines/DESeq2/A431_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

M231_padj<-read.table("other_cell_lines/DESeq2/M231_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

SJSA_GSE89807_padj<-read.table("other_cell_lines/DESeq2/SJSA_GSE89807_DMSO_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

MRC5_padj<-read.table("other_cell_lines/DESeq2/MRC5_control_treatment_padj_0.05.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )


#merge together different FC option with padj or less results into one file
#one row per cell line.
#padj

name<-"RE_occurances/padj_0.05_FC"
i<-0
while (i<4.5)
{
  temp_HCT116_Nut<-HCT116_Nut_padj[abs(HCT116_Nut_padj$log2FoldChange) >= i,]
  
  temp_HCT116_FU<-HCT116_FU_padj[abs(HCT116_FU_padj$con_log2FoldChange) >= i,]
  temp_HCT116_FU<-temp_HCT116_FU[abs(temp_HCT116_FU$DMSO_log2FoldChange) >= i,]
  
  temp_IMR90_Nut<-IMR90_Nut_padj[abs(IMR90_Nut_padj$con_log2FoldChange) >= i,]
  temp_IMR90_Nut<-temp_IMR90_Nut[abs(temp_IMR90_Nut$DMSO_log2FoldChange) >= i,]
  
  temp_IMR90_Nut<-IMR90_Nut_padj[abs(IMR90_Nut_padj$con_log2FoldChange) >= i,]
  temp_IMR90_Nut<-temp_IMR90_Nut[abs(temp_IMR90_Nut$DMSO_log2FoldChange) >= i,]
  
  temp_IMR90_FU<-IMR90_FU_padj[abs(IMR90_FU_padj$con_log2FoldChange) >= i,]
  temp_IMR90_FU<-temp_IMR90_FU[abs(temp_IMR90_FU$DMSO_log2FoldChange) >= i,]
  
  temp_Saos_2<-Saos_2_padj[abs(Saos_2_padj$log2FoldChange) >= i,]
  
  temp_GM06170<-GM06170_padj[abs(GM06170_padj$log2FoldChange) >= i,]
  
  temp_GM00011<-GM00011_padj[abs(GM00011_padj$log2FoldChange) >= i,]
  
  temp_MCF7_GSE47042<-MCF7_GSE47042_padj[abs(MCF7_GSE47042_padj$log2FoldChange) >= i,]
  
  temp_A549<-A549_padj[abs(A549_padj$log2FoldChange) >= i,]
  
  temp_Fibro<-Fibro_padj[abs(Fibro_padj$log2FoldChange) >= i,]
  
  temp_MCF10A<-MCF10A_padj[abs(MCF10A_padj$log2FoldChange) >= i,]
  
  temp_MCF7_GSE86221<-MCF7_GSE86221_padj[abs(MCF7_GSE86221_padj$log2FoldChange) >= i,]
  
  temp_SJSA<-SJSA_padj[abs(SJSA_padj$log2FoldChange) >= i,]
  
  temp_MV4_11<-MV4_11_padj[abs(MV4_11_padj$log2FoldChange) >= i,]
  
  temp_PAEC<-PAEC_padj[abs(PAEC_padj$log2FoldChange) >= i,]
  
  temp_IMR90_GSE139003<-IMR90_GSE139003_padj[abs(IMR90_GSE139003_padj$log2FoldChange) >= i,]
  
  temp_MCF10A_hr4<-MCF10A_hr4_padj[abs(MCF10A_hr4_padj$log2FoldChange) >= i,]
  
  temp_hNCCs<-hNCCs_padj[abs(hNCCs_padj$log2FoldChange) >= i,]
  
  temp_hSMCs<-hSMCs_padj[abs(hSMCs_padj$log2FoldChange) >= i,]
  
  temp_A431<-A431_padj[abs(A431_padj$log2FoldChange) >= i,]
  
  temp_M231<-M231_padj[abs(M231_padj$log2FoldChange) >= i,]
  
  temp_SJSA_GSE89807<-SJSA_GSE89807_padj[abs(SJSA_GSE89807_padj$log2FoldChange) >= i,]
  
  temp_MRC5<-MRC5_padj[abs(MRC5_padj$log2FoldChange) >= i,]
  
  temp<-merge(temp_HCT116_Nut[,c("RE","padj")], temp_HCT116_FU[,c("RE","con_padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU")
  
  temp<-merge(temp, temp_IMR90_Nut[,c("RE","con_padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin")
  
  temp<-merge(temp, temp_IMR90_FU[,c("RE","con_padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU")
  
  temp<-merge(temp, temp_Saos_2[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin")
  
  temp<-merge(temp, temp_GM06170[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin")
  
  temp<-merge(temp, temp_GM00011[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin")
  
  temp<-merge(temp, temp_MCF7_GSE47042[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin")
  
  temp<-merge(temp, temp_A549[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin")
  
  temp<-merge(temp, temp_Fibro[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin")
  
  temp<-merge(temp, temp_MCF10A[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin")
  
  temp<-merge(temp, temp_MCF7_GSE86221[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin")
  
  temp<-merge(temp, temp_SJSA[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin")
  
  temp<-merge(temp, temp_MV4_11[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin")
  
  temp<-merge(temp, temp_PAEC[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin")
  
  temp<-merge(temp, temp_IMR90_GSE139003[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin")
  
  temp<-merge(temp, temp_MCF10A_hr4[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin")
  
  temp<-merge(temp, temp_hNCCs[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin",
                 "hNCCs_Nutlin")
  
  temp<-merge(temp, temp_hSMCs[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin",
                 "hNCCs_Nutlin","hSMCs_Nutlin")
  
  temp<-merge(temp, temp_A431[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin",
                 "hNCCs_Nutlin","hSMCs_Nutlin","A431_Doxorubicin")
  
  temp<-merge(temp, temp_M231[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin",
                 "hNCCs_Nutlin","hSMCs_Nutlin","A431_Doxorubicin","M231_5_FU")
  
  temp<-merge(temp, temp_SJSA_GSE89807[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin",
                 "hNCCs_Nutlin","hSMCs_Nutlin","A431_Doxorubicin","M231_5_FU","SJSA_GSE89807_5_FU")
  
  temp<-merge(temp, temp_MRC5[,c("RE","padj")],by="RE", all=TRUE)
  names(temp)<-c("RE","HCT116_Nutlin","HCT116_5_FU", "IMR90_Nutlin", "IMR90_5_FU","Saos_2_doxycyclin","GM06170_Doxorubicin",
                 "GM00011_Doxorubicin","MCF7_GSE47042_Nutlin","A549_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","MCF7_GSE86221_Nutlin",
                 "SJSA_Nutlin", "MV4_11_Nutlin", "PAEC_Nutlin", "IMR90_GSE139003_hr12_Nutlin", "MCF10A_hr4_Nutlin",
                 "hNCCs_Nutlin","hSMCs_Nutlin","A431_Doxorubicin","M231_5_FU","SJSA_GSE89807_5_FU","MRC5_UV")
  
  #label TE with differential anyalisis 1 and without 0 for each cell line
  temp$HCT116_Nutlin[!is.na(temp$HCT116_Nutlin)]<-1
  temp$HCT116_Nutlin[is.na(temp$HCT116_Nutlin)]<-0
  
  temp$HCT116_5_FU[!is.na(temp$HCT116_5_FU)]<-1
  temp$HCT116_5_FU[is.na(temp$HCT116_5_FU)]<-0
  
  temp$IMR90_Nutlin[!is.na(temp$IMR90_Nutlin)]<-1
  temp$IMR90_Nutlin[is.na(temp$IMR90_Nutlin)]<-0
  
  temp$IMR90_5_FU[!is.na(temp$IMR90_5_FU)]<-1
  temp$IMR90_5_FU[is.na(temp$IMR90_5_FU)]<-0
  
  temp$Saos_2_doxycyclin[!is.na(temp$Saos_2_doxycyclin)]<-1
  temp$Saos_2_doxycyclin[is.na(temp$Saos_2_doxycyclin)]<-0
  
  temp$GM06170_Doxorubicin[!is.na(temp$GM06170_Doxorubicin)]<-1
  temp$GM06170_Doxorubicin[is.na(temp$GM06170_Doxorubicin)]<-0
  
  temp$GM00011_Doxorubicin[!is.na(temp$GM00011_Doxorubicin)]<-1
  temp$GM00011_Doxorubicin[is.na(temp$GM00011_Doxorubicin)]<-0
  
  temp$MCF7_GSE47042_Nutlin[!is.na(temp$MCF7_GSE47042_Nutlin)]<-1
  temp$MCF7_GSE47042_Nutlin[is.na(temp$MCF7_GSE47042_Nutlin)]<-0
  
  temp$A549_Nutlin[!is.na(temp$A549_Nutlin)]<-1
  temp$A549_Nutlin[is.na(temp$A549_Nutlin)]<-0
  
  temp$Fibro_Nutlin[!is.na(temp$Fibro_Nutlin)]<-1
  temp$Fibro_Nutlin[is.na(temp$Fibro_Nutlin)]<-0
  
  temp$MCF10A_Nutlin[!is.na(temp$MCF10A_Nutlin)]<-1
  temp$MCF10A_Nutlin[is.na(temp$MCF10A_Nutlin)]<-0
  
  temp$MCF7_GSE86221_Nutlin[!is.na(temp$MCF7_GSE86221_Nutlin)]<-1
  temp$MCF7_GSE86221_Nutlin[is.na(temp$MCF7_GSE86221_Nutlin)]<-0
  
  temp$SJSA_Nutlin[!is.na(temp$SJSA_Nutlin)]<-1
  temp$SJSA_Nutlin[is.na(temp$SJSA_Nutlin)]<-0
  
  temp$MV4_11_Nutlin[!is.na(temp$MV4_11_Nutlin)]<-1
  temp$MV4_11_Nutlin[is.na(temp$MV4_11_Nutlin)]<-0
  
  temp$PAEC_Nutlin[!is.na(temp$PAEC_Nutlin)]<-1
  temp$PAEC_Nutlin[is.na(temp$PAEC_Nutlin)]<-0
  
  temp$IMR90_GSE139003_hr12_Nutlin[!is.na(temp$IMR90_GSE139003_hr12_Nutlin)]<-1
  temp$IMR90_GSE139003_hr12_Nutlin[is.na(temp$IMR90_GSE139003_hr12_Nutlin)]<-0
  
  temp$MCF10A_hr4_Nutlin[!is.na(temp$MCF10A_hr4_Nutlin)]<-1
  temp$MCF10A_hr4_Nutlin[is.na(temp$MCF10A_hr4_Nutlin)]<-0
  
  temp$hNCCs_Nutlin[!is.na(temp$hNCCs_Nutlin)]<-1
  temp$hNCCs_Nutlin[is.na(temp$hNCCs_Nutlin)]<-0
  
  temp$hSMCs_Nutlin[!is.na(temp$hSMCs_Nutlin)]<-1
  temp$hSMCs_Nutlin[is.na(temp$hSMCs_Nutlin)]<-0
  
  temp$A431_Doxorubicin[!is.na(temp$A431_Doxorubicin)]<-1
  temp$A431_Doxorubicin[is.na(temp$A431_Doxorubicin)]<-0
  
  temp$M231_5_FU[!is.na(temp$M231_5_FU)]<-1
  temp$M231_5_FU[is.na(temp$M231_5_FU)]<-0
  
  temp$SJSA_GSE89807_5_FU[!is.na(temp$SJSA_GSE89807_5_FU)]<-1
  temp$SJSA_GSE89807_5_FU[is.na(temp$SJSA_GSE89807_5_FU)]<-0
  
  temp$MRC5_UV[!is.na(temp$MRC5_UV)]<-1
  temp$MRC5_UV[is.na(temp$MRC5_UV)]<-0
  
  name_temp<-paste(name,i, sep="_")
  name_temp<-paste(name_temp,"csv",sep=".")
  
  write.table(temp, name_temp, sep=',', row.names = F,col.names = T, quote = F)
  
  
  
  i<-i+0.5
  
}
#summ row was added for some of the tables for further anyalisis.