#Overlaps MAC2 results (ChIP-seq) with the telescope DESeq2 result for for HCT116 5-FU

library("GenomicRanges") #overlaps via the genomic location


setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope/")


HCT116_FC_0<-read.table("comparisons/FC_0_padj_5-FU_HCT116_only.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
HCT116_FC_1<-read.table("comparisons/FC_1_padj_5-FU_HCT116_only.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

IMR90_FC_0<-read.table("comparisons/FC_0_padj_5-FU_IMR90_only.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
IMR90_FC_1<-read.table("comparisons/FC_1_padj_5-FU_IMR90_only.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

HCT116_peaks<-read.table("RE_peaks/pilot_study_peaks/Espinosa_csem_sampling_peaks.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #MAC2 ChIP-Seq peak file
HCT116_peaks$neg_LOG10_pvalue<-10^(-HCT116_peaks$neg_LOG10_pvalue) #convert neg_LOG10_pvalue to regular pvalue
HCT116_peaks$neg_LOG10_qvalue<-10^(-HCT116_peaks$neg_LOG10_qvalue) #convert neg_LOG10_qvalue to regular qvalue

HCT116_peaks<-subset(HCT116_peaks, neg_LOG10_qvalue<0.05)
HCT116_peaks$name<-paste("peak",row.names(HCT116_peaks), sep="_")

IMR90_peaks<-read.table("RE_peaks/pilot_study_peaks/IMR90_csem_sampling_peaks.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #MAC2 ChIP-Seq peak file
IMR90_peaks$neg_LOG10_pvalue<-10^(-IMR90_peaks$neg_LOG10_pvalue) #convert neg_LOG10_pvalue to regular pvalue
IMR90_peaks$neg_LOG10_qvalue<-10^(-IMR90_peaks$neg_LOG10_qvalue) #convert neg_LOG10_qvalue to regular qvalue

IMR90_peaks<-subset(IMR90_peaks, neg_LOG10_qvalue<0.05)
IMR90_peaks$name<-paste("peak",row.names(IMR90_peaks), sep="_")

#create GenomicRanges objects from the various tables that have been read into script.
HCT116_peaks_gr<-GRanges(seqnames = Rle(HCT116_peaks$chr), ranges=IRanges(start=(HCT116_peaks$abs_summit-2000), end = (HCT116_peaks$abs_summit+2000), names = HCT116_peaks$name), strand = NULL)
IMR90_peaks_gr<-GRanges(seqnames = Rle(IMR90_peaks$chr), ranges=IRanges(start=(IMR90_peaks$abs_summit-2000), end = (IMR90_peaks$abs_summit+2000), names = IMR90_peaks$name), strand = NULL)

#overlap with REs and then create GRanges objects.

#HCT116 FC 0
i<-1
REs<-c()
chrs<-c()
starts<-c()
ends<-c()

while (nrow(HCT116_FC_0)>i-1)
{
  RE<-HCT116_FC_0[i,1]
  REs<-c(REs,RE)
  
  temp<-strsplit(RE,"_")[[1]]
  
  chr<-temp[grep("chr",temp):(length(temp)-2)]
  chr<-paste(chr,collapse = "_")
  chrs<-c(chrs,chr)
  
  starts<-c(starts,as.numeric(temp[(length(temp)-1)]))
  ends<-c(ends,as.numeric(temp[length(temp)]))
  
  
  
  i<-i+1
}
temp<-data.frame(RE=REs,chr=chrs,start=starts,end=ends)
HCT116_FC_0<-merge(temp,HCT116_FC_0, by="RE")

HCT116_FC_0_gr<-GRanges(seqnames = Rle(HCT116_FC_0$chr), ranges=IRanges(start=HCT116_FC_0$start, end = HCT116_FC_0$end, names = HCT116_FC_0$RE), strand = NULL) #create GenomicRanges objects from the various tables that have been read into script.

#HCT116 Fc 1
i<-1
REs<-c()
chrs<-c()
starts<-c()
ends<-c()

while (nrow(HCT116_FC_1)>i-1)
{
  RE<-HCT116_FC_1[i,1]
  REs<-c(REs,RE)
  
  temp<-strsplit(RE,"_")[[1]]
  
  chr<-temp[grep("chr",temp):(length(temp)-2)]
  chr<-paste(chr,collapse = "_")
  chrs<-c(chrs,chr)
  
  starts<-c(starts,as.numeric(temp[(length(temp)-1)]))
  ends<-c(ends,as.numeric(temp[length(temp)]))
  
  
  
  i<-i+1
}
temp<-data.frame(RE=REs,chr=chrs,start=starts,end=ends)
HCT116_FC_1<-merge(temp,HCT116_FC_1, by="RE")

HCT116_FC_1_gr<-GRanges(seqnames = Rle(HCT116_FC_1$chr), ranges=IRanges(start=HCT116_FC_1$start, end = HCT116_FC_1$end, names = HCT116_FC_1$RE), strand = NULL) #create GenomicRanges objects from the various tables that have been read into script.

#IMR90 FC 0

i<-1
REs<-c()
chrs<-c()
starts<-c()
ends<-c()

while (nrow(IMR90_FC_0)>i-1)
{
  RE<-IMR90_FC_0[i,1]
  REs<-c(REs,RE)
  
  temp<-strsplit(RE,"_")[[1]]
  
  chr<-temp[grep("chr",temp):(length(temp)-2)]
  chr<-paste(chr,collapse = "_")
  chrs<-c(chrs,chr)
  
  starts<-c(starts,as.numeric(temp[(length(temp)-1)]))
  ends<-c(ends,as.numeric(temp[length(temp)]))
  
  
  
  i<-i+1
}
temp<-data.frame(RE=REs,chr=chrs,start=starts,end=ends)
IMR90_FC_0<-merge(temp,IMR90_FC_0, by="RE")

IMR90_FC_0_gr<-GRanges(seqnames = Rle(IMR90_FC_0$chr), ranges=IRanges(start=IMR90_FC_0$start, end = IMR90_FC_0$end, names = IMR90_FC_0$RE), strand = NULL) #create GenomicRanges objects from the various tables that have been read into script.

#IMR90 FC 1

i<-1
REs<-c()
chrs<-c()
starts<-c()
ends<-c()

while (nrow(IMR90_FC_1)>i-1)
{
  RE<-IMR90_FC_1[i,1]
  REs<-c(REs,RE)
  
  temp<-strsplit(RE,"_")[[1]]
  
  chr<-temp[grep("chr",temp):(length(temp)-2)]
  chr<-paste(chr,collapse = "_")
  chrs<-c(chrs,chr)
  
  starts<-c(starts,as.numeric(temp[(length(temp)-1)]))
  ends<-c(ends,as.numeric(temp[length(temp)]))
  
  
  
  i<-i+1
}
temp<-data.frame(RE=REs,chr=chrs,start=starts,end=ends)
IMR90_FC_1<-merge(temp,IMR90_FC_1, by="RE")

IMR90_FC_1_gr<-GRanges(seqnames = Rle(IMR90_FC_1$chr), ranges=IRanges(start=IMR90_FC_1$start, end = IMR90_FC_1$end, names = IMR90_FC_1$RE), strand = NULL) #create GenomicRanges objects from the various tables that have been read into script.



#overlaps
#HCT116 FC 0
hits <- findOverlaps(HCT116_FC_0_gr, HCT116_peaks_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <-HCT116_FC_0_gr[queryHits(hits)] #where the query hits overlap 
overlaps_chip<-HCT116_peaks_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("chip"=overlaps_chip) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_chip))

temp1<-cbind(temp1, names(overlaps_RE))
temp1<-temp1[,c("names(overlaps_chip)","names(overlaps_RE)")]
names(temp1)<-c("chipID", "RE")
temp2<-merge(HCT116_peaks,temp1,by.x="name", by.y="chipID")
temp2<-temp2[,c("name", "RE", "chr","start","end","abs_summit")]
names(temp2)<-c("peak_ID", "RE", "chr","start","end","abs_summit")
temp2<-merge(temp2,HCT116_FC_0,by.x="RE", by.y="RE")
temp2<-temp2[,c("RE","peak_ID","chr.x","start.x","end.x","abs_summit","start.y","end.y","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")]
names(temp2)<-c("RE","peak_ID","chr","peak_start","peak_end","abs_summit","RE_start","RE_end","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")

dupes<-temp2[duplicated(temp2$RE),]

temp2<-subset(temp2, !(temp2$RE %in% dupes$RE))

write.table(temp2, "comparisons/unique_RE_HcT116_5-FU_only_FC_0.csv", sep=',', row.names = F, col.names = T, quote = F) #write out 2kb TE overlaps

#HCT116 FC 1
hits <- findOverlaps(HCT116_FC_1_gr, HCT116_peaks_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <-HCT116_FC_1_gr[queryHits(hits)] #where the query hits overlap 
overlaps_chip<-HCT116_peaks_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("chip"=overlaps_chip) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_chip))

temp1<-cbind(temp1, names(overlaps_RE))
temp1<-temp1[,c("names(overlaps_chip)","names(overlaps_RE)")]
names(temp1)<-c("chipID", "RE")
temp2<-merge(HCT116_peaks,temp1,by.x="name", by.y="chipID")
temp2<-temp2[,c("name", "RE", "chr","start","end","abs_summit")]
names(temp2)<-c("peak_ID", "RE", "chr","start","end","abs_summit")
temp2<-merge(temp2,HCT116_FC_1,by.x="RE", by.y="RE")
temp2<-temp2[,c("RE","peak_ID","chr.x","start.x","end.x","abs_summit","start.y","end.y","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")]
names(temp2)<-c("RE","peak_ID","chr","peak_start","peak_end","abs_summit","RE_start","RE_end","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")

dupes<-temp2[duplicated(temp2$RE),]

temp2<-subset(temp2, !(temp2$RE %in% dupes$RE))

write.table(temp2, "comparisons/unique_RE_HcT116_5-FU_only_FC_1.csv", sep=',', row.names = F, col.names = T, quote = F) #write out 2kb TE overlaps

#IMR90 FC 0
hits <- findOverlaps(IMR90_FC_0_gr, IMR90_peaks_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <-IMR90_FC_0_gr[queryHits(hits)] #where the query hits overlap 
overlaps_chip<-IMR90_peaks_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("chip"=overlaps_chip) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_chip))

temp1<-cbind(temp1, names(overlaps_RE))
temp1<-temp1[,c("names(overlaps_chip)","names(overlaps_RE)")]
names(temp1)<-c("chipID", "RE")
temp2<-merge(IMR90_peaks,temp1,by.x="name", by.y="chipID")
temp2<-temp2[,c("name", "RE", "chr","start","end","abs_summit")]
names(temp2)<-c("peak_ID", "RE", "chr","start","end","abs_summit")
temp2<-merge(temp2,IMR90_FC_0,by.x="RE", by.y="RE")
temp2<-temp2[,c("RE","peak_ID","chr.x","start.x","end.x","abs_summit","start.y","end.y","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")]
names(temp2)<-c("RE","peak_ID","chr","peak_start","peak_end","abs_summit","RE_start","RE_end","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")

dupes<-temp2[duplicated(temp2$RE),]

temp2<-subset(temp2, !(temp2$RE %in% dupes$RE))

write.table(temp2, "comparisons/unique_RE_IMR90_5-FU_only_FC_0.csv", sep=',', row.names = F, col.names = T, quote = F) #write out 2kb TE overlaps

#IMR90 FC 1

hits <- findOverlaps(IMR90_FC_1_gr, IMR90_peaks_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <-IMR90_FC_1_gr[queryHits(hits)] #where the query hits overlap 
overlaps_chip<-IMR90_peaks_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("chip"=overlaps_chip) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_chip))

temp1<-cbind(temp1, names(overlaps_RE))
temp1<-temp1[,c("names(overlaps_chip)","names(overlaps_RE)")]
names(temp1)<-c("chipID", "RE")
temp2<-merge(IMR90_peaks,temp1,by.x="name", by.y="chipID")
temp2<-temp2[,c("name", "RE", "chr","start","end","abs_summit")]
names(temp2)<-c("peak_ID", "RE", "chr","start","end","abs_summit")
temp2<-merge(temp2,IMR90_FC_1,by.x="RE", by.y="RE")
temp2<-temp2[,c("RE","peak_ID","chr.x","start.x","end.x","abs_summit","start.y","end.y","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")]
names(temp2)<-c("RE","peak_ID","chr","peak_start","peak_end","abs_summit","RE_start","RE_end","con_1",
                "con_2","hr12_DMSO_1","hr12_DMSO_2","hr12_5_FU_1","hr12_5_FU_2","con_padj","DMSO_padj")

dupes<-temp2[duplicated(temp2$RE),]

temp2<-subset(temp2, !(temp2$RE %in% dupes$RE))

write.table(temp2, "comparisons/unique_RE_IMR90_5-FU_only_FC_1.csv", sep=',', row.names = F, col.names = T, quote = F) #write out 2kb TE overlaps