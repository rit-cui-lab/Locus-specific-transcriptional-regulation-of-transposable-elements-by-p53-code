#Multiple number of scripts resulting in the FC_0_occur_2_RE_overlaps_classes.tsv file.

#Associate UTR file and full_gene file with various occur files using genomic means.
library("GenomicRanges") #overlaps via the genomic location

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope")

UTR<-read.table("overlaps/UCSC_UTR.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #UTR data from USC genome browser
full_gene<-read.table("overlaps/UCSC_exons_introns.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #Genomic data from USC genome browser

# occurrence 2

FC_0_occur_2<-read.table("data_sets_locs/padj_0.05_FC_0_occur_2.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_1_occur_2<-read.table("data_sets_locs/padj_0.05_FC_1_occur_2.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_2_occur_2<-read.table("data_sets_locs/padj_0.05_FC_2_occur_2.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

# occurrence 3
FC_0_occur_3<-read.table("data_sets_locs/padj_0.05_FC_0_occur_3.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_1_occur_3<-read.table("data_sets_locs/padj_0.05_FC_1_occur_3.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_2_occur_3<-read.table("data_sets_locs/padj_0.05_FC_2_occur_3.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

#create GenomicRanges objects from the various tables that have been read into script.
UTR_gr<-GRanges(seqnames = Rle(UTR$chr), ranges=IRanges(start=UTR$start, end = UTR$end, names = UTR$UTR_name), strand = UTR$strand)
gene_gr<-GRanges(seqnames = Rle(full_gene$chr), ranges=IRanges(start=full_gene$start, end = full_gene$end, names = full_gene$ID), strand = full_gene$strand)


FC_0_occur_2_gr<-GRanges(seqnames = Rle(FC_0_occur_2$chr), ranges=IRanges(start=FC_0_occur_2$start, end = FC_0_occur_2$end, names = FC_0_occur_2$RE), strand = NULL)
FC_1_occur_2_gr<-GRanges(seqnames = Rle(FC_1_occur_2$chr), ranges=IRanges(start=FC_1_occur_2$start, end = FC_1_occur_2$end, names = FC_1_occur_2$RE), strand = NULL)
FC_2_occur_2_gr<-GRanges(seqnames = Rle(FC_2_occur_2$chr), ranges=IRanges(start=FC_2_occur_2$start, end = FC_2_occur_2$end, names = FC_2_occur_2$RE), strand = NULL)

FC_0_occur_3_gr<-GRanges(seqnames = Rle(FC_0_occur_3$chr), ranges=IRanges(start=FC_0_occur_3$start, end = FC_0_occur_3$end, names = FC_0_occur_3$RE), strand = NULL)
FC_1_occur_3_gr<-GRanges(seqnames = Rle(FC_1_occur_3$chr), ranges=IRanges(start=FC_1_occur_3$start, end = FC_1_occur_3$end, names = FC_1_occur_3$RE), strand = NULL)
FC_2_occur_3_gr<-GRanges(seqnames = Rle(FC_2_occur_3$chr), ranges=IRanges(start=FC_2_occur_3$start, end = FC_2_occur_3$end, names = FC_2_occur_3$RE), strand = NULL)


#FC 0 occurrence 2
hits <- findOverlaps(FC_0_occur_2_gr, UTR_gr) #Find the overlaps between two GenomicRanges objects
overlaps_FC_0_occur_2 <- FC_0_occur_2_gr[queryHits(hits)] #where the query hits overlap 
overlaps_UTR<-UTR_gr[subjectHits(hits)] #where the subject hits overlap


test_list<-GRangesList("mine"=overlaps_FC_0_occur_2) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_FC_0_occur_2))
temp1<-cbind(temp1, names(overlaps_UTR))
temp1<-temp1[,c("names(overlaps_FC_0_occur_2)","names(overlaps_UTR)")]
names(temp1)<-c("RE", "UTR_ID")
temp1<-merge(temp1,FC_0_occur_2, by="RE")
temp1<-merge(temp1,UTR, by.x="UTR_ID", by.y="UTR_name")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances")]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_0_occur_2_UTR_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_UTR<-temp1[,c("UTR_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_UTR)<-c("UTR_ID","UCSC_ID","gene_name","RE","chr","start_UTR","end_UTR","strand"  )
write.table(temp_UTR, "overlaps/FC_0_occur_2_UTR_RES_UTR_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#exon and intron
RE_temp<-subset(FC_0_occur_2, !(FC_0_occur_2$RE %in% temp_FC$RE))

#create GenomicRanges objects from the various tables that have been read into script.
RE_gr<-GRanges(seqnames = Rle(RE_temp$chr), ranges=IRanges(start=RE_temp$start, end = RE_temp$end, names = RE_temp$RE), strand = RE_temp$strand)

hits <- findOverlaps(RE_gr, gene_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <- RE_gr[queryHits(hits)] #where the query hits overlap
overlaps_gene<-gene_gr[subjectHits(hits)] #where the query hits overlap

test_list<-GRangesList("mine"=overlaps_RE) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_RE))
temp1<-cbind(temp1, names(overlaps_gene))
temp1<-temp1[,c("names(overlaps_RE)","names(overlaps_gene)")]
names(temp1)<-c("RE", "gene_ID")
temp1<-merge(temp1,RE_temp, by="RE")
temp1<-merge(temp1,full_gene, by.x="gene_ID", by.y="ID")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances" )]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_0_occur_2_intron_exon_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_gene<-temp1[,c("gene_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_gene)<-c("gene_ID","UCSC_ID","gene_name","RE","chr","start_gene","end_gene","strand"  )
write.table(temp_gene, "overlaps/FC_0_occur_2_genes_RES_intron_exon_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#no overlap

no_overlap<-subset(RE_temp, !(RE_temp$RE %in% temp_FC$RE))
no_overlap<-cbind(no_overlap, "gene_name"=NA)
no_overlap<-cbind(no_overlap, "strand"=NA)

write.table(no_overlap, "overlaps/FC_0_occur_2_no_overlap_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

#FC 0 occurrence 3
hits <- findOverlaps(FC_0_occur_3_gr, UTR_gr) #Find the overlaps between two GenomicRanges objects
overlaps_FC_0_occur_3 <- FC_0_occur_3_gr[queryHits(hits)] #where the query hits overlap
overlaps_UTR<-UTR_gr[subjectHits(hits)] #where the subject hits overlap


test_list<-GRangesList("mine"=overlaps_FC_0_occur_3) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_FC_0_occur_3))
temp1<-cbind(temp1, names(overlaps_UTR))
temp1<-temp1[,c("names(overlaps_FC_0_occur_3)","names(overlaps_UTR)")]
names(temp1)<-c("RE", "UTR_ID")
temp1<-merge(temp1,FC_0_occur_3, by="RE")
temp1<-merge(temp1,UTR, by.x="UTR_ID", by.y="UTR_name")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances")]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_0_occur_3_UTR_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_UTR<-temp1[,c("UTR_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_UTR)<-c("UTR_ID","UCSC_ID","gene_name","RE","chr","start_UTR","end_UTR","strand"  )
write.table(temp_UTR, "overlaps/FC_0_occur_3_UTR_RES_UTR_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#exon and intron
RE_temp<-subset(FC_0_occur_3, !(FC_0_occur_3$RE %in% temp_FC$RE))

RE_gr<-GRanges(seqnames = Rle(RE_temp$chr), ranges=IRanges(start=RE_temp$start, end = RE_temp$end, names = RE_temp$RE), strand = RE_temp$strand) #create GenomicRanges objects from the various tables that have been read into script.

hits <- findOverlaps(RE_gr, gene_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <- RE_gr[queryHits(hits)] #where the query hits overlap
overlaps_gene<-gene_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("mine"=overlaps_RE) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_RE))
temp1<-cbind(temp1, names(overlaps_gene))
temp1<-temp1[,c("names(overlaps_RE)","names(overlaps_gene)")]
names(temp1)<-c("RE", "gene_ID")
temp1<-merge(temp1,RE_temp, by="RE")
temp1<-merge(temp1,full_gene, by.x="gene_ID", by.y="ID")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances" )]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_0_occur_3_intron_exon_REs.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

temp_gene<-temp1[,c("gene_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_gene)<-c("gene_ID","UCSC_ID","gene_name","RE","chr","start_gene","end_gene","strand"  )
write.table(temp_gene, "overlaps/FC_0_occur_3_genes_RES_intron_exon_info.csv", sep=',', row.names = F, col.names = T, quote = F)

#no overlap

no_overlap<-subset(RE_temp, !(RE_temp$RE %in% temp_FC$RE))
no_overlap<-cbind(no_overlap, "gene_name"=NA)
no_overlap<-cbind(no_overlap, "strand"=NA)

write.table(no_overlap, "overlaps/FC_0_occur_3_no_overlap_REs.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#FC 1 occurence 2
hits <- findOverlaps(FC_1_occur_2_gr, UTR_gr) #Find the overlaps between two GenomicRanges objects
overlaps_FC_1_occur_2 <- FC_1_occur_2_gr[queryHits(hits)] #where the query hits overlap
overlaps_UTR<-UTR_gr[subjectHits(hits)] #where the subject hits overlap


test_list<-GRangesList("mine"=overlaps_FC_1_occur_2) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_FC_1_occur_2))
temp1<-cbind(temp1, names(overlaps_UTR))
temp1<-temp1[,c("names(overlaps_FC_1_occur_2)","names(overlaps_UTR)")]
names(temp1)<-c("RE", "UTR_ID")
temp1<-merge(temp1,FC_1_occur_2, by="RE")
temp1<-merge(temp1,UTR, by.x="UTR_ID", by.y="UTR_name")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances")]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_1_occur_2_UTR_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_UTR<-temp1[,c("UTR_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_UTR)<-c("UTR_ID","UCSC_ID","gene_name","RE","chr","start_UTR","end_UTR","strand"  )
write.table(temp_UTR, "overlaps/FC_1_occur_2_UTR_RES_UTR_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#exon and intron
RE_temp<-subset(FC_1_occur_2, !(FC_1_occur_2$RE %in% temp_FC$RE))

#create GenomicRanges objects from the various tables that have been read into script.
RE_gr<-GRanges(seqnames = Rle(RE_temp$chr), ranges=IRanges(start=RE_temp$start, end = RE_temp$end, names = RE_temp$RE), strand = RE_temp$strand)

hits <- findOverlaps(RE_gr, gene_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <- RE_gr[queryHits(hits)] #where the query hits overlap
overlaps_gene<-gene_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("mine"=overlaps_RE) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_RE))
temp1<-cbind(temp1, names(overlaps_gene))
temp1<-temp1[,c("names(overlaps_RE)","names(overlaps_gene)")]
names(temp1)<-c("RE", "gene_ID")
temp1<-merge(temp1,RE_temp, by="RE")
temp1<-merge(temp1,full_gene, by.x="gene_ID", by.y="ID")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances" )]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_1_occur_2_intron_exon_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_gene<-temp1[,c("gene_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_gene)<-c("gene_ID","UCSC_ID","gene_name","RE","chr","start_gene","end_gene","strand"  )
write.table(temp_gene, "overlaps/FC_1_occur_2_genes_RES_intron_exon_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#no overlap

no_overlap<-subset(RE_temp, !(RE_temp$RE %in% temp_FC$RE))
no_overlap<-cbind(no_overlap, "gene_name"=NA)
no_overlap<-cbind(no_overlap, "strand"=NA)

write.table(no_overlap, "overlaps/FC_1_occur_2_no_overlap_REs.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#FC 1 occurence 3
hits <- findOverlaps(FC_1_occur_3_gr, UTR_gr) #Find the overlaps between two GenomicRanges objects
overlaps_FC_1_occur_3 <- FC_1_occur_3_gr[queryHits(hits)] #where the query hits overlap
overlaps_UTR<-UTR_gr[subjectHits(hits)] #where the subject hits overlap


test_list<-GRangesList("mine"=overlaps_FC_1_occur_3) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_FC_1_occur_3))
temp1<-cbind(temp1, names(overlaps_UTR))
temp1<-temp1[,c("names(overlaps_FC_1_occur_3)","names(overlaps_UTR)")]
names(temp1)<-c("RE", "UTR_ID")
temp1<-merge(temp1,FC_1_occur_3, by="RE")
temp1<-merge(temp1,UTR, by.x="UTR_ID", by.y="UTR_name")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances")]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_1_occur_3_UTR_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_UTR<-temp1[,c("UTR_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_UTR)<-c("UTR_ID","UCSC_ID","gene_name","RE","chr","start_UTR","end_UTR","strand"  )
write.table(temp_UTR, "overlaps/FC_1_occur_3_UTR_RES_UTR_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#exon and intron
RE_temp<-subset(FC_1_occur_3, !(FC_1_occur_3$RE %in% temp_FC$RE))

RE_gr<-GRanges(seqnames = Rle(RE_temp$chr), ranges=IRanges(start=RE_temp$start, end = RE_temp$end, names = RE_temp$RE), strand = RE_temp$strand) #create GenomicRanges objects from the various tables that have been read into script.

hits <- findOverlaps(RE_gr, gene_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <- RE_gr[queryHits(hits)] #where the query hits overlap
overlaps_gene<-gene_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("mine"=overlaps_RE) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_RE))
temp1<-cbind(temp1, names(overlaps_gene))
temp1<-temp1[,c("names(overlaps_RE)","names(overlaps_gene)")]
names(temp1)<-c("RE", "gene_ID")
temp1<-merge(temp1,RE_temp, by="RE")
temp1<-merge(temp1,full_gene, by.x="gene_ID", by.y="ID")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances" )]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_1_occur_3_intron_exon_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_gene<-temp1[,c("gene_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_gene)<-c("gene_ID","UCSC_ID","gene_name","RE","chr","start_gene","end_gene","strand"  )
write.table(temp_gene, "overlaps/FC_1_occur_3_genes_RES_intron_exon_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#no overlap

no_overlap<-subset(RE_temp, !(RE_temp$RE %in% temp_FC$RE))
no_overlap<-cbind(no_overlap, "gene_name"=NA)
no_overlap<-cbind(no_overlap, "strand"=NA)

write.table(no_overlap, "overlaps/FC_1_occur_3_no_overlap_REs.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#FC 2 occurence 2
hits <- findOverlaps(FC_2_occur_2_gr, UTR_gr) #Find the overlaps between two GenomicRanges objects
overlaps_FC_2_occur_2 <- FC_2_occur_2_gr[queryHits(hits)] #where the query hits overlap
overlaps_UTR<-UTR_gr[subjectHits(hits)] #where the subject hits overlap


test_list<-GRangesList("mine"=overlaps_FC_2_occur_2) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_FC_2_occur_2))
temp1<-cbind(temp1, names(overlaps_UTR))
temp1<-temp1[,c("names(overlaps_FC_2_occur_2)","names(overlaps_UTR)")]
names(temp1)<-c("RE", "UTR_ID")
temp1<-merge(temp1,FC_2_occur_2, by="RE")
temp1<-merge(temp1,UTR, by.x="UTR_ID", by.y="UTR_name")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances")]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_2_occur_2_UTR_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_UTR<-temp1[,c("UTR_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_UTR)<-c("UTR_ID","UCSC_ID","gene_name","RE","chr","start_UTR","end_UTR","strand"  )
write.table(temp_UTR, "overlaps/FC_2_occur_2_UTR_RES_UTR_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#exon and intron
RE_temp<-subset(FC_2_occur_2, !(FC_2_occur_2$RE %in% temp_FC$RE))

RE_gr<-GRanges(seqnames = Rle(RE_temp$chr), ranges=IRanges(start=RE_temp$start, end = RE_temp$end, names = RE_temp$RE), strand = RE_temp$strand) #create GenomicRanges objects from the various tables that have been read into script.

hits <- findOverlaps(RE_gr, gene_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <- RE_gr[queryHits(hits)] #where the query hits overlap
overlaps_gene<-gene_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("mine"=overlaps_RE) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_RE))
temp1<-cbind(temp1, names(overlaps_gene))
temp1<-temp1[,c("names(overlaps_RE)","names(overlaps_gene)")]
names(temp1)<-c("RE", "gene_ID")
temp1<-merge(temp1,RE_temp, by="RE")
temp1<-merge(temp1,full_gene, by.x="gene_ID", by.y="ID")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances" )]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_2_occur_2_intron_exon_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_gene<-temp1[,c("gene_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_gene)<-c("gene_ID","UCSC_ID","gene_name","RE","chr","start_gene","end_gene","strand"  )
write.table(temp_gene, "overlaps/FC_2_occur_2_genes_RES_intron_exon_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#no overlap

no_overlap<-subset(RE_temp, !(RE_temp$RE %in% temp_FC$RE))
no_overlap<-cbind(no_overlap, "gene_name"=NA)
no_overlap<-cbind(no_overlap, "strand"=NA)

write.table(no_overlap, "overlaps/FC_2_occur_2_no_overlap_REs.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#FC 2 occurrence 3
hits <- findOverlaps(FC_2_occur_3_gr, UTR_gr) #Find the overlaps between two GenomicRanges objects
overlaps_FC_2_occur_3 <- FC_2_occur_3_gr[queryHits(hits)] #where the query hits overlap
overlaps_UTR<-UTR_gr[subjectHits(hits)] #where the subject hits overlap


test_list<-GRangesList("mine"=overlaps_FC_2_occur_3) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_FC_2_occur_3))
temp1<-cbind(temp1, names(overlaps_UTR))
temp1<-temp1[,c("names(overlaps_FC_2_occur_3)","names(overlaps_UTR)")]
names(temp1)<-c("RE", "UTR_ID")
temp1<-merge(temp1,FC_2_occur_3, by="RE")
temp1<-merge(temp1,UTR, by.x="UTR_ID", by.y="UTR_name")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances")]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_2_occur_3_UTR_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_UTR<-temp1[,c("UTR_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_UTR)<-c("UTR_ID","UCSC_ID","gene_name","RE","chr","start_UTR","end_UTR","strand"  )
write.table(temp_UTR, "overlaps/FC_2_occur_3_UTR_RES_UTR_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data


#exon and intron
RE_temp<-subset(FC_2_occur_3, !(FC_2_occur_3$RE %in% temp_FC$RE))

RE_gr<-GRanges(seqnames = Rle(RE_temp$chr), ranges=IRanges(start=RE_temp$start, end = RE_temp$end, names = RE_temp$RE), strand = RE_temp$strand) #create GenomicRanges objects from the various tables that have been read into script.

hits <- findOverlaps(RE_gr, gene_gr) #Find the overlaps between two GenomicRanges objects
overlaps_RE <- RE_gr[queryHits(hits)] #where the query hits overlap
overlaps_gene<-gene_gr[subjectHits(hits)] #where the subject hits overlap

test_list<-GRangesList("mine"=overlaps_RE) #get a list of the locations of overlaps to turn into  data.frame.
temp1<-as.data.frame(test_list)
temp1<-cbind(temp1,names(overlaps_RE))
temp1<-cbind(temp1, names(overlaps_gene))
temp1<-temp1[,c("names(overlaps_RE)","names(overlaps_gene)")]
names(temp1)<-c("RE", "gene_ID")
temp1<-merge(temp1,RE_temp, by="RE")
temp1<-merge(temp1,full_gene, by.x="gene_ID", by.y="ID")

temp_FC<-temp1[,c("RE","gene_name","chr.x","start.x","end.x","strand", "RE_occurances" )]
temp_FC<-unique(temp_FC)
names(temp_FC)<-c("RE","gene_name","chr","start","end","strand", "RE_occurances")

write.table(temp_FC, "overlaps/FC_2_occur_3_intron_exon_REs.csv", sep=',', row.names = F, col.names = T, quote = F)

temp_gene<-temp1[,c("gene_ID","UCSC_ID", "gene_name","RE","chr.y","start.y","end.y","strand")]
names(temp_gene)<-c("gene_ID","UCSC_ID","gene_name","RE","chr","start_gene","end_gene","strand"  )
write.table(temp_gene, "overlaps/FC_2_occur_3_genes_RES_intron_exon_info.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#no overlap

no_overlap<-subset(RE_temp, !(RE_temp$RE %in% temp_FC$RE))
no_overlap<-cbind(no_overlap, "gene_name"=NA)
no_overlap<-cbind(no_overlap, "strand"=NA)

write.table(no_overlap, "overlaps/FC_2_occur_3_no_overlap_REs.csv", sep=',', row.names = F, col.names = T, quote = F) #write overlap data

#updated pie chart data script
library("GenomicRanges")
library(tidyr)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/overlaps")
files<-read.table("paired_files.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

UTR<-files[grep("UTR_REs",files$data),]
gene<-files[grep("intron_exon",files$data),]
intergenic<-files[grep("no_overlap",files$data),]

i<-1
file_names<-c()
while (nrow(UTR)>i-1)
{
  temp<-strsplit(UTR$data[i],"_")[[1]]
  num<-grep("occur",temp)
  temp_name<-paste(temp[1:(num+1)],collapse="_")
  temp_name<-paste(temp_name,"overlap_percentages",sep="_")
  
  if (grepl("up|down",temp[length(temp)]))
  {
    temp_name<-paste(temp_name,temp[length(temp)], sep="_")
    
    
  }
  
  else
  {
    temp_name<-paste(temp_name,"csv",sep=".")
    
  }
  
  file_names<-c(file_names,temp_name)
  
  i<-i+1
}


UTR<-cbind(names=file_names,UTR,stringsAsFactors=FALSE)
gene<-cbind(names=file_names,gene,stringsAsFactors=FALSE)
intergenic<-cbind(names=file_names,intergenic,stringsAsFactors=FALSE)

file_names<-merge(UTR,gene, by="names")
names(file_names)<-c("names","data_UTR","info_UTR","data_gene","info_gene")
file_names<-merge(file_names,intergenic[,c("names","data")], by="names")
names(file_names)<-c("names","data_UTR","info_UTR","data_gene","info_gene", "data_intergenic")


i<-1
while (nrow(file_names)>i-1)
{
  
  UTR_data<-read.table(file_names$data_UTR[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  UTR_info<-read.table(file_names$info_UTR[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  gene_data<-read.table(file_names$data_gene[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  gene_info<-read.table(file_names$info_gene[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  intergenic_data<-read.table(file_names$data_intergenic[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  #UTR
  if (nrow(UTR_data)!=0)
  {
    
    
    UTR_info<-separate(data = UTR_info, col = UTR_ID, into = c("ID","type_prime","prime", "UTR","UTR_num"), sep = "_")
    UTR_info<-cbind(UTR_info, "UTR_type"= paste(UTR_info$type_prime, "prime","UTR", sep="_"))
    UTR_info<-UTR_info[,c("RE","UTR_type")]
    UTR_info<-unique(UTR_info)
    
    REs<-unique(UTR_data$RE)
    
    prime_RE<-c()
    primes<-c()
    
    for (temp_RE in REs)
    {
      
      temp_info<-subset(UTR_info, RE==temp_RE)
      
      temp_prime<-paste(temp_info$UTR_type, collapse=", ")
      prime_RE<-c(prime_RE,temp_RE)
      primes<-c(primes,temp_prime)
      
    }
    
    temp<-data.frame(RE=prime_RE,overlap=primes)
    
    UTR_data<-merge(UTR_data,temp, by="RE")
    UTR_data<-UTR_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    
    
  }
  
 
  #gene
  
  if (nrow(gene_data)!=0)
  {
    
    gene_info<-separate(data = gene_info, col = gene_ID, into = c("ID","gene_type","num_type_gene"), sep = "_")
    gene_info<-gene_info[,c("RE","gene_type")]
    gene_info<-unique(gene_info)
    
    
    
    #FC_0
    REs<-unique(gene_data$RE)
    
    prime_RE<-c()
    primes<-c()
    
    for (temp_RE in REs)
    {
      
      temp_info<-subset(gene_info, RE==temp_RE)
      
      temp_prime<-paste(temp_info$gene_type, collapse=", ")
      prime_RE<-c(prime_RE,temp_RE)
      primes<-c(primes,temp_prime)
      
    }
    
    temp<-data.frame(RE=prime_RE,overlap=primes)
    
    gene_data<-merge(gene_data,temp, by="RE")
    gene_data<-gene_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    
  }
  
 
  
  
  
  #no overlap
  if (nrow(intergenic_data)!=0)
  {
    intergenic_data<-cbind(intergenic_data, "overlap"="Intergenic")
    intergenic_data<-intergenic_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
  }
  
  
  
  
  #combine
  
  #check for 
  temp_data<-rbind(UTR_data, gene_data)
  temp_data<-rbind(temp_data,intergenic_data)
  temp_data<-temp_data[order(temp_data$RE_occurances, decreasing=TRUE),]
  
  
  overlaps<-c("UTR", "intron","exon","Intergenic")
  
  
  overlap_cat<-c()
  num<-c()
  
  j<-1
  while (length(overlaps)>j-1)
  {
    overlap<-overlaps[j]
    overlap_cat<-c(overlap_cat,overlap)
    
    temp<-temp_data[grep(overlap,temp_data$overlap),]
    num<-c(num,nrow(temp))
    
   
    j<-j+1
  }
  
  
  per<-(num/nrow(temp_data))*100
  per<-round(per,1)

  
  temp<-data.frame(overlaps=overlap_cat, counts=num, percents=per)
  
  print(file_names$names[i])
  write.table(temp, file_names$names[i], sep=',', row.names = F, col.names = T, quote = F)
  
  
  
  i<-i+1
}


#Create and update the class and genomic data and write out.
library("GenomicRanges")
library(tidyr)

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/overlaps")
files<-read.table("paired_files.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #File that show files that are paired together for the anyalisis.

UTR<-files[grep("UTR_REs",files$data),]
gene<-files[grep("intron_exon",files$data),]
intergenic<-files[grep("no_overlap",files$data),]

i<-1
file_names<-c()
while (nrow(UTR)>i-1)
{
  temp<-strsplit(UTR$data[i],"_")[[1]]
  num<-grep("occur",temp)
  temp_name<-paste(temp[1:(num+1)],collapse="_")
  temp_name<-paste(temp_name,"overlap_percentages",sep="_")
  
  if (grepl("up|down",temp[length(temp)]))
  {
    temp_name<-paste(temp_name,temp[length(temp)], sep="_")
    
    
  }
  
  else
  {
    temp_name<-paste(temp_name,"csv",sep=".")
    
  }
  
  file_names<-c(file_names,temp_name)
  
  i<-i+1
}


UTR<-cbind(names=file_names,UTR,stringsAsFactors=FALSE)
gene<-cbind(names=file_names,gene,stringsAsFactors=FALSE)
intergenic<-cbind(names=file_names,intergenic,stringsAsFactors=FALSE)

file_names<-merge(UTR,gene, by="names")
names(file_names)<-c("names","data_UTR","info_UTR","data_gene","info_gene")
file_names<-merge(file_names,intergenic[,c("names","data")], by="names")
names(file_names)<-c("names","data_UTR","info_UTR","data_gene","info_gene", "data_intergenic")


i<-1
while (nrow(file_names)>i-1)
{
  
  UTR_data<-read.table(file_names$data_UTR[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  UTR_info<-read.table(file_names$info_UTR[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  gene_data<-read.table(file_names$data_gene[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  gene_info<-read.table(file_names$info_gene[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  intergenic_data<-read.table(file_names$data_intergenic[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  #UTR
  if (nrow(UTR_data)!=0)
  {
    
    
    UTR_info<-separate(data = UTR_info, col = UTR_ID, into = c("ID","type_prime","prime", "UTR","UTR_num"), sep = "_")
    UTR_info<-cbind(UTR_info, "UTR_type"= paste(UTR_info$type_prime, "prime","UTR", sep="_"))
    UTR_info<-UTR_info[,c("RE","UTR_type")]
    UTR_info<-unique(UTR_info)
    
    REs<-unique(UTR_data$RE)
    
    prime_RE<-c()
    primes<-c()
    
    for (temp_RE in REs)
    {
      
      temp_info<-subset(UTR_info, RE==temp_RE)
      
      temp_prime<-paste(temp_info$UTR_type, collapse=", ")
      prime_RE<-c(prime_RE,temp_RE)
      primes<-c(primes,temp_prime)
      
    }
    
    temp<-data.frame(RE=prime_RE,overlap=primes)
    
    UTR_data<-merge(UTR_data,temp, by="RE")
    UTR_data<-UTR_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    
    
  }
  
  
  #gene
  
  if (nrow(gene_data)!=0)
  {
    
    gene_info<-separate(data = gene_info, col = gene_ID, into = c("ID","gene_type","num_type_gene"), sep = "_")
    gene_info<-gene_info[,c("RE","gene_type")]
    gene_info<-unique(gene_info)
    
    
    
    #FC_0
    REs<-unique(gene_data$RE)
    
    prime_RE<-c()
    primes<-c()
    
    for (temp_RE in REs)
    {
      
      temp_info<-subset(gene_info, RE==temp_RE)
      
      temp_prime<-paste(temp_info$gene_type, collapse=", ")
      prime_RE<-c(prime_RE,temp_RE)
      primes<-c(primes,temp_prime)
      
    }
    
    temp<-data.frame(RE=prime_RE,overlap=primes)
    
    gene_data<-merge(gene_data,temp, by="RE")
    gene_data<-gene_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    
  }
  
  
  
  
  
  #no overlap
  if (nrow(intergenic_data)!=0)
  {
    intergenic_data<-cbind(intergenic_data, "overlap"="Intergenic")
    intergenic_data<-intergenic_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
  }
  
  
  
  
  #combine
  
  #check for 
  temp_data<-rbind(UTR_data, gene_data)
  temp_data<-rbind(temp_data,intergenic_data)
  temp_data<-temp_data[order(temp_data$RE_occurances, decreasing=TRUE),]
  
  
  overlaps<-c("UTR","3_prime_UTR" ,"5_prime_UTR","intron","exon","Intergenic")
  
  
  overlap_cat<-c()
  num<-c()

    
    overlap_cat<-c("all_UTR")
    temp<-temp_data[grep("UTR",temp_data$overlap),]
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"3_prime")
    temp<-subset(temp_data, overlap=="3_prime_UTR")
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"5_prime")
    temp<-subset(temp_data, overlap=="5_prime_UTR")
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"3_prime_and_5_prime")
    temp<-subset(temp_data, overlap=="3_prime_UTR, 5_prime_UTR" | overlap=="5_prime_UTR, 3_prime_UTR")
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"all_intron")
    temp<-temp_data[grep("intron",temp_data$overlap),]
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"all_exon")
    temp<-temp_data[grep("exon",temp_data$overlap),]
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"only_intron")
    temp<-subset(temp_data, overlap=="intron")
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"only_exon")
    temp<-subset(temp_data, overlap=="exon")
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"intron_and_exon")
    temp<-subset(temp_data, overlap=="exon, intron" | overlap=="intron, exon")
    num<-c(num,nrow(temp))
    
    overlap_cat<-c(overlap_cat,"Intergenic")
    temp<-temp_data[grep("Intergenic",temp_data$overlap),]
    num<-c(num,nrow(temp))
    
  
  per<-(num/nrow(temp_data))
  per<-round(per,3)
  
  
  temp<-data.frame(overlaps=overlap_cat, counts=num, percents=per)
  
  temp_name<-paste("prime_percentages", file_names$names[i], sep="/")
  print(temp_name)
  write.table(temp, temp_name, sep=',', row.names = F, col.names = T, quote = F) #write out the percentages to different files.
  
  
  
  i<-i+1
}
#Use heatmap script to get file used further down in script.
library(gplots)
library(DESeq2)
library("edgeR")

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope")

#counts
#produce count table for TE RE for each telescope 
HCT116_Nutlin<-read.table("HCT116/HCT116_total_Nutlin_my_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) 
HCT116_Nutlin<-HCT116_Nutlin+1
HCT116_Nutlin <- as.data.frame(t(scale(t(HCT116_Nutlin), center=TRUE, scale=TRUE)))
HCT116_Nutlin<-HCT116_Nutlin[,c("Nutlin_1","Nutlin_2")]
names(HCT116_Nutlin)<-c("HCT116_Nutlin_1","HCT116_Nutlin_2")
HCT116_Nutlin_mean<-as.data.frame(rowMeans(HCT116_Nutlin))
names(HCT116_Nutlin_mean)<-c("HCT116_Nutlin")

HCT116_GSE137297_wt_Nutlin<-read.table("other_cell_lines/HCT116_GSE137297_wt_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
HCT116_GSE137297_wt_Nutlin<-HCT116_GSE137297_wt_Nutlin+1
HCT116_GSE137297_wt_Nutlin <- as.data.frame(t(scale(t(HCT116_GSE137297_wt_Nutlin), center=TRUE, scale=TRUE)))
HCT116_GSE137297_wt_Nutlin<-HCT116_GSE137297_wt_Nutlin[,c("Nutlin_1","Nutlin_2","Nutlin_3")]
names(HCT116_GSE137297_wt_Nutlin)<-c("HCT116_GSE137297_wt_Nutlin_1","HCT116_GSE137297_wt_Nutlin_2","HCT116_GSE137297_wt_Nutlin_3")
HCT116_GSE137297_wt_Nutlin_mean<-as.data.frame(rowMeans(HCT116_GSE137297_wt_Nutlin))
names(HCT116_GSE137297_wt_Nutlin_mean)<-c("HCT116_GSE137297_wt_Nutlin")

HCT116_5_FU<-read.table("time_course/HCT116/HCT116_hr12_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
HCT116_5_FU<-HCT116_5_FU+1
HCT116_5_FU<-as.data.frame(t(scale(t(HCT116_5_FU), center=TRUE, scale=TRUE)))
HCT116_5_FU<-HCT116_5_FU[,c("hr12_5_FU_1","hr12_5_FU_2")]
names(HCT116_5_FU)<-c("HCT116_5_FU_1","HCT116_5_FU_2")
HCT116_5_FU_mean<-as.data.frame(rowMeans(HCT116_5_FU))
names(HCT116_5_FU_mean)<-c("HCT116_5_FU")

IMR90_Nutlin<-read.table("IMR90/IMR90_total_Nutlin_my_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
IMR90_Nutlin<-IMR90_Nutlin+1
IMR90_Nutlin<-as.data.frame(t(scale(t(IMR90_Nutlin), center=TRUE, scale=TRUE)))
IMR90_Nutlin<-IMR90_Nutlin[,c("Nutlin_1","Nutlin_2")]
names(IMR90_Nutlin)<-c("IMR90_Nutlin_1","IMR90_Nutlin_2")
IMR90_Nutlin_mean<-as.data.frame(rowMeans(IMR90_Nutlin))
names(IMR90_Nutlin_mean)<-c("IMR90_Nutlin")

IMR90_5_FU<-read.table("time_course/IMR90/IMR90_hr12_5_FU_time_course_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
IMR90_5_FU<-IMR90_5_FU+1
IMR90_5_FU<-as.data.frame(t(scale(t(IMR90_5_FU), center=TRUE, scale=TRUE)))
IMR90_5_FU<-IMR90_5_FU[,c("hr12_5_FU_1","hr12_5_FU_2")]
names(IMR90_5_FU)<-c("IMR90_5_FU_1","IMR90_5_FU_2")
IMR90_5_FU_mean<-as.data.frame(rowMeans(IMR90_5_FU))
names(IMR90_5_FU_mean)<-c("IMR90_5_FU")

Saos_2_doxycyclin<-read.table("other_cell_lines/Saos-2_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
Saos_2_doxycyclin<-Saos_2_doxycyclin+1
Saos_2_doxycyclin<-as.data.frame(t(scale(t(Saos_2_doxycyclin), center=TRUE, scale=TRUE)))
Saos_2_doxycyclin<-Saos_2_doxycyclin[,c("dox_1","dox_2")]
names(Saos_2_doxycyclin)<-c("Saos_2_doxycyclin_1","Saos_2_doxycyclin_2")
Saos_2_doxycyclin_mean<-as.data.frame(rowMeans(Saos_2_doxycyclin))
names(Saos_2_doxycyclin_mean)<-c("Saos_2_doxycyclin")

GM06170_Doxorubicin<-read.table("other_cell_lines/GM06170_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
GM06170_Doxorubicin<-GM06170_Doxorubicin+1
GM06170_Doxorubicin<-as.data.frame(t(scale(t(GM06170_Doxorubicin), center=TRUE, scale=TRUE)))
GM06170_Doxorubicin<-GM06170_Doxorubicin[,c("dox_1","dox_2")]
names(GM06170_Doxorubicin)<-c("GM06170_Doxorubicin_1","GM06170_Doxorubicin_2")
GM06170_Doxorubicin_mean<-as.data.frame(rowMeans(GM06170_Doxorubicin))
names(GM06170_Doxorubicin_mean)<-c("GM06170_Doxorubicin")

GM00011_Doxorubicin<-read.table("other_cell_lines/GM00011_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
GM00011_Doxorubicin<-GM00011_Doxorubicin+1
GM00011_Doxorubicin<-as.data.frame(t(scale(t(GM00011_Doxorubicin), center=TRUE, scale=TRUE)))
GM00011_Doxorubicin<-GM00011_Doxorubicin[,c("dox_1","dox_2")]
names(GM00011_Doxorubicin)<-c("GM00011_Doxorubicin_1","GM00011_Doxorubicin_2")
GM00011_Doxorubicin_mean<-as.data.frame(rowMeans(GM00011_Doxorubicin))
names(GM00011_Doxorubicin_mean)<-c("GM00011_Doxorubicin")

MCF7_GSE47042_Nutlin<-read.table("other_cell_lines/MCF7_GSE47042_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF7_GSE47042_Nutlin<-MCF7_GSE47042_Nutlin+1
MCF7_GSE47042_Nutlin<-as.data.frame(t(scale(t(MCF7_GSE47042_Nutlin), center=TRUE, scale=TRUE)))
MCF7_GSE47042_Nutlin<-MCF7_GSE47042_Nutlin[,c("Nutlin_1","Nutlin_2")]
names(MCF7_GSE47042_Nutlin)<-c("MCF7_GSE47042_Nutlin_1","MCF7_GSE47042_Nutlin_2")
MCF7_GSE47042_Nutlin_mean<-as.data.frame(rowMeans(MCF7_GSE47042_Nutlin))
names(MCF7_GSE47042_Nutlin_mean)<-c("MCF7_GSE47042_Nutlin")

Fibro_Nutlin<-read.table("other_cell_lines/Fibro_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
Fibro_Nutlin<-Fibro_Nutlin+1
Fibro_Nutlin<-as.data.frame(t(scale(t(Fibro_Nutlin), center=TRUE, scale=TRUE)))
Fibro_Nutlin<-Fibro_Nutlin[,c("Nutlin_1","Nutlin_2", "Nutlin_3")]
names(Fibro_Nutlin)<-c("Fibro_Nutlin_1","Fibro_Nutlin_2", "Fibro_Nutlin_3")
Fibro_Nutlin_mean<-as.data.frame(rowMeans(Fibro_Nutlin))
names(Fibro_Nutlin_mean)<-c("Fibro_Nutlin")


MCF10A_Nutlin<-read.table("other_cell_lines/MCF10A_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF10A_Nutlin<-MCF10A_Nutlin+1
MCF10A_Nutlin<-as.data.frame(t(scale(t(MCF10A_Nutlin), center=TRUE, scale=TRUE)))
MCF10A_Nutlin<-MCF10A_Nutlin[,c("Nut_1","Nut_2", "Nut_3")]
names(MCF10A_Nutlin)<-c("MCF10A_Nutlin_1","MCF10A_Nutlin_2", "MCF10A_Nutlin_3")
MCF10A_Nutlin_mean<-as.data.frame(rowMeans(MCF10A_Nutlin))
names(MCF10A_Nutlin_mean)<-c("MCF10A_Nutlin")

SJSA_Nutlin<-read.table("other_cell_lines/SJSA_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
SJSA_Nutlin<-SJSA_Nutlin+1
SJSA_Nutlin<-as.data.frame(t(scale(t(SJSA_Nutlin), center=TRUE, scale=TRUE)))
SJSA_Nutlin<-SJSA_Nutlin[,c("Nutlin_1","Nutlin_2")]
names(SJSA_Nutlin)<-c("SJSA_Nutlin_1","SJSA_Nutlin_2")
SJSA_Nutlin_mean<-as.data.frame(rowMeans(SJSA_Nutlin))
names(SJSA_Nutlin_mean)<-c("SJSA_Nutlin")


MV4_11_Nutlin<-read.table("other_cell_lines/MV4_11_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MV4_11_Nutlin<-MV4_11_Nutlin+1
MV4_11_Nutlin<-as.data.frame(t(scale(t(MV4_11_Nutlin), center=TRUE, scale=TRUE)))
MV4_11_Nutlin<-MV4_11_Nutlin[,c("Nut_1","Nut_2","Nut_3")]
names(MV4_11_Nutlin)<-c("MV4_11_Nutlin_1","MV4_11_Nutlin_2","MV4_11_Nutlin_3")
MV4_11_Nutlin_mean<-as.data.frame(rowMeans(MV4_11_Nutlin))
names(MV4_11_Nutlin_mean)<-c("MV4_11_Nutlin")

PAEC_Nutlin<-read.table("other_cell_lines/PAEC_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
PAEC_Nutlin<-PAEC_Nutlin+1
PAEC_Nutlin<-as.data.frame(t(scale(t(PAEC_Nutlin), center=TRUE, scale=TRUE)))
PAEC_Nutlin<-PAEC_Nutlin[,c("Nut_1","Nut_2","Nut_3")]
names(PAEC_Nutlin)<-c("PAEC_Nutlin_1","PAEC_Nutlin_2","PAEC_Nutlin_3")
PAEC_Nutlin_mean<-as.data.frame(rowMeans(PAEC_Nutlin))
names(PAEC_Nutlin_mean)<-c("PAEC_Nutlin")

IMR90_GSE139003_hr12_Nutlin<-read.table("other_cell_lines/IMR90_GSE139003_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
IMR90_GSE139003_hr12_Nutlin<-IMR90_GSE139003_hr12_Nutlin+1
IMR90_GSE139003_hr12_Nutlin<-as.data.frame(t(scale(t(IMR90_GSE139003_hr12_Nutlin), center=TRUE, scale=TRUE)))
IMR90_GSE139003_hr12_Nutlin<-IMR90_GSE139003_hr12_Nutlin[,c("hr12_Nut_1","hr12_Nut_2")]
names(IMR90_GSE139003_hr12_Nutlin)<-c("IMR90_GSE139003_hr12_Nutlin_1","IMR90_GSE139003_hr12_Nutlin_2")
IMR90_GSE139003_hr12_Nutlin_mean<-as.data.frame(rowMeans(IMR90_GSE139003_hr12_Nutlin))
names(IMR90_GSE139003_hr12_Nutlin_mean)<-c("IMR90_GSE139003_hr12_Nutlin")

MCF10A_hr4_Nutlin<-read.table("other_cell_lines/MCF10A_hr4_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF10A_hr4_Nutlin<-MCF10A_hr4_Nutlin+1
MCF10A_hr4_Nutlin<-as.data.frame(t(scale(t(MCF10A_hr4_Nutlin), center=TRUE, scale=TRUE)))
MCF10A_hr4_Nutlin<-MCF10A_hr4_Nutlin[,c("Nut_1","Nut_2","Nut_3")]
names(MCF10A_hr4_Nutlin)<-c("MCF10A_hr4_Nutlin_1","MCF10A_hr4_Nutlin_2","MCF10A_hr4_Nutlin_3")
MCF10A_hr4_Nutlin_mean<-as.data.frame(rowMeans(MCF10A_hr4_Nutlin))
names(MCF10A_hr4_Nutlin_mean)<-c("MCF10A_hr4_Nutlin")

hNCCs_Nutlin<-read.table("other_cell_lines/hNCCs_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
hNCCs_Nutlin<-hNCCs_Nutlin+1
hNCCs_Nutlin<-as.data.frame(t(scale(t(hNCCs_Nutlin), center=TRUE, scale=TRUE)))
hNCCs_Nutlin<-hNCCs_Nutlin[,c("Nut_1","Nut_2","Nut_3")]
names(hNCCs_Nutlin)<-c("hNCCs_Nutlin_1","hNCCs_Nutlin_2","hNCCs_Nutlin_3")
hNCCs_Nutlin_mean<-as.data.frame(rowMeans(hNCCs_Nutlin))
names(hNCCs_Nutlin_mean)<-c("hNCCs_Nutlin")

hSMCs_Nutlin<-read.table("other_cell_lines/hSMCs_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
hSMCs_Nutlin<-hSMCs_Nutlin+1
hSMCs_Nutlin<-as.data.frame(t(scale(t(hSMCs_Nutlin), center=TRUE, scale=TRUE)))
hSMCs_Nutlin<-hSMCs_Nutlin[,c("Nut_1","Nut_2","Nut_3")]
names(hSMCs_Nutlin)<-c("hSMCs_Nutlin_1","hSMCs_Nutlin_2","hSMCs_Nutlin_3")
hSMCs_Nutlin_mean<-as.data.frame(rowMeans(hSMCs_Nutlin))
names(hSMCs_Nutlin_mean)<-c("hSMCs_Nutlin")

A431_Doxorubicin<-read.table("other_cell_lines/A431_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
A431_Doxorubicin<-A431_Doxorubicin+1
A431_Doxorubicin<-as.data.frame(t(scale(t(A431_Doxorubicin), center=TRUE, scale=TRUE)))
A431_Doxorubicin<-A431_Doxorubicin[,c("Dox_1","Dox_2","Dox_3")]
names(A431_Doxorubicin)<-c("A431_Doxorubicin_1","A431_Doxorubicin_2","A431_Doxorubicin_3")
A431_Doxorubicin_mean<-as.data.frame(rowMeans(A431_Doxorubicin))
names(A431_Doxorubicin_mean)<-c("A431_Doxorubicin")

M231_5_FU<-read.table("other_cell_lines/M231_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
M231_5_FU<-M231_5_FU+1
M231_5_FU<-as.data.frame(t(scale(t(M231_5_FU), center=TRUE, scale=TRUE)))
M231_5_FU<-M231_5_FU[,c("X5_FU_1","X5_FU_2","X5_FU_3")]
names(M231_5_FU)<-c("M231_5_FU_1","M231_5_FU_2","M231_5_FU_3")
M231_5_FU_mean<-as.data.frame(rowMeans(M231_5_FU))
names(M231_5_FU_mean)<-c("M231_5_FU")

SJSA_GSE89807_5_FU<-read.table("other_cell_lines/SJSA_GSE89807_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
SJSA_GSE89807_5_FU<-SJSA_GSE89807_5_FU+1
SJSA_GSE89807_5_FU<-as.data.frame(t(scale(t(SJSA_GSE89807_5_FU), center=TRUE, scale=TRUE)))
SJSA_GSE89807_5_FU<-SJSA_GSE89807_5_FU[,c("X5_FU_1","X5_FU_2")]
names(SJSA_GSE89807_5_FU)<-c("SJSA_GSE89807_5_FU_1","SJSA_GSE89807_5_FU_2")
SJSA_GSE89807_5_FU_mean<-as.data.frame(rowMeans(SJSA_GSE89807_5_FU))
names(SJSA_GSE89807_5_FU_mean)<-c("SJSA_GSE89807_5_FU")

A549_Nutlin<-read.table("other_cell_lines/A549_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
A549_Nutlin<-A549_Nutlin+1
A549_Nutlin<-as.data.frame(t(scale(t(A549_Nutlin), center=TRUE, scale=TRUE)))
A549_Nutlin<-A549_Nutlin[,c("Nutlin_1","Nutlin_2","Nutlin_3")]
names(A549_Nutlin)<-c("A549_Nutlin_1","A549_Nutlin_2","A549_Nutlin_3")
A549_Nutlin_mean<-as.data.frame(rowMeans(A549_Nutlin))
names(A549_Nutlin_mean)<-c("A549_Nutlin")

MCF7_GSE86221_Nutlin<-read.table("other_cell_lines/MCF7_GSE86221_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MCF7_GSE86221_Nutlin<-MCF7_GSE86221_Nutlin+1
MCF7_GSE86221_Nutlin<-as.data.frame(t(scale(t(MCF7_GSE86221_Nutlin), center=TRUE, scale=TRUE)))
MCF7_GSE86221_Nutlin<-MCF7_GSE86221_Nutlin[,c("Nutlin_1","Nutlin_2")]
names(MCF7_GSE86221_Nutlin)<-c("MCF7_GSE86221_Nutlin_1","MCF7_GSE86221_Nutlin_2")
MCF7_GSE86221_Nutlin_mean<-as.data.frame(rowMeans(MCF7_GSE86221_Nutlin))
names(MCF7_GSE86221_Nutlin_mean)<-c("MCF7_GSE86221_Nutlin")

MRC5_UV<-read.table("other_cell_lines/MRC5_filter.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
MRC5_UV<-MRC5_UV+1
MRC5_UV<-as.data.frame(t(scale(t(MRC5_UV), center=TRUE, scale=TRUE)))
MRC5_UV<-MRC5_UV[,c("UV_hr3_1","UV_hr3_2")]
names(MRC5_UV)<-c("MRC5_UV_1","MRC5_UV_2")
MRC5_UV_mean<-as.data.frame(rowMeans(MRC5_UV))
names(MRC5_UV_mean)<-c("MRC5_UV")

data_sets<-c("IMR90_Nutlin","Fibro_Nutlin","MCF10A_Nutlin","PAEC_Nutlin","IMR90_GSE139003_hr12_Nutlin","MCF10A_hr4_Nutlin",
             "hNCCs_Nutlin","hSMCs_Nutlin","IMR90_5_FU","GM06170_Doxorubicin","GM00011_Doxorubicin","MRC5_UV",
             "HCT116_Nutlin","HCT116_GSE137297_wt_Nutlin","MCF7_GSE47042_Nutlin","A549_Nutlin","MCF7_GSE86221_Nutlin",
             "SJSA_Nutlin","MV4_11_Nutlin","HCT116_5_FU","M231_5_FU","SJSA_GSE89807_5_FU","A431_Doxorubicin","Saos_2_doxycyclin")

#Produce the means for each dataset
data_sets_mean<-paste(data_sets,"mean",sep="_")

count_heat<-get(data_sets[1])
count_heat_mean<-get(data_sets_mean[1])

count_heat<-cbind("RE"=row.names(count_heat), count_heat)
count_heat_mean<-cbind("RE"=row.names(count_heat_mean), count_heat_mean)

i<-2
while (length(data_sets)>i-1)
{
  temp<-get(data_sets[i])
  count_heat<-merge(count_heat,temp, by.x="RE",by.y="row.names", all=TRUE)
  
  
  i<-i+1
  
}


row.names(count_heat)<-count_heat$RE
count_heat<-count_heat[,2:ncol(count_heat)]
count_heat[is.na(count_heat)] <- -2.0


i<-2
while (length(data_sets_mean)>i-1)
{
  temp<-get(data_sets_mean[i])
  count_heat_mean<-merge(count_heat_mean,temp, by.x="RE",by.y="row.names", all=TRUE)
  
  
  i<-i+1
  
}


row.names(count_heat_mean)<-count_heat_mean$RE
count_heat_mean<-count_heat_mean[,2:ncol(count_heat_mean)]
count_heat_mean[is.na(count_heat_mean)] <- -2.0

#TE occurance.
FC_0<-read.table("RE_occurances/padj_0.05_FC_0.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_0<-FC_0[,c("RE","RE_occurances")]
FC_0_2<-subset(FC_0, RE_occurances>=2)
FC_0_3<-subset(FC_0, RE_occurances>=3)

FC_1<-read.table("RE_occurances/padj_0.05_FC_1.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_1<-FC_1[,c("RE","RE_occurances")]
FC_1_2<-subset(FC_1, RE_occurances>=2)
FC_1_3<-subset(FC_1, RE_occurances>=3)

FC_2<-read.table("RE_occurances/padj_0.05_FC_2.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )
FC_2<-FC_2[,c("RE","RE_occurances")]
FC_2_2<-subset(FC_2, RE_occurances>=2)
FC_2_3<-subset(FC_2, RE_occurances>=3)


#count_heat subset
count_heat_FC_0_2<-subset(count_heat, row.names(count_heat) %in% FC_0_2$RE)
count_heat_FC_0_3<-subset(count_heat, row.names(count_heat) %in% FC_0_3$RE)


count_heat_FC_1_2<-subset(count_heat, row.names(count_heat) %in% FC_1_2$RE)
count_heat_FC_1_3<-subset(count_heat, row.names(count_heat) %in% FC_1_3$RE)

count_heat_FC_2_2<-subset(count_heat, row.names(count_heat) %in% FC_2_2$RE)
count_heat_FC_2_3<-subset(count_heat, row.names(count_heat) %in% FC_2_3$RE)

#count_heat_mean
count_heat_mean_FC_0_2<-subset(count_heat_mean, row.names(count_heat_mean) %in% FC_0_2$RE)
count_heat_mean_FC_0_3<-subset(count_heat_mean, row.names(count_heat_mean) %in% FC_0_3$RE)


count_heat_mean_FC_1_2<-subset(count_heat_mean, row.names(count_heat_mean) %in% FC_1_2$RE)
count_heat_mean_FC_1_3<-subset(count_heat_mean, row.names(count_heat_mean) %in% FC_1_3$RE)

count_heat_mean_FC_2_2<-subset(count_heat_mean, row.names(count_heat_mean) %in% FC_2_2$RE)
count_heat_mean_FC_2_3<-subset(count_heat_mean, row.names(count_heat_mean) %in% FC_2_3$RE)

#Sort by family name, this is the important part of this script.
family<-read.table("RE_family.csv", header=TRUE, sep=',', stringsAsFactors=FALSE )

family$repClass<-gsub("?","",family$repClass, fixed=TRUE)
family$repFamily<-gsub("?","",family$repFamily, fixed=TRUE)

RE_ID<-row.names(count_heat)

i<-1
REs<-c()

while (length(RE_ID)>i-1)
{
  temp_ID<-RE_ID[i]
  temp<-strsplit(temp_ID," ")[[1]]
  temp_RE<-temp[1]
  REs<-c(REs,temp_RE)
  

  i<-i+1
}

REs<-data.frame("RE_ID"=RE_ID,"RE"=REs, stringsAsFactors = FALSE)

family_names<-unique(family$repFamily)
class_names<-unique(family$repClass)

i<-1
families<-c()
nums<-c()

while (length(family_names)>i-1)
{
  
  temp<-subset(family, repFamily==family_names[i])
  families<-c(families,family_names[i])
  nums<-c(nums,nrow(temp))
  
  i<-i+1
  
}

family_num<-data.frame("family"=families,"num"=nums, stringsAsFactors = FALSE)
family_num<-family_num[order(family_num$num, decreasing=TRUE),]
write.table(family_num, "Heatmaps/occurance_families.csv", sep=',', row.names = F,col.names = T, quote = F)


i<-1
classes<-c()
nums<-c()

while (length(class_names)>i-1)
{
  
  temp<-subset(family, repClass==class_names[i])
  classes<-c(classes,class_names[i])
  nums<-c(nums,nrow(temp))
  
  i<-i+1
  
}

class_num<-data.frame("class"=classes,"num"=nums, stringsAsFactors = FALSE)
class_num<-class_num[order(class_num$num, decreasing=TRUE),]
write.table(class_num, "Heatmaps/occurance_classes.csv", sep=',', row.names = F,col.names = T, quote = F)
#create IDs
family<-cbind(family, "ID"=paste("ID",row.names(family), sep=""),stringsAsFactors=FALSE  )

#remove any duplicate entries
dupes<-family[duplicated(family$repName),]
dupes<-subset(family, family$repName %in% dupes$repName)
dupes<-dupes[order(dupes$repName, decreasing=TRUE),]

dupes<-merge(dupes, family_num, by.x="repFamily", by.y="family")
names(dupes)<-c("repFamily","repName","repClass","ID","family_num")

dupes<-merge(dupes, class_num, by.x="repClass", by.y="class")
names(dupes)<-c("repClass","repFamily","repName","ID","family_num", "class_num")
dupes<-dupes[order(dupes$repName, decreasing=TRUE),]

RE_names<-unique(dupes$repName)
to_keep<-c()

i<-1
while (length(RE_names)>i-1)
{
  temp<-subset(dupes, repName== RE_names[i])
  temp_loc<-temp$ID[1]
  print(temp_loc)
  
  j<-2
  
  while(nrow(temp)>j-1)
  {
    if (temp[j,]$class_num==temp[j-1,]$class_num)
    {
      if (temp[j,]$family_num>temp[j-1,]$family_num)
      {
        
        temp_loc<-temp$ID[j]
        
      }
      
    }
    
    else if (temp[j,]$class_num>temp[j-1,]$class_num)
    {
      
      temp_loc<-temp$ID[j]
      
      
    }
    print(temp_loc)
    
    j<-j+1
  }
  
  
  to_keep<-c(to_keep,temp_loc)
  print("")
  
  i<-i+1
}
#remove the duplicate entries.
dupes<-subset(dupes, !(dupes$ID %in% to_keep))

family<-subset(family, !(family$ID %in% dupes$ID))


REs_merge<-merge(REs,family, by.x="RE", by.y="repName")
missing<-subset(REs, !(REs$RE %in% REs_merge$RE))

missing<-cbind(missing, repClass="LTR",repFamily="ERVL")
REs_merge<-rbind(REs_merge[,c("RE","RE_ID","repClass","repFamily")],missing, stringsAsFactors=FALSE)

REs_merge<-cbind(REs_merge,color_cat="other", stringsAsFactors=FALSE)
REs_merge$color_cat[REs_merge$repClass=="Simple_repeat"] = "Simple_repeat"
REs_merge$color_cat[REs_merge$repClass=="LINE"] = "LINE"
REs_merge$color_cat[REs_merge$repClass=="SINE"] = "SINE"
REs_merge$color_cat[REs_merge$repClass=="LTR"] = "LTR"
write.table(REs_merge, "Heatmaps/RE_family_class.csv", sep=',', row.names = F,col.names = T, quote = F) #this will be used in the next part of the script.


library(tidyr)*

setwd("C:/Users/[USER NAME]/Documents/[LOCATION ON DRIVE]/telescope")

classes<-read.table("Heatmaps/RE_family_class_update.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #made in previous entry
classes<-classes[,c("RE","RE_ID","repClass","color_cat")]
names(classes)<-c("RE","RE_ID","repClass","major_class")

all<-unique(classes$repClass) #make sure that entries are unique.
pie<-c("SINE","LINE","LTR","other") #a subject of the most important classes.

ht<-"RE_occurances"
k<-list.files(ht, full.names=TRUE, pattern="padj")
k<-k[grep("0.csv|1.csv|2.csv",k)]
k<-k[grep("ChIP-seq",k,invert = TRUE)]


i<-1
while(length(k)>i-1)
{
  
  occur<-read.table(k[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  occur<-occur[,c("RE","RE_occurances")]
  occur_2<-subset(occur, RE_occurances>=2)
  occur_3<-subset(occur, RE_occurances>=3)
  
  name<-strsplit(k[i],"/")[[1]]
  name<-strsplit(name[2],".csv")[[1]]

  
  occur_2_class<-merge(classes,occur_2, by.y="RE",by.x="RE_ID")
  occur_3_class<-merge(classes,occur_3, by.y="RE",by.x="RE_ID")

  
  temp_name<-paste(name,"occur_2",sep="_")
  temp_name<-paste(temp_name,"RE_classes", sep="_")
  temp_name<-paste(temp_name,"csv",sep=".")
 
  temp_name<-paste("classes","RE_classes",temp_name,sep="/")
  
  print(temp_name)
  write.table(occur_2_class, temp_name, sep=',', row.names = F, col.names = T, quote = F)
  
  
  temp_name<-paste(name,"occur_3",sep="_")
  temp_name<-paste(temp_name,"RE_classes", sep="_")
  temp_name<-paste(temp_name,"csv",sep=".")
  
  temp_name<-paste("classes","RE_classes",temp_name,sep="/")
  
  print(temp_name)
  write.table(occur_3_class, temp_name, sep=',', row.names = F, col.names = T, quote = F)
  
  
  
  i<-i+1
  
}


#RE genomic location

files<-read.table("overlaps/paired_files.csv", header=TRUE, sep=',', stringsAsFactors=FALSE ) #File that show files that are paired together for the anyalisis.

UTR<-files[grep("UTR_REs",files$data),]
gene<-files[grep("intron_exon",files$data),]
intergenic<-files[grep("no_overlap",files$data),]

i<-1
file_names<-c()

while (nrow(UTR)>i-1)
{
  temp<-strsplit(UTR$data[i],"_")[[1]]
  num<-grep("occur",temp)
  temp_name<-paste(temp[1:(num+1)],collapse="_")
  temp_name<-paste("overlaps","RE_overlaps", temp_name, sep="/")
  
  
  if (grepl("up|down",temp[length(temp)]))
  {
    reg<-strsplit(temp[length(temp)],".csv")[[1]]
    temp_name<-paste(temp_name,reg, sep="_")
    temp_name<-paste(temp_name,"RE_overlaps",sep="_")
    temp_name<-paste(temp_name,"tsv",sep=".")
    
    
  }
  
  else
  {
    temp_name<-paste(temp_name,"RE_overlaps",sep="_")
    temp_name<-paste(temp_name,"tsv",sep=".")
    
    
  }
  
  file_names<-c(file_names,temp_name)
  
  i<-i+1
}



UTR$data<-paste("overlaps",UTR$data, sep="/")
UTR$info<-paste("overlaps",UTR$info, sep="/")

gene$data<-paste("overlaps",gene$data, sep="/")
gene$info<-paste("overlaps",gene$info, sep="/")

intergenic$data<-paste("overlaps",intergenic$data, sep="/")


UTR<-cbind(names=file_names,UTR,stringsAsFactors=FALSE)
gene<-cbind(names=file_names,gene,stringsAsFactors=FALSE)
intergenic<-cbind(names=file_names,intergenic,stringsAsFactors=FALSE)

file_names<-merge(UTR,gene, by="names")
names(file_names)<-c("names","data_UTR","info_UTR","data_gene","info_gene")
file_names<-merge(file_names,intergenic[,c("names","data")], by="names")
names(file_names)<-c("names","data_UTR","info_UTR","data_gene","info_gene", "data_intergenic")


i<-1
while (nrow(file_names)>i-1)
{
  
  UTR_data<-read.table(file_names$data_UTR[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  UTR_info<-read.table(file_names$info_UTR[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  gene_data<-read.table(file_names$data_gene[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  gene_info<-read.table(file_names$info_gene[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  intergenic_data<-read.table(file_names$data_intergenic[i], header=TRUE, sep=',', stringsAsFactors=FALSE )
  
  #UTR
  if (nrow(UTR_data)!=0)
  {
    
    
    UTR_info<-separate(data = UTR_info, col = UTR_ID, into = c("ID","type_prime","prime", "UTR","UTR_num"), sep = "_")
    UTR_info<-cbind(UTR_info, "UTR_type"= paste(UTR_info$type_prime, "prime","UTR", sep="_"))
    UTR_info<-UTR_info[,c("RE","UTR_type")]
    UTR_info<-unique(UTR_info)
    
    REs<-unique(UTR_data$RE)
    
    prime_RE<-c()
    primes<-c()
    
    for (temp_RE in REs)
    {
      
      temp_info<-subset(UTR_info, RE==temp_RE)
      
      temp_prime<-paste(temp_info$UTR_type, collapse=", ")
      prime_RE<-c(prime_RE,temp_RE)
      primes<-c(primes,temp_prime)
      
    }
    
    temp<-data.frame(RE=prime_RE,overlap=primes)
    
    UTR_data<-merge(UTR_data,temp, by="RE")
    UTR_data<-UTR_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    UTR_data<-cbind(UTR_data,overlap_cat="UTR")
    UTR_data<-UTR_data[,c("RE", "gene_name","overlap","overlap_cat","chr","start","end","strand","RE_occurances")]
    
    
    
  }
  

  #gene
  
  if (nrow(gene_data)!=0)
  {
    
    gene_info<-separate(data = gene_info, col = gene_ID, into = c("ID","gene_type","num_type_gene"), sep = "_")
    gene_info<-gene_info[,c("RE","gene_type")]
    gene_info<-unique(gene_info)
    
    
    
    #FC_0
    REs<-unique(gene_data$RE)
    
    prime_RE<-c()
    primes<-c()
    
    for (temp_RE in REs)
    {
      
      temp_info<-subset(gene_info, RE==temp_RE)
      
      temp_prime<-paste(temp_info$gene_type, collapse=", ")
      prime_RE<-c(prime_RE,temp_RE)
      primes<-c(primes,temp_prime)
      
    }
    
    temp<-data.frame(RE=prime_RE,overlap=primes)
    
    gene_data<-merge(gene_data,temp, by="RE")
    gene_data<-gene_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    gene_data<-cbind(gene_data,overlap_cat=gene_data$overlap)
    gene_data<-gene_data[,c("RE", "gene_name","overlap","overlap_cat","chr","start","end","strand","RE_occurances")]
    
    
    
  }
  
  
  
  
  
  #no overlap
  if (nrow(intergenic_data)!=0)
  {
    intergenic_data<-cbind(intergenic_data, "overlap"="Intergenic")
    intergenic_data<-intergenic_data[,c("RE", "gene_name","overlap","chr","start","end","strand","RE_occurances")]
    intergenic_data<-cbind(intergenic_data,overlap_cat=intergenic_data$overlap)
    intergenic_data<-intergenic_data[,c("RE", "gene_name","overlap","overlap_cat","chr","start","end","strand","RE_occurances")]
    
  }
  
  
  
  
  #combine
  
  #check for 
  temp_data<-rbind(UTR_data, gene_data)
  temp_data<-rbind(temp_data,intergenic_data)
  temp_data<-temp_data[order(temp_data$RE_occurances, decreasing=TRUE),]
  
  
  
  temp_name<-paste("prime_percentages", file_names$names[i], sep="/")
  print(temp_name)
  write.table(temp_data, file_names$names[i], sep='\t', row.names = F, col.names = T, quote = F)
  
  
  
  i<-i+1
}


#combine


ht<-"classes/RE_classes"
classes<-list.files(ht, full.names=TRUE, pattern=".csv")
all_classes<-classes[grep("normal|cancer",classes, invert=TRUE)]
cancer_classes<-classes[grep("cancer",classes)]
normal_classes<-classes[grep("normal",classes)]


ht<-"overlaps/RE_overlaps"
overlaps<-list.files(ht, full.names=TRUE, pattern=".tsv")
all_overlaps<-overlaps[grep("normal|cancer",overlaps, invert=TRUE)]
cancer_overlaps<-overlaps[grep("cancer",overlaps)]
normal_overlaps<-overlaps[grep("normal",overlaps)]

subs<-c("FC_0_occur_2","FC_0_occur_3","FC_1_occur_2","FC_1_occur_2","FC_2_occur_2","FC_2_occur_3")

#all
i<-1

while (length(subs)>i-1)
{
temp_class<-all_classes[grep(subs[i],all_classes)]
temp_class_up<-temp_class[grep("up",temp_class)][1]
temp_class_down<-temp_class[grep("down",temp_class)][1]
temp_class_all<-temp_class[grep("down|up",temp_class, invert = TRUE)][1]

temp_overlap<-all_overlaps[grep(subs[i],all_overlaps)]
temp_overlap_up<-temp_overlap[grep("up",temp_overlap)][1]
temp_overlap_down<-temp_overlap[grep("down",temp_overlap)][1]
temp_overlap_all<-temp_overlap[grep("down|up",temp_overlap, invert = TRUE)][1]

#all
temp_class_all<-read.table(temp_class_all, header=TRUE, sep=',', stringsAsFactors=FALSE )
temp_overlap_all<-read.table(temp_overlap_all, header=TRUE, sep='\t', stringsAsFactors=FALSE )
all<-merge(temp_class_all,temp_overlap_all, by.x="RE_ID", by.y="RE")
all<-all[,c("RE_ID","RE","repClass","major_class","gene_name","overlap","overlap_cat","RE_occurances.x",
            "chr","start","end","strand")]
names(all)<-c("RE_ID","RE","repClass","major_class","gene_name","overlap","overlap_cat","RE_occurances",
              "chr","start","end","strand")

name<-paste("RE_peaks/RE_files",subs[i],sep="/")
name<-paste(name,"RE_overlaps_classes",sep="_")
name<-paste(name,"tsv",sep=".")
write.table(all, name, sep='\t', row.names = F, col.names = T, quote = F)

#up
temp_class_up<-read.table(temp_class_up, header=TRUE, sep=',', stringsAsFactors=FALSE )
temp_overlap_up<-read.table(temp_overlap_up, header=TRUE, sep='\t', stringsAsFactors=FALSE )
up<-merge(temp_class_up,temp_overlap_up, by.x="RE_ID", by.y="RE")
up<-up[,c("RE_ID","RE","repClass","major_class","gene_name","overlap","overlap_cat","RE_occurances.x",
            "chr","start","end","strand")]
names(up)<-c("RE_ID","RE","repClass","major_class","gene_name","overlap","overlap_cat","RE_occurances",
              "chr","start","end","strand")

name<-paste("RE_peaks/RE_files",subs[i],sep="/")
name<-paste(name,"up_RE_overlaps_classes",sep="_")
name<-paste(name,"tsv",sep=".")
write.table(up, name, sep='\t', row.names = F, col.names = T, quote = F)

#down
temp_class_down<-read.table(temp_class_down, header=TRUE, sep=',', stringsAsFactors=FALSE )
temp_overlap_down<-read.table(temp_overlap_down, header=TRUE, sep='\t', stringsAsFactors=FALSE )
down<-merge(temp_class_down,temp_overlap_down, by.x="RE_ID", by.y="RE")
down<-down[,c("RE_ID","RE","repClass","major_class","gene_name","overlap","overlap_cat","RE_occurances.x",
          "chr","start","end","strand")]
names(down)<-c("RE_ID","RE","repClass","major_class","gene_name","overlap","overlap_cat","RE_occurances",
             "chr","start","end","strand")

name<-paste("RE_peaks/RE_files",subs[i],sep="/")
name<-paste(name,"down_RE_overlaps_classes",sep="_")
name<-paste(name,"tsv",sep=".")
write.table(down, name, sep='\t', row.names = F, col.names = T, quote = F)

i<-i+1
}
