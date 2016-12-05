library(gplots)
EVE_names=c("AT3G46640","AT2G25930", "AT2G40080")
names(EVE_names)=c("lux", "elf3", "elf4")


a=load("rna_seq_filtered.RData")
lux22=rna_seq_filtered[["lux_22"]]
col0=rna_seq_filtered[["col0_22"]]
elf322=rna_seq_filtered[["elf3_22"]]

phy22=rna_seq_filtered[["phy_22"]]
phy27=rna_seq_filtered[["phy_27"]]
ler022=rna_seq_filtered[["ler0_22"]]
ler027=rna_seq_filtered[["ler0_27"]]

rownames(col0)= col0[,1]
rownames(lux22)= lux22[,1]
rownames(elf322)= elf322[,1]

#for some reason a few genes are listed twice-- just delete them for now: none are crucial for EC story
phy22=phy22[which(as.character(phy22[,1]) %in% as.character(names(table(phy22[,1])))[which(table(phy22[,1])==1)]),]
rownames(phy22)=substring(phy22[,1], 6, 14)

phy27=phy27[which(as.character(phy27[,1]) %in% as.character(names(table(phy27[,1])))[which(table(phy27[,1])==1)]),]
rownames(phy27)=substring(phy27[,1], 6, 14)

ler022=ler022[which(as.character(ler022[,1]) %in% as.character(names(table(ler022[,1])))[which(table(ler022[,1])==1)]),]
rownames(ler022)=substring(ler022[,1], 6, 14)

ler027=ler027[which(as.character(ler027[,1]) %in% as.character(names(table(ler027[,1])))[which(table(ler027[,1])==1)]),]
rownames(ler027)=substring(ler027[,1], 6, 14)

lux27=rna_seq_filtered[["lux_27"]]
col027=rna_seq_filtered[["col0_27"]]
elf327=rna_seq_filtered[["elf3_27"]]
rownames(col027)= col027[,1]
rownames(lux27)= lux27[,1]
rownames(elf327)= elf327[,1]

#read in list of G-box motifs
#elf4_fl2=read.table("beds/ELF4beds_FL2_matches.txt", header=TRUE, sep="\t") ########################################### DO
#read in list of LBS motifs
#elf4_gboxes=read.table("beds/ELF4beds_Gbox_matches.txt", header=TRUE, sep="\t") ########################################### DO

#load gene targets
luxTargets=read.table("geneTargets/LUXgenes2.txt")
elf3Targets=read.table("geneTargets/ELF3genes2.txt")
elf4Targets=read.table("geneTargets/ELF4genes2.txt")

########################

#identify which potential gene targets have differential expression in at least one time point (which I'm going to define as 1.5 )

#load htseq counts
library(edgeR)

#rna-seq_20140814_ler0_22C_ZT-2_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct
tp=c("-2", "-4", "-8", "-12", "0", "1", "4", "8")
col22_files=c("htCounts/rna-seq_20140814_col_22C_ZT-2_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct",
              paste("htCounts/rna-seq_2201404_col_ZT", tp[2:5], "-22c_R1_raw_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct", sep=""),
              "htCounts/rna-seq_20140814_col_22C_ZT1_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct",
              paste("htCounts/rna-seq_2201404_col_ZT", tp[7:8], "-22c_R1_raw_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct", sep=""))
col22_counts=lapply(col22_files, function(i){
  a=read.table(i)
  if(nchar(as.character(a[1,1]))>9){
    rownames(a)=substring(as.character(a[,1]), 6, 14)
  }else{
    rownames(a)=a[,1] #substring(as.character(a[,1]), 6, 14)
  }
  a
})


#read.table("htCounts/rna-seq_20140814_lux-4_22C_ZT-2_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_unstranded_sorted_rmdup_picard.bam.ct")
lux22_files=paste("htCounts/rna-seq_20140814_lux-4_22C_ZT", tp, "_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_unstranded_sorted_rmdup_picard.bam.ct", sep="")
lux22_counts=lapply(lux22_files, function(i){
  a=read.table(i)
  rownames(a)=a[,1]
  a
})

#elf3: rna-seq_20140814_elf3-1_27C_ZT-12_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct
elf3_22_files=paste("htCounts/rna-seq_20140814_elf3-1_27C_ZT", tp, "_trimmo_paired_2_10_5_1_tophat_TAIR10_ensembl_nomixed_sorted_rmdup_picard.bam.ct", sep="")
elf3_22_counts=lapply(elf3_22_files, function(i){
  a=read.table(i)
  rownames(a)=a[,1]
  a
})


top5percent_de=lapply(c(1:length(tp)), function(i){
  sharedGenes=rownames(col22_counts[[i]])[which(rownames(col22_counts[[i]]) %in% rownames(lux22_counts[[i]]))]
  extractShared=cbind(col22_counts[[i]][sharedGenes,2], lux22_counts[[i]][sharedGenes,2])
  colnames(extractShared)=c(paste("col22_ZT", tp[i], sep=""), paste("lux22_ZT", tp[i], sep=""))
  rownames(extractShared)=sharedGenes
  group <- factor(c(1,2))
  
  y <- DGEList(counts=as.matrix(extractShared),group=group)
  et <- exactTest(y, pair=1:2, dispersion=0.1) #should be 0.1
  tt=topTags(et, n=floor(0.05*length(sharedGenes))) 
  write.table(tt, file="testRankDispersion.1.txt")
  tt
}) 

top5percent_elf3_de=lapply(c(1:length(tp)), function(i){
  sharedGenes=rownames(col22_counts[[i]])[which(rownames(col22_counts[[i]]) %in% rownames(elf3_22_counts[[i]]))]
  extractShared=cbind(col22_counts[[i]][sharedGenes,2], elf3_22_counts[[i]][sharedGenes,2])
  colnames(extractShared)=c(paste("col22_ZT", tp[i], sep=""), paste("elf3_22_ZT", tp[i], sep=""))
  rownames(extractShared)=sharedGenes
  group <- factor(c(1,2))
  
  y <- DGEList(counts=as.matrix(extractShared),group=group)
  et <- exactTest(y, pair=1:2, dispersion=0.1)
  tt=topTags(et, n=floor(0.05*length(sharedGenes)))
  tt
}) 

compileTopGenes=unique(c(unlist(lapply(top5percent_de, function(i){rownames(i)})),
                         unlist(lapply(top5percent_elf3_de, function(i){rownames(i)}))))  #####Now in the following graphs, I can verify that the genes are in this vector
elf3Targets_de=as.character(elf3Targets[which(as.character(elf3Targets[,1]) %in% compileTopGenes),1])
elf4Targets_de=as.character(elf4Targets[which(as.character(elf4Targets[,1]) %in% compileTopGenes),1])
luxTargets_de=as.character(luxTargets[which(as.character(luxTargets[,1]) %in% compileTopGenes),1])
lux17CTargets=as.character(read.table("geneTargets/LUX_targetgenes_17C_TP3.txt")[,1])
lux17CTargets_de=as.character(lux17CTargets)[which(as.character(lux17CTargets) %in% compileTopGenes)]

####input for venn diagram:
#do the binary table trick
uniqueGenes=unique(c(lux17CTargets_de,  luxTargets_de, elf4Targets_de,  elf3Targets_de))
mat=sapply(uniqueGenes, function(i){
  sapply(list(lux17CTargets_de,  luxTargets_de, elf4Targets_de,  elf3Targets_de), function(j){
    if(as.character(i) %in% as.character(j)){1}else{0}
  })
})
binary=apply(mat, 2, function(i){paste(i[1], i[2], i[3], i[4], sep="")})
table(binary)

write.table(elf3Targets_de, file="goldStandardTargets_elf3.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(elf4Targets_de, file="goldStandardTargets_elf4.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(luxTargets_de, file="goldStandardTargets_lux.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)


