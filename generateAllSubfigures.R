


##########This file will generate most of the subfigures used in the paper
#####Please run "findTargetGenes.R" first to load all necessary input files and place them in the correct data structures


source("zscore.R")

rna_to_combine=list(col0, col027, elf322, elf327, lux22, lux27)
names(rna_to_combine)=c("col0", "col0", "elf3", "elf3", "lux", "lux")
zscore_rnaSeq_de=lapply(c(1:length(rna_to_combine)), function(id){
  i=rna_to_combine[[id]]
  print("hi")
  a=i[compileTopGenes,2:9]
  print(dim(a))
  a=t(apply(a, 1, function(j){zscore(j)}))
  print(dim(a))
  colnames(a)=paste(names(rna_to_combine)[id], substring(colnames(i)[2:9], 4, nchar(colnames(i)[2:9])))
  rownames(a)=compileTopGenes
  a
})



##########Make figure of temperature responsiveness for new Figure 4

allGenesCorr=sapply(c(2:9), function(i){
  a=sapply(compileTopGenes, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(elf322[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[compileTopGenes[which(a)],i]/col0[compileTopGenes[which(a)],i]), 
      log(elf322[compileTopGenes[which(a)],i]/col0[compileTopGenes[which(a)],i]))
})
names(allGenesCorr)=colnames(col0)[2:9]

ELF4targetsCorr=sapply(c(2:9), function(i){
  a=sapply(elf4Targets_de, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(elf322[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[elf4Targets_de[which(a)],i]/col0[elf4Targets_de[which(a)],i]), 
      log(elf322[elf4Targets_de[which(a)],i]/col0[elf4Targets_de[which(a)],i]))
})


ELF3targetsCorr=sapply(c(2:9), function(i){
  a=sapply(elf3Targets_de, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(elf322[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[elf3Targets_de[which(a)],i]/col0[elf3Targets_de[which(a)],i]), 
      log(elf322[elf3Targets_de[which(a)],i]/col0[elf3Targets_de[which(a)],i]))
})

LUXtargetsCorr=sapply(c(2:9), function(i){
  a=sapply(luxTargets_de, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(elf322[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[luxTargets_de[which(a)],i]/col0[luxTargets_de[which(a)],i]), 
      log(elf322[luxTargets_de[which(a)],i]/col0[luxTargets_de[which(a)],i]))
})





allGenesCorr_lux=sapply(c(2:9), function(i){
  a=sapply(compileTopGenes, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(lux22[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[compileTopGenes[which(a)],i]/col0[compileTopGenes[which(a)],i]), 
      log(lux22[compileTopGenes[which(a)],i]/col0[compileTopGenes[which(a)],i]))
})
names(allGenesCorr_lux)=colnames(col0)[2:9]

ELF4targetsCorr_lux=sapply(c(2:9), function(i){
  a=sapply(elf4Targets_de, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(lux22[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[elf4Targets_de[which(a)],i]/col0[elf4Targets_de[which(a)],i]), 
      log(lux22[elf4Targets_de[which(a)],i]/col0[elf4Targets_de[which(a)],i]))
})


ELF3targetsCorr_lux=sapply(c(2:9), function(i){
  a=sapply(elf3Targets_de, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(lux22[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[elf3Targets_de[which(a)],i]/col0[elf3Targets_de[which(a)],i]), 
      log(lux22[elf3Targets_de[which(a)],i]/col0[elf3Targets_de[which(a)],i]))
})

LUXtargetsCorr_lux=sapply(c(2:9), function(i){
  a=sapply(luxTargets_de, function(j){
    b=log(col027[j,i]/col0[j,i])
    c=log(lux22[j,i]/col0[j,i])
    if(is.na(b) | is.na(c) | is.nan(b) | is.nan(c) | is.infinite(b) | is.infinite(c)){
      FALSE
    }else{TRUE}
  })
  
  cor(log(col027[luxTargets_de[which(a)],i]/col0[luxTargets_de[which(a)],i]), 
      log(lux22[luxTargets_de[which(a)],i]/col0[luxTargets_de[which(a)],i]))
})


pdf(paste("Figure_supplement_correlationWithTempTranscriptome.pdf", sep=""), width=8, height=4,pointsize = 10) 
par(mfrow=c(1,2));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

timepoints=c(-4, -2, 0, 1, 4, 8, 12, 16)
plot(timepoints, allGenesCorr, type="l", ylim=c(-0.1,1), ylab="Pearson's R", xlab="time (h)", lwd=2)
lines(timepoints, ELF3targetsCorr, col="blue", lwd=2)
lines(timepoints, ELF4targetsCorr, col="deeppink3", lwd=2)
lines(timepoints, LUXtargetsCorr, col="goldenrod", lwd=2)
abline(v=c(0, 8), col="grey")
legend(10, 1.05, c("all", "ELF3", "ELF4", "LUX"), col=c("black", "blue", "deeppink3", "goldenrod"), lwd=2, bty='n')

timepoints=c(-4, -2, 0, 1, 4, 8, 12, 16)
plot(timepoints, allGenesCorr_lux, type="l", ylim=c(-0.1,1), ylab="Pearson's R", xlab="time (h)", lwd=2)
lines(timepoints, ELF3targetsCorr_lux, col="blue", lwd=2)
lines(timepoints, ELF4targetsCorr_lux, col="deeppink3", lwd=2)
lines(timepoints, LUXtargetsCorr_lux, col="goldenrod", lwd=2)
abline(v=c(0, 8), col="grey")
#legend(10, 1.1, c("all", "ELF3", "ELF4", "LUX"), col=c("black", "blue", "deeppink3", "goldenrod"), lwd=2)

dev.off()


pdf(paste("Figure_temperatureResponseElf3LUX_v2.pdf", sep=""), width=8, height=8,pointsize = 10) 
par(mfrow=c(2,2));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)
plot(log(col027[compileTopGenes,"TPM_ZT0.27c"]/col0[compileTopGenes,"TPM_ZT0.22c"]), 
     log(elf322[compileTopGenes,"TPM_ZT0.22c"]/col0[compileTopGenes,"TPM_ZT0.22c"]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.1), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("27oC")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     main="ZT0"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,"TPM_ZT0.27c"]/col0[elf3Targets_de,"TPM_ZT0.22c"]), 
       log(elf322[elf3Targets_de,"TPM_ZT0.22c"]/col0[elf3Targets_de,"TPM_ZT0.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,"TPM_ZT0.27c"]/col0[elf4Targets_de,"TPM_ZT0.22c"]), 
       log(elf322[elf4Targets_de,"TPM_ZT0.22c"]/col0[elf4Targets_de,"TPM_ZT0.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,"TPM_ZT0.27c"]/col0[luxTargets_de,"TPM_ZT0.22c"]), 
       log(elf322[luxTargets_de,"TPM_ZT0.22c"]/col0[luxTargets_de,"TPM_ZT0.22c"]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]), 
       log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]),
       col="forestgreen", pch=19)

legend(-5.2, 5.2, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), "blue3", "mediumvioletred", "goldenrod3", "forestgreen"), bty="n", cex=0.8, pch=19)


####

plot(log(col027[compileTopGenes,"TPM_ZT.4.27c"]/col0[compileTopGenes,"TPM_ZT.4.22c"]), 
     log(elf322[compileTopGenes,"TPM_ZT.4.22c"]/col0[compileTopGenes,"TPM_ZT.4.22c"]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.1), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("27oC")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     main="ZT20"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,"TPM_ZT.4.27c"]/col0[elf3Targets_de,"TPM_ZT.4.22c"]), 
       log(elf322[elf3Targets_de,"TPM_ZT.4.22c"]/col0[elf3Targets_de,"TPM_ZT.4.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,"TPM_ZT.4.27c"]/col0[elf4Targets_de,"TPM_ZT.4.22c"]), 
       log(elf322[elf4Targets_de,"TPM_ZT.4.22c"]/col0[elf4Targets_de,"TPM_ZT.4.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,"TPM_ZT.4.27c"]/col0[luxTargets_de,"TPM_ZT.4.22c"]), 
       log(elf322[luxTargets_de,"TPM_ZT.4.22c"]/col0[luxTargets_de,"TPM_ZT.4.22c"]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.27c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.22c"]), 
       log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.22c"]),
       col="forestgreen", pch=19)


plot(log(col027[compileTopGenes,"TPM_ZT0.27c"]/col0[compileTopGenes,"TPM_ZT0.22c"]), 
     log(lux22[compileTopGenes,"TPM_ZT0.22c"]/col0[compileTopGenes,"TPM_ZT0.22c"]),
     xlim=c(-4, 4), ylim=c(-4, 4), col=rgb(0.3, 0.3, 0.3, 0.1), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("27oC")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     main="ZT0"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,"TPM_ZT0.27c"]/col0[elf3Targets_de,"TPM_ZT0.22c"]), 
       log(lux22[elf3Targets_de,"TPM_ZT0.22c"]/col0[elf3Targets_de,"TPM_ZT0.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,"TPM_ZT0.27c"]/col0[elf4Targets_de,"TPM_ZT0.22c"]), 
       log(lux22[elf4Targets_de,"TPM_ZT0.22c"]/col0[elf4Targets_de,"TPM_ZT0.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,"TPM_ZT0.27c"]/col0[luxTargets_de,"TPM_ZT0.22c"]), 
       log(lux22[luxTargets_de,"TPM_ZT0.22c"]/col0[luxTargets_de,"TPM_ZT0.22c"]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]), 
       log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT0.22c"]),
       col="forestgreen", pch=19)

# legend(-5.2, 5.2, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), "blue3", "mediumvioletred", "goldenrod3", "forestgreen"), bty="n", cex=0.8, pch=19)


####

plot(log(col027[compileTopGenes,"TPM_ZT.4.27c"]/col0[compileTopGenes,"TPM_ZT.4.22c"]), 
     log(lux22[compileTopGenes,"TPM_ZT.4.22c"]/col0[compileTopGenes,"TPM_ZT.4.22c"]),
     xlim=c(-4, 4), ylim=c(-4, 4), col=rgb(0.3, 0.3, 0.3, 0.1), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("27oC")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     main="ZT20"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,"TPM_ZT.4.27c"]/col0[elf3Targets_de,"TPM_ZT.4.22c"]), 
       log(lux22[elf3Targets_de,"TPM_ZT.4.22c"]/col0[elf3Targets_de,"TPM_ZT.4.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,"TPM_ZT.4.27c"]/col0[elf4Targets_de,"TPM_ZT.4.22c"]), 
       log(lux22[elf4Targets_de,"TPM_ZT.4.22c"]/col0[elf4Targets_de,"TPM_ZT.4.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,"TPM_ZT.4.27c"]/col0[luxTargets_de,"TPM_ZT.4.22c"]), 
       log(lux22[luxTargets_de,"TPM_ZT.4.22c"]/col0[luxTargets_de,"TPM_ZT.4.22c"]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.27c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.22c"]), 
       log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.4.22c"]),
       col="forestgreen", pch=19)






dev.off()
















combined_zscores=do.call("cbind", zscore_rnaSeq_de)
ids=apply(combined_zscores, 1, function(i){length(which(is.na(i) | is.nan(i) | is.infinite(i)))==0})
geneIDs=ids
ahm=heatmap.2(combined_zscores[ids,], trace="none", Colv=FALSE, dendrogram="row", margins=c(12,8))

#flip order of zscoring
rna_to_combine=list(col0, col027, elf322, elf327, lux22, lux27)
names(rna_to_combine)=c("col0", "col0", "elf3", "elf3", "lux", "lux")

#reorder columns too: to make all night together
nozscore_rnaSeq_de=lapply(c(1:length(rna_to_combine)), function(id){
  i=rna_to_combine[[id]]
  print("hi")
  a=i[compileTopGenes,c(4:9, 2, 3)]
  print(dim(a))
  a=t(apply(a, 1, function(j){j}))
  print(dim(a))
  colnames(a)=paste(names(rna_to_combine)[id], substring(colnames(i)[c(4:9, 2, 3)], 4, nchar(colnames(i)[c(4:9, 2, 3)])))
  rownames(a)=compileTopGenes
  a
})

combined_Nozscores=do.call("cbind", nozscore_rnaSeq_de)
#do z-score now
combined_zscore_rowwise=t(apply(combined_Nozscores, 1, function(i){zscore(i)}))
ids=apply(combined_zscore_rowwise, 1, function(i){length(which(is.na(i) | is.nan(i) | is.infinite(i)))==0})
my_palette <- colorRampPalette(c("darkmagenta", "white", "darkgreen"))(n = 299)

image(t(combined_zscore_rowwise[names(geneIDs)[which(geneIDs)][ahm$rowInd],]))#, col=my_palette)
geneOrder=names(geneIDs)[which(geneIDs)][ahm$rowInd]
colCodeGenes=rep("white", length(geneOrder))
# colCodeGenes[which(geneOrder %in% elf3Targets_de)]="blue"
colCodeGenes[which(geneOrder %in% luxTargets_de)]="red"
#heatmap(combined_zscore_rowwise[names(geneIDs)[which(geneIDs)][ahm$rowInd],], Rowv=NA, Colv=NA, scale="none", col=my_palette, labRow = NA, RowSideColors = colCodeGenes)

# heatmap.2(combined_zscore_rowwise[names(geneIDs)[which(geneIDs)][ahm$rowInd],], dendrogram="none", Colv=NA, Rowv=trace="none", labRow = NA, col=my_palette, margins=c(12,8), scale="none")


a= heatmap.2(combined_zscore_rowwise[ids,1:16], Colv=NA, trace="none", labRow = NA, col=my_palette, margins=c(12,8), scale="none")


jpeg(paste("Figure_hugeClustering3_v3.jpg", sep=""), width=5200, height=5200, units="px", pointsize = 100) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=30)
# a=t(combined_zscore_rowwise[ids,])

#heatmap.2(combined_zscore_rowwise[ids[a$rowInd],], dendrogram="none", Colv=NA, Rowv=NA, trace="none", labRow = NA, col=my_palette, margins=c(12,8), scale="none")

#day and night color bars
daynightCol=rep("lightyellow",length(colnames(combined_zscore_rowwise)))
daynightCol[grep("ZT.", colnames(combined_zscore_rowwise), fixed=TRUE)]="black"
#col-code for text by temperature
tempCol=rep("blue", length(colnames(combined_zscore_rowwise)))
tempCol[grep("27", colnames(combined_zscore_rowwise))]="red"
#make labels that are ZT only
#add text labels for col-0, elf3, and lux separately
colLabels=rep(c("ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16", "ZT20", "ZT22"),6)
hc.rows<- hclust(dist(combined_zscore_rowwise[ids,]))
ct<- cutree(hc.rows, k=7)
colsRow=c("cornflowerblue", "darkblue", "orange", "lightgreen", "yellow", "red", "grey")[ct[rownames(combined_zscore_rowwise[ids,])]]

heatmap.2(combined_zscore_rowwise[ids,], dendrogram = "row", Colv=NA, trace="none", labRow = NA, col=my_palette, margins=c(12,8), scale="none", ColSideColors = daynightCol, RowSideColors=colsRow, colCol=tempCol, labCol=colLabels, key.title="", key.xlab="", key.ylab="")

# heatmap.2(a, trace="none", labCol = NA, dendrogram="row", Colv=NA, col=my_palette, margins=c(12,8), scale="none")

dev.off() 

#####Also do just ratio of elf3/col0 and lux/col0 at 22 and 27
#flip order of zscoring
rna_to_combine=list(col0, col027, elf322, elf327, lux22, lux27)
names(rna_to_combine)=c("col0", "col0", "elf3", "elf3", "lux", "lux")

allElf3=log(cbind(elf322[compileTopGenes, 2:9], elf327[compileTopGenes, 2:9]))
allCol=log(cbind(col0[compileTopGenes, 2:9], col027[compileTopGenes, 2:9]))
allLux=log(cbind(lux22[compileTopGenes, 2:9], lux27[compileTopGenes, 2:9]))
mat=cbind((allElf3-allCol), (allLux-allCol))

ids=apply(mat, 1, function(i){length(which(is.na(i) | is.nan(i) | is.infinite(i) | min(i)<(-6)))==0})
my_palette <- colorRampPalette(c("darkmagenta", "white", "darkgreen"))(n = 299)





#################Change color code
pdf(paste("Figure_transcriptomeElf3vsLux_v4.pdf", sep=""), width=7, height=3.5,pointsize = 10) 
par(mfrow=c(1,2));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

plot(log(lux22[compileTopGenes,"TPM_ZT.12.22c"]/col0[compileTopGenes,"TPM_ZT.12.22c"]), 
     log(elf322[compileTopGenes,"TPM_ZT.12.22c"]/col0[compileTopGenes,"TPM_ZT.12.22c"]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     main="ZT12 (4h dark)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(lux22[elf3Targets_de,"TPM_ZT.12.22c"]/col0[elf3Targets_de,"TPM_ZT.12.22c"]), 
       log(elf322[elf3Targets_de,"TPM_ZT.12.22c"]/col0[elf3Targets_de,"TPM_ZT.12.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(lux22[elf4Targets_de,"TPM_ZT.12.22c"]/col0[elf4Targets_de,"TPM_ZT.12.22c"]), 
       log(elf322[elf4Targets_de,"TPM_ZT.12.22c"]/col0[elf4Targets_de,"TPM_ZT.12.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(lux22[luxTargets_de,"TPM_ZT.12.22c"]/col0[luxTargets_de,"TPM_ZT.12.22c"]), 
       log(elf322[luxTargets_de,"TPM_ZT.12.22c"]/col0[luxTargets_de,"TPM_ZT.12.22c"]),
       col="goldenrod3", pch=19)

points(log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.12.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.12.22c"]), 
       log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.12.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],"TPM_ZT.12.22c"]),
       col="forestgreen", pch=19)

legend(-5.2, 5.2, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), "blue3", "mediumvioletred", "goldenrod3", "forestgreen"), bty="n", cex=0.8, pch=19)

plot(log(lux22[compileTopGenes,6]/col0[compileTopGenes,6]), 
     log(elf322[compileTopGenes,6]/col0[compileTopGenes,6]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     main="ZT4 (4h light)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(lux22[elf3Targets_de,6]/col0[elf3Targets_de,6]), 
       log(elf322[elf3Targets_de,6]/col0[elf3Targets_de,6]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(lux22[elf4Targets_de,6]/col0[elf4Targets_de,6]), 
       log(elf322[elf4Targets_de,6]/col0[elf4Targets_de,6]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(lux22[luxTargets_de,6]/col0[luxTargets_de,6]), 
       log(elf322[luxTargets_de,6]/col0[luxTargets_de,6]),
       col="goldenrod3", pch=19)

points(log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]), 
       log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]),
       col="forestgreen", pch=19)

#legend(-5, 5, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), rgb(0, 0.6, 0, 1), rgb(0, 0, 0.6, 1), rgb(0.9, 0, 0, 1), "orange"), bty="n", cex=0.8, pch=19)

dev.off()

#############
#phyB_17=as.character(read.table("phyb_17C_targets_bound.txt", sep="\t", quote="", comment.char = "", header=TRUE)[,1])
#phyB_27=as.character(read.table("phyb_27C_targets_bound.txt", sep="\t", quote="", comment.char = "", header=TRUE)[,1])


names=c("bELF4", "bLUX17", "bLUX", "bphyB_17", "bphyB_27") #, "PIF1", "PIF3", "PIF4", "PIF5")

beds=lapply(names, function(i){
  read.table(paste("intersections/ELF3_", i, "_loj_intersect2.bed", sep=""))
})

noMatchRowsOnly=lapply(beds, function(i){
  i[which(i[,5]==-1),]
})

allUniquePeaks=unique(paste(beds[[1]][,1], beds[[1]][,2], beds[[1]][,3]))

mat=matrix(rep(1,length(allUniquePeaks)*length(names)),  nrow=length(allUniquePeaks), ncol=length(names))
rownames(mat)=allUniquePeaks
colnames(mat)=names

for(i in c(1:length(noMatchRowsOnly))){
  m=noMatchRowsOnly[[i]]
  n=paste(m[,1], m[,2], m[,3])
  mat[which(rownames(mat) %in% n),i]=0
}

pdf(paste("Figure_bedOverlaps_all_v3.pdf", sep=""), width=3, height=3,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

mat_ordered=mat[order(mat[,1], mat[,2], mat[,3], mat[,4]),]
heatmap(-mat_ordered, scale="none", labRow=NA, Colv=NA, Rowv=NA)
mat_subset=mat[which(apply(mat, 1, function(i){sum(i)>0})),]
mat_subset_ordered=mat_subset[order(mat_subset[,1], mat_subset[,2], mat_subset[,3], mat_subset[,4]),]


dev.off()

pdf(paste("Figure_bedOverlaps_all_v4.pdf", sep=""), width=3, height=3,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

mat_ordered=mat[order(mat[,1], mat[,2], mat[,3], mat[,4]),]
mat_ordered=mat_ordered[,c("bELF4", "bLUX", "bphyB_17", "bphyB_27")]
heatmap(-mat_ordered, scale="none", labRow=NA, Colv=NA, Rowv=NA)
mat_subset=mat[which(apply(mat, 1, function(i){sum(i)>0})),]
mat_subset_ordered=mat_subset[order(mat_subset[,1], mat_subset[,2], mat_subset[,3], mat_subset[,4]),]


dev.off()



#############
#show how these genes are affected by phyq
phy22_v2=read.table("PhyABCDEsample_22C", header=TRUE)
phy27_v2=read.table("PhyABCDEsample_27C", header=TRUE)
#remove duplicate rows
geneNames=substring(as.character(phy22_v2[,1]), 6, 14)
phy22_v2=phy22_v2[which(sapply(geneNames, function(i){(length(which(geneNames==i))==1)})),]
rownames(phy22_v2)=substring(as.character(phy22_v2[,1]), 6, 14)
phy27_v2=phy27_v2[which(sapply(geneNames, function(i){(length(which(geneNames==i))==1)})),]
rownames(phy27_v2)=substring(as.character(phy27_v2[,1]), 6, 14)
phy22_v2=phy22_v2[,-2]
phy27_v2=phy27_v2[,-2]
pdf(paste("Figure_transcriptomePhyQvsLux_v2_color.pdf", sep=""), width=7, height=3.5,pointsize = 10) 
par(mfrow=c(1,2));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

plot(log(lux22[compileTopGenes, "TPM_ZT.12.22c"]/col0[compileTopGenes, "TPM_ZT.12.22c"]), 
     log(phy22_v2[compileTopGenes, "TPM_ZT.12.22c"]/ler022[compileTopGenes, "TPM_ZT.12.22c"]),
     xlim=c(-6.3, 6.3), ylim=c(-6.3, 6.3), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("phyABCDE")),
     main="ZT12 (4h dark)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(lux22[elf3Targets_de, "TPM_ZT.12.22c"]/col0[elf3Targets_de, "TPM_ZT.12.22c"]), 
       log(phy22_v2[elf3Targets_de, "TPM_ZT.12.22c"]/ler022[elf3Targets_de, "TPM_ZT.12.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(lux22[elf4Targets_de, "TPM_ZT.12.22c"]/col0[elf4Targets_de, "TPM_ZT.12.22c"]), 
       log(phy22_v2[elf4Targets_de, "TPM_ZT.12.22c"]/ler022[elf4Targets_de, "TPM_ZT.12.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(lux22[luxTargets_de, "TPM_ZT.12.22c"]/col0[luxTargets_de, "TPM_ZT.12.22c"]), 
       log(phy22_v2[luxTargets_de, "TPM_ZT.12.22c"]/ler022[luxTargets_de, "TPM_ZT.12.22c"]),
       col="goldenrod3", pch=19)

points(log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.12.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.12.22c"]), 
       log(phy22_v2[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.12.22c"]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.12.22c"]),
       col="forestgreen", pch=19)

legend(-6.3, 6.3, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), "blue3", "mediumvioletred", "goldenrod3", "forestgreen"), bty="n", cex=0.8, pch=19)

plot(log(lux22[compileTopGenes,6]/col0[compileTopGenes,6]), 
     log(phy22_v2[compileTopGenes,6]/ler022[compileTopGenes,6]),
     xlim=c(-6.3, 6.3), ylim=c(-6.3, 6.3), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("phyABCDE")),
     main="ZT4 (4h light)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(lux22[elf3Targets_de,6]/col0[elf3Targets_de,6]), 
       log(phy22_v2[elf3Targets_de,6]/ler022[elf3Targets_de,6]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(lux22[elf4Targets_de,6]/col0[elf4Targets_de,6]), 
       log(phy22_v2[elf4Targets_de,6]/ler022[elf4Targets_de,6]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(lux22[luxTargets_de,6]/col0[luxTargets_de,6]), 
       log(phy22_v2[luxTargets_de,6]/ler022[luxTargets_de,6]),
       col="goldenrod3", pch=19)

points(log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]), 
       log(phy22_v2[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],6]),
       col="forestgreen", pch=19)

#legend(-4.5, 4.5, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), rgb(0, 0.6, 0, 1), rgb(0, 0, 0.6, 1), rgb(0.9, 0, 0, 1), "orange"), bty="n", cex=0.8, pch=19)

dev.off()


#show how these genes are affected by phyq
pdf(paste("Figure_transcriptomeCauseOfTemperatureResponsiveness_v2_color.pdf", sep=""), width=7, height=3.5,pointsize = 10) 
par(mfrow=c(1,2));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

plot(log(ler027[compileTopGenes,"TPM_ZT.8.27c"]/ler022[compileTopGenes,"TPM_ZT.8.22c"]), 
     log(phy22_v2[compileTopGenes,"TPM_ZT.8.22c"]/ler022[compileTopGenes,"TPM_ZT.8.22c"]),
     xlim=c(-6.3, 6.3), ylim=c(-6.3, 6.3), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("ler")~"27"*degree*"C"),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("phyABCDE")),
     main="ZT16 (8h dark)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(ler027[elf3Targets_de, "TPM_ZT.8.27c"]/ler022[elf3Targets_de, "TPM_ZT.8.22c"]), 
       log(phy22_v2[elf3Targets_de, "TPM_ZT.8.22c"]/ler022[elf3Targets_de, "TPM_ZT.8.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(ler027[elf4Targets_de, "TPM_ZT.8.27c"]/ler022[elf4Targets_de, "TPM_ZT.8.22c"]), 
       log(phy22_v2[elf4Targets_de, "TPM_ZT.8.22c"]/ler022[elf4Targets_de, "TPM_ZT.8.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(ler027[luxTargets_de, "TPM_ZT.8.27c"]/ler022[luxTargets_de, "TPM_ZT.8.22c"]), 
       log(phy22_v2[luxTargets_de, "TPM_ZT.8.22c"]/ler022[luxTargets_de, "TPM_ZT.8.22c"]),
       col="goldenrod3", pch=19)

points(log(ler027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.27c"]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.22c"]), 
       log(phy22_v2[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.22c"]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.22c"]),
       col="forestgreen", pch=19)

#legend(-6, 6, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), rgb(0, 0.6, 0, 1), rgb(0, 0, 0.6, 1), rgb(0.9, 0, 0, 1), "orange"), bty="n", cex=0.8, pch=19)

plot(log(col027[compileTopGenes, "TPM_ZT.8.27c"]/col0[compileTopGenes, "TPM_ZT.8.22c"]), 
     log(lux22[compileTopGenes, "TPM_ZT.8.22c"]/col0[compileTopGenes, "TPM_ZT.8.22c"]),
     xlim=c(-6.3, 6.3), ylim=c(-6.3, 6.3), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("col0")~"27"*degree*"C"),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     main="ZT16 (8h dark)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de, "TPM_ZT.8.27c"]/col0[elf3Targets_de, "TPM_ZT.8.22c"]), 
       log(lux22[elf3Targets_de, "TPM_ZT.8.22c"]/col0[elf3Targets_de, "TPM_ZT.8.22c"]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de, "TPM_ZT.8.27c"]/col0[elf4Targets_de, "TPM_ZT.8.22c"]), 
       log(lux22[elf4Targets_de, "TPM_ZT.8.22c"]/col0[elf4Targets_de, "TPM_ZT.8.22c"]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de, "TPM_ZT.8.27c"]/col0[luxTargets_de, "TPM_ZT.8.22c"]), 
       log(lux22[luxTargets_de, "TPM_ZT.8.22c"]/col0[luxTargets_de, "TPM_ZT.8.22c"]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.27c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.22c"]), 
       log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.22c"]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))], "TPM_ZT.8.22c"]),
       col="forestgreen", pch=19)

#legend(-6.3, 6.3, c("all genes", "ELF3 only", "ELF3/ELF4 only", "ELF3/ELF4/LUX", "ELF3/LUX only"), col=c(rgb(0.3, 0.3, 0.3, 0.3), rgb(0, 0.6, 0, 1), rgb(0, 0, 0.6, 1), rgb(0.9, 0, 0, 1), "orange"), bty="n", cex=0.8, pch=19)

dev.off()

pdf(paste("Figure_S3.pdf", sep=""), width=10.4, height=5.2,pointsize = 10) 
par(mfrow=c(2,4));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)
labels=c("", "ZT20 (Night)", "ZT22 (Night)", "ZT0 (Day)", "ZT1 (Day)", "ZT4 (Day)", "ZT8 (Day)", "ZT12 (Night)", "ZT16 (Night)")
for(colID in c(2:9)){
  
  plot(log(lux22[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
       log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]),
       xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
       xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
       ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
       main=labels[colID]
  )
  abline(v=0)
  abline(h=0)
  abline(c(0,1), col="grey")
  
  points(log(lux22[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
         log(elf322[elf3Targets_de,colID]/col0[elf3Targets_de,colID]),
         col="blue3", pch=19)
  #subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
  points(log(lux22[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
         log(elf322[elf4Targets_de,colID]/col0[elf4Targets_de,colID]),
         col="mediumvioletred", pch=19)
  #subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
  points(log(lux22[luxTargets_de,colID]/col0[luxTargets_de,colID]), 
         log(elf322[luxTargets_de,colID]/col0[luxTargets_de,colID]),
         col="goldenrod3", pch=19)
  
  points(log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
         log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
         col="forestgreen", pch=19)
  
}

dev.off()

#lux17CTargets_de
####Now I need to make scatterplots that only highlight genes mentioned in the main text
circadian=c("AT1G22770",
            "AT3G46640",
            "AT5G02810",
            "AT2G46790",
            "AT2G46830",
            "AT5G62430",
            "AT5G39660")
names(circadian)=c("GI", "LUX", "PRR7", "PRR9", "CCA1", "CDF1", "CDF2")

growth=c("AT2G43010",
         "AT5G45830",
         "AT5G15160",
         "AT5G17300")
names(growth)=c("PIF4", "DOG1", "BNQ2", "RVE1")

temperature=c("AT3G54900",
              "AT4G25490",
              "AT4G25470",
              "AT2G40340",
              "AT4G36040")
names(temperature)=c("CXIP1", "CBF1", "CBF2", "DREB2C", "DNAJ11")

photosystem=c("AT3G54890",
              "AT2G46820",
              "AT4G02770",
              "AT1G03130",
              "AT2G20260",
              "AT1G31330",
              "AT1G55670",
              "AT1G52230",
              "AT4G12800",
              "AT1G08380",
              "AT1G67740",
              "AT2G05070",
              "AT3G27690",
              "AT2G34420",
              "AT2G30570",
              "AT4G05180",
              "AT1G79040",
              "AT4G10340",
              "AT5G01530",
              "AT3G08940",
              "AT1G15820")




pdf(paste("Figure_S3_byFunction.pdf", sep=""), width=10.4, height=5.2,pointsize = 10) 
par(mfrow=c(2,4));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)
labels=c("", "ZT20 (Night)", "ZT22 (Night)", "ZT0 (Day)", "ZT1 (Day)", "ZT4 (Day)", "ZT8 (Day)", "ZT12 (Night)", "ZT16 (Night)")
for(colID in c(2:9)){
  
  plot(log(lux22[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
       log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]),
       xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.7, 0.7, 0.7, 0.3), pch=19,
       xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
       ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
       main=labels[colID]
  )
  abline(v=0)
  abline(h=0)
  abline(c(0,1), col="grey")
  
  
  points(log(lux22[photosystem,colID]/col0[photosystem,colID]), 
         log(elf322[photosystem,colID]/col0[photosystem,colID]),
         col="lightgreen", pch=19)
  
  
  points(log(lux22[circadian,colID]/col0[circadian,colID]), 
         log(elf322[circadian,colID]/col0[circadian,colID]),
         col="blue3", pch=19)
  #text(log(lux22[circadian,colID]/col0[circadian,colID])-0.2, 
  #      log(elf322[circadian,colID]/col0[circadian,colID])-0.2, 
  #      labels=names(circadian), cex=0.6)
  
  #subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
  points(log(lux22[temperature,colID]/col0[temperature,colID]), 
         log(elf322[temperature,colID]/col0[temperature,colID]),
         col="mediumvioletred", pch=19)
  #text(log(lux22[temperature,colID]/col0[temperature,colID])-0.2, 
  #     log(elf322[temperature,colID]/col0[temperature,colID])-0.2, 
  #     labels=names(temperature), cex=0.6)
  
  #subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
  points(log(lux22[growth,colID]/col0[growth,colID]), 
         log(elf322[growth,colID]/col0[growth,colID]),
         col="goldenrod3", pch=19)
  # text(log(lux22[growth,colID]/col0[growth,colID])-0.2, 
  #       log(elf322[growth,colID]/col0[growth,colID])-0.2, 
  #      labels=names(growth), cex=0.6)
  
  
  if(colID==2){
    legend(-5.3, 5.3, c("photosystem I/II", 
                        "circadian", 
                        "temperature", 
                        "growth"), 
           col=c("lightgreen", "blue3", "mediumvioletred", "goldenrod3"), bty="n", cex=0.7, pch=19)
  }
  
  
  
}

dev.off()



pdf(paste("Figure_S3_byFunction_labeled.pdf", sep=""), width=5.2, height=5.2,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)
labels=c("", "ZT20 (Night)", "ZT22 (Night)", "ZT0 (Day)", "ZT1 (Day)", "ZT4 (Day)", "ZT8 (Day)", "ZT12 (Night)", "ZT16 (Night)")
for(colID in c(9)){
  
  plot(log(lux22[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
       log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]),
       xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.7, 0.7, 0.7, 0.3), pch=19,
       xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
       ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
       main=labels[colID]
  )
  abline(v=0)
  abline(h=0)
  abline(c(0,1), col="grey")
  
  
  points(log(lux22[photosystem,colID]/col0[photosystem,colID]), 
         log(elf322[photosystem,colID]/col0[photosystem,colID]),
         col="lightgreen", pch=19)
  
  
  points(log(lux22[circadian,colID]/col0[circadian,colID]), 
         log(elf322[circadian,colID]/col0[circadian,colID]),
         col="blue3", pch=19)
  text(log(lux22[circadian,colID]/col0[circadian,colID])-0.2, 
       log(elf322[circadian,colID]/col0[circadian,colID])-0.2, 
       labels=names(circadian), cex=0.6)
  
  #subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
  points(log(lux22[temperature,colID]/col0[temperature,colID]), 
         log(elf322[temperature,colID]/col0[temperature,colID]),
         col="mediumvioletred", pch=19)
  text(log(lux22[temperature,colID]/col0[temperature,colID])-0.2, 
       log(elf322[temperature,colID]/col0[temperature,colID])-0.2, 
       labels=names(temperature), cex=0.6)
  
  #subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
  points(log(lux22[growth,colID]/col0[growth,colID]), 
         log(elf322[growth,colID]/col0[growth,colID]),
         col="goldenrod3", pch=19)
  text(log(lux22[growth,colID]/col0[growth,colID])-0.2, 
       log(elf322[growth,colID]/col0[growth,colID])-0.2, 
       labels=names(growth), cex=0.6)
  
}

dev.off()










pdf(paste("Figure_S3ii_v2.pdf", sep=""), width=10.4, height=5.2,pointsize = 10) 
par(mfrow=c(2,4));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)
labels=c("", "ZT20 (Night)", "ZT22 (Night)", "ZT0 (Day)", "ZT1 (Day)", "ZT4 (Day)", "ZT8 (Day)", "ZT12 (Night)", "ZT16 (Night)")
for(colID in c(2:9)){
  
  plot(log(lux22[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
       log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]),
       xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
       xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
       ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
       main=labels[colID]
  )
  abline(v=0)
  abline(h=0)
  abline(c(0,1), col="grey")
  
  points(log(lux22[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
         log(elf322[elf3Targets_de,colID]/col0[elf3Targets_de,colID]),
         col="yellow", pch=19)
  #subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
  points(log(lux22[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
         log(elf322[elf4Targets_de,colID]/col0[elf4Targets_de,colID]),
         col=rgb(0, 0, 0.6, 1), pch=19)
  #subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
  points(log(lux22[lux17CTargets_de,colID]/col0[lux17CTargets_de,colID]), 
         log(elf322[lux17CTargets_de,colID]/col0[lux17CTargets_de,colID]),
         col=rgb(0.9, 0, 0, 1), pch=19)
  ####doubles and triple
  #elf3 and elf4 only: lightgreen
  #lux and elf4: doesn't exist
  intersection=elf3Targets_de[which(elf3Targets_de %in% elf4Targets_de)]
  points(log(lux22[intersection,colID]/col0[intersection,colID]), 
         log(elf322[intersection,colID]/col0[intersection,colID]),
         col="lightgreen", pch=19)
  intersection=lux17CTargets_de[which(lux17CTargets_de %in% elf4Targets_de)]
  points(log(lux22[intersection,colID]/col0[intersection,colID]), 
         log(elf322[intersection,colID]/col0[intersection,colID]),
         col="darkmagenta", pch=19)
  intersection=lux17CTargets_de[which(lux17CTargets_de %in% elf3Targets_de)]
  points(log(lux22[intersection,colID]/col0[intersection,colID]), 
         log(elf322[intersection,colID]/col0[intersection,colID]),
         col="orange", pch=19)
  #triple: dark magenta
  intersection=lux17CTargets_de[which((lux17CTargets_de %in% elf3Targets_de) & (lux17CTargets_de %in% elf4Targets_de))]
  points(log(lux22[intersection,colID]/col0[intersection,colID]), 
         log(elf322[intersection,colID]/col0[intersection,colID]),
         col="darkmagenta", pch=19)
  
  if(colID==2){
    legend(-5.3, 5.3, c("all genes", 
                        "ELF3 only", 
                        "LUX 17C only", 
                        "ELF3/ELF4/LUX 17C", 
                        "ELF3/LUX 17C only",
                        "ELF3/ELF4 only"), 
           col=c(rgb(0.3, 0.3, 0.3, 0.3), "yellow", "red", "darkmagenta", "orange", "lightgreen"), bty="n", cex=0.4, pch=19)
  }
}
dev.off()



############Draw relevant supplements for extra phytochrome analysis

pdf(paste("Figure_S6.pdf", sep=""), width=10.4, height=5.2,pointsize = 10) 
par(mfrow=c(2,3));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

colID=9
###########change in elf3 vs. col-0 at 27vs22:
plot(log(col027[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
     log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("col0")~"at "*27~degree*C),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     main="ZT16 (8 hrs darkness)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
       log(elf322[elf3Targets_de,colID]/col0[elf3Targets_de,colID]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
       log(elf322[elf4Targets_de,colID]/col0[elf4Targets_de,colID]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,colID]/col0[luxTargets_de,colID]), 
       log(elf322[luxTargets_de,colID]/col0[luxTargets_de,colID]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
       log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
       col="forestgreen", pch=19)

phy22_v2=read.table("alternativePhyABCDEsample_22C", header=TRUE)
phy27_v2=read.table("alternativePhyABCDEsample_27C", header=TRUE)
#remove duplicate rows
geneNames=substring(as.character(phy22_v2[,1]), 6, 14)
phy22_v2=phy22_v2[which(sapply(geneNames, function(i){(length(which(geneNames==i))==1)})),]
rownames(phy22_v2)=substring(as.character(phy22_v2[,1]), 6, 14)
phy27_v2=phy27_v2[which(sapply(geneNames, function(i){(length(which(geneNames==i))==1)})),]
rownames(phy27_v2)=substring(as.character(phy27_v2[,1]), 6, 14)
phy22_v2=phy22_v2[,-2]
phy27_v2=phy27_v2[,-2]
#name rows based on gene name (remember to parse the Gene:)
#redo everything
##########
#Now I need phyABCDE vs. elf3 at night and day
colID=9
###########change in elf3 vs. col-0 at 27vs22:
plot(log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
     log(phy22_v2[compileTopGenes,colID]/ler022[compileTopGenes,colID]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("phyABCDE")),
     main="ZT16 (8 hrs darkness)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(elf322[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
       log(phy22_v2[elf3Targets_de,colID]/ler022[elf3Targets_de,colID]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(elf322[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
       log(phy22_v2[elf4Targets_de,colID]/ler022[elf4Targets_de,colID]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(elf322[luxTargets_de,colID]/col0[luxTargets_de,colID]), 
       log(phy22_v2[luxTargets_de,colID]/ler022[luxTargets_de,colID]),
       col="goldenrod3", pch=19)

points(log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
       log(phy22_v2[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
       col="forestgreen", pch=19)



colID=6
###########change in elf3 vs. col-0 at 27vs22:
plot(log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
     log(phy22_v2[compileTopGenes,colID]/ler022[compileTopGenes,colID]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("phyABCDE")),
     main="ZT4 (4 hrs light)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(elf322[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
       log(phy22_v2[elf3Targets_de,colID]/ler022[elf3Targets_de,colID]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(elf322[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
       log(phy22_v2[elf4Targets_de,colID]/ler022[elf4Targets_de,colID]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(elf322[luxTargets_de,colID]/col0[luxTargets_de,colID]), 
       log(phy22_v2[luxTargets_de,colID]/ler022[luxTargets_de,colID]),
       col="goldenrod3", pch=19)

points(log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
       log(phy22_v2[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
       col="forestgreen", pch=19)


#######Now look at phyABCDE at 27 vs 22C in the day (temperature shouldn't influence it in the day)  
colID=6
###########change in elf3 vs. col-0 at 27vs22:
plot(log(ler027[compileTopGenes,colID]/ler022[compileTopGenes,colID]), 
     log(phy22_v2[compileTopGenes,colID]/ler022[compileTopGenes,colID]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("ler")~"at "*27~degree*C),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("phyABCDE")),
     main="ZT4 (4 hrs light)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(ler027[elf3Targets_de,colID]/ler022[elf3Targets_de,colID]), 
       log(phy22_v2[elf3Targets_de,colID]/ler022[elf3Targets_de,colID]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(ler027[elf4Targets_de,colID]/ler022[elf4Targets_de,colID]), 
       log(phy22_v2[elf4Targets_de,colID]/ler022[elf4Targets_de,colID]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(ler027[luxTargets_de,colID]/ler022[luxTargets_de,colID]), 
       log(phy22_v2[luxTargets_de,colID]/ler022[luxTargets_de,colID]),
       col="goldenrod3", pch=19)

points(log(ler027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
       log(phy22_v2[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/ler022[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
       col="forestgreen", pch=19)



#######Now look at lux at 27 vs 22C in the day (temperature shouldn't influence it in the day)  
colID=6
###########change in lux vs. col-0 at 27vs22:
plot(log(col027[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
     log(lux22[compileTopGenes,colID]/col0[compileTopGenes,colID]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("col0")~"at "*27~degree*C),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("lux-4")),
     main="ZT4 (4 hrs light)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
       log(lux22[elf3Targets_de,colID]/col0[elf3Targets_de,colID]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
       log(lux22[elf4Targets_de,colID]/col0[elf4Targets_de,colID]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,colID]/col0[luxTargets_de,colID]), 
       log(lux22[luxTargets_de,colID]/col0[luxTargets_de,colID]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
       log(lux22[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
       col="forestgreen", pch=19)



colID=6
###########change in elf3 vs. col-0 at 27vs22:
plot(log(col027[compileTopGenes,colID]/col0[compileTopGenes,colID]), 
     log(elf322[compileTopGenes,colID]/col0[compileTopGenes,colID]),
     xlim=c(-5, 5), ylim=c(-5, 5), col=rgb(0.3, 0.3, 0.3, 0.3), pch=19,
     xlab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("col0")~"at "*27~degree*C),
     ylab=expression("log-fold"~"change"~"in"~"expression"~"in"~italic("elf3-1")),
     main="ZT4 (4 hrs light)"
)
abline(v=0)
abline(h=0)
abline(c(0,1), col="grey")

points(log(col027[elf3Targets_de,colID]/col0[elf3Targets_de,colID]), 
       log(elf322[elf3Targets_de,colID]/col0[elf3Targets_de,colID]),
       col="blue3", pch=19)
#subset_elf3_ids= which(as.character(colTPM[,1]) %in% as.character(elf3_genes[,1]))
points(log(col027[elf4Targets_de,colID]/col0[elf4Targets_de,colID]), 
       log(elf322[elf4Targets_de,colID]/col0[elf4Targets_de,colID]),
       col="mediumvioletred", pch=19)
#subset_elf4_ids= which(as.character(colTPM[,1]) %in% as.character(elf4_genes[,1]))
points(log(col027[luxTargets_de,colID]/col0[luxTargets_de,colID]), 
       log(elf322[luxTargets_de,colID]/col0[luxTargets_de,colID]),
       col="goldenrod3", pch=19)

points(log(col027[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]), 
       log(elf322[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]/col0[luxTargets_de[which(!(luxTargets_de %in% elf4Targets_de))],colID]),
       col="forestgreen", pch=19)





dev.off()



pdf('Figure_EC_complex_effect_v5.pdf', width=7.4, height=7.4,pointsize = 12) 
par(mfrow=c(1,1));
par(mar=c(5, 5, 10, 10)+0.1);
par(cex=0.9);

# #22C
# logdiff=sapply(colnames(elf322)[2:length(colnames(elf322))], function(i){
#   log((elf322[as.character(elf4Targets[,1]),i])/(col0[as.character(elf4Targets[,1]),i]))})
# rownames(logdiff)=rownames(elf322[as.character(elf4Targets[,1]),])
# logdiff=logdiff[which(!is.na(logdiff[,2])),]
# logdiff1=logdiff[which(apply(logdiff, 1, function(i){!is.infinite(sum(i)) & !is.na(sum(i))})),]
# 
# #27C
# logdiff=sapply(colnames(elf327)[2:length(colnames(elf327))], function(i){
#   log((elf327[as.character(elf4Targets[,1]),i])/(col027[as.character(elf4Targets[,1]),i]))})
# rownames(logdiff)=rownames(elf327[as.character(elf4Targets[,1]),])
# logdiff2=logdiff[rownames(logdiff1),]
# 
# special_cycle=c("Dof-type TF", "BBX28", "PRR9", "PIF4", "PRR7", "GI") #, "Lux")
# names(special_cycle)=c("AT1G69570", "AT4G27310", "AT2G46790", "AT2G43010", "AT5G02810", "AT1G22770") #, "AT3G46640")
# 
# special_photo=c("LHCB6", "LHCB3", "LHCA4", "LHCB1B1", "LHCA3", "LHB1B2")
# names(special_photo)=c("AT1G15820", "AT5G54270", "AT3G47470", "AT2G34430", "AT1G61520", "AT2G34420")
# 
# logdiff_total=cbind(logdiff1, logdiff2)
# 
# #NOTE: I AM SPECIFICALLY REMOVING A SMALL uORF WHICH HAS A REALLY RANDOM PATTERN AND RUINS THE CLUSTERING
# logdiff_total=logdiff_total[which(rownames(logdiff_total)!="AT5G15948"),]
# cols=sapply(colnames(logdiff_total), function(i){sub(".", "-", substr(i, 5, nchar(i)-4), fixed=TRUE)})
# rows=sapply(rownames(logdiff_total), function(i){
#   if(as.character(i) %in% as.character(names(special_photo))){
#     special_photo[i]
#   }else{
#     if(as.character(i) %in% as.character(names(special_cycle))){
#       special_cycle[i]
#     }else{
#       i
#     }
#   }
# })
# 
# rownames(logdiff_total)=rows
# temp=c(rep("blue", 8), rep("red", 8))
# timeCol=c(rep("black", 2), rep("lightyellow", 4), rep("black", 4), rep("lightyellow", 4), rep("black", 2))
# heatmap.2(logdiff_total, Colv=NA, dendrogram="row", labCol=cols, colCol=temp, ColSideColors=timeCol, trace="none", key.title="", key.ylab="frequency", key.xlab=expression("expression"~"change"~"in"~italic("elf3")))

#22C
logdiff=sapply(colnames(elf322)[2:length(colnames(elf322))], function(i){
  log((elf322[elf4Targets_de,i])/(col0[elf4Targets_de,i]))})
rownames(logdiff)=rownames(elf322[elf4Targets_de,])
logdiff=logdiff[which(!is.na(logdiff[,2])),]
logdiff1=logdiff[which(apply(logdiff, 1, function(i){!is.infinite(sum(i)) & !is.na(sum(i))})),]

#27C
logdiff=sapply(colnames(elf327)[2:length(colnames(elf327))], function(i){
  log((elf327[elf4Targets_de,i])/(col027[elf4Targets_de,i]))})
rownames(logdiff)=rownames(elf327[elf4Targets_de,])
logdiff2=logdiff[rownames(logdiff1),]

special_cycle=c("Dof-type TF", "BBX28", "PRR9", "PIF4", "PRR7", "GI", "PIF5", "Lux")
names(special_cycle)=c("AT1G69570", "AT4G27310", "AT2G46790", "AT2G43010", "AT5G02810", "AT1G22770", "AT3G59060", "AT3G46640")

special_photo=c("LHCB6", "LHCB3", "LHCA4", "LHCB1B1", "LHCA3", "LHB1B2")
names(special_photo)=c("AT1G15820", "AT5G54270", "AT3G47470", "AT2G34430", "AT1G61520", "AT2G34420")

logdiff_total=cbind(logdiff1, logdiff2)

#NOTE: I AM SPECIFICALLY REMOVING A SMALL uORF WHICH HAS A REALLY RANDOM PATTERN AND RUINS THE CLUSTERING
logdiff_total=logdiff_total[which(! (rownames(logdiff_total) %in% c("AT5G15948", "AT5G64190"))),]
cols=sapply(colnames(logdiff_total), function(i){sub(".", "-", substr(i, 5, nchar(i)-4), fixed=TRUE)})
rows=sapply(rownames(logdiff_total), function(i){
  if(as.character(i) %in% as.character(names(special_photo))){
    special_photo[i]
  }else{
    if(as.character(i) %in% as.character(names(special_cycle))){
      special_cycle[i]
    }else{
      i
    }
  }
})

rownames(logdiff_total)=rows
temp=c(rep("blue", 8), rep("red", 8))
rowCols=c("cornflowerblue", "darkmagenta", "grey")[cutree(hclust(dist(logdiff_total)), k=3)]
timeCol=c(rep("black", 2), rep("lightyellow", 4), rep("black", 4), rep("lightyellow", 4), rep("black", 2))
cols=c("ZT20", "ZT22", "ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16","ZT20", "ZT22", "ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16")
heatmap.2(logdiff_total, Colv=NA, dendrogram="row", labCol=cols, margins=c(12,8), keysize=1.9, colCol=temp, ColSideColors=timeCol, RowSideColors=rowCols, trace="none", key.title="", key.ylab="frequency", key.xlab=expression("expression"~"change"~"in"~italic("elf3-1")))

dev.off()

##########Supplement7: temperature change in elf3-1 and phyABCDE
pdf(paste("Figure_S7a.pdf", sep=""), width=5, height=5,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

arr=as.matrix(log(lux27[elf4Targets_de,2:9])-log(lux22[elf4Targets_de,2:9]))
arr=arr[which(apply(arr, 1, function(i){!is.na(sum(as.numeric(i))) & !is.infinite(sum(as.numeric(i)))})), ]
rownames(arr)[which(rownames(arr) %in% names(special_photo))]=special_photo[rownames(arr)[which(rownames(arr) %in% names(special_photo))]]
rownames(arr)[which(rownames(arr) %in% names(special_cycle))]=special_cycle[rownames(arr)[which(rownames(arr) %in% names(special_cycle))]]


timeCol=c(rep("black", 2), rep("lightyellow", 4), rep("black", 4), rep("lightyellow", 4), rep("black", 2))
cols=c("ZT20", "ZT22", "ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16","ZT20", "ZT22", "ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16")
heatmap.2(arr, Colv=NA, dendrogram="row", labCol=cols[1:8], margins=c(12,8), ColSideColors=timeCol[1:8], trace="none", key.title="", key.ylab="frequency", key.xlab=expression("expression"~"change"~"in"~italic("lux-4")~"at "*27~degree*C))

dev.off()

pdf(paste("Figure_S7b.pdf", sep=""), width=5, height=5,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 2, 1)+0.1);
par(cex=0.9)

arr=as.matrix(log(phy27_v2[elf4Targets_de,2:9])-log(phy22_v2[elf4Targets_de,2:9]))
arr=arr[which(apply(arr, 1, function(i){!is.na(sum(as.numeric(i))) & 
    !is.infinite(sum(as.numeric(i)))})), ]
rownames(arr)[which(rownames(arr) %in% names(special_photo))]=special_photo[rownames(arr)[which(rownames(arr) %in% names(special_photo))]]
rownames(arr)[which(rownames(arr) %in% names(special_cycle))]=special_cycle[rownames(arr)[which(rownames(arr) %in% names(special_cycle))]]


timeCol=c(rep("black", 2), rep("lightyellow", 4), rep("black", 4), rep("lightyellow", 4), rep("black", 2))
cols=c("ZT20", "ZT22", "ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16","ZT20", "ZT22", "ZT0", "ZT1", "ZT4", "ZT8", "ZT12", "ZT16")
heatmap.2(arr, Colv=NA, dendrogram="row", labCol=cols[1:8], margins=c(12,8), ColSideColors=timeCol[1:8], trace="none", key.title="", key.ylab="frequency", key.xlab=expression("expression"~"change"~"in"~italic("phyABCDE")~"at "*27~degree*C))

dev.off()



