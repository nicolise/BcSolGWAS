#Make plotting variables
Top50SNP$Index = NA
ticks = NULL
lastbase = 0
Top50SNPhold <- Top50SNP
#subsample DF to work on script
Top50SNPpr <- Top50SNP[sample(nrow(Top50SNP), 20),]

Top50SNPpr$Contig <- paste(Top50SNPpr$Chrom, Top50SNPpr$Segment, sep=".")
Top50SNPpr$Contig <-recode(Top50SNPpr$Contig, "1.0=1;2.0=2;3.1=3")

#Redo the positions to make them sequential		-- accurate position indexing
#go through each Chromosome sequentially (i from 1 to 16)
for (i in unique(Top50SNPpr$Chrom)) {
  print(i)
  #special rules for Chromosome 1 because counting from zero
  lastbase=lastbase+max(subset(Top50SNPpr,Top50SNPpr$Chrom==i-1)$Pos, 1)
  Top50SNPpr[Top50SNPpr$Chrom==i, ]$Index=Top50SNPpr[Top50SNPpr$Chrom==i, ]$Pos+lastbase
}
ticks=c(ticks, Top50SNPpr[Top50SNPpr$Chrom==i, ]$Index[floor(length(Top50SNPpr[Top50SNPpr$Chrom==i, ]$Index)/2)+1])
}
ticklim=c(min(Top50SNPpr$Index),max(Top50SNPpr$Index))