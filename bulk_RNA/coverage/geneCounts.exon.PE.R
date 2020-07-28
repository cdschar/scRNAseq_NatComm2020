library("Rsamtools")
library("ShortRead")
library("GenomicRanges")
library("GenomicFeatures")

#load bistool functions
source("bisTools.R");

experimentName = "BcKO.Div.RNAseq"

#set working directory
homeDir = "/BcKO_RNA/coverage/";
setwd(homeDir);

#read in manifest file
fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[files$include,]

#set genome specific options
genome = "mm9"
if ( genome == "mm9"){
	library("org.Mm.eg.db")
	library("TxDb.Mmusculus.UCSC.mm9.knownGene")
	tx = TxDb.Mmusculus.UCSC.mm9.knownGene
	org = org.Mm.eg.db
} else if ( genome == "hg19"){
	library("org.Hs.eg.db")
	library("TxDb.Hsapiens.UCSC.hg19.knownGene")
	tx = TxDb.Hsapiens.UCSC.hg19.knownGene
	org = org.Hs.eg.db
} else if ( genome == "mm10"){
	library("org.Mm.eg.db")
	library("TxDb.Mmusculus.UCSC.mm10.knownGene")
	tx = TxDb.Mmusculus.UCSC.mm10.knownGene
	org = org.Mm.eg.db
} else if ( genome == "hg38"){
	library("org.Hs.eg.db")
	library("TxDb.Hsapiens.UCSC.hg38.knownGene")
	tx = TxDb.Hsapiens.UCSC.hg38.knownGene
	org = org.Hs.eg.db
}

#set up exon and gene tables for annotation of reads
genes = genes(tx)
exons = exonsBy(tx, by = "gene")

geneCts = matrix(0, ncol = dim(files)[1], nrow = length(exons), dimnames = list(names(exons), files$sample))

for (i in 1:dim(files)[1]) {

	print(paste(files$sample[i], Sys.time()))

	si = seqinfo(BamFile(paste0(files$dir[i], files$bamFile[i])))

	sbp = ScanBamParam(flag = scanBamFlag(isDuplicate = F))

	#read in paired alignment
	reads = readGAlignmentPairs(paste0(files$dir[i], files$bamFile[i]), param = sbp)

	#Get counts
	so = summarizeOverlaps(exons, reads, mode = "IntersectionNotEmpty", singleEnd = F, ignore.strand = T)
	cts = assays(so)$counts

	geneCts[rownames(cts), files$sample[i]] = cts
	
	rm(reads)
	gc()
}


#Get gene symbol from entrez id row names
geneSym = select(org, columns = "SYMBOL", keytype = "ENTREZID", keys = row.names(geneCts))

geneCts = cbind(geneSym, geneCts)
write.table(geneCts, file = "geneCts.exon.csv", sep = ",", quote = F, row.names = F)

#add in gene length info
width = sum(width(exons))
geneCts$length = width

#make RPM normalized gene counts file
geneRpm = geneCts
geneRpm = geneRpm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))]
colnames(geneRpm) = paste0(files$sample, ".rpm")
geneRpm = t(t(geneRpm) * 1e6 / colSums(geneRpm))
geneRpm = round(geneRpm, 3)
geneRpm = cbind(geneCts[, !grepl(paste(files$sample, collapse = "|"), colnames(geneCts))], geneRpm)

#make RPKM normalized gene counts file
geneRpkm = geneRpm
geneRpkm = geneRpkm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpkm))]
colnames(geneRpkm) = gsub(".rpm", ".rpkm", names(geneRpkm))
geneRpkm = (geneRpkm * 1000)/ geneCts$length
geneRpkm = round(geneRpkm, 3)
geneRpkm = cbind(geneCts[, !grepl(paste(files$sample, collapse = "|"), colnames(geneCts))], geneRpkm)

#plot distribution of RPM
plot = geneRpm[,grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))]

#cutoffs
rpm10 = log2(10+0.01)
rpm5 = log2(5+0.01)
rpm3 = log2(3+0.01)

cairo_pdf(file = paste0(experimentName, ".Transcript.Density.plot.pdf"), height = 5, width = 5);
	par(mai = c(1.5, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
	plot(NA, ylim = c(0, 0.2), xlim = c(-10,15), xlab = "log2 RPM+0.01", ylab = "density")
	#annotate where 3 cut offs are
	lines(c(rpm10, rpm10), c(0, 0.2), lty = 2, col = "blue", lwd = 1)
	lines(c(rpm5, rpm5), c(0, 0.2), lty = 2, col = "blue", lwd = 1)
	lines(c(rpm3, rpm3), c(0, 0.2), lty = 2, col = "blue", lwd = 1)
	text(rpm10+0.8, 0.15, labels = "10", col = "blue", cex = 0.8)
	text(rpm5+0.5, 0.12, labels = "5", col = "blue", cex = 0.8)
	text(rpm3-0.5, 0.11, labels = "3", col = "blue", cex = 0.8)
	#plot data
	for(i in 1:dim(plot)[2]){
		lines(density(log2(plot[,i]+0.01)), lwd = 0.2)
		}
dev.off();

#use cutoff of 3, filter for detected gene based on groups
cutoff = 3
detected = apply(geneRpm[,grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))], 1, function(x) filter_detected(x, files$group, cutoff) )
geneCts_detected = geneCts[detected,]
write.table(geneCts_detected, file = paste0(experimentName, ".geneCts.detected.RPM.", cutoff, ".exon.csv"), sep = ",", quote = F, row.names = F)

geneRpm_detected = geneRpm[detected,]
write.table(geneRpm_detected, file = paste0(experimentName, ".geneRpm.detected.RPM.", cutoff, ".exon.csv"), sep = ",", quote = F, row.names = F);

geneRpkm_detected = geneRpkm[detected,]
write.table(geneRpkm_detected, file = paste0(experimentName, ".geneRpkm.detected.RPM.", cutoff, ".exon.csv"), sep = ",", quote = F, row.names = F);

