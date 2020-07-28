#calculate ERCC read coverage

source("bisTools.R");

homeDir = "/BcKO_RNA/coverage/";
setwd(homeDir);

fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);

#extract ERCCs from bamfile header
bamFile = paste0(files$dir[1], files$bamFile[1])	
si = seqinfo(BamFile(bamFile))

erccChrs = si@seqnames[grepl("ERCC", si@seqnames)]
erccCts = matrix(0, ncol = dim(files)[1], nrow = length(erccChrs), dimnames = list(erccChrs, files$sample))

for (i in 1:dim(files)[1]) {

	bamFile = paste0(files$dir[i], files$bamFile[i])	
	si = seqinfo(BamFile(bamFile))

	for (chr in erccChrs) {

		#read in paired alignment
		sbp = ScanBamParam(which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si@seqlengths[si@seqnames == chr])))

		reads = readGAlignments(bamFile, param = sbp)

		erccCts[chr, files$sample[i]] = length(reads)

	}

	print(paste(files$sample[i], Sys.time()))
}

write.table(erccCts, file = "erccCts.csv", sep = ",");
