source("bisTools.R");

#set working directory
homeDir = "/BcKO_RNA/coverage/";
setwd(homeDir);

#read in manifest file
fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[files$include, ]

#ERCC concentration information in each mix 
erccInfoDir = "/BcKO_RNA/coverage/";
erccInfoFile = "ERCC.concentrations.txt";
erccInfo = read.table(paste0(erccInfoDir, erccInfoFile), sep = "\t", header = T)

#read in coverage files
geneCtsFile = "geneCts.exon.csv"
geneCts = read.table(geneCtsFile, sep = ",", header = T)

geneRpkmFile = "geneRpkm.exon.csv"
geneRpkm = read.table(geneRpkmFile, sep = ",", header = T)

# Sum read counts in each sample
readCts = data.frame(reads = colSums(geneCts[grepl(paste0(files$sample, collapse = "|"), names(geneCts))]))


#####
#normalize ERCC coverage

#read in ERCC coverage file
erccCtsFile = "erccCts.csv";
erccCts = read.table(erccCtsFile, sep = ",", header = T)

#fix this in previous ERCC coverage script!!!
#this file really needs an ID with each ERCC to make sure they stay annotated right
#add unique column with ERCC id
erccCts = cbind(ERCC.ID = erccInfo$ERCC.ID, erccCts)

# Sum of of ercc counts in each sample
erccPercRna = data.frame(reads = colSums(erccCts[grepl(paste0(files$sample, collapse = "|"), names(erccCts))]))

# ercc table has counts of ercc in each rnaseq sample and the ercc conc. info
setkey(as.data.table(erccInfo), ERCC.ID)
setkey(as.data.table(erccCts), ERCC.ID)
ercc = merge(erccInfo, erccCts) 

#rpm normalize ercc
totalReads = erccPercRna$reads + readCts$reads
erccRpm = erccCts[, grepl(paste(files$sample, collapse = "|"), colnames(erccCts))]
erccRpm = data.frame(t(t(erccRpm) * 1e6 / totalReads))
colnames(erccRpm) = paste0(files$sample, ".rpm")
erccRpm = cbind(erccCts$ERCC.ID, erccRpm)
colnames(erccRpm)[1] = "ERCC.ID"
#ercc = merge(ercc, erccRpm, by= "ERCC.ID")

#fpkm normalize ercc (normalizing rpm for length(per kilobase) of the transcript
#erccFpkm = erccRpm * 1e3 / ercc$length[match(row.names(erccRpm), ercc$ERCC.ID)]
erccFpkm = erccRpm[, grepl(paste(files$sample, collapse = "|"), colnames(erccRpm))]
erccFpkm = data.frame(erccFpkm * 1e3 / ercc$length)
colnames(erccFpkm) = gsub(".rpm", ".fpkm", colnames(erccFpkm))
erccFpkm = cbind(erccCts$ERCC.ID, erccFpkm)
colnames(erccFpkm)[1] = "ERCC.ID"
#ercc = merge(ercc, erccFpkm, by= "ERCC.ID")

write.table(ercc, gsub(".csv$", ".rpm.fpkm.csv", erccCtsFile), sep = ",", row.names = F, quote = F)

#######
#mpc normalize gene coverage file

erccPercRna = cbind(erccPercRna, readCts)


#set up MPC matrix and loop through to calculate for each sample
mpc = geneRpkm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpkm))]

erccMpc = ercc[, !grepl(".rpm", names(ercc)) & !(names(ercc) %in% files$sample)]
names(erccMpc) = gsub(".fpkm", "", names(erccMpc))

for (sample in files$sample) {
	print(sample)
	
	#identify sample specific SpikeIn info
	cells = files$cells[which(files$sample == sample)]
	spikeinDilution = files$spikeinVolume[which(files$sample == sample)] / files$spikeinDilution[which(files$sample == sample)]
	cellsSorted = files$cellsSorted[which(files$sample == sample)]
	
	#calculate ercc molecules added
	#ercc$spikein.molecules.per.cell = erccInfo$mix1.molecules.per.ul * spikeinDilution / cells 
	ercc$spikein.molecules.per.cell = (erccInfo$mix1.molecules.per.ul * spikeinDilution / cellsSorted) * (cells/cellsSorted)
	
	#convert sample fpkm to mRNAs per cell (mpc)
	mpc[[sample]] = mpc[[sample]] * sum(ercc$spikein.molecules.per.cell) / sum(erccFpkm[[paste0(sample, ".fpkm")]]); ## converts based on total fpkm / molecules.
	erccMpc[[sample]] = erccMpc[[sample]] * sum(ercc$spikein.molecules.per.cell) / sum(erccFpkm[[paste0(sample, ".fpkm")]]); ## converts based on total fpkm / molecules.

}


#make QC plots
noOfErccDetected = NULL
noOfGenesDetected = NULL

######### plot ERCC controls
for (sample in files$sample) {
	print(sample)

    # to plot ercc that did not get detected
	tempErccFpkm = ercc[[paste0(sample, ".fpkm")]]
	tempErccFpkm = log(tempErccFpkm,2)
	tempErccMpc = ercc$spikein.molecules.per.cell
	tempErccMpc = log(tempErccMpc,2)
	mpcInf = tempErccMpc[which(tempErccFpkm == -Inf)]
	fpkmInf = tempErccFpkm[which(tempErccFpkm == -Inf)]
	fpkmInf = replace(fpkmInf, fpkmInf == -Inf, 0)
	tempErccMpc = tempErccMpc[which(tempErccFpkm != -Inf)]
	tempErccFpkm = tempErccFpkm[which(tempErccFpkm != -Inf)]
	
	noOfErccDetected[sample] = length(tempErccFpkm)
 	noOfGenesDetected[sample] = length(grep("TRUE", geneCts[,sample] > 0))
 	
	#plot ercc ctrls
	cairo_pdf(paste0("plots/ercc.", sample, ".pdf"), height = 5, width = 5)
	#png(file=paste0("plots/ercc.", sample, ".png"), height = 5, width = 5, units = "in", res=300);
	par(mai = c(0.75, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
	xlim <- c(0, round(max(tempErccFpkm), 0))
	plot(tempErccFpkm , tempErccMpc, xlab = "fpkm (log2)", ylab = "Ercc Spike in Mols/cell (log2)", pch = 19, cex = 0.5, cex.lab = 0.6, cex.axis = 0.6, font.lab = 2, xlim = xlim)
	points(fpkmInf, mpcInf, col = "red", pch = 19, cex = 0.5)
	lm = lm(tempErccMpc[is.finite(tempErccMpc) & is.finite(tempErccFpkm)] ~ tempErccFpkm[is.finite(tempErccMpc) & is.finite(tempErccFpkm)])
	abline(lm)
	legend("bottomleft", paste("y =", round(lm$coefficients[2], 3), "x + ", round(lm$coefficients[1], 3)), cex = 0.6, bty = "n")
	dev.off()
	
}


write.table(mpc, gsub(".csv$", ".mRNAsPerCell.csv", geneCtsFile), sep = ",", row.names = F, quote = F)
write.table(erccMpc, "ercc.mRNAsPerCell.csv", sep = ",", row.names = F, quote = F)

cairo_pdf(paste0("plots/No.of.92ERCCs.detected.pdf"), height = 5, width = 5)
par(mar=c(7,4,4,2))
barplot(noOfErccDetected, names.arg = names(noOfErccDetected), las = 2, ylim = c(0,92), ylab = "No. of 92-ERCCs Detected")
dev.off()

cairo_pdf(paste0("plots/No.of.Genes.detected.pdf"), height = 5, width = 5)
par(mar=c(7,4,4,2))
barplot(noOfGenesDetected, names.arg = names(noOfGenesDetected), las = 2, ylab = "No. of Genes Detected")
dev.off()


cairo_pdf(paste0("plots/rnaseq.molecules.per.cell.pdf"), height = 5, width = 5)
par(mar=c(7,4,4,2))
barplot(colSums(mpc), names.arg = colnames(mpc), las = 2)
dev.off()

cairo_pdf(paste0("plots/ercc.molecules.per.cell.pdf"), height = 5, width = 5)
par(mar=c(7,4,4,2))
barplot(colSums(erccMpc[c(12:dim(erccMpc)[2])]), names.arg = colnames(erccMpc[c(12:dim(erccMpc)[2])]), las = 2)
dev.off()

dfErccPerc = as.data.frame(erccPercRna)
dfErccPerc[,3] = dfErccPerc[,1]/(dfErccPerc[,1]+dfErccPerc[,2])*100
cairo_pdf(paste0("plots/PercentOfErccInReads.pdf"), height = 5, width = 5)
par(mar=c(7,4,4,2))
barplot(dfErccPerc[,3], names.arg = rownames(dfErccPerc), las = 2, main = "Percent Of Ercc In reads", ylab = "Percentage")
dev.off()

############
#to delete

#Calculating the exon size
tx = TxDb.Mmusculus.UCSC.mm9.knownGene
exons = exonsBy(tx, by = "gene")
geneExonSize = lapply(lapply(lapply(exons, reduce), width), sum)
