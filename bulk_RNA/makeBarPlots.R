#cycle through groups and summarize as barplots in a "plots/" subfolder

#load basicPlots functions
source("/BcKO_RNA/analysis/barPlots/basicPlots.R")

#set working directory
homeDir = "/BcKO_RNA/analysis/barPlots/";
setwd(homeDir);

#load sample manifest
fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
samples = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
samples = samples[samples$include, ]

#load gene expression data to plot
genesDir = "/BcKO_RNA/analysis/diff/";
genesFile = "BcKOvWt.Div.RNAseq.diff.significant.glm.Bcko_RNA.fdr.txt";
genes = read.table(paste0(genesDir, genesFile), sep = "\t", header = T)

#load and use MPC expression info
mpcFile = "/BcKO_RNA/coverage/geneCts.exon.mRNAsPerCell.csv"
mpc = read.table(mpcFile, sep = ",", header = T)

#color scheme
#set to be the same as the group order
wt.col = rgb(255,51,51, maxColorValue = 255)
ko.col = rgb(0,204,255, maxColorValue = 255)

cols = rep(c(wt.col, ko.col), 6);

#set group order for plotting, this should match above order
grpOrder = c("WT_div0", "KO_div0", "WT_div1", "KO_div1", "WT_div3", "KO_div3", "WT_div5", "KO_div5", "WT_div8neg", "KO_div8neg", "WT_div8pos")

type = "both"
cex = 0.5

#cycle through all genes and plot
for (i in 1:dim(genes)[1]) {

	#print gene name in progress
	print(paste(genes$SYMBOL[i]));
	
	#select gene specific data
	eGene = genes[i, match(paste0(samples$sample, ".rpkm"), colnames(genes))]
	
	#set group order
	grps = samples$group[match(paste0(samples$sample, ".rpkm"), names(eGene))]
	grps = factor(grps, levels = grpOrder, ordered = T)
	
	#call barPlot function
	bbPlot(as.numeric(eGene), grps, las = 3, mai = c(1.2, 0.6, 0.2, 0.1), mgp = c(1.5, 0.5, 0.2), main = genes$SYMBOL[i], 
		ylab = "FPKM", font.lab = 2, outFile = paste0("plots/", genes$SYMBOL[i], ".bbPlot.pdf"), cols = cols, cex.dots = cex, names = grpOrder)
	
}

##################################
# plot MPC info

#cycle through all genes and plot 
for (i in 1:dim(genes)[1]) {

	#print gene name in progress
	print(paste(genes$SYMBOL[i]));
	
	#select gene specific data
	eGene = mpc[mpc$ENTREZID == genes$ENTREZID[i], match(samples$sample, colnames(mpc))]
	
	#set group order
	grps = samples$group[match(samples$sample, names(eGene))]
	grps = factor(grps, levels = grpOrder, ordered = T)
	
	#call barPlot function
	barPlot(as.numeric(eGene), grps, las = 3, mai = c(1.2, 0.6, 0.2, 0.1), mgp = c(1.5, 0.5, 0.2), main = genes$SYMBOL[i], 
		ylab = "MPC", font.lab = 2, outFile = paste0("plots/", genes$SYMBOL[i], ".barplot.MPC.pdf"), cols = cols)
	
}
