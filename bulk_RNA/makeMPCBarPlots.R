#load basicPlots functions
source("/BcKO_RNA/analysis/barPlots/basicPlots.R")

#set working directory
homeDir = "/BcKO_RNA/analysis/MPC.BarPlots/";
setwd(homeDir);

#load sample manifest
fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
samples = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
samples = samples[samples$include, ]

#load gene expression data to plot
genesFile = "Mouse.total.mpc.txt";
genes = read.table(paste0(homeDir, genesFile), sep = "\t", header = T)

#color scheme
#set to be the same as the group order
ko.col = rgb(255,51,51, maxColorValue = 255)
wt.col = rgb(0,204,255, maxColorValue = 255)

cols = rep(c(wt.col, ko.col), 6);

#set group order for plotting, this should match above order
grpOrder = c("WT_div0", "KO_div0", "WT_div1", "KO_div1", "WT_div3", "KO_div3", "WT_div5", "KO_div5",  "WT_div8neg", "KO_div8neg", "WT_div8pos")

#choose plot type
type = "both" #both, bee, bar

#plot options
ylab = "MPC" #y-axis label
cex.dots = 0.5 #beeswarm point size

i = 1

#print gene name in progress
print(paste(genes$Mouse[i]));

#select gene specific data
eGene = genes[i, match(samples$sample, colnames(genes))]

#set group order
grps = samples$group[match(samples$sample, names(eGene))]
grps = factor(grps, levels = grpOrder, ordered = T)

#call barPlot function
bbPlot(as.numeric(eGene), grps, type = type, main = "MPC/cell", outFile = "MPC.per.cell.bbplot.pdf", cols = cols, 
       cex.dots = cex.dots, ylab = ylab, ylim = c(0,6000000000000), names = grpOrder)
