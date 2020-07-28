library("som")
library("data.table")
#install.packages("vegan")
library("vegan");


setwd("/BcKO_RNA/analysis/pca/")

#sample manifest
fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste(fqDir, fqFile, sep = ""), sep = "\t", header = TRUE, as.is = TRUE);
files = files[files$include, ]

#sort by group
files = files[order(files$group), ]

#read in RNA file
dataDir = "/BcKO_RNA/analysis/diff/";
dataFile = "BcKOvWt.Div.RNAseq.diff.significant.glm.Bcko_RNA.fdr.txt";
data = read.table(file = paste0(dataDir, dataFile), sep = "\t", header = T, quote = "", comment.char = "");

#plot colors
D0_WT= "red"
D1_WT= "green"
D3_WT= "blue"
D5_WT= "yellow"
D8_neg_WT= "pink"
D8_pos_WT= "dark green"
D0_KO= "orange"
D1_KO= "magenta"
D3_KO= "light blue"
D5_KO= "black"
D8_neg_KO= "gray"
wt.col = rgb(255,51,51, maxColorValue = 255)
ko.col = rgb(0,204,255, maxColorValue = 255)

cols = c(rep(ko.col, 19), rep(wt.col, 11))
grpcol = c(ko.col, wt.col)

grpcols = function(colnms) {if(grepl("WT_div0", colnms)) wt.col else if (grepl("WT_div1", colnms)) wt.col
							else if (grepl("WT_div3", colnms)) wt.col else if (grepl("WT_div5", colnms)) wt.col
							else if (grepl("WT_div8neg", colnms)) wt.col else if(grepl("WT_div8pos", colnms)) wt.col
							else if(grepl("KO_div0", colnms)) ko.col else if (grepl("KO_div1", colnms)) ko.col
							else if (grepl("KO_div3", colnms)) ko.col else if (grepl("KO_div5", colnms)) ko.col
							else if (grepl("KO_div8neg", colnms)) ko.col }
#apply color function
cols = unlist(lapply(files$group, grpcols))

#set colors for confidence colors
ordiel.col = c(D0_KO, D1_KO, D3_KO, D5_KO, D8_neg_KO, D0_WT, D1_WT, D3_WT, D5_WT, D8_neg_WT, D8_pos_WT)
ordiel.col = c(rep(ko.col, 5), rep(wt.col, 6))

#select data rpm data
data.pca = data[ ,grep(paste(files$sample, collapse = "|"), names(data))]

#remove extensions for headers
colnames(data.pca) = gsub(".rppm|.rpm|.fpkm|.rpkm", "", names(data.pca))
#make the same order as the manifest
data.pca = data.pca[,as.vector(files$sample)]

#z-score normalize data by row
range(data.pca);
data.norm <- som::normalize(data.pca, byrow=TRUE);
range(data.norm);

#perform PCA analysis
pca = rda(data.norm, scale = F)

#extract variance
pca.var = eigenvals(pca)/sum(eigenvals(pca))

#select components to plot
p1 = 1;
p2 = 2;

#axis labels with variance
xlab = paste0("PC", p1, " ", signif(pca.var[p1] * 100, digits = 3), "%")
ylab = paste0("PC", p2, " ", signif(pca.var[p2] * 100, digits = 3), "%")

#access ordination values for select pc
pcas = scores(pca, choices = c(p1, p2))

#plot
cairo_pdf(paste0("Blimp1cKO.PCA.fig3.pc", p1, ".v.pc", p2, ".pdf"), height = 5, width = 5)
par(mai = c(0.8, 0.8, 0.5, 0.5), mgp = c(2, 0.8, 0), family = "Arial")

#pca plot
#plot(pcas$species[, 1], pcas$species[, 2], pch = 19, xlab = xlab, ylab = ylab, cex.main = 2, cex.lab = 0.8, font.lab = 2, cex.axis = 0.6, xlim = c(-6,6), ylim = c(-5,7), col = cols, main = paste("PCA", dim(data)[1]), cex = 1.5)
plot(pcas$species[, 1], pcas$species[, 2], pch = 19, xlab = xlab, ylab = ylab,
     cex.main = 2, cex.lab = 0.8, font.lab = 2, cex.axis = 0.6, col = cols,
     main = paste("PCA", dim(data)[1]), cex = 1.5)
#add in 99% confidence intervals
ordiellipse(pca, as.character(files$group), conf = 0.99, kind = "se", display = "species", choices = c(p1, p2), col = ordiel.col)
dev.off()
