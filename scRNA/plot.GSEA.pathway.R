#set working directory
homeDir = "~/GSEA/";
setwd(homeDir);

#read in collapsed PreRanked Gene list used for GSEA (resides in the "edb" folder in the output directory)(
geneListDir = "~/GSEA/";
geneListFile = "PreRanked.geneList.rnk";
geneList = read.table(paste0(geneListDir, geneListFile), sep = "\t")
geneList = geneList[order(geneList[, 2], decreasing = T), ]

#set colors for each group if desired
upCol = NA
dwnCol = NA

#set group names if desired
upPheno = NA
dwnPheno = NA

###
#point to .xls file for the gene set to plot
#example here is for the ISHIDA_E2F_TARGETS
dataDir = "~/GSEA/";
dataFile = "REACTOME_INTERFERON_SIGNALING.xls";
data = read.table(paste0(dataDir, dataFile), sep = "\t", header = T, comment.char="", quote="")

#set columns names to parse file
esCol = "RUNNING.ES"
rankCol = "RANK.IN.GENE.LIST"

#set axis range based on length of PreRanked gene list length
xlim = c(0, dim(geneList)[1])

#Edit ylim depending on whether they are positively (0, 1) or negatively (-1, 0) enriched
ylim = c(-1, 1)

#set color ramp for ranks ranked gene list plot
if(is.na(upCol)) 
	{upCol = rgb(31, 73, 125, maxColorValue = 255)}
if(is.na(dwnCol)) {dwnCol = rgb(192, 80, 77, maxColorValue = 255)}

rampCols = c(upCol, rgb(0.5, 0.5,0.5), rgb(1,1,1), rgb(0.5, 0.5, 0.5), dwnCol)    # 5 colors scale
#rampCols = c(rgb(0.5, 0.5,0.5))      # all gray color scale
ramp = colorRampPalette(rampCols)

#ploting parameters
height = 2.5; 
width  = 2.5;

#set the header for the plot by simplifying the file name
main = tolower(gsub("_", " ", gsub("\\..*", "", dataFile)))

#set the output file name or output to screen
#X11(height = height, width = width)
pdf(paste0(homeDir, gsub(".xls", ".pdf", dataFile)), height = height, width = width); 

#set plot layout for 3 sections
layout(matrix(1:3, nrow = 3), widths = 1, heights = c(1.5, 0.5, 0.75))

#section 1: enrichment plot
#section 1 layout options
par(mai = c(0, 0.5, 0.25, 0.25), mgp = c(1.2, 0.4, 0), bty = "n", xaxs = "i", yaxs = "i")
#plot axis
plot(NA, xlim = xlim, ylim = ylim, xaxt = "n", xlab = NA, ylab = "Enrichment Score", main = main, cex.main = 1, cex.lab = 0.8, cex.axis = 0.8, font.lab = 2);
lines(c(0, dim(geneList)[1]), c(0,0), lwd = 1, col = rgb(0.5, 0.5, 0.5))
#plot the green enrichment line
lines(c(0, data[, rankCol], dim(geneList)[1]), c(0, data[, esCol], 0), lwd = 3, col = rgb(0, 0.75, 0))

#section 2: tick section showing where all the genes are in your gene set
#section 2 layout options
par(mai = c(0.25, 0.5, 0, 0.25), mgp = c(0.2, 0, 0), bty = "o", xaxs = "i", yaxs = "i", las = 1)
#plot box
plot(NA, xlim = xlim, ylim = c(0, 1), ylab = NA, yaxt = "n", xlab = "Gene Set Rank", cex.lab = 0.8, cex.axis = 0.8, font.lab = 2, xaxt = "n");
#loop through and plot a tick for each gene
for (i in 1:dim(data)[1]) lines(rep(data[i, rankCol], 2), c(0, 1), lwd = 0.2, col = rgb(0, 0, 0))

#section 3: barplot section showing the distribution of rank values
#section 3 layout options
par(mai = c(0.25, 0.5, 0, 0.25), mgp = c(0.2, 1, 0), bty = "o", xaxs = "i", yaxt = "n", las = 1)
#loop through and plot a tick for each gene
barplot(geneList[,2], ylab = NA, xlab = "Preranked Genes", cex.lab = 0.8, cex.axis = 0.8, font.lab = 2, border=ramp(length(geneList[,2])))
med = median(which(geneList[,2]==0))
lines(c(med,med), range(geneList[,2]), lwd = 0.5, col = rgb(0, 0, 0), lty=2)
#add phenotype labels
if(!is.na(dwnPheno)) {text(dim(geneList)[1]-(dim(geneList)[1] * 0.001), geneList[dim(geneList)[1],2] * 0.5, dwnPheno, col = dwnCol)}
if(!is.na(upPheno)) {text(1+(dim(geneList)[1] * 0.2), geneList[1,2] * 0.5, upPheno, col = upCol)}

dev.off()

