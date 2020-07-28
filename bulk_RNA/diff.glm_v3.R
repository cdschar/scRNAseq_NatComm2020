library("edgeR");
library("limma")

homeDir = "/BcKO_RNA/analysis/diff/"
setwd(homeDir);

#set significance threshold
sigFDR = 0.05
sigFC = 1

experimentName = "BcKOvWt.Div.RNAseq"

#read in manifest file
fqDir = "/BcKO_RNA/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste(fqDir, fqFile, sep = ""), sep = "\t", header = TRUE, as.is = TRUE);
files = files[files$include, ]

#read in raw and rpm normalized files
covDir = "/BcKO_RNA/coverage/";
cts = read.table(paste0(covDir, "BcKO.Div.RNAseq.geneCts.detected.RPM.3.exon.csv"), sep = ",", header = T)
rpkm = read.table(paste0(covDir, "BcKO.Div.RNAseq.geneRpkm.detected.RPM.3.exon.csv"), sep = ",", header = T);

#establish groups for comparisons
grps= unique(files$group)
grpPairs = combn(grps, 2)
#filter for specific comparisons
grpPairs= grpPairs[,c(1:6, 16,25,33,35,40, 45:49)]

##########################

#establish edgeR DGE object
dge = DGEList(as.matrix(cts[, grepl(paste(files$sample, collapse = "|"), colnames(cts))]), group = files$group[match(colnames(cts[, grepl(paste(files$sample, collapse = "|"), colnames(cts))]), files$sample)])

#calculate normalization factor
dge = calcNormFactors(dge)

#plotMeanVar(dge)
cairo_pdf(file = paste0(experimentName, ".MDSplot.pdf"), height = 5, width = 5);
par(mai = c(1.5, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
plotMDS.DGEList(dge, main = "MDS Plot for Count Data", labels = gsub("Sample\\_[0-9]+\\_", "", colnames(dge$counts)), cex = 0.3)
dev.off();

#######################################################################
## GLM Function

#set up differential matrix to append data to
diff = rpkm

#create design matrix and estimate dispersion
design = model.matrix(~0 + files$group + files$mouse)

#change design headers to simple strings
#this may have to be customized
colnames(design) = gsub("files\\$group", "", colnames(design))
colnames(design) = gsub("files\\$", "", colnames(design))

#estimate dispersion
dge = estimateDisp(dge, design)

#call glmfit function
fit = glmFit(dge, design)
rownames(design) = colnames(dge)

#loop through each comparison and append stats to diff matrix

for (i in 1:dim(grpPairs)[2]) {
    
    #set comparisons
    comp = paste0(grpPairs[2, i], ".v.", grpPairs[1, i]);
    print(paste(comp,  Sys.time()))
    
    #establish contrast matrix
    cont = makeContrasts(contrasts = paste0(grpPairs[2, i], "-", grpPairs[1, i]), levels = design)
    
    #create matrix for this comparisons and use glmLRT function for differential analysis
    lrt = glmLRT(fit, contrast = cont)
    
    #FDR correcting Pvalues
    lrt$table$fdr = p.adjust(lrt$table$PValue, method = "fdr")
    
    #mark genes that are significant by FDR criteria
    #lrt$table$sig = lrt$table$fdr < sigFDR
    lrt$table$sig = (lrt$table$fdr < sigFDR & abs(lrt$table$logFC) >= sigFC)
 
    #add header to the lrt$table slot that describes the comparison performed
    colnames(lrt$table) = paste0(comp, ".", names(lrt$table))
    
    #cbind diff test onto the diff matrix
    diff[names(lrt$table)] = NA
    diff[rownames(lrt$table), names(lrt$table)] = lrt$table
}

#add column that tells if a gene is significant in any comparison
diff$sigAny <- apply(diff[,(grep(".sig" == "TRUE", diff))], 1, any)

#remove LR columns and write to file
diff = diff[,-grep(".LR",names(diff))]
write.table(diff, file = paste0(homeDir, experimentName, ".diff.glm.Bcko_RNA.txt"), sep = "\t", quote = F, row.names = F)

#filter only the significant genes
sig.diff = diff[which(diff$sigAny == "TRUE"), ]
#sig.diff= sig.diff[complete.cases(sig.diff),]
write.table(sig.diff, file = paste0(experimentName, ".diff", ".significant.glm.Bcko_RNA.txt"), sep = "\t",, quote = F, row.names = F)

#################################################################
#write table that tells how many genes are significant for each comparison
cols = colnames(diff)[grepl(".sig", colnames(diff))]
stats = data.frame(Comp = gsub(".sig", "", cols))

for(i in 1:dim(stats)[1]){
  
  print(stats$Comp[i])
  print(table(diff[, cols[i]]))
  
  stats$notSig[i] = as.numeric(table(diff[, cols[i]])[1])
  stats$Sig[i] = as.numeric(table(diff[, cols[i]])[2])
  
}

write.table(stats, file = paste0("Sig.stats.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F);

