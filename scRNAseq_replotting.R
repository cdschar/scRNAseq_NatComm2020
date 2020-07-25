library(data.table)
source("scRNAseq_lib.R")

#set working directory
setwd("/count-IRF4null/")

#loading monocle object from Rdata file and set working directory
load(".RData")

#load magic normalized data and sort to be in the same order as monocle object
magic = data.frame(fread("../MAGIC/wt_magic_out.txt"), row.names = 1)
magic = magic[order(row.names(magic)),]

########
#plot tsne clusters or gene expression overlay
########
#Note:  make sure magic matrix ordered by row.names as above

#Additional parameters for plot_cell_clusters
#colors = c("white", "grey", "dark blue") #supply custom color scheme
#breaks = c(0,1) #set custom expression scales

#plot t-sne plot using phenotype table column includes Cluster, State, Pseudotime etc.
p = plot_cell_clusters(my_cds_subset, color_by="Cluster")
cairo_pdf("Irf4KO.tsne.clusters.pdf", width=5, height=6)
print(p)
dev.off()

#plot t-sne plot using raw gene expression, can use multiple markers at once
gname = "Adam17"
p = plot_cell_clusters(my_cds_subset, markers=c(gname))
cairo_pdf(paste0(gname, ".umi.tsne.pdf"), width=5, height=6)
print(p)
dev.off()

#plot t-sne plot using magic estimated gene expression
#gname = "Klf6"
p = plot_cell_clusters(my_cds_subset, color_by=magic[, gname], breaks = c(0, 2.5))
cairo_pdf(paste0(gname, ".magic.tsne.pdf"), width=5, height=6)
print(p)
dev.off()

########
#plot pseudotime trajectory
########

#plot expression along linear pseudotime trajectory

#some alternate options are listed below
#genes can be one, a short list or a vector
#max_norm = TRUE/FALSE is useful for a list of genes and they will be normalized to each gene's maximal value to ensure they are all on the same scale
#export = "Filename.txt" #put a file name if you want to save the data for other analyses
#lwd controls the line width

p = plot_genes_time(my_cds_subset, genes = as.vector(batf$gene), max_norm=TRUE, lwd=0.5, export = "Filename.txt")
#p = plot_genes_time(my_cds_subset, genes = c("Irf4", "Irf8"), max_norm=TRUE, lwd=0.5)
cairo_pdf(paste0("Irf4.Irf8.norm.pseudotime.pdf"), width=4, height=4)
print(p)
dev.off()

#plot single cell trajectory using phenotype table columns such as Cluster, State, Pseudotime, etc.
p = plot_cell_trajectory(my_cds_subset, color_by="Cluster")
cairo_pdf(paste0("Irf4KO.sct.pdf"), width=4, height=4)
print(p)
dev.off()

#plot single cell trajectory using magic estimated expression
gname = "Irf4"
p = plot_cell_trajectory(my_cds_subset, markers = gname, magic_array=magic[, gname], use_color_gradient=T)
cairo_pdf(paste0(gname, ".magic.expression.sct.pdf"), width=5, height=6)
print(p)
dev.off()

########
#violin plots
########

#reorder magic matrix
magic = magic[as.vector(pData(my_cds_subset)$barcode),]

#plot violin plot based on clusters information
gname = "Xbp1"
id = fData(my_cds_subset)$id[which(fData(my_cds_subset)$gene_short_name == gname)]
p = plot_genes_violin(my_cds_subset[id,], color_by="Cluster", grouping="ClusterAdj")
cairo_pdf(paste0("violin.plots/", gname, ".pdf"), width=5, height=6)
print(p)
dev.off()

#make violin plot using Magic expression
gname = "Lgals1"
df = data.frame(value=magic[,gname], cluster=pData(my_cds_subset)$Cluster, 
                clusterAdj=pData(my_cds_subset)$ClusterAdj)
p = ggplot(df, aes(clusterAdj, value, fill=cluster)) + geom_violin(scale = "width")+
    stat_summary(fun.data=data_summary)+ylab("MAGIC Expression")+xlab("ClusterAdj")+
    monocle_theme_opts()
cairo_pdf(paste0("violin.plots/", gname, ".pdf"), width=5, height=6)
print(p)
dev.off()


######
#plot tSNE using gene signature or collection of genes
######

#batf = read.table("/Volumes/GRAID2/scRNAseq_IRF4/SCENIC/Batf/BATF.targets.info.txt", header = T, sep = "\t")

#colors = c("#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#993404", "#662506") #orange/brown
colors = c("#ffffcc", "#ffeda0", "#ffeda0", "#fed976", "#feb24c", "#fc4e2a", "#bd0026", "#800026") #white, pink, red

scenic = read.table("/Volumes/GRAID2/scRNAseq_IRF4/SCENIC/SCENIC_run/output/Step2_regulonTargetsInfo.tsv", header = T, sep = "\t")
tf = "Batf"
#filter for tf
filter = scenic[which(scenic$TF == tf), ]

#subset magic matrix by gene names
score = magic[, colnames(magic) %in% filter$gene]
score$mean = rowMeans(score)
range((score$mean)^3) #scaling factor

p = plot_cell_clusters(my_cds_subset, color_by=(score$mean)^3, breaks = c(0, 100),  colors = colors)

cairo_pdf(paste0("/Volumes/GRAID2/scRNAseq_IRF4/SCENIC/", tf, ".IRF4ko.score.scaled.red.pdf"), width=5, height=6)
print(p)
dev.off()


######
#adjust cluster names to be in a user defined order
######
pData(my_cds_subset)$ClusterAdj <- plyr::revalue(as.character(pData(my_cds_subset)$Cluster),
                                        c("1" = "5",
                                        "2" = "1",
                                        "3" = "3",
                                        "4" = "4",
                                        "5" = "6",
                                        "6" = "2"))
