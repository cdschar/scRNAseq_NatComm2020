library(cellrangerRkit)
library(monocle)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tibble)
library(viridis)

source("scRNAseq_lib.R")

#set cellranger output directory for data import
dirs = c("/count-IRF4null/")

#set perplexity for t-SNE plots
perplexity = 30

#provide certain mark genes for plots
marker_list = c("Cd19", "Sdc1", "Irf4", "Irf8", "Tbx21")

#choose to save the RData image in folder or not
save_image = TRUE

# load the outs directory of cellranger
gbm <- load_cellranger_matrix(dir)
setwd(dir)
dir.create("DEG")

# rename gene symbol column and create cell data set
my_feat <- fData(gbm)
names(my_feat) <- c('id', 'gene_short_name')
my_cds <- newCellDataSet(exprs(gbm), phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = my_feat), 
                         lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# normalization and variance estimation
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

# gene detection QC plot
my_cds <- detectGenes(my_cds, min_expr = 0.1)
x <- pData(my_cds)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
df <- data.frame(x = x_1)
cairo_pdf("qc_num_genes_expressed.pdf", width=4, height=4)
p=ggplot(df, aes(x)) + geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
print(p)
dev.off()
expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))

# UMI Read count QC plot
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
cairo_pdf("qc_UMI_num_genes_expressed.pdf", width=4, height=4)
p=ggplot(pData(my_cds), aes(num_genes_expressed, UMI))+geom_point()+xlim(1000,8000)+ylim(0,125000)
print(p)
dev.off()

######################################################
# unsupervised clustering using PCA reduced data
disp_table <- dispersionTable(my_cds)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)
cairo_pdf("qc_expressed_genes.pdf", width=4, height=4)
p=plot_ordering_genes(my_cds)
print(p)
dev.off()
cairo_pdf("pca_variance.pdf", width=4, height=4)
p=plot_pc_variance_explained(my_cds, return_all = FALSE)
print(p)
dev.off()
my_cds <- reduceDimension(my_cds, max_components = 2, num_dim = 10,
                          reduction_method = 'tSNE', verbose = TRUE)
my_cds <- clusterCells(my_cds, num_clusters = 10)
my_cluster_dim_10 <- pData(my_cds)$Cluster

cairo_pdf("clusters_10.pdf", width=4, height=4)
p=plot_cell_clusters(my_cds)
print(p)
dev.off()

for(gene in marker_list){
  cairo_pdf( paste0("clusters_10_marker_", gene, ".pdf"), width=4, height=5)
  p=plot_cell_clusters(my_cds, marker=c(gene))
  print(p)
  dev.off()
}

######################################################
# DEG based on clusters
de_cluster_one <- differentialGeneTest(my_cds[unsup_clustering_genes$gene_id,],
                                       fullModelFormulaStr = '~Cluster',
                                       cores = 8)
# highlight top 15 DEG
deg <- de_cluster_one[order(de_cluster_one$qval),]
for (i in 1:15){
  gene_id <- deg$id[i]
  gene_name <- deg$gene_short_name[i]
  cairo_pdf(paste0("./DEG/clusters_10_DEG_", gene_name, ".pdf"), width=6, height=3)
  p=plot_genes_violin(my_cds[gene_id,], grouping = "Cluster", color_by = "Cluster")
  print(p)
  dev.off()
}

for (gene_name in marker_list){
  gene_id = fData(my_cds)$id[which(fData(my_cds)$gene_short_name==gene_name)]
  cairo_pdf(paste0("./DEG/clusters_10_Genes_", gene_name, ".pdf"), width=4, height=3)
  p=plot_genes_violin(my_cds[gene_id,], grouping = "Cluster", color_by = "Cluster")
  print(p)
  dev.off()
}

######################################################
#using dpFeature to select subset of genes for single cell trajectory
expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))
if(!"ENSMUSG00000020592" %in% expressed_genes){
  expressed_genes <- c(expressed_genes, "ENSMUSG00000020592")
}
my_cds_subset <- my_cds[expressed_genes, ]
fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)

my_cds_subset <- reduceDimension(my_cds_subset,
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 10,
                                 reduction_method = 'tSNE',
                                 verbose = TRUE,
                                 perplexity = perplexity)
my_cds_subset <- clusterCells(my_cds_subset, verbose = FALSE)
cairo_pdf("dpFeature_rho_delta.pdf", width=4, height=4)
p=plot_rho_delta(my_cds_subset, rho_threshold = 2, delta_threshold = 10)
print(p)
dev.off()

#make dpFeature cluster plots and marker plots
my_cds_subset <- clusterCells(my_cds_subset,
                              rho_threshold = 5 ,
                              delta_threshold = 15,
                              skip_rho_sigma = T,
                              verbose = FALSE)
cairo_pdf("dpFeature_clusters.pdf", width=4, height=4)
p=plot_cell_clusters(my_cds_subset, color_by = "Cluster")
print(p)
dev.off()

for(gene in marker_list){
  cairo_pdf( paste0("dpFeature_clusters_marker.", gene, ".pdf"), width=4, height=5)
  p=plot_cell_clusters(my_cds_subset, marker=c(gene))
  print(p)
  dev.off()
}

#using top 1000 DEG for SCT
clustering_DEG_genes <- differentialGeneTest(my_cds_subset, fullModelFormulaStr = '~Cluster', cores = 8)
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')

my_cds_subset <- orderCells(my_cds_subset, reverse = T)
cairo_pdf("sct_pseudotime.pdf", width=4, height=4)
p=plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime")
print(p)
dev.off()
cairo_pdf("sct_clusters.pdf", width=4, height=4)
p=plot_cell_trajectory(my_cds_subset, color_by = "Cluster")
print(p)
dev.off()

#plot DEG expression based on dpFeature clusters and marker list
deg <- clustering_DEG_genes[order(clustering_DEG_genes$qval),]
for (i in 1:15){
  gene_id <- deg$id[i]
  gene_name <- deg$gene_short_name[i]
  cairo_pdf(paste0("./DEG/dpFeature_DEG_", gene_name, ".pdf"), width=4, height=3)
  p=plot_genes_jitter(my_cds_subset[gene_id,], grouping = "Cluster")
  print(p)
  dev.off()
}
gene_names<-marker_list
gene_names<-gene_names[gene_names %in% fData(my_cds_subset)$gene_short_name]
for (gene_name in gene_names){
  gene_id = fData(my_cds_subset)$id[which(fData(my_cds_subset)$gene_short_name==gene_name)]
  cairo_pdf(paste0("./DEG/dpFeature_Genes_", gene_name, ".pdf"), width=4, height=3)
  p=plot_genes_jitter(my_cds_subset[gene_id,], grouping = "Cluster")
  print(p)
  dev.off()
}

#plot expression along pseudotime for selected genes
my_genes <- row.names(subset(fData(my_cds_subset), gene_short_name %in% marker_list))
cds_plotted <- my_cds_subset[my_genes,]
cairo_pdf("sct_expr_markers.pdf", width=6, height=6)
p=plot_genes_in_pseudotime(cds_plotted, color_by="Cluster")
print(p)
dev.off()

cairo_pdf("sct_combined_markers.pdf", width=4, height=4)
p=plot_genes_time(my_cds_subset, marker_list)
print(p)
dev.off()

if(save_image) {save.image(".RData")}
