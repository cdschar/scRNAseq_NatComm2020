library(SCENIC)
library(Biobase)
library(AUCell)

#set the colors based on t-SNE clusters
colVars <- list(CellType=setNames(c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"), 
    c(1,2,3,4,5,6,7,8)))
saveRDS(colVars, file="int/colVars.Rds")

#read in MAGIC transformed data matrix
exprMat <- read.table(file = "/Volumes/GRAID2/scRNAseq_IRF4/MAGIC/wt_expr.csv", header=T, sep = ",", row.names=1)
exprMat <- t(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=10,
                           minSamples=9)

#set scenic run parameters
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 8
scenicOptions@settings$seed <- 123

#gene filtering and data normalization
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered, scenicOptions)

#three basic steps of SCENIC run
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat, skipTsne = TRUE, skipHeatmap = TRUE)

#function to map the cluster colors
map_cluster <- function(x) { if(x==1) "#F8766D" else if(x==2) "#CD9600" else if(x==3) "#7CAE00" else if(x==4) "#00BE67" 
                  else if(x==5) "#00BFC4" else if(x==6) "#00A9FF" else if(x==7) "#C77CFF" else if(x==8) "#FF61CC"}
                  
map_state <- function(x) { if(x==1) "#F8766D" else if(x==2) "#C49A00" else if(x==3) "#53B400" else if(x==4) "#00C094" 
                  else if(x==5) "#00B6Eb" else if(x==6) "#A58AFF" else if(x==7) "#FB61D7" }

Colcols <- unlist(lapply(pdata$Cluster, map_cluster))

#make a heatmap of the results scenic scores of different TFs
heatmap3(as.matrix(data1), outFile = heatFile, rowClust = T, colClust = T, col = ramp(100), annotCol = Colcols, 
heatDim = c(3, 3), annotWid = 0.25, clusterWid = 0.4, key = T, keyWid = 0.4, flatten = T, res = 600, labDim = c(0.7, 0.4), height = 3, width = 3, 
outMai = c(0.2, 0.1, 0.1, 0.1), flatColClust = F, clustMeth = "ward.D", scale = c(-1.5, 1.5), rowLab = row.names(mat), labCex = 0.05)

#AUC heatmaps and density plots
auc<-readRDS("~/Desktop/SCENIC_test/int/3.4_regulonAUC.Rds")
auc<-getAUC(auc)
auc_cutoff<-read.table("~/Desktop/SCENIC_test/int/3.5_AUCellThresholds_Info.tsv", header=T, sep="\t")
wt_aucell<-as.data.frame(t(auc))
x1<-strsplit(colnames(wt_aucell), split="/.")
tf_names<-sapply(x1, function(x) x[1])
colnames(wt_aucell)<-tf_names

on_off_mat = NA
for(i in 1:dim(auc)[1]){
  tf_name = tf_names[i]

  #p<-ggplot(wt_aucell, aes_string(tf_name))+geom_density()
  #ggsave(paste0("./plots/",tf_name,"_density.png"), p, width=5, height=5)
  
  cutoff = auc_cutoff$threshold[i]
  #p<-plot_cell_clusters(my_cds_subset, color_by=as.numeric(wt_aucell[,tf_name]>cutoff))
  on_off_mat = cbind(on_off_mat, as.numeric(wt_aucell[,tf_name]>cutoff))
  
  #ggsave(paste0("./plots/",tf_name,"_AUCell.png"), p, width=5, height=5)
}

#write aucell scores for TFs
pdata = read.table(file = "/Volumes/GRAID2/scRNAseq_IRF4/count-Wt/Phenotype.data.txt", header = T, sep = "\t")
wt_aucell = read.table(file = "../aucell_matrix.txt", header = T, sep = "\t", row.names = 1, check.names=F)
wt_aucell = t(wt_aucell)
x1<-strsplit(colnames(wt_aucell), split=" ", fixed=TRUE)
tf_names<-sapply(x1, function(x) x[1])
expr<-read.table("../../MAGIC/wt_magic_out.txt", header=T, row.names=1)
wt_aucell = wt_aucell[row.names(expr),]
x1<-strsplit(tf_names, split="_", fixed=TRUE)
tf_genes<-sapply(x1, function(x) x[1])

#AUC differential test between clusters
comb = combn(1:8, 2)
for(j in 1:28){
  p_list = c()
  fc_list = c()
  expr_fc = c()
  first = comb[1, j]
  second = comb[2, j]
  for(i in 1:dim(wt_aucell)[2]){
    set1 = wt_aucell[pdata$Cluster == first, i]
    set2 = wt_aucell[pdata$Cluster == second, i]
    fit = wilcox.test(set1, set2)
    p_list = c(p_list, fit$p.value)
    fc = log2(mean(set1)/mean(set2))
    fc_list = c(fc_list, fc)
    expr1 = expr[pdata$Cluster == first, tf_genes[i]]
    expr2 = expr[pdata$Cluster == second, tf_genes[i]]
    fc2 = log2(mean(expr1)/mean(expr2))
    expr_fc = c(expr_fc, fc2)
  }
  out = data.frame(TF=tf_names, p_value=p_list, logFC=fc_list, expr_FC=expr_fc)
  out = out[order(out$p_value),]
  fname = paste0("Wilcox_test_out_Cluster_", first, ".v.", second, ".txt")
  write.table(out, fname, quote=F, sep="\t", row.names=F)
}



