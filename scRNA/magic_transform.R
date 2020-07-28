library(Rmagic)
library(ggplot2)
library(readr)
library(viridis)

#read in raw UMI matrix and select expression values
data <- read_csv("wt_expr.csv")
cells_wt <- data[,1]
data <- data[,2:ncol(data)]
keep_cols <- colSums(data > 0) > 10
wt_expr<-data[,keep_cols]

#filtering of detected cells and genes
keep_rows <- rowSums(wt_expr) > 1000
wt_expr <- wt_expr[keep_rows,]

#data normalization
wt_expr <- library.size.normalize(wt_expr)
wt_expr <- sqrt(wt_expr)

#magic transformation of all genes
wt_magic <- magic(wt_expr, genes="all_genes")
wt_magic <- wt_magic$result
row.names(wt_magic) <- cells_wt$cell
wt_magic <- wt_magic[sort(row.names(wt_magic)),]

#various plotting functions
ggplot(wt_expr)+geom_point(aes(Sdc1, Ptprc))+scale_colour_viridis(option="B")
  
ggplot(wt_magic) + geom_point(aes(Sdc1, Ptprc)) + scale_colour_viridis(option="B")+xlim(0,2.5)+ylim(0,10)
  
wt_magic <- magic(wt_expr, genes=c("Sdc1", "Irf4", "Uggt1"), t=4, init=wt_magic)

wt_magic_pca <- magic(wt_expr, genes="pca_only")
ggplot(wt_magic_pca) +
  geom_point(aes(x=PC1, y=PC2, color=wt_magic$result$Uggt1)) +
  scale_color_viridis(option="B") +
  labs(color="Uggt1")

#the same process ofr IRF4 KO data set
data<-read_csv("./IRF4null_expr.csv")
cells_ko <- data[,1]
data <- data[,2:ncol(data)]
keep_cols <- colSums(data > 0) > 10
ko_expr<-data[,keep_cols]

keep_rows <- rowSums(ko_expr) > 1000
ko_expr <- ko_expr[keep_rows,]
ko_expr <- library.size.normalize(ko_expr)
ko_expr <- sqrt(ko_expr)

ko_magic <- magic(ko_expr, genes="all_genes")
ko_magic <- ko_magic$result
row.names(ko_magic) <- cells_ko$cell
ko_magic <- ko_magic[sort(row.names(ko_magic)),]

ggplot(ko_magic) + geom_point(aes(Sdc1, Ptprc)) + scale_colour_viridis(option="B")+xlim(-0.05,0.25)+ylim(0,10)
