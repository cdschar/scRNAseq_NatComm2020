library(preprocessCore)
library(FNN)

#read in magic normalized expression matrix
cds<-read.table("../MAGIC/wt_expr.csv", sep=",", header=T, row.names=1)
cds<-t(cds)
cds<-as.data.frame(cds)
cds$SYMBOL=row.names(cds)

diff<-read.table("../KNN_prediction/diff.significant.glm.Bcell.GSK343.txt", sep="\t", header=T)
#sig_names<-names(diff)[grep(".sig", names(diff))]
#sig_names<-sig_names[!grepl("KO", sig_names)]
#sig<-apply(diff[,sig_names], 1, any)
#diff<-diff[sig,]

combined<-merge(cds, diff)
comb_expr<-combined[, grepl("rpkm|-2|-1", names(combined))]
comb_expr<-comb_expr[, !grepl("GSK", names(comb_expr))]
#comb_expr<-comb_expr[, !grepl("2011|2012|2014|2037", names(comb_expr))]

comb_norm<-normalize.quantiles(as.matrix(comb_expr))
testing<-comb_norm[,1:8368]
training<-comb_norm[,8369:8377]
cl<-c("Naive","Naive","Naive","actB","Pb","actB","Pb","actB","Pb")
#cl<-rep(c("D0","D1","D3","D5","D8neg","D8pos"), 2)
#cl<-c("Nav_IgM","Nav_IgM","Nav_IgM","Mem_IgM","Mem_IgM","Mem_IgG","Mem_IgM","Mem_IgG")

fit<-knn(t(training), t(testing), cl=cl, k=2)
class_Scharer<-as.vector(fit)
pData(my_cds_subset)$class_Scharer<-class_Scharer
plot_cell_clusters(my_cds_subset, color_by="class_Scharer")
