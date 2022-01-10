library(FactoMineR)
library(factoextra)
library(ggplot2)

setwd("E:/pca_csv")# 
x <- read.csv("all.csv")####
res.pca=PCA(x,scale.unit = T,ncp=2,graph = F)
PCA(x)
var=get_pca_var(res.pca)
var$coord
fviz_pca_var(res.pca,col.var="black")
fviz_pca_var(res.pca,col.var = "contrib",
             repel = TRUE, 
             alpha.var = 10,
          gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
         # title = "PCA",
             geom = c("arrow","text"))


