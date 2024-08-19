###################fix error so that even when 1 row, it can map p_mat#####
remove(list=ls())
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
selection <-args[3]
filename="ButtocksRash.txt"
outname="GO_ButtocksRash"
#filename="PerianalRash.txt"
#outname="GO_PerianalRash"
selection ="CC"
if(selection == "MF"){
  categoryS="molecular_function"
  outname=paste0(outname,".MF")
}else if (selection == "CC"){
  categoryS="cellular_component"
  outname=paste0(outname,".CC")
} else{
  outname=paste0(outname,".BP")
  categoryS="biological_process"
}
create_pk_mat <- function(pk, adjpk) {
  nrow_pk <- nrow(pk)
  ncol_pk <- ncol(pk)
  
  pk_mat <- matrix("", nrow = nrow_pk, ncol = ncol_pk)
  
  # Set the row and column names
  rownames(pk_mat) <- rownames(pk)
  colnames(pk_mat) <- colnames(pk)
  
  pk_mat[!is.na(pk) & pk <= 0.1] <- "*"
  pk_mat[!is.na(pk) & pk <= 0.05] <- "**"
  pk_mat[!is.na(adjpk) & adjpk <= 0.1] <- "***"
  pk_mat[!is.na(adjpk) & adjpk <= 0.05] <- "****"
  
  return(pk_mat)
}

A1<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)

rownames(A1) <-paste( A1[,1], A1[,2])

d<-dim(A1)
rownames(A1) <-paste( A1[,1], A1[,2])
A<-A1[A1$GOCategory==categoryS, ]
d<-dim(A)
print(d)
Hpixel=800+26*d[1]

create_selected_matrix <- function(selected_columns, pmat) {
  selected_mat <- matrix("", nrow = nrow(pmat), ncol = length(selected_columns))
  colnames(selected_mat) <- colnames(pmat)[selected_columns]
  
  for (i in 1:length(selected_columns)) {
    selected_mat[, i] <- pmat[, selected_columns[i]]
  }
  
  return(selected_mat)
}

X <- data.frame(A[, c(16:29, 71:83)])
Y <- data.frame(A[, c(43:56, 97:109)])
pmat <- create_pk_mat(X, Y)

dim(pmat)

selected_columns <- c(1, 15, 5, 19, 11, 25)
SigPK <- create_selected_matrix(selected_columns, pmat)

selected_columns <- c(6, 7, 8, 20, 21, 22)
pfcRash_mat <- create_selected_matrix(selected_columns, pmat)

selected_columns <- c(9, 10, 24, 25)
pcorRash_mat <- create_selected_matrix(selected_columns, pmat)

raRash <- data.frame(A[, c(8,9,10,62,63,64)])
ZraRash <-t(apply(raRash, 1, function(x) (x - mean(x)) / sd(x)))

fcRash<-data.frame(A[,c(33,34,35,88,89,90)])
corRash<-data.frame(A[,c(36,37,91,92)])


#GOsplit<-A[,3] 
  

library(circlize)
ha1 = HeatmapAnnotation(
  group=GOsplit,
  annotation_name_side = "left"
)
col_letters = c( "*" = "#FFCCFF", "**" = "#FF99FF", "***"="#FF66FF", "****"="#FF00FF")
#ht0 <- Heatmap(GOsplit, name="Sig")
ht1 <- Heatmap(ZraRash,name="relative% Zscore", column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE,row_names_side = "right",
                     row_names_gp = gpar(fontsize = 6, fontface = "italic"), cluster_rows = FALSE,
               row_title = NULL,show_row_names = FALSE
               )
ht2 <- Heatmap(SigPK,name="pkruskal-wallis", column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE,row_names_side = "right",
               row_names_gp = gpar(fontsize = 6, fontface = "italic"), cluster_rows = FALSE,
               row_title = NULL,show_row_names = FALSE, 
               col = col_letters,
               )




ht3 <- Heatmap(fcRash, name="Fold Change", column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE,row_names_side = "right",
               row_names_gp = gpar(fontsize = 6, fontface = "italic"), cluster_rows = FALSE, 
               col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, w, h, col) {
                 grid.text(pfcRash_mat[i, j], x, y)
               },
               row_title = NULL,show_row_names = FALSE)
ht4 <- Heatmap(corRash,name="Correlation", column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE, 
               row_names_side = "right",
               row_names_gp = gpar(fontsize = 6, fontface = "italic"), cluster_rows = FALSE, 
               cell_fun = function(j, i, x, y, w, h, col) {
                 grid.text(pcorRash_mat[i, j], x, y)
               },
               col = colorRamp2(c(-0.3, 0, 0.3), c("blue", "white", "red")))


ht_list = ht1 +ht2 +ht3+ht4
jpeg(paste0(outname, ".jpg"), width = 2800, height = Hpixel, res=300 )
draw(ht_list, heatmap_legend_side = "left", #row_split = GOsplit, 
     column_title = paste( filename, "\n", categoryS))
dev.off()
