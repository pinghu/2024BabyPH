remove(list=ls())
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
filename="combineRash.txt"
outname="Pathway_Rash"

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

A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d<-dim(A)
rownames(A) <- A[,1]
X<-data.frame(A[,2:7])
Y <-data.frame(A[,8:13])
SigPK <-create_pk_mat(X, Y)
#raPH <- data.frame(A[, 15:18])
raRash<- data.frame(A[, 20:25])
#ZraPH <-t(apply(raPH, 1, function(x) (x - mean(x)) / sd(x)))
pRash<-data.frame(A[,68:73])
fdrRash<-data.frame(A[,92:97])
pfcRash_mat<-create_pk_mat(pRash, fdrRash) 
#pPH<- data.frame(A[,42:43])
#fdrPH<-data.frame(A[,58:59])
#pfcPH_mat <-create_pk_mat(pPH, fdrPH)

#fcPH<-data.frame(A[,c(50,51)])
fcRash<-data.frame(A[,80:85])
#################################

#corPH<-data.frame(A[,53:56])
corRash<-data.frame(A[,87:90])

# X<-data.frame(A[,45:48])
# Y<-data.frame(A[,61:64])
# pcorPH_mat <-create_pk_mat(X, Y)

X<-data.frame(A[,75:78])
Y<-data.frame(A[,99:102])
pcorRash_mat <-create_pk_mat(X, Y)

#PHsplit<-A[,65] 

Rashsplit <- A[,125] 



# raRash<-A[,32:37]
# pRash<-A[,39:50]
# fcRash<-A[,52:61]
# fdrRash<-A[,63:74]
# 
# raRashPH <-A[,76:86]
# pRashPH<-A[,88:91]
# fcRashPH <-A[,93:94]
# fdrRashPH <- A[,96:99]

# cal_z_score <- function(x){
#   (x - mean(x)) / sd(x)
# }
# data_subset_norm <- t(apply(B, 1, cal_z_score))
# average_Z<-t(apply(C, 1, cal_z_score))
#data <- as.matrix(raRash)

#data.frame(data)
#normalized_data <- t(apply(data, 1, function(x) (x - mean(x)) / sd(x)))
#data.fram(normalized_data)
library(circlize)
# ha1 = HeatmapAnnotation(
#   group=PHsplit,
#   annotation_name_side = "left"
# )
col_letters = c( "*" = "#FFCCFF", "**" = "#FF99FF", "***"="#FF66FF", "****"="#FF00FF")
ht0 <- Heatmap(Rashsplit, name="Sig")
ht1 <- Heatmap(raRash,name="relative%", column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE,row_names_side = "right",
                     row_names_gp = gpar(fontsize = 6, fontface = "italic"), cluster_rows = FALSE,
               row_title = NULL,show_row_names = FALSE
               )
ht2 <- Heatmap(SigPK,name="pkruskal-wallis", column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE,row_names_side = "right",
               row_names_gp = gpar(fontsize = 6, fontface = "italic"), cluster_rows = FALSE,
               row_title = NULL,show_row_names = FALSE, 
               col = col_letters,
               )
# jpeg("ht1.jpg", width = 2800, height = 2800, res=300 )
# draw(ht1, ,heatmap_legend_side = "left")
# dev.off()



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


ht_list = ht0+ht1 +ht2 +ht3+ht4
jpeg(paste0(outname, ".jpg"), width = 2800, height = 2300, res=300 )
draw(ht_list, ,heatmap_legend_side = "left", row_split = Rashsplit)
dev.off()




# 
# 
# jpeg(paste0(outname, ".ht1.png"), width = 1200, height = 2800,  res = 300)
# draw(ht1, heatmap_legend_side = "left")
# dev.off()
# 
# row_ha <- rowAnnotation(LFC = anno_barplot(row_avg))
# column_ha <- HeatmapAnnotation(LFC = anno_barplot(col_avg))
# 
# 
# 
# rownames(data)=NULL
# ht_list <- Heatmap(data, name = "Log2FoldChange", top_annotation = column_ha,
#                    left_annotation = row_ha,
#                    column_names_gp = gpar(fontsize = 9),
#                    row_split = bacteria,
#                    row_names_gp = gpar(fontsize = 6, fontface = "italic"),  # Reduce row names font size
#                    column_names_gp = gpar(fontsize = 6),  # Reduce column names font size
#                    heatmap_legend_param = list(
#                      title_gp = gpar(fontsize = 6),  # Reduce legend title font size
#                      labels_gp = gpar(fontsize = 6)  # Reduce legend labels font size
#                    )
#                    )
# 



