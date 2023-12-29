library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(NMF)
library(viridis)

scRNA_list<-readRDS("data/scRNA_list.rds")
res_list<-list(c(1:5))
fs_list<-list(c(1:5))
DL<-list(c(1:5))
rank<-c(2,2,5,4,2)
for(id in c(1:5)){
  library(readr)
  cnv_scores <- read_csv(paste0("data/ref_new/cnv_scores_",id,".csv"))
  #cnv_scores$barcode<-Replace(cnv_scores$barcode,from="_1",to="_1")
  library(NMF)
  sce<-scRNA_list[[id]]
  c<-colnames(sce)%in%cnv_scores$barcode
  rownames(cnv_scores)<-cnv_scores$barcode
  b<-colnames(sce)
  D<-data.frame(barcode=b,cancer=c(1:length(b)))
  for(i in 1:length(c)){
    if(c[i]==TRUE){
      D[i,"cancer"]<-cnv_scores[b[i],"cancer"]
    }else{
      D[i,"cancer"]<-0.005
    }
  }
  rownames(D)<-b
  D$barcode<-NULL
  #redefine tumor cells
  sce<-AddMetaData(sce,D,col="cancer")
  sce@meta.data$cancer[sce@meta.data$cancer>=0.3]<-"tumor"
  sce@meta.data$cancer[sce@meta.data$cancer<0.3]<-"normal"
  table(sce@meta.data$cell,sce@meta.data$cancer)
  sce@meta.data$cell[sce@meta.data$cancer=="tumor"&sce@meta.data$cell=="Mesothelial cells"]<-"Tumour cells"
  sce@meta.data$cell[sce@meta.data$cancer=="tumor"&sce@meta.data$cell=="Endothelial cells"]<-"Tumour cells"
  table(sce@meta.data$cell,sce@meta.data$cancer)
  Idents(sce)<-"cell"
  sub_sce<-subset(sce,ident="Tumour cells" )
  sub_sce=CreateSeuratObject(
    counts = sub_sce@assays$RNA@counts,
    meta.data = sub_sce@meta.data
  )
  sub_sce = NormalizeData(sub_sce) %>% FindVariableFeatures() %>% ScaleData(do.center = F)
  #perfom NMF analysis
  vm <- sub_sce@assays$RNA@scale.data
  res <- nmf(vm,rank=rank[id],method = "snmf/r") 
  coefmap(object=res, Colv = TRUE,annRow = NA, annCol = NA,file=paste0("plot/",id,"p0.pdf"))
  basismap(object=res, Colv = TRUE,annRow = NA, annCol = NA,file=paste0("plot/",id,"p2.pdf"))
  consensusmap(object=res, annRow = NA,Colv = TRUE,
               annCol = NA,color = "YlOrRd:50",
               main = "Consensus matrix", info = FALSE,file=paste0("plot/",id,"p1.pdf"))
  #extract the top 50 features
  fs <- extractFeatures(res,50L)
  fs <- lapply(fs,function(x)rownames(res)[x])
  write.csv(fs,file=paste0("data/",id,"_b_fs.csv"))
  res_list[[id]]<-res
  #add the score of program to each sample
  Subset<-sub_sce
  data<- as.matrix(Subset@assays$RNA@scale.data)
  rk=rank[id]
  s<-fs
  data1 <- c()
  for (i in 1:as.numeric(rk)) {
    data1 <- rbind(data1, as.matrix(Subset@assays$RNA@scale.data)[s[[i]],])
    Subset <- AddModuleScore(object = Subset, features = list(s[[i]]), name = paste("Pg", i, sep = "_"))
  }
  Pg_mat <- Subset@meta.data[ ,paste('Pg',"_",seq(1:rk), rep('1',rk), sep="")]
  SD <- apply(Pg_mat, 2, sd)
  data2 <- c()
  for (i in which(SD>0.2)) {
    data2 <- rbind(data2, as.matrix(Subset@assays$RNA@scale.data)[s[[i]],])
  }
  data2[data2>3] <- 3
  data2[data2< -3] <- -3
  #heatmap of programs in each sample
  x <- pheatmap(data2, cluster_cols=T, cluster_rows=F, scale="none",color=colorRampPalette(rev(brewer.pal(n = 7, name ="PdBu")))(100),#pals::kovesi.diverging_bwr_40_95_c42(100),
                fontsize=4, #treeheight_row=0, treeheight_col=0, cellheight = 0.8,cellwidth = 0.1,
                show_rownames=T, show_colnames=F, clustering_method = "ward.D",
                border_color=NA, width=10,height=10,
                filename = paste0("plot/",id,"_b_Target_tumor.pdf"))
  
  #write.table(data2, file=paste0("plot/",id,"_top50_gene.xls"), quote=F, sep="\t")
  DL[[id]]<-data2
}
#### Meta-programme
#determine similarity between different Programs
jac <- function(x, y) {
  inter <- intersect(x, y)
  total <- union(x, y)
  similarity <- length(inter)/length(total)
  return(similarity)
}

#Pgs <- read.table("/home/jovyan/farm/tm/RCC_final/RCC/Program/program_geneset_all.xls", header=T)
library(readr)
Pgs<-read.csv("data/fs_new.csv")
Pgs$gene<-NULL
Mat <- matrix(0, ncol = ncol(Pgs), nrow = ncol(Pgs))
rownames(Mat) <- colnames(Pgs)
colnames(Mat) <- colnames(Pgs)

for (i in 1:ncol(Pgs)) {
  for (j in 1:ncol(Pgs)) {
    ss <- jac(as.vector(Pgs[,i]), as.vector(Pgs[,j]))
    Mat[i,j] <- ss*100
  }
}

#custom_magma <- c(colorRampPalette(c("yellow", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
custom_magma<-colorRampPalette(c("yellow","purple"))(n=299)
#plot the heatmap of similarity between programs
Mat[Mat>25] <-25
x <- pheatmap(as.matrix(Mat), cluster_cols=T, cluster_rows=T, 
              clustering_distance_rows="euclidean", color=colorRampPalette(c("#f0eedf","#fdbb2b","#E4336F","#000F5D"))(1000),
              fontsize=2,treeheight_row=0,treeheight_col=30, 
              cellheight = 7,cellwidth = 7,show_rownames=T, 
              show_colnames=F,clustering_method = "ward.D2", border_color = "NA",width=20,height=20,
              filename = paste0("plot/","1_Target_tumor.pdf"))

#identify the marker genes for each meta program
pg<-list(c(1:3))
pg[[1]]<-c("X2.2","X1.2","X4.3")
pg[[2]]<-c("X2.1","X1.1","X5.2")
pg[[3]]<-c("X3.3","X5.1")
#pg[[4]]<-c("X3.1","X4.2")
name<-c("Pg1","Pg2","Pg3")
t<-c(1,1,1)
sig<-list()
for(k in 1:3){
  l<-pg[[k]]
  d<-c()
  for(i in l){
    d<-append(d,Pgs[,i])
  }
  co<-c()
  for(j in d){
    if(sum(d==j)>t[k]){
      co<-append(co,j)
    }
  }
  co<-unique(co)
  sig[[name[k]]]<-co
  write.csv(co,file=paste0("data/",k,"_b_co.csv"))
}

#check the overlap with the markers of each tumor cluster
name<-c("Tumor0","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
Gene<-list(1:6)
Portion<-list(1:6)
for(i in 1:6){
  gene<-c(1:6)
  portion<-c(1:6)
  for(j in 1:6){
    df<-rownames(sub_sce.markers)[sub_sce.markers$cluster==name[j]]
    gene_df<-df[df%in%as.character(sig[[i]])]
    if(length(gene_df)!=0){
      gene[j]<-gene_df
      portion[j]<-length(df)/length(as.character(sig[[i]]))
    }else{
      gene[j]<-"Nothing"
      portion[j]<-0
    }
  }
  Gene[[i]]<-gene
  Portion[[i]]<-portion
}
for(i in 1:6){
  write.csv(sig[[i]],file=paste0("data/Pgs_sig_",i,".csv"))
}

