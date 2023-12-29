####add_cell_score####
library(readr)
library(Matrix)
#st_list <- readRDS("~/spatial/data/st_list.rds")
cell_Sig<- read_csv("~/spatial/data/cells.csv")
a<-cell_Sig$...1
cell_Sig$...1<-NULL
rownames(cell_Sig)<-a
cell<-colnames(cell_Sig)
exp<-st@assays$Spatial@scale.data
marker<-c()
score<-data.frame(colnames(exp))
for(i in 1:length(cell)){
  cell_vector<-cell_Sig[,cell[i]]
  for(j in 1:283){
    if(cell_vector[j,1]==1){
      marker<-append(marker,a[j])
    }
  }
  marker<-marker[marker%in%rownames(exp)]
  cell_exp<-exp[marker,]
  cell_score<-colMeans(cell_exp)
  score[i]<-cell_score
  marker<-c()
}
rownames(score)<-colnames(exp)
colnames(score)<-cell
Score_df<-data.frame(barcode=rownames(score))
for(i in colnames(score)){
  Score<-score[,i]
  normalized<- (Score - min(Score)) / (max(Score) - min(Score))
  Score_df[,i]<-normalized
}

a<-Score_df$barcode
rownames(Score_df)<-a
Score_df$barcode<-NULL
st<-AddMetaData(st,Score_df,col=colnames(Score_df))
point_size<-c(1.5,1.9,2,1.5,1.9,2,1.5,1.4)
image<-c("image","image_2","image_3","image_4","image_5","image_6","image_7","image_8")
for(id in 1:8){
  SpatialPlot(st,features=colnames(st@meta.data[14:24]),images = image[id],
              pt.size.factor= point_size[id],ncol=6)
  ggsave(file=paste0("plot_3/",id,"cell_spatial.pdf"),width=30,height=13)
}

####addinferCNV score####
# Scoring
cnv_table <- read.table("data/cnv_ts/infercnv.observations.txt", header=T)
# Score cells based on their CNV scores 
# Replicate the table 
cnv_score_table <- as.matrix(cnv_table)
cnv_score_mat <- as.matrix(cnv_table)
cnv_score_table_pts<-abs(cnv_score_mat-1)

cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
head(cell_scores_CNV)
write.csv(x = cell_scores_CNV, file = "data/cnv_ts/cnv_scores.csv")
cell_scores_CNV <- read.csv("data/cnv_ts/cnv_scores.csv")
Idents(st)<-"sample"
library(do)
cell_scores_CNV$X<-Replace(data=cell_scores_CNV$X,from="X",to="")
cell_scores_CNV$X<-Replace(data=cell_scores_CNV$X,from=".1",to="-1")
a<-cell_scores_CNV$X
cell_scores_CNV$X<-NULL
rownames(cell_scores_CNV)<-a

cnv_score<-data.frame(barcode=rownames(st@meta.data))
a<-rownames(st@meta.data)
b<-rownames(cell_scores_CNV)
rownames(cnv_score)<-a
name_1<-b[b%in%a]
library(Hmisc)
name_2<-a[a%nin%b]
cnv_score[name_1,"cnv_score"]<-cell_scores_CNV[name_1,"cnv_score"]
cnv_score[name_2,"cnv_score"]<-5
a<-cnv_score$barcode
rownames(cnv_score)<-a
cnv_score$barcode <- NULL
st<-AddMetaData(st,cnv_score,col="cnv_score")

point_size<-c(1.5,1.9,2,1.5,1.9,2,1.5,1.4)
image<-c("image","image_2","image_3","image_4","image_5","image_6","image_7","image_8")
for(id in 1:8){
  df<-st@meta.data[st$sample == paste0("patient",id),"cnv_score"]
  SpatialFeaturePlot(st,features="cnv_score",images = image[id],
                     pt.size.factor= point_size[id],image.alpha = 0)+
    scale_fill_gradientn(colors = pals::kovesi.linear_blue_5_95_c73(20),limits= range(df))
  ggsave(file=paste0("plot_3/",id,"p_cnv_score_spatial.pdf"),width=7.5,height=6)
}


normalized<-data.frame(normalize=cnv_score$cnv_score/max(cnv_score$cnv_score))
rownames(normalized)<-rownames(cnv_score)
st<-AddMetaData(st,normalized,col="n_cnv_score")
point_size<-c(1.5,1.9,2,1.5,1.9,2,1.5,1.4)
image<-c("image","image_2","image_3","image_4","image_5","image_6","image_7","image_8")
for(id in 1:8){
  SpatialPlot(st,features="n_cnv_score",images = image[id],
              pt.size.factor= point_size[id],image.alpha = 0)+
    scale_fill_gradientn(colors = pals::kovesi.linear_blue_5_95_c73(20),limits =range(st@meta.data$n_cnv_score))
  ggsave(file=paste0("plot_3/",id,"n_nv_score_spatial.pdf"),width=7.5,height=6)
}

####Tumor Score####
df<-st@meta.data[,c("Tumour.cells","n_cnv_score")]
score<-c()
for(i in rownames(df)){
  key1<-df[i,"n_cnv_score"]
  key2<-df[i,"Tumour.cells"]
  score<-append(score,as.numeric(sqrt(key1*key2)))
}
df[,"tumor_score"]<-score
df$Tumour.cells<-NULL
df$n_cnv_score<-NULL
rownames(df)<-rownames(st@meta.data)
st<-AddMetaData(st,df,col="tumor_score")

point_size<-c(1.5,1.9,2,1.5,1.9,2,1.5,1.4)
image<-c("image","image_2","image_3","image_4","image_5","image_6","image_7","image_8")
for(id in 1:8){
  df<-st@meta.data[st$sample == paste0("patient",id),"tumor_score"]
  SpatialPlot(st,features="tumor_score",images = image[id],
              pt.size.factor= point_size[id],image.alpha = 0)+
    scale_fill_gradientn(colors = pals::kovesi.linear_blue_5_95_c73(20),limits =range(df))
  ggsave(file=paste0("plot_3/",id,"p_tumor_score_spatial.pdf"),width=7.5,height=6)
  
  SpatialPlot(st,features="tumor_score",images = image[id],
              pt.size.factor= point_size[id],image.alpha = 0)+
    scale_fill_gradientn(colors = pals::kovesi.linear_blue_5_95_c73(20),limits =range(st@meta.data$tumor_score))
  ggsave(file=paste0("plot_3/",id,"a_tumor_score_spatial.pdf"),width=7.5,height=6)
}
