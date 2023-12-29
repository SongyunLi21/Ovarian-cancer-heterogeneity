####load package####
library(dplyr)
library(rhdf5)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringi)
library(ROCit)
library(png)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(reshape2)
library(SeuratObject)
library(DESeq2)
library(clusterProfiler)
library(pheatmap)
library(ggnewscale)
library(enrichplot)
library(sceasy)
library(reticulate)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(msigdbr)  
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(Spaniel)
library(RColorBrewer)
####load data####
st_index_list<-c(1:8)
st_list <- sapply(st_index_list,function(id){
  # Load the expression data
  expr.url <- paste0("data/st1/", id)
  expr.data <- Read10X(data.dir =  expr.url )
  stref <- Seurat::CreateSeuratObject(counts = expr.data, project = 'anterior1', assay = 'Spatial')
  # Load the image data
  path <- paste0('data\\st1\\',id,'/spatial')
  img <- Seurat::Read10X_Image(image.dir = path)
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = stref)]
  stref[['image']] <- img
  #calculate the percentage of MT genes。
  stref[["percent.mt"]] <- PercentageFeatureSet(stref, pattern = "^MT-")
  #plot the portion of MT and feature scores
  #VlnPlot(stref, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3)
  #plot <- SpatialFeaturePlot(stref, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"))
  #plot 
  #plot1 <- FeatureScatter(stref, feature1 = "nCount_Spatial", feature2 = "percent.mt") 
  #plot2 <- FeatureScatter(stref, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") 
  #plot1 + plot2
  stref <- subset(stref, subset = nFeature_Spatial > 400)#fitering
  stref$sample<-paste0("patient",id)#patient id 
  return(stref)
})
saveRDS(st_list,"data/st_list_h.rds")#save

####Integrated data and clustering####
#merge data 
library(harmony)
st <- merge(st_list[[1]], y = st_list[c(2:8)],add.cell.ids = c(1:8), merge.data = TRUE)
unique(sapply(X = strsplit(colnames(st), split = "_"), FUN = "[", 1))
#normalization
st <- NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 1e4)
#dim reduction
st <- FindVariableFeatures(st,selection.method = 'vst', nfeatures = 2000)
st <- ScaleData(st)
st <- RunPCA(st, features = VariableFeatures(object = st))
#use harmony for integration
st <- RunHarmony(st, group.by.vars = "sample")
saveRDS(st,"data/st_h.rds")
#save(st,file="data/st_h.RData")
#clustering
st <- RunUMAP(st, dims = 1:20,reduction = "harmony")
st <- FindNeighbors(st, dims = 1:20,reduction = "harmony")
st <- FindClusters(st, resolution = 0.8)
Idents(st)<-"state"
DimPlot(st, reduction = 'umap',label = FALSE, cols=color)
ggsave(file=paste0("plot_1/","all_patient_state","cluster_dim_reduction",".pdf"),width=7,height=7)


####Identify the markers for clusters####
#identify the markers
st.markers <- FindAllMarkers(st, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(st.markers,"data/st_h_markers.rds")
top10markers<-st.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
features<- st.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
features <- features$gene
DotPlot(st, features = unique(features)) + RotatedAxis()
ggsave(file=paste0("plot_1/","all_","dot_plot.pdf"),height=10,width=10,dpi=600)
#visulization of the markers
a<-brewer.pal(name="PuBu", n = 9)
FeaturePlot(object = st,label.size = 2,cols=c(a[1], a[6]),
            features = c("WFDC2", "EPCAM", "CLDN4", "MUC16", "KRT8", "KRT18", 
                         "KRT19","KRT7", "MUC1","DSC2"))
ggsave(file=paste0("tumor","_feature_plot",".pdf"),width=12,height=10)
FeaturePlot(object = st,label.size = 0.5,cols=c(a[1], a[6]),ncol = 1,
            features = c("WFDC2", "EPCAM", "KRT8"))
ggsave(file=paste0("plot_1/","tumor","_feature_plot",".pdf"),width=3,height=10)
b<-brewer.pal(name="RdPu", n = 9)
FeaturePlot(object = st,label.size = 1,cols=c(b[1], b[5]),
            features = c("COL1A1", "COL1A2","DCN", "VIM","FN1","ACTA2",
                         "CD36"))
ggsave(file=paste0("stromal","_feature_plot",".pdf"),width=9,height=10)
FeaturePlot(object = st,label.size = 0.5,cols=c(b[1], b[5]),ncol = 1,
            features = c("COL1A1","DCN","FN1"))
ggsave(file=paste0("plot_1/","stromal_feature_plot",".pdf"),width=3,height=10)

####Annotation####
#annotate each cluster
Idents(st)<-"seurat_clusters"
#new.cluster.ids <- c("tumor_1","stromal_1","tumor_2","hybrid_1","stromal_2","stromal_3",
#                     "stromal_4","stromal_5","tumor_3","tumor_4","tumor_5","hybrid_2",
#                     "stromal_6","tumor_6","stromal_7","stromal_8","unknown","stromal_9",
#                     "tumor_7","tumor_8")
new.cluster.ids <- c("tumor_1","stromal_1","tumor_2","hybrid_1","stromal_2","stromal_3","unknown_1","stromal_4","tumor_3","unknown_2",
                     "unknown_3","unknown_4","unknown_5","tumor_4","unknown_6","stromal_5","unknown_7","stromal_6","tumor_5","unknown_8")
names(new.cluster.ids) <- levels(st)
st <- RenameIdents(st, new.cluster.ids)
st$state<-st@active.ident
#annotate tumor, stromal and unknown region
Idents(st)<-"seurat_clusters"
#new.cluster.ids <- c("tumor","stromal","tumor","hybrid","stromal","stromal",
#                     "stromal","stromal","tumor","tumor","tumor","hybrid",
#                     "stromal","tumor","stromal","stromal","unknown","stromal",
#                     "tumor","tumor")
new.cluster.ids <- c("tumor","stromal","tumor","hybrid","stromal","stromal","unknown","stromal","tumor","unknown",
                     "unknown","unknown","unknown","tumor","unknown","stromal","unknown","stromal","tumor","unknown")
names(new.cluster.ids) <- levels(st)
st <- RenameIdents(st, new.cluster.ids)
st$class<-st@active.ident
#Annotation using sample + tumor, stromal and unknown region
df<-st@meta.data[,c("class","sample")]
for(i in 1:8){
  df_sub<-df[df$sample==paste0("patient",i),]
  e<-c()
  for(j in rownames(df_sub)){
    e<-append(e,paste0(df_sub[j,"class"],"_",df_sub[j,"sample"]))
  }
  df_sub$identclass<-e
  colnames(df_sub)<-c("class","sample","identclass")
  if(i==1){
    D<-df_sub
  }else{
    D<-rbind(D,df_sub)
  }
}
D$class<-NULL
D$sample<-NULL
st<-AddMetaData(st,D,col="identclass")
Idents(st)<-"identclass"
#Annotation using sample + cluster id
df<-st@meta.data[,c("state","sample")]
for(i in 1:8){
  df_sub<-df[df$sample==paste0("patient",i),]
  e<-c()
  for(j in rownames(df_sub)){
    e<-append(e,paste0(df_sub[j,"state"],"_",df_sub[j,"sample"]))
  }
  df_sub$identstate<-e
  colnames(df_sub)<-c("state","sample","identstate")
  if(i==1){
    D<-df_sub
  }else{
    D<-rbind(D,df_sub)
  }
}
D$state<-NULL
D$sample<-NULL
st<-AddMetaData(st,D,col="identstate")
Idents(st)<-"identstate"

Idents(st)<-"class"
st_sub<-subset(st,ident="tumor")
Idents(st_sub)<-"identstate"

####color for spatial plot and dimensional reduction plot####
color<-c()
order<-c()
tumor<-c()
stromal<-c()
other<-c()
color_list<-list()
n1<-0
n2<-0
n3<-0
color_list[[1]]<-brewer.pal(name="RdPu", n = 9)
color_list[[2]]<-brewer.pal(name="YlGn", n = 9)
color_list[[3]]<-brewer.pal(name="YlOrRd", n = 9)
ident<- c("tumor","stromal","tumor","hybrid","stromal","stromal","unknown","stromal","tumor","unknown",
          "unknown","unknown","unknown","tumor","unknown","stromal","unknown","stromal","tumor","unknown")
for(i in 1:length(ident)){
  if(ident[i]=="tumor"){
    tumor<-append(tumor,i)
    n1<-n1+1
    color<-append(color,color_list[[1]][7-n1])
  }else if(ident[i]=="stromal"){
    stromal<-append(stromal,i)
    n2<-n2+1
    color<-append(color,color_list[[2]][7-n2])
  }else{
    other<-append(other,i)
    n3<-n3+1
    color<-append(color,color_list[[3]][n3])
  }
}

#color for dim_reduction plot
new.cluster.ids <- c("tumor_1","stromal_1","tumor_2","hybrid_1","stromal_2","stromal_3","unknown_1","stromal_4","tumor_3","unknown_2",
                     "unknown_3","unknown_4","unknown_5","tumor_4","unknown_6","stromal_5","unknown_7","stromal_6","tumor_5","unknown_8")
order<-new.cluster.ids[c(tumor,stromal,other)]
Idents(st)<-"state"
Idents(st) <- factor(Idents(st), levels= order)
st$state<-st@active.ident

#color for spatial plot
color<-c(color_list[[1]][c(6,5,4,3,2)],color_list[[2]][c(6,5,4,3,2,1)],color_list[[3]][c(1:9)])
color_d<-list()
for(i in 1: length(levels(st@active.ident))){
  key<-levels(st@active.ident)[i]
  color_d[[key]]<-color[i]
}

####spatial plot####
#spatial plot of tumor, stromal, unknown region
#Idents(st)<-"state"
Idents(st)<-"class"
SpatialDimPlot(st, label = F, label.size = 3,ncol=3, image.alpha=0,pt.size.facto=2)
ggsave(file=paste0("plot_1/","all_patient_class","_spatial_plot",".pdf"),width=15,height=15)

#spatial plot of annotated clusters(all clusters in one plot)
Idents(st)<-"state"
#Idents(st)<-"seurat_clusters"
point_size<-c(1.6,2,2,1.6,2,2,1.6,1.6)
image<-c("image","image_2","image_3","image_4","image_5","image_6","image_7","image_8")
df<-as.data.frame(table(st$state,st$sample))
for(i in 1:8){
  key<-as.character(df[df$Var2==paste0("patient",i)&df$Freq!=0,"Var1"])
  color_factor<-c()
  for(j in key){
    color_factor<-append(color_factor,color_d[[j]])
  }
  SpatialDimPlot(st, images = image[i], label = FALSE, label.size = 1,ncol=1,pt.size.factor =point_size[i],image.alpha=0)+scale_fill_manual(values=color_factor) 
  ggsave(file=paste0("plot_1/","all","_state_cluster_spatial_h",i,".pdf"),width=5,height=5)
}

#spatial plot of tumor and stromal clusters(one cluster one plot)
Idents(st)<-"seurat_clusters"
key<-c("tumor","stromal")
ident<-list()
ident[[1]]<-c(0,2,8,13,18)
ident[[2]]<-c(1,4,5,7,15,17)
color<-list()
color[[1]]<-c("gold","black")
color[[2]]<-c("pink","black")
for(i in 1:8){
  for(j in 1:2){
    SpatialDimPlot(st, images = image[[i]] ,cells.highlight = CellsByIdentities(object = st, idents = ident[[j]]), 
                   facet.highlight = TRUE,label = FALSE, label.size = 3,ncol=3,
                   image.alpha=0,cols.highlight = color[[j]],
                   pt.size.factor =point_size[[i]])
    ggsave(file=paste0("plot_1/",key[j],"cluster_spatial_h",i,".pdf"),width=20,height=20)
  }
}

####Percentage chart####
library(scales)
cluster<-list()
cluster[[1]]<-"tumor"
cluster[[2]]<-"stromal"
cluster[[3]]<-c("hybrid","unknown")
cluster[[4]]<-c("tumor","stromal","hybrid","unknown")
key<-c("Tumor","Stromal","Other","all")
color<-list()
color[[1]]<-color_list[[1]][c(6,5,4,3,2)]
color[[2]]<-color_list[[2]][c(6,5,4,3,2,1)]
color[[3]]<-color_list[[3]][c(1:9)]
color[[4]]<-c(color_list[[1]][c(6,5,4,3,2)],color_list[[2]][c(6,5,4,3,2,1)],color_list[[3]][c(1:9)])
for(i in 1:4){
  df<-st@meta.data[st@meta.data$class==cluster[[i]],]
  df<-df[order(df$state,decreasing = TRUE),]
  ggplot(df, aes(x =sample , fill = state)) +
    geom_bar(width = 0.3, position = "fill") + 
    scale_fill_manual(values = color[[i]]) +     
    scale_y_continuous(labels = percent) +     
    guides(fill=guide_legend(title = "cluster")) +
    labs(title = paste0(key[i],"_portion"),
         x = "patients",
         y = "percentage") +
    theme_minimal() +
    coord_flip() 
  ggsave(
    filename = paste0("plot_1/bar",key[i],"_bar.pdf"), 
    width = 10,            
    height = 10,           
    units = "in",          
    dpi = 600              
  )
}
ggplot(df, aes(x =sample , fill = class)) +
  geom_bar(width = 0.3, position = "fill") + 
  #scale_fill_brewer(palette = kovesi.rainbow_bgyrm_35_85_c71(20)) +     # 调色板{RColorBrewer}
  scale_y_continuous(labels = percent) +     
  guides(fill=guide_legend(title = "cluster")) +
  labs(title = "all_portion",
       x = "patients",
       y = "percentage") +
  theme_minimal() +
  coord_flip() # 倒转x与y轴
ggsave(
  filename = paste0("plot/bar","all_","_bar.pdf"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 10,             # 宽
  height = 10,            # 高
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)


####GSVA####
library(GSEABase)
library(GSVA)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
#load reference datasets
library(msigdbr)
#geneset_1 = msigdbr(species = "Homo sapiens",category='C8')
#length(unique(table(geneset_1$gs_name)))
#tail(table(geneset_1$gs_name))
geneset_2 = msigdbr(species = "Homo sapiens",category='H')
length(unique(table(geneset_2$gs_name)))
tail(table(geneset_2$gs_name))
#all_gene_sets<-rbind(geneset_1,geneset_2)
all_gene_sets<-geneset_2
length(unique(table(all_gene_sets$gs_name)))
tail(table(all_gene_sets$gs_name))
gs=split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs = lapply(gs, unique)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
geneset<-gsc
#expression matrix
df<-table(st$identstate)
df<-as.data.frame(df)
name<-as.character(df$Var1[df$Freq>100])
Idents(st)<-"identstate"
st_sub<-subset(st,ident=name)
Idents(st_sub)<-"class"
#st_sub<-subset(st_sub,ident= c("stromal"))#or "tumor"
exp <-AverageExpression(st_sub,
                        group.by = "identstate",
                        assays = "Spatial",
                        slot = "data")
exp<-exp$Spatial
cg=names(tail(sort(apply(exp,1, sd)),50))
pheatmap::pheatmap(cor(exp[cg,]),
                   filename = "Plot_1/all_deg_cor.pdf",width = 7, height = 7)

#GSVA
es.max <- gsva(exp, geneset)
#筛选+画图
library(RColorBrewer)
cg = names(tail(sort(apply(es.max, 1, sd) ),50))
pheatmap(es.max[cg,],show_rownames = T,fontsize_row = 7,fontsize_col=7,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         filename = "Plot_1/all_deg_GSVA.pdf",width = 5, height = 10)
