library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringi)
library(scales)
####scRNA seq Preporocessing####
ref_list<-c(1:5)
scRNA_list <- sapply(ref_list,function(id){
  path <- paste0("data/ref_new/",id,".txt")
  #path <- paste0("ref_new/",id,".txt")
  #ref <- Read10X(data.dir = path)
  ref <- read.table(file=path,head=TRUE,sep="\t")
  ref <- CreateSeuratObject(counts = ref, project = "pbmc3k", min.cells = 3, min.features = 500)
  #计算线粒体基因占总基因数目的百分比，pattern后可以是正则表达式。
  ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
  #小提琴图展示meta.data中的数据，ncol 表示显示多个图时的列数
  VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(ref, feature1 = "nCount_RNA", feature2 = "percent.mt") 
  plot2 <- FeatureScatter(ref, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  plot1 + plot2
  ref <- subset(ref, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15)
  #normalization
  ref <- SCTransform(ref, assay = "RNA", verbose = T)
  return(ref)
})

####Add Annotation####
library(gridExtra)
library(stringi)
library(ROCit)
library(png)
library(readxl)
annot<- read_excel("data/annot.xlsx")
annot_list <- c("MJ10","MJ11","Y2","Y3","Y5")
ref_list<-c(1:5)
#ADD annotation
for(scRNA in ref_list){
  #correct ident
  scRNA_list[[scRNA]]$orig.ident <- scRNA
  annot_table <-annot[annot$Sample==annot_list[scRNA],]
  colnames(annot_table)<-c("Barcode","Sample","Cells")
  a<-annot_table$Barcode
  annot_table$Barcode <- NULL
  annot_table$Sample<- NULL
  rownames(annot_table) <- a
  scRNA_list[[scRNA]] <- AddMetaData(object= scRNA_list[[scRNA]],
                                     metadata = annot_table,
                                     col.name = "cell")
}
saveRDS(scRNA_list,file="data/scRNA_list.rds")
save(ref, file="ref_new/ref.RData")
####Batch Effect####
integ_features <- SelectIntegrationFeatures(object.list = scRNA_list, 
                                            nfeatures = 3000) 
scRNA_list <- PrepSCTIntegration(object.list = scRNA_list, 
                                 anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = scRNA_list, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
ref <- IntegrateData(anchorset = integ_anchors,normalization.method = "SCT",k.weight=76)
saveRDS(ref,file="data/ref.rds")

####scRNA Seq for MIA####
Idents(ref) <- 'cell'
ref.markers <- FindAllMarkers(ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ref.markers<- read_csv("data/ref_mia_marker.csv")
ref.markers['cluster'] %>% summary(maxsum=50) 
saveRDS(ref.markers,file="data/ref_markers_mia.rds")

####ST data for MIA####
load("data/st.RData")
Idents(st)<-"class"
pbmc<-subset(st,ident=c("tumor","stromal"))
#when perform analysis for between stromal groups, choose ident = "stromal"
Idents(pbmc)<-"identstate"
#divide stromal cluster in to two group based on their GSVA result
#new.cluster.ids <- c("group1","group1","group1","group1","group1","group1","group1","group1",
#                     "group2","group2","group2", "group2","group2","group2",
#                     "group2","group2","group2","group2","group2","group2","group2")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$stromal_deg<-pbmc@active.ident
Idents(pbmc)<-"stromal_deg"
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers <-pbmc.markers[pbmc.markers$p_val_adj<0.05&pbmc.markers$avg_log2FC>=0.5,]
write.csv(pbmc.markers,file="data/st_mia_identitycluster_markers.csv")
pbmc.markers['cluster'] %>% summary()
#saveRDS(pbmc.markers,file="data/mia_CNV_markers_st.rds")

#### Create a list object containing the marker genes for each ST region:####
pbmc.clusts <- Idents(pbmc)%>% levels()
N <- length(pbmc.clusts)
pbmc.marker.list <- vector(mode = 'list', length = N)
names(pbmc.marker.list) <- pbmc.clusts
for(i in pbmc.clusts) {
  pbmc.marker.list[[i]] <- pbmc.markers[pbmc.markers$cluster == i,'gene']
}

# Create a list object containing the marker genes for each cell type:
ref.clusts <- Idents(ref) %>% levels()
M <- length(ref.clusts)
ref.marker.list <- vector(mode = 'list', length = M)
names(ref.marker.list) <- ref.clusts
for (i in ref.clusts) {
  ref.marker.list[[i]] <- ref.markers[ref.markers$cluster == i,'gene']
}

####MIA####
# Initialize a dataframe for us to store values in:
N <- length(pbmc.clusts) ; M <- length(ref.clusts)
MIA.results <- matrix(0,nrow = M, ncol = N)
row.names(MIA.results) <- ref.clusts
colnames(MIA.results) <- pbmc.clusts
# Gene universe
gene.universe <- length(rownames(pbmc))
# Loop over ST clusters
for (i in 1:N) {
  # Then loop over SC clusters
  for (j in 1:M) {
    genes1 <- pbmc.marker.list[[pbmc.clusts[i]]]
    genes2 <- ref.marker.list[[ref.clusts[j]]]
    
    # Hypergeometric    
    A <- length(intersect(genes1,genes2))
    B <- length(genes1)
    C <- length(genes2)
    enr <- -log10(phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    dep <- -log10(1-phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    if (enr < dep) {
      MIA.results[j,i] = -dep
    } else {
      MIA.results[j,i] = enr
    }
  }
}
# Some results were -Inf...check why this is the case...
MIA.results[is.infinite(MIA.results)] <- 0
# Visualize as heatmap
heatmap_df <- data.frame('Cell types' = melt(MIA.results)[,1],
                         'Tissue regions' = melt(MIA.results)[,2],
                         enrichment = melt(MIA.results)[,3])
heatmap_df$Tissue.regions <- factor(heatmap_df$Tissue.regions, levels=order)
p<-ggplot(data = heatmap_df, aes(x = Tissue.regions, y = Cell.types, fill = enrichment)) +
  geom_tile() + 
  scale_fill_gradient2(low = "#7EC0EE", high = "pink", mid = "white", #low = "#104E8B", high = "#B4EEB4", mid = "#00CED1"
                       midpoint = 10, limit = c(0,20),space = "Lab", 
                       oob=squish, name="Enrichment \n -log10(p)") +
  ylim(heatmap_df$Cell.types %>% levels() %>% sort() %>% rev())
ggsave(
  filename = "mia_stromal_DEG.pdf", 
  width = 5,             
  height = 7,           
  units = "in",          
  dpi = 600              
)
