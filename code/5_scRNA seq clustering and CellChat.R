library(CellChat)
####Load data####
load("data/ref.RData")
Idents(ref)<-"cell"
####Macrophages####
#clustering
ref_sub<-subset(ref,ident="Macrophages")
ref_sub <- FindVariableFeatures(ref_sub, selection.method = 'vst', nfeatures = 2000)
ref_sub <- RunPCA(ref_sub, features = VariableFeatures(ref_sub)) 
ref_sub <- FindNeighbors(ref_sub, dims = 1:10)
ref_sub<- FindClusters(ref_sub, resolution = 0.5, graph.name = "integrated_nn")
ref_sub <- RunUMAP(ref_sub, dims = 1:10)
DimPlot(ref_sub, reduction = 'umap',label = F,cols=brewer.pal(name="PiYG",n=4),label.size = 1, pt.size=2.5)
ggsave(file=paste0("plot_5/clustering/","Macrophage_","cluster_dim_reduction",".pdf"),width=7,height=5)
#identify the markers for each cluster
ref.markers <- FindAllMarkers(ref_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ref.markers<-ref.markers[ref.markers$p_val_adj<0.05,]
write.csv(ref.markers,file="data/macrophages.csv")
ref.markers<-read.csv("data/macrophages.csv")
top10markers<-ref.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#features<- ref.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#annotation
FeaturePlot(object = ref_sub,
            features = c("C1QA", "CQQB", "C1QC","CD83","HLA-DQA1","HLA-DQB1","HLA-DEB5",
                         "S100A8","S100A9","S100A12","VCAN","MARCO","FCN1","SIGLEC1","CX3CR1"))
FeaturePlot(object = ref_sub,
            features =c("ID3", "MITF", "RUNX2", "MAF",
                        "BCL3","NR4A1", "RXRA", "TCF25"))
FeaturePlot(object = ref_sub,
            features = c("CD14","CD16","CD1C","CLEC10A"))
new.cluster.ids <- c("HLA TAM1","HLA TAM2","PKM TAM","S100A8 monocytes")
Idents(ref_sub)<-"seurat_clusters"
names(new.cluster.ids) <- levels(ref_sub)
ref_sub <- RenameIdents(ref_sub, new.cluster.ids)
ref_sub@meta.data[,"new_cluster"]<-ref_sub@active.ident
macrophages<-data.frame(barcode=rownames(ref_sub@meta.data))
macrophages$cluster<-ref_sub@meta.data$new_cluster
features<-c("HLA-DQA1","HLA-DQB1","HLA-DEB5","APOE","S100A8","S100A9","CD14","LITAF","CTSH","LAP3")
DotPlot(ref_sub, features = unique(features),cols = "RdYlBu") #+ RotatedAxis()
#pals::kovesi.diverging_isoluminant_cjo_70_c25(2)
ggsave(file=paste0("plot/","Macrophage_","dot_plot.pdf"),height=5,width=15,dpi=600)
DoHeatmap(ref_sub, features = top10markers$gene, size = 3,group.colors = brewer.pal(4,"Set1"))+
  scale_fill_gradientn(colors =pals::kovesi.linear_bmy_10_95_c71(20)) 
ggsave(file=paste0("plot/","Macrophage_","heat_plot.pdf"),height=3,width=3,dpi=600)
#GSVA
library(GSEABase)
library(GSVA)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
#加载MsigDb数据集
library(msigdbr)
geneset_1 = msigdbr(species = "Homo sapiens",category='C8')
length(unique(table(geneset_1$gs_name)))
tail(table(geneset_1$gs_name))
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
ref_sub$state<-ref@active.ident
exp <-AverageExpression(ref_sub ,
                        group.by = "new_cluster",
                        assays = "integrated",
                        slot = "data")
exp<-exp$integrated
cg=names(tail(sort(apply(exp,1, sd)),50))
pheatmap::pheatmap(cor(exp[cg,]))

#GSVA
es.max <- gsva(exp, geneset)
#筛选+画图
library(RColorBrewer)
cg = names(tail(sort(apply(es.max, 1, sd) ),50))
pheatmap(es.max[cg,],show_rownames = T,fontsize_row = 7,fontsize_col=7,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))

####Tumor cells####
#redifine tumor cell using CNV score and tumor markers
sce<-ref
remove(ref)
library(readr)
cnv_scores <- read_csv(paste0("data/ref_new/cnv_scores",".csv"))
#cnv_scores$barcode<-Replace(cnv_scores$barcode,from="_1",to="_1")
colnames(cnv_scores)<-c("barcode","cnv_scores","cancer")
library(NMF)
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
sce<-AddMetaData(sce,D,col="cancer")
sce@meta.data$cancer[sce@meta.data$cancer>=0.3]<-"tumor"
sce@meta.data$cancer[sce@meta.data$cancer<0.3]<-"normal"
table(sce@meta.data$cell,sce@meta.data$cancer)
sce@meta.data$cell[sce@meta.data$cancer=="tumor"&sce@meta.data$cell=="Mesothelial cells"]<-"Tumour cells"
sce@meta.data$cell[sce@meta.data$cancer=="tumor"&sce@meta.data$cell=="Endothelial cells"]<-"Tumour cells"
table(sce@meta.data$cell,sce@meta.data$cancer)
Idents(sce)<-"cell"
sub_sce<-subset(sce,ident="Tumour cells" )
#clustering
sub_sce <- FindVariableFeatures(sub_sce, selection.method = 'vst', nfeatures = 2000)
sub_sce <- RunPCA(sub_sce, features = VariableFeatures(sub_sce)) 
sub_sce <- FindNeighbors(sub_sce, dims = 1:10)
sub_sce<- FindClusters(sub_sce, resolution = 0.5, graph.name = "integrated_nn")
sub_sce <- RunUMAP(sub_sce, dims = 1:10)
DimPlot(sub_sce, reduction = 'umap',label = TRUE,cols="Paired")
ggsave(file=paste0("plot_5/clustering/","tumor_","cluster_dim_reduction",".pdf"),width=7,height=7)
ref<-sce
remove(sce)

new.cluster.ids <- c("Tumor0","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
names(new.cluster.ids) <- levels(sub_sce)
sub_sce <- RenameIdents(sub_sce, new.cluster.ids)
sub_sce@meta.data[,"new_cluster"]<-sub_sce@active.ident

ref$new_cluster<-c(sub_sce$new_cluster,ref_sub$new_cluster)
for(i in rownames(ref@meta.data[ref@meta.data$cell=="Tumour cells"|ref@meta.data$cell=="Macrophages",])){
  ref@meta.data[i,"cell"]<-as.character(ref@meta.data[i,"new_cluster"])
}
Idents(ref)<-"new_cluster"
ref.markers <- FindAllMarkers(ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ref.markers,"data/ref_subcluster.csv")
ref.markers['cluster'] %>% summary(maxsum=50) 

#calculate the score of tumor and macrophages in ST data
exp<-st@assays$Spatial@counts
name<-c("HLA TAM1","HLA TAM2","PKM TAM","S100A8 monocytes","Tumor0","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
marker<-c()
score<-data.frame(colnames(exp))
for(i in 1:length(name)){
  marker<-ref.markers[ref.markers$cluster==name[i],]$gene
  marker<-marker[marker%in%rownames(exp)]
  cell_exp<-exp[marker,]
  cell_score<-colMeans(cell_exp)
  score[i]<-cell_score
  marker<-c()
}
rownames(score)<-colnames(exp)
colnames(score)<-name
st<-AddMetaData(st,score,col=colnames(score))
SpatialPlot(st,features=name,image="image_8",pt.size.factor =2,ncol=4,image.alpha = 0)#+ggsave(file=paste0("data/",id,"spatial.pdf"),width=15,height=15)

marker_list<-list(1:length(name))
for(i in 1:length(name)){
  marker_list[[name[i]]]<-ref.markers[ref.markers$cluster==name[i],]$gene
}
marker_list[[1]]<-NULL
#heatmap of the enrchiment of each cell types in each sample of ST data
exp <-AverageExpression(st ,
                        group.by = "sample",
                        assays = "Spatial",
                        slot = "scale.data")
exp<-exp$Spatial

library(fgsea)
for(i in 1:length(colnames(exp))){
  result <- fgsea(stats =exp[,i], pathways=sig)#marker_list
  if(i == 1){
    Result<-data.frame(result$ES)
    rownames(Result)<-result$pathway
  }else{
    Result<-cbind(Result,result$ES)
  }
}

colnames(Result)<-colnames(exp)
library(RColorBrewer)
pheatmap(Result,show_rownames = T,fontsize_row = 14,fontsize_col=14,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),file="plot_5/gsea_pg.pdf",width=8,height=4)#

####cell chat####
sub_sce@meta.data[,"seurat_clusters"]<-sub_sce@active.ident
tumor<-data.frame(barcode=rownames(sub_sce@meta.data))
tumor$cluster<-sub_sce@meta.data$new_cluster
cluster<-rbind(tumor,macrophages)
a<-cluster$barcode
cluster$barcode<-NULL
rownames(cluster)<-a

name<-c("HLA TAM1","HLA TAM2","PKM TAM","S100A8 monocytes","Tumor0","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
Idents(ref)<-"new_cluster"
ref_sub<-subset(ref,ident=name)
Idents(ref_sub)<-"new_cluster"
#ref_sub@meta.data$TAM<-ref_sub@active.ident
data.input = ref_sub@assays$SCT@data
meta = ref_sub@meta.data # a dataframe with rownames containing cell mata data
cell.use = rownames(meta) # extract the cell names from disease data

# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "new_cluster")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "new_cluster") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

#Preprocessing
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)
#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
#general circle plot
#number of interactions or the total interaction strength (weights) between any two cell groups
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#can compare edge weights between differet networks.
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#circle plot for specific genes
#access by cellchat@netP$pathways
pathways.show <- c("GAS","TNF","VEGF","ANNEXIN","CXCL","CCL","MK","MIF")
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq() # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
svglite(paste0("plot_5/CellChat/specific/",pathways.show,"Circle plot",".svg"),width=10,height=10)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
#> Plot the aggregated cell-cell communication network at the signaling pathway level
# Heatmap
library(svglite)
svglite(paste0("plot_5/CellChat/","pathways.show","all","heatmap",".svg"),width=7,height=7)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
#> Do heatmap based on a single object
svglite(paste0("plot/CellChat","pathways.show","all","barplot",".svg"),width=10,height=10)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

#chord for specific L-R
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
for(i in 1:17){
  LR.show <- pairLR.CXCL[i,] # show one ligand-receptor pair
  # Chord diagram
  svglite(paste0("plot_5/CellChat/specific/",LR.show,"Circle plot",".svg"),width=8,height=8)
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  dev.off()
}
# Hierarchy plot
vertex.receiver = seq() # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
svglite(paste0("plot/","pathways.show","1","all","Circle plot",".svg"),width=10,height=10)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
dev.off()
#> [[1]]

#save
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
#general result
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
name<-c("TAM0","TAM1","TAM2","TAM3")
key<-c(7,8,9,10)
for(i in 1:4){
  netVisual_bubble(cellchat, sources.use =key[i], targets.use = c(1:6), remove.isolate = FALSE)
  ggsave(paste0("plot_5/CellChat/",name[i],"dot_Plot",".pdf"),width=10,height=10,dpi=600)
}
netVisual_bubble(cellchat, sources.use =c(1:6), targets.use = c(7:10), remove.isolate = FALSE)
ggsave(paste0("plot_5/CellChat/overall/","Tumor_all","dot_Plot",".pdf"),width=10,height=10,dpi=600)
name<-c("Tumor0","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
key<-c(1,2,3,4,5,6)
for(i in 1:6){
  netVisual_bubble(cellchat, sources.use =key[i], targets.use = c(7:10), remove.isolate = FALSE)
  ggsave(paste0("plot_5/CellChat/",name[i],"dot_Plot",".pdf"),width=10,height=10,dpi=600)
}
netVisual_bubble(cellchat, sources.use =c(7:10), targets.use = c(7:10), remove.isolate = FALSE)
ggsave(paste0("plot_5/CellChat/overall/","TAM_to_TAM_all","dot_Plot",".pdf"),width=10,height=10,dpi=600)
#> Comparing communications on a single object
#chord plot
name<-c("HLA TAM1","HLA TAM2","PKM TAM","S100A8 monocytes")
key<-c(7,8,9,10)
for(i in 1:4){
  svglite(paste0("plot_5/CellChat/","TAM_TO_Tumor",name[i],"circle_Plot",".svg"),width=7,height=7)
  netVisual_chord_gene(cellchat, sources.use = as.numeric(key[i]), targets.use = c(1:6), slot.name = "netP", legend.pos.x = 10)
  dev.off()
}
svglite(paste0("plot_5/CellChat/","TAM_TO_Tumor","all","circle_Plot",".svg"),width=10,height=10)
netVisual_chord_gene(cellchat, sources.use = c(7:10), targets.use = c(1:6), slot.name = "netP", legend.pos.x = 10)
dev.off()

name<-c("Tumor0","Tumor1","Tumor2","Tumor3","Tumor4","Tumor5")
key<-list(1,2,3,4,5,6)
for(i in 1:6){
  svglite(paste0("plot_5/CellChat/","Tumor_TO_TAM",name[i],"circle_Plot",".svg"),width=7,height=7)
  netVisual_chord_gene(cellchat, sources.use = as.numeric(key[i]), targets.use = c(7:10), slot.name = "netP", legend.pos.x = 10)
  dev.off()
}
svglite(paste0("plot_5/CellChat/","Tumor_TO_TAM","all","circle_Plot",".svg"),width=12,height=12)
netVisual_chord_gene(cellchat, sources.use = c(1:6), targets.use = c(7:10), slot.name = "netP", legend.pos.x = 10)
dev.off()
#vilin plot
plotGeneExpression(cellchat, signaling = cellchat@netP$pathways)
ggsave(paste0("plot_5/CellChat/","all","violin_Plot",".pdf"),width=10,height=20,dpi=600)
#Spatial Plot of LR enrichment
load("data/st.RData")
gene<-c("CXCL1","CXCL2","CXCL3","CXCL8","CXCL12","ACKR1","CXCR4")
point_size<-c(1.6,2,2,1.6,2,2,1.6,1.6,2)
image<-c("image","image.1","image.2","image.3","image.4","image.5","image.6","image.7","slice1")
for(i in 1:9){
  SpatialPlot(st,features=gene, images = image[i], label = FALSE, label.size = 1,ncol=3,pt.size.factor =point_size[i],image.alpha=0,cols=) 
  ggsave(file=paste0("plot_5/CellChat/","all","CXCL_spatial_",i,".pdf"),width=15,height=15)
}
FeaturePlot(object = st,features = gene)
ggsave(file=paste0("plot_5/CellChat/","all","feature_map","CXCL",".pdf"),width=15,height=15)
gene<-c("CCL2","CCL5","ACKR1","GAS6","AXL","TNF","TNFRSF1B","VEGFA","VEGFR1","ANXA1","FPR2","FPR3")
for(i in 1:9){
  SpatialPlot(st,features=gene, images = image[i], label = FALSE, label.size = 1,ncol=3,pt.size.factor =point_size[i],image.alpha=0,cols=) 
  ggsave(file=paste0("plot_5/CellChat/","all","Other_spatial_",i,".pdf"),width=15,height=20)
}
FeaturePlot(object = st,features = gene)
ggsave(file=paste0("plot_5/CellChat/","all","feature_map","other",".pdf"),width=15,height=15)
#Identify signaling rolesof cell groups
# Compute the network centrality scores
pathways.show<- c("CXCL", "CCL","GAS","TNF","VEGF","ANNEXIN")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",signaling = c("CXCL", "CCL","GAS","TNF","VEGF","ANNEXIN"))

