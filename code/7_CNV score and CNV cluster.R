####Normalize CNV Score and plot i in spatial####
#load CNV score
cnv_table <- read.table("data/pbmc/infercnv.observations_4.txt", header=T)
cnv_score_table <- as.matrix(cnv_table)
#normalization
cnv_score_table <-abs(cnv_table-1)
cell_scores_CNV <- as.data.frame(colSums(cnv_score_table))
colnames(cell_scores_CNV) <- "cnv_score"
head(cell_scores_CNV)
cell_scores_CNV$cancer<-cell_scores_CNV$cnv_score/max(cell_scores_CNV$cnv_score)
#write.csv(x = cell_scores_CNV, file = "cnv_score_pbmc.csv")
#load position information
library(readr)  
# Load the expression data
expr.url <- paste0("data/st1/", 4,"/spatial/tissue_positions_list.csv")
position_file <- read_csv(file =expr.url,col_names = F)
positions <- position_file[,c("X3","X4")]
colnames(positions)<-c("X","Y")
bar<-position_file$X1
#for(i in 1:length(bar)){
#  bar[i] <- paste0(bar[i],"_",id)
#}
positions <- as.data.frame(positions)
rownames(positions)<-bar
#path<-paste0("data/CNV/cnv_score_",id,".xlsx")
#cnv_score_1 <- read_excel(path)
#b<-cnv_score_1$barcode
#cnv_score_1$barcode<-NULL
#rownames(cnv_score_1)<-b
#link CNV information to position information
library(do)
rownames(cell_scores_CNV)<-Replace(rownames(cell_scores_CNV),from=".1",to="-1")
c<-rownames(positions)%in%rownames(cell_scores_CNV)
d<-rownames(positions)
for(i in 1:4992){
  if(c[i]==TRUE){
    positions[d[i],"cnv_scores"]  =  cell_scores_CNV[d[i],"cancer"]###4
  }
  else{
    positions[i,"cnv_scores"]  =  0.005
  }
}
for(i in 1:4992){
  if(is.na(positions[i,"cnv_scores"])){
    positions[i,"cnv_scores"]<-0.005
  }
}
#position_list[[id]] <- positions

#positions<-position_list[[k]]
library(ggplot2)
p=ggplot(positions, aes(x = X, y = Y, col = cnv_scores)) + geom_point()+
  scale_colour_gradientn(colors = pals::kovesi.linear_blue_5_95_c73(20),limits = c(0,1))
p

####Spatial Plot of CNV cluster####
#spatial plot of cnv cluster
st_list <- readRDS("~/spatial/data/st_list_h.rds")
library(readr)
cell_grouping <- read_csv("data/pbmc/cell_grouping_4.csv")
a<-cell_grouping$cell
cell_grouping$cell<-NULL
rownames(cell_grouping)<-a
pbmc<-st_list[[4]]
pbmc<-AddMetaData(pbmc,cell_grouping,col="cell_group")
Idents(pbmc) <- 'cell_group'
SpatialPlot(pbmc,pt.size.factor =1.5,image.alpha = 0,cols="Paired")
ggsave(file=paste0("plot_3/","spatial_8_","all_plot.pdf"),height=15,width=15,dpi=600)
SpatialDimPlot(pbmc,pt.size.factor =1.45,
               image.alpha = 0,
               ncol=4,
               cells.highlight = CellsByIdentities(object = pbmc, idents = levels(pbmc)),
               cols.highlight = c("gold","black"),facet.highlight = TRUE)
ggsave(file=paste0("plot_3/","spatial_4_","_plot.pdf"),height=10,width=20,dpi=600)
marker_cnv<-FindAllMarkers(pbmc,assay = "SCT",slot="counts")

#identify the marker of CNV cluster
gene<-FindMarkers(pbmc,assay = "SCT",slot="counts",ident.1 ="M",ident.2 = "L")
SpatialFeaturePlot(object=pbmc,image.alpha = 0, pt.size.factor = 1.4,ncol = 3,
                   features = c("IGFBP4", "GRB7", "EZH1", "ERBB2", "EIF1","CAVIN1","ATP6V0A1","AOC3","ITFM","IRF7"))
ggsave(file=paste0("plot_3/","spatial_feature_4_","_plot.pdf"),height=14,width=10,dpi=600)
k=1
for(i in levels(pbmc)){
  gene<-cnv_genes[cnv_genes$cell_group_name==i,]
  marker<-marker_cnv[marker_cnv$cluster==i,]
  name<-rownames(marker)[rownames(marker)%in%rownames(gene)]
  marker<-marker[name,]
  gene<-gene[name,]
  if(k==1){
    df_1<-marker
    df_2<-gene
    
  }else{
    df_1<-rbind(df_1,marker)
    df_2<-rbind(df_2,gene)
  }
  k=0
}
write.csv(df_2,"data/cnv_chromosome_genes_8.csv")
marker_cnv<-marker_cnv[marker_cnv$p_val_adj<0.05,]#
features<- marker_cnv %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
features <- features$gene
DotPlot(pbmc, features = unique(features),cols = pals::kovesi.diverging_isoluminant_cjo_70_c25(2)) + RotatedAxis()
ggsave(file=paste0("plot_3/","cnv_8_","dot_plot.pdf"),height=7,width=7,dpi=600)
write.csv(marker_cnv,file="data/marker_8.csv")
top10markers<-marker_cnv %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
library(RColorBrewer)
DoHeatmap(pbmc, features = geneset$gene_symbol, size = 3,group.colors = brewer.pal(8,"Pastel2"))+
  scale_fill_gradientn(colors =pals::kovesi.diverging_cwm_80_100_c22(20)) 
ggsave(file=paste0("plot_3/","heatmap_8_","heat_plot.pdf"),height=20,width=20,dpi=600)

#gsva for CNV cluster
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
exp <-AverageExpression(pbmc ,
                        group.by = "cell_group",
                        assays = "SCT",
                        slot = "counts")
exp<-exp$SCT
cg=names(tail(sort(apply(exp,1, sd)),50))
pheatmap::pheatmap(cor(exp[cg,]),filename = paste0("plot_3/","GSVA_4_cor",".pdf"))

#GSVA
es.max <- gsva(exp, geneset)
#筛选+画图
library(RColorBrewer)
cg = names(tail(sort(apply(es.max, 1, sd) ),50))
pheatmap(es.max[cg,],show_rownames = T,fontsize_row = 7,fontsize_col=7,color = colorRampPalette(rev(brewer.pal(n = 7, name ="BrBG")))(100),
         filename = paste0("plot_3/","GSVA_4",".pdf"),width = 6,height=7)
