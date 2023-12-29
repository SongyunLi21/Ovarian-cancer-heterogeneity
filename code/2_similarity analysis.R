load("data/st.RData")
df<-table(st$identstate)
df<-as.data.frame(df)
name<-as.character(df$Var1[df$Freq>100])
Idents(st)<-"identstate"
st_sub<-subset(st,ident=name)
Idents(st_sub)<-"class"
#only perform the analysis on tumor and stromal region
st_sub<-subset(st_sub,ident= c("stromal","tumor"))
Idents(st_sub)<-"identstate"
exp_cluster<-as.data.frame(rownames(st_sub@assays$Spatial@scale.data))
word<-levels(st_sub@active.ident)
for(i in 1:length(word)){
  table<-st_sub@meta.data[st_sub@meta.data$identstate==word[[i]],]
  exp<-st_sub@assays$Spatial@scale.data[,rownames(table)]
  exp_cluster[,i]<-rowMeans(exp)
}
rownames(exp_cluster)<-rownames(st_sub@assays$Spatial@scale.data)
colnames(exp_cluster)<-levels(st_sub@active.ident)

# PCA
library(factoextra)
pca <- princomp(exp_cluster, cor = TRUE)
load<-as.data.frame(c(1:length(word)))
for(i in 1:5){
  load[,i]<-pca$loadings[,i]
}
rownames(load)<-levels(st_sub@active.ident)
colnames(load)<-c("PC1","PC2","PC3","PC4","PC5")
stromal<-sort(unique(c(brewer.pal(name="Greens", n = 9),brewer.pal(name="BuGn", n = 9),brewer.pal(name="YlGn", n = 9))))
tumor<-sort(unique(c(brewer.pal(name="OrRd", n = 9),brewer.pal(name="Reds", n = 9),brewer.pal(name="Oranges", n = 9))))
#plot dimensional reduction plot using PC1, PC2
ggplot(load ,aes(x=PC1,y=PC2,col=rownames(load)),label=F)+theme_classic()+scale_color_manual(values=c(stromal[c(1:21)],tumor[c(7:21)]))+ 
  geom_point(shape=16,size=4)+geom_text(label=rownames(load),size=2,color="black",check_overlap = F)
ggsave(file=paste0("plot_1/","all_patient_","dim_rd_label",".pdf"),width=11,height=10)
