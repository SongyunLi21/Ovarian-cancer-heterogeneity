load("data/st_h.RData")
annotation <- data.frame(st@active.ident)
colnames(annotation) <- NULL
write.table(annotation, "data/annotation_file.txt",sep = "\t")
#1
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=pbmc@assays$Spatial@data,
                                    annotations_file="data/annotation_file.txt",
                                    delim="\t",
                                    gene_order_file="data/gene_ordering_file.txt",
                                    ref_group_names=NULL)

#perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,#use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",# dir is auto-created for storing outputs
                             cluster_by_groups=F,,
                             #tumor_subcluster_partition_method = "qnorm",
                             denoise=T,
                             HMM=T,
                             num_threads=1,
                             output_format="pdf")