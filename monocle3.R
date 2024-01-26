library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(readr)
library(tidyr)
library(ggforce)
library(pals)
library(pheatmap)
library(scales)
library(ggthemes)
library(clustree)
library(ComplexHeatmap)
library(circlize)
library(ggrastr)
library(ggrepel)
library(purrr)
library(reshape2)
library(scatterpie)

setwd("/mnt/public3/chaohy/pancancer/fig3/")

rds <- readRDS("/mnt/public3/chaohy/pancancer/fig3/Ductal0115_umap.rds")
rds

library(monocle3)

##创建CDS对象并预处理数据
data <- GetAssayData(rds, assay = 'RNA', slot = 'counts')
cell_metadata <- rds@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="major_cluster") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(rds, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="major_cluster") + ggtitle('int.umap')
p2

## Monocle3聚类分区
cds <- cluster_cells(cds)
plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")


## 识别轨迹
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) -> p

##细胞按拟时排序
# # 使用辅助线选择root细胞
# p + geom_vline(xintercept = seq(-9.5,-9,0.25)) + geom_hline(yintercept = seq(0,0.25,0.25))
# embed <- data.frame(Embeddings(scRNAsub, reduction = "umap"))
# embed <- subset(embed, UMAP_1 > -8 & UMAP_1 < -7.5 & UMAP_2 > 0.24 & UMAP_2 < 0.25)
# root.cell <- rownames(embed)
# cds <- order_cells(cds, root_cells = root.cell)

# 使用在线工具选择
cds <- order_cells(cds)

# pdf("/mnt/public3/chaohy/pancancer/fig3/pseudotime.pdf",width=6,height=5)
png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotime.png",width=6,height=5, units = "cm",bg = "white",res = 300)
plot_cells(cds, color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, 
           label_leaves = FALSE,  
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           show_trajectory_graph = FALSE,
           rasterize=T) + NoAxes() + NoLegend()
dev.off()


png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotimebyCourse.png",width=13,height=10, units = "cm",bg = "white",res = 300)
plot_cells(cds, color_cells_by = "Course", 
           label_cell_groups = FALSE, 
           label_leaves = FALSE,  
           label_branch_points = FALSE, 
           label_roots = FALSE, 
           show_trajectory_graph = TRUE,
           rasterize=T)
dev.off()


plot_cells(cds, color_cells_by = "major_cluster", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,rasterize=TRUE)

plot_cells(cds, color_cells_by = "Course", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,rasterize=TRUE)

png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotimebyGeneExp.png",width=18,height=9, units = "cm",bg = "white",res = 300)
plot_genes_in_pseudotime(cds[c("HSALNG0050644","HSALNG0136338",
                               "HSALNG0051132","HSALNG0093627",
                               "HSALNG0103881","HSALNG0092235"),], 
                         color_cells_by="Course",
                         min_expr=1, nrow = 3,ncol=3) 
dev.off()

png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotimebyGeneExp_low2high.png",width=9,height=4.5, units = "cm",bg = "transparent",res = 300)
plot_genes_in_pseudotime(cds[c("FXYD2",),], 
                         color_cells_by="pseudotime",
                         min_expr=0.5,nrow=3) 
dev.off()

png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotimebyGeneExp_high2low.png",width=18,height=9, units = "cm",bg = "transparent",res = 300)
plot_genes_in_pseudotime(cds[c("TFF1","LINC01133","HSALNG0092235"),], 
                         color_cells_by="pseudotime",
                         min_expr=0.5,nrow=3)  -> p1
dev.off()


overlap_normal <- read.table("./overlap_normal.txt",header = F,sep = "\t")
read.table("/mnt/public3/chaohy/pancancer/fig2/gene.bed", sep="\t", header=F) -> gene_bed
gene_bed[gene_bed$V5=="lncRNA",] -> gene_bed_lncRNA
Track_genes_overlap_normal <- Track_genes[Track_genes$gene_short_name %in% overlap_normal$V1, ] %>% arrange(desc(morans_I))
Track_genes_overlap_normal_lnc <- Track_genes_overlap_normal[Track_genes_overlap_normal$gene_short_name %in% gene_bed_lncRNA$V4, ] %>% arrange(desc(morans_I))




plot_genes_in_pseudotime(cds[c("TFF1","HSALNG0051132","HSALNG0092235","FXYD2","HSALNG0111708","HSALNG0052479"),], 
                         color_cells_by="pseudotime",
                         min_expr=0.5,nrow=3,ncol=2)  -> p2

png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotimebyGeneExp_ALL.png",width=15,height=10, units = "cm",bg = "transparent",res = 300)
p2
dev.off()


##寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=50)
Track_genes %>% dplyr::filter(p_value < 0.001 & morans_I > 0.2) -> Track_genes_filter
write.table(Track_genes_filter, "/mnt/public3/chaohy/pancancer/fig3/Track_genes_by_monolce3.txt",col.names = T,sep = "\t", quote = F, row.names = F)

saveRDS(Track_genes,"./Track_genes.rds")

# 读取overlap的基因

overlap <- read.table("/mnt/public3/chaohy/pancancer/fig3/overlap.txt",header = F,sep="\t")

Track_genes_overlap <-  Track_genes[Track_genes$gene_short_name %in% overlap$V1, ] %>% arrange(morans_I)

read.table("/mnt/public3/chaohy/pancancer/fig2/gene.bed", sep="\t", header=F) -> gene_bed
gene_bed[gene_bed$V5=="lncRNA",] -> gene_bed_lncRNA

Track_genes_overlap_lnc <- Track_genes_overlap[Track_genes_overlap$gene_short_name %in% gene_bed_lncRNA$V4, ] %>% arrange(desc(morans_I))

Track_genes_lnc <- Track_genes[Track_genes$gene_short_name %in% gene_bed_lncRNA$V4, ] %>% arrange(desc(morans_I))
  

png(filename = "/mnt/public3/chaohy/pancancer/fig3/pseudotimebyGeneExp11Lnc.png",width=18,height=18, units = "cm",bg = "white",res = 300)
plot_genes_in_pseudotime(cds[Track_genes_overlap_lnc$gene_short_name,], 
                         color_cells_by="Course",
                         min_expr=0.5,nrow=6,ncol=2) 
dev.off()



#挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图
plot_genes_in_pseudotime(cds[c()], color_cells_by="Course", min_expr=0.5)


ggsave("/mnt/public3/chaohy/pancancer/fig3/fig3e.pdf",p,width =7, height =4)

plot_genes_violin(cds["HSALNG0051132"], group_cells_by = "Course") + theme(axis.text.x = element_text(angle = 45,hjust = 1))

#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

##寻找共表达模块
genelist <- Track_genes %>% filter(q_value < 0.001 & morans_I > 0.3 ) %>% pull(gene_short_name) %>% as.character() 

gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)

cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)

agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")


plot_cell_trajectory(cds,color_by = "Course")


# 教程：https://www.jianshu.com/p/c402b6588e17




Track_genes %>% filter(q_value < 0.001 & morans_I > 0.1) %>% grepl("^HSALNG", gene) %>% sum()

Track_genes %>% filter(q_value < 0.001 & morans_I > 0.1 & grepl("^HSALNG", gene_short_name)) %>% dim()