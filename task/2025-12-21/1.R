library(Seurat)
library(ggplot2)
library(dplyr)
library(scutilsR)
library(tidyverse)
library(celda)
library(sceasy)
library(enrichR)
library(Nebulosa)
library(magrittr)
library(SeuratWrappers)
library(glue)
library(data.table)
source("function/read_loom.R")
source("function/preprocess.R")

#### 导入数据 ####
samples <- list.dirs("data/filter/", full.names = F, recursive = F)
samples <- samples[grepl("^A", samples)]
seu.list <- pbapply::pblapply(samples, function(sn) {
  counts <- Read10X(file.path("data/filter/", sn))
  sn <- gsub("_", "-", sn) # 注意"_"在`CreateSeuratObject()`里有特殊的意义
  colnames(counts) <- paste(sn, colnames(counts), sep = "_")
  seu <- CreateSeuratObject(counts = counts)
  return(seu)
})
## 合并样本
seu <- base::Reduce(f = merge, x = seu.list)

##QC
scobj <- seu
rm(seu)
rm(seu.list)

##QC
##线粒体
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")

## 其他文章给的FFPE阈值
scobj$QC <- "low"
scobj$QC[scobj$nFeature_RNA > 500 & scobj$nFeature_RNA < 5000 & scobj$percent.mt < 50 & scobj$nCount_RNA > 400 & scobj$nCount_RNA < 25000] <- "high"
table(scobj$QC)

## 细胞周期
scobj <- NormalizeData(scobj)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
scobj <- CellCycleScoring(scobj, s.features = s.genes, g2m.features = g2m.genes)

## 双细胞
library(scDblFinder)
counts <- scobj@assays$RNA@counts
sce <- SingleCellExperiment(list(counts = counts))

sce <- scDblFinder(sce)

scobj$scDblFinder.class <- sce$scDblFinder.class
scobj$scDblFinder.score <- sce$scDblFinder.score

rm(sce)

## QC score
QC_genes <- c("ACADL", "AP2S1", "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC3", 
              "ATP5PF", "ATP5MF", "ATP5ME", "ATP5MG", "ATP6V1F", 
              "CHCHD10", "COA3", "COX5B", "COX6A1", "COX6B1", 
              "COX6C", "COX7A2", "COX7A2L", "COX7B", "CYCS", 
              "EDF1", "EEF1B2", "EIF5A", "FAU", "FKBP3", "FTL", 
              "GUK1", "HEPH", "HRAS", "MIF", "MRAP", "NACA", 
              "NDUFA1", "NDUFA2", "NDUFA4", "NDUFA5", "NDUFB7", 
              "NDUFC1", "NDUFS7", "NDUFV3", "NECAP1", "NLRP4", 
              "PDXP", "PFN2", "POLR2M", "RAB3A", "RTL8A", 
              "SLC16A2", "SNRPD2", "SNU13", "TAF1C", "TIMM8B", 
              "TPT1", "UBB", "UQCR11", "UQCRB", "UQCRQ", "USP50")
expr_matrix <- GetAssayData(scobj, assay = "RNA", slot = "data")
log_expr_values <- log1p(expr_matrix[QC_genes, ])
qc_scores <- colSums(log_expr_values, na.rm = TRUE)
scobj$qc_score <- qc_scores

## 应激
hot_shock <- c("FOS", "CXCL2", "ZFP36", "FOSB", "DUSP1", "ATF3", "CXCL8", 
               "NR4A1", "CXCL3", "PPP1R15A", "JUNB", "EGR1", "HSPA1A", "HSPA1B", 
               "SOCS3", "KLF6", "JUN", "IER2", "CXCL1", "NFKBIA", "HSPA6", "DNAJB1", 
               "IER3", "CCNL1", "MTRNR2L2", "IER5", "ID1", "CEBPD", "KRT6A", 
               "CYR61", "DEPP1", "CLDN4", "IRF1", "DUSP2", "BTG2", "PLAUR", 
               "MAFF", "KLF4", "PHLDA2", "TNFAIP3", "ACTG1", "BTG1", "DNAJB4", 
               "ERRFI1", "H3F3B", "HSPB1", "PCF11", "PXDC1", "SDC4", "SRF", 
               "TPM3", "USP2", "GADD45G", "ANKRD1", "FAM132B", "HIPK3", "HSPH1", 
               "IRF8", "KLF9", "NFKBIZ", "PDE4B", "RAP1B", "SERPINE1", "TPPP3", 
               "WAC", "HSPE1", "ARID5A", "DCN", "DUSP8", "HSP90AA1", "ID3", 
               "ITPKC", "LITAF", "NOP58", "PER1", "RASSF1", "SKIL", "SRSF7", 
               "TRA2A", "ZC3H12A", "CCRN4L", "DDX3X", "HSP90AB1", "IDI1", "LMNA", 
               "MYADM", "NPPC", "PHLDA1", "RHOB", "SLC10A6", "STAT3", "TRA2B", 
               "ZFAND5", "KCNE4", "ATF4", "CEBPB", "DDX5", "EGR2", "FOSL2", 
               "MYC", "PNP", "RHOH", "SLC38A2", "TAGLN2", "TRIB1", "BAG3", "DES", 
               "GADD45A", "JUND", "MAFK", "MYD88", "ODC1", "PNRC1", "RIPK1", 
               "SLC41A1", "TIPARP", "TUBB4B", "ZFP36L1", "BHLHE40", "CEBPG", 
               "DNAJA1", "EIF5", "GCC1", "HSPA5", "IFRD1", "KLF2", "MCL1", "NCKAP5L", 
               "OSGIN1", "SAT1", "TUBB6", "ZFP36L2", "BRD2", "CSRNP1", "ERF", 
               "GEM", "HSPA8", "IL6", "MIDN", "NCOA7", "OXNAD1", "SBNO2", "SQSTM1", 
               "TNFAIP6", "UBC", "ZYX", "MIR22HG", "MT1A", "SRSF5", "MT2A", 
               "EIF1", "PPP1CC", "ACTB", "ADAMTS1", "ADAMTS9", "AHNAK", "ANKRD11", 
               "ARF4", "AZIN1", "BAIAP2", "BAZ1A", "CAMK1D", "CCDC138", "CDKN1A", 
               "CHD4", "CHKA", "CLIC4", "CMSS1", "COL1A1", "CTNNB1", "CX3CR1", 
               "ELF2", "EP400", "ERN1", "ETF1", "FBXL18", "FLT1", "GADD45B", 
               "GLS", "GNAS", "GSK3A", "GSN", "HIVEP2", "INTS6", "JAK1", "JDP2", 
               "KDM6B", "KPNA1", "LSMEM1", "LUZP1", "MAGI3", "MAN1A1", "MAPKAPK2", 
               "MAPRE1", "MED13", "MSN", "MYLIP", "NABP1", "NASP", "NUFIP2", 
               "NUP210L", "PEAK1", "PECAM1", "POLG2", "PPP1CB", "PRKCG", "RNF19B", 
               "RTN4", "SERTAD2", "SGPL1", "SIK3", "SPAG9", "TAF4B", "TEX14", 
               "TOB2", "TOP1", "DIAPH1", "NEAT1", "PTMA", "ARIH1")
hot_shock_genes <- intersect(hot_shock, rownames(scobj))
scobj[["percent.stress"]] <- PercentageFeatureSet(scobj, features = hot_shock_genes)

## RNA污染
QuickCluster <- function(object) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object)
  return(object)
}
seu.list <- SplitObject(scobj, split.by = "orig.ident")
seu.list <- lapply(seu.list, QuickCluster)
clusters <- lapply(seu.list, function(xx) xx$seurat_clusters) %>% base::Reduce(c, .)
scobj$quick_clusters <- clusters[rownames(scobj@meta.data)]
rm(seu.list)
gc()
scobj <- scutilsR::RemoveAmbientRNAs(scobj, split.by = "orig.ident", cluster.name = "quick_clusters")
scobj[["decontX"]] <- NULL 

##红细胞基因
erythrocyte_genes <- c(
  "GATA1",    # 红细胞转录因子
  "EPOR",     # 红细胞生成素受体
  "HBB",      # 血红蛋白β链
  "HBA1",     # 血红蛋白α链1
  "HBA2",     # 血红蛋白α链2
  "KLF1",     # 红细胞发育关键转录因子
  "SLC4A1",   # 红细胞膜蛋白
  "EKLF",     # 红细胞特异性转录因子
  "ALAS2",    # 红细胞δ氨基酮戊酸合成酶
  "TAL1"      # T细胞白血病/淋巴瘤蛋白1
)
erythrocyte_genes <- intersect(erythrocyte_genes, rownames(scobj))
scobj[["percent.ery"]] <- PercentageFeatureSet(scobj, features = erythrocyte_genes)

##计算核糖体基因
ribo.genes <- ProjectSVR::ribo.genes
ribo.genes <- intersect(ribo.genes, rownames(scobj))
scobj[["percent.ribo"]] <- PercentageFeatureSet(scobj, features = ribo.genes)

##计算内含子占比
## 读取loom文件
samples <- list.files("data/loom/")
ldat <- pbapply::pblapply(samples, function(fn) {
  message(glue::glue("Loading {fn} ..."))
  read.loom.matrices(file.path("data/loom/", fn))
})

## 重命名
names(ldat) <- sub(".loom", "", gsub("_", "-", samples))
sapply(ldat[[1]], dim)

## 合并loom文件
matrix.name <- names(ldat[[1]])
matrix.name

ldat.merged <- lapply(matrix.name, function(mn){
  mat.list <- lapply(ldat, function(xx) xx[[mn]])
  do.call(cbind, mat.list) ## merge by columns (cell.ID)
})
names(ldat.merged) <- matrix.name
sapply(ldat.merged, dim)

## Note: the cell.IDs in loom are inconsistent with those in `Seurat` object.
head(colnames(ldat.merged$spliced))
head(colnames(scobj))

## fix the cell.IDs in loom
## 不一定一样，目的是把ldat.merged（继承FASTQ名字）改成seu的
fix_id <- function(x) {
  # 提取中间部分和后面的部分
  main_part <- sub(".*-(A[0-9]+).*", "\\1", x)  # 提取以 A 开头的部分（中间 ID）
  suffix_part <- sub(".*:(.*?)x$", "\\1", x)    # 提取末尾部分，去掉最后的 -1
  
  # 合并并返回结果
  return(paste0(main_part, "_", suffix_part))
}

fix_id(colnames(ldat.merged$spliced)) %>% head()

for (i in seq_along(ldat.merged)) {
  colnames(ldat.merged[[i]]) %<>% fix_id()
}

## check the fixed cellID，输出为T
colnames(scobj) %in% colnames(ldat.merged$spliced) %>% all() 

## filter cells，只保留seu有的细胞
for (i in seq_along(ldat.merged)) {
  ldat.merged[[i]] <- ldat.merged[[i]][, colnames(scobj)]
}
sapply(ldat.merged, dim)

## 计算内含子占比
emat <- ldat.merged$spliced   # exonic read (spliced) expression matrix
nmat <- ldat.merged$unspliced # intronic read (unspliced) expression matrix
percent.intron <- colSums(nmat) / (colSums(nmat) + colSums(emat))
scobj[["percent.intron"]] <- percent.intron[colnames(scobj)]

##确认有没有NA
table(is.na(scobj@meta.data$percent.intron))

## 下游分析
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, nfeatures = 2000)
scobj <- ScaleData(scobj)
scobj <- RunPCA(scobj)

DimPlot(scobj, reduction = "pca", group.by = "Phase",raster = F)

## 矫正细胞周期的影响
scobj<- ScaleData(scobj, vars.to.regress = c("S.Score", "G2M.Score"))
scobj <- RunPCA(scobj)

DimPlot(scobj, reduction = "pca", group.by = "Phase",raster = F)
DimPlot(scobj, reduction = "pca", group.by = "orig.ident", label = T,raster = F)

## 检查批次效应
scobj <- RunUMAP(scobj, reduction = "pca", dims = 1:30)
DimPlot(scobj, reduction = "umap", group.by = "orig.ident",raster = F)
DimPlot(scobj, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident",raster = F) + ggsci::scale_color_d3()

## 去除批次效应
library(harmony)
scobj <- RunHarmony(scobj, group.by.vars = c("orig.ident"), reduction.use = "pca", dims.use = 1:50)

scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30)
DimPlot(scobj, reduction = "umap", group.by = "orig.ident",raster = F) + ggsci::scale_color_d3()
DimPlot(scobj, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident",raster = F) + ggsci::scale_color_d3("category20")

## 聚类
resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1)
scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:30, k.param = 20)
scobj <- FindClusters(scobj, resolution = resolutions)

DimPlot(scobj, reduction = "umap", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T,raster = F) & NoLegend()
Idents(scobj) <- "RNA_snn_res.1"
DimPlot(scobj,label = T,reduction = "umap",raster = F)

## QC
DimPlot(scobj,group.by = "QC",reduction = "umap",raster = F)
DimPlot(scobj,group.by = "scDblFinder.class",reduction = "umap",raster = F)
FeaturePlot(scobj, features = "percent.stress",reduction = "umap",raster = F)
FeaturePlot(scobj,features = "decontX_contamination",reduction = "umap",raster = F)
FeaturePlot(scobj,features = "qc_score",reduction = "umap",raster = F)
FeaturePlot(scobj,features = "nFeature_RNA",reduction = "umap",raster = F)
FeaturePlot(scobj,features = "percent.ery",reduction = "umap",raster = F)
FeaturePlot(scobj,features = "percent.ribo",reduction = "umap",raster = F)
FeaturePlot(scobj,features = "percent.mt",reduction = "umap",raster = F)

############################################################################################## QC矩阵，遍历每个样本（数据集）
##基于样本QC
##饼图展示
library(ggplot2)
library(dplyr)
library(ggsci)
# 计算数据
count_data <- table(scobj@meta.data$orig.ident)
df <- data.frame(
  group = names(count_data),
  value = as.numeric(count_data)
)

# 计算百分比
df <- df %>%
  mutate(
    percentage = round(value / sum(value) * 100, 1),
    label = paste0(group, " (", percentage, "%)")
  )

ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_nejm() +  # 使用Nature杂志配色
  geom_text(
    aes(label = label), 
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold"
  ) +
  labs(
    title = "Sample Distribution",
    fill = "Sample"
  )+ ggsci::scale_fill_d3("category20") 

## 构建QC矩阵
table(scobj@meta.data$orig.ident)
sample = "A3291"
seu1=subset(scobj, subset = orig.ident==sample)

metadata=seu1@meta.data
qc.names = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.intron",
             "percent.ery", "scDblFinder.score",
             "decontX_contamination", "qc_score","percent.stress")

qc.mat = metadata[, qc.names, drop = FALSE]

## MinMax缩放
qc.mat <- as.data.frame(lapply(qc.mat, as.numeric))
rownames(qc.mat) <- rownames(metadata)
qc.minmax <- as.data.frame(
  apply(qc.mat, 2, function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  })
)
summary(qc.minmax)

## 把矩阵添加到Seurat
new_embeddings <- as.matrix(qc.minmax)
colnames(new_embeddings) <- paste0("QC_", 1:ncol(new_embeddings))

seu1[["qc"]] <- CreateDimReducObject(
  embeddings = new_embeddings,
  key = "qc_", 
  assay = DefaultAssay(seu1),
  
)

##可视化
seu1 <- RunUMAP(seu1, reduction = "qc", dims = 1:9,reduction.name = "umap_qc")
seu1 <- FindNeighbors(seu1, reduction = "qc", dims = 1:9, k.param = 20)

## 最好保持一个resolution，要么就删了这列，不然后面merge不起来
seu1 <- FindClusters(seu1, resolution = 0.6,cluster.name = "qc_cluster")
VlnPlot(seu1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.intron",
                           "percent.ery",  "scDblFinder.score",
                           "decontX_contamination", "qc_score","percent.stress"))
FeaturePlot(seu1, features =  c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.intron",
                                "percent.ery",  "scDblFinder.score",
                                "decontX_contamination", "qc_score","percent.stress"), order = TRUE,ncol=3,reduction="umap_qc")
DimPlot(seu1, reduction = "umap_qc",label = T) +ggsci::scale_color_d3("category20", na.value = "grey80")

## label
Idents(seu1)
for(i in 1:nrow(seu1@meta.data)) {
  if(seu1@meta.data$seurat_clusters[i] == 0) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }  
  if(seu1@meta.data$seurat_clusters[i] == 1) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
  if(seu1@meta.data$seurat_clusters[i] == 2) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
  if(seu1@meta.data$seurat_clusters[i] == 3) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
  if(seu1@meta.data$seurat_clusters[i] == 4) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
  if(seu1@meta.data$seurat_clusters[i] == 5) {
    seu1@meta.data$QC_Matrix[i] <- "High"
  }
  if(seu1@meta.data$seurat_clusters[i] == 6) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
  if(seu1@meta.data$seurat_clusters[i] == 7) {
    seu1@meta.data$QC_Matrix[i] <- "High"
  }
  if(seu1@meta.data$seurat_clusters[i] == 8) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
  if(seu1@meta.data$seurat_clusters[i] == 9) {
    seu1@meta.data$QC_Matrix[i] <- "Low"
  }
}

DimPlot(seu1, reduction = "umap_qc",label = T,group.by = "QC_Matrix") + ggsci::scale_color_d3("category20", na.value = "grey80")
metadata=seu1@meta.data
write.csv(metadata,glue("E:/task/PNI/data_mining/output/metadata/{sample}.csv"))

##合并数据
# 初始化一个空的列表以存储读取的数据框
meta_list <- list()

# 获取目录下的所有CSV文件
csv_files <- list.files("E:/task/PNI/data_mining/output/metadata", pattern = "*.csv", full.names = TRUE)

# 循环读取每个CSV文件并加入到列表中
for (file in csv_files) {
  tmp <- fread(file, data.table = FALSE)  # 读取CSV文件
  rownames(tmp) <- tmp$V1                  # 设置行名
  tmp <- tmp[, -1]                         # 删除原始行名列
  meta_list[[file]] <- tmp                 # 将数据框添加到列表
}

# 将列表中的所有数据框合并为一个大的数据框，缺失的列用NA填充
meta_list <- lapply(meta_list, function(df) {
  df[,"QC_Matrix", drop = FALSE]  
})

qc <- bind_rows(meta_list)
metadata=scobj@meta.data
qc=qc[rownames(metadata),,drop = F]
metadata=cbind(metadata,qc)
scobj@meta.data=metadata
table(scobj@meta.data$QC_Matrix,useNA="ifany")

p1 = DimPlot(scobj, 
        label.size = 4,        # 调整标签大小
        repel = TRUE,          # 启用标签自动防重叠
        pt.size = 0.5,       # 调整点的大小
        group.by="QC_Matrix"
)  + ggsci::scale_color_d3("category20", na.value = "grey80")

df <- data.frame(
  group = c("High QC", "Low QC"),
  count = c(22422, 104266 - 22422)  # 81844
)

p2 <- ggplot(df, aes(x = group, y = count, fill = group)) +
  geom_col(width = 0.55, color = "#2d3436", alpha = 0.9) +
  geom_text(aes(label = format(count, big.mark = ",")),
            vjust = -0.4, size = 4, color = "#1f2933") +
  scale_fill_manual(values = c("High QC" = "#6c63ff", "Low QC" = "#00bcd4")) +
  labs(title = "QC", x = NULL, y = "Cell Count") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    axis.text = element_text(color = "#37474f"),
    axis.title = element_text(color = "#263238"),
    plot.title = element_text(color = "#1f2933", face = "bold", hjust = 0.5)
  )+ggsci::scale_fill_d3("category20") 

p1+p2

seu_filter <- subset(scobj, subset = QC_Matrix %in% "High")
seu_filter <- NormalizeData(seu_filter)
seu_filter <- FindVariableFeatures(seu_filter, nfeatures = 2000)
seu_filter <- ScaleData(seu_filter)
seu_filter <- RunPCA(seu_filter)

library(harmony)
seu_filter <- RunHarmony(seu_filter, group.by.vars = c("orig.ident"), reduction.use = "pca", dims.use = 1:50)

seu_filter <- RunUMAP(seu_filter, reduction = "harmony", dims = 1:30)
DimPlot(seu_filter, reduction = "umap", group.by = "orig.ident") + ggsci::scale_color_d3()

resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1)
seu_filter <- FindNeighbors(seu_filter, reduction = "pca", dims = 1:30, k.param = 20)
seu_filter <- FindClusters(seu_filter, resolution = resolutions)
DimPlot(seu_filter, reduction = "umap", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T) & NoLegend()

Idents(seu_filter) <- "RNA_snn_res.1"
DimPlot(seu_filter,label = T)

## T细胞 22 
FeaturePlot(seu_filter, features = c("CD3E", "CD8A","IL7R","CD4"),reduction = "umap", order = TRUE,ncol=2)

## B细胞
FeaturePlot(seu_filter, features = c("MS4A1", "CD79A"),reduction = "umap", order = TRUE,ncol=2)

## NK细胞
FeaturePlot(seu_filter, features = c("GNLY", "NKG7"),reduction = "umap", order = TRUE,ncol=2)

## 单核细胞
FeaturePlot(seu_filter, features = c("CD14", "FCGR3A", "LYZ","FCN1","S100A8"),reduction = "umap", order = TRUE,ncol=2)

## 巨噬细胞
FeaturePlot(seu_filter, features = c("CD68",  "CSF1R","C1QA"),reduction = "umap", order = TRUE,ncol=2)

## DC
FeaturePlot(seu_filter, features = c("CD1C", "CD1E", "FCER1A", "CLEC10A", "ITGAX", "CD80", "CD86"),reduction = "umap", order = TRUE,ncol=2)

##髓系 14
FeaturePlot(seu_filter, features = c("CSF1R", "ITGAX","CD80","CD86"),reduction = "umap", order = TRUE,ncol=2)

## 肥大细胞
FeaturePlot(seu_filter, features =c("MS4A2", "HDC", "CTSG", "KIT", "FCER1A", "TPSB2", "CMA1"),reduction = "umap", order = TRUE,ncol=2)

## 内皮细胞 6 
FeaturePlot(seu_filter, features = c("PECAM1", "CDH5", "VWF"),reduction = "umap", order = TRUE,ncol=2)

## 纤维细胞 8
FeaturePlot(seu_filter, features = c("COL3A1", "DCN", "LUM"),reduction = "umap", order = TRUE,ncol=2)

## 平滑肌细胞 16
FeaturePlot(seu_filter, features = c("ACTA2", "TAGLN", "CNN1"),reduction = "umap", order = TRUE,ncol=2)

## 肌成纤维细胞
FeaturePlot(seu_filter, features = c("ACTA2", "COL1A1", "COL3A1", "CCN2", "TAGLN", "VIM"),reduction = "umap", order = TRUE,ncol=2)

## 周细胞 23
FeaturePlot(seu_filter, features = c("PDGFRB","CSPG4","RGS5","ANPEP"),reduction = "umap", order = TRUE)

## 施万细胞 24
FeaturePlot(seu_filter, features = c("MPZ","MBP","CADM2","NCAM1"),reduction = "umap", order = TRUE,ncol=2)

## 基底细胞 13
FeaturePlot(seu_filter, features = c("KRT5", "KRT14", "TP63"),reduction = "umap", order = TRUE,ncol=2)

## 腔上皮细胞
FeaturePlot(seu_filter, features = c("KRT18", "MSMB",  "CD38", "KLK3", "ACP3", "NKX3-1"),reduction = "umap", order = TRUE,ncol=2)

## 肿瘤细胞
FeaturePlot(seu_filter, features = c("AMACR", "PCA3", "PCAT14", "TACSTD2", "FOLH1", "AR", "KLK3"),reduction = "umap", order = TRUE,ncol=2)

## club cell 18
FeaturePlot(seu_filter, features = c("LTF", "MMP7", "PIGR", "CP"),reduction = "umap", order = TRUE,ncol=2)

## Hillock cells
FeaturePlot(seu_filter, features =c("KRT13", "CLDN4", "LY6D"),reduction = "umap", order = TRUE,ncol=2)

DimPlot(seu_filter,reduction ="umap",cells.highlight = colnames(subset(seu_filter,subset= RNA_snn_res.1%in% "27")) )

Idents(seu_filter) <- "RNA_snn_res.1"
seu_filter <- FindSubCluster(
  seu_filter,
  cluster="18",
  graph.name = "RNA_snn",
  subcluster.name = "RNA_snn_res_S1",
  resolution = 0.5
)
DimPlot(seu_filter,label = T,group.by = "RNA_snn_res_S1")

Idents(seu_filter) <- "RNA_snn_res_S1"
seu_filter <- FindSubCluster(
  seu_filter,
  cluster="5",
  graph.name = "RNA_snn",
  subcluster.name = "RNA_snn_res_S2",
  resolution = 0.5
)
DimPlot(seu_filter,label = T,group.by = "RNA_snn_res_S2")

Idents(seu_filter) <- "RNA_snn_res_S2"
seu_filter <- FindSubCluster(
  seu_filter,
  cluster="20",
  graph.name = "RNA_snn",
  subcluster.name = "RNA_snn_res_S3",
  resolution = 0.5
)
DimPlot(seu_filter,label = T,group.by = "RNA_snn_res_S3")

Idents(seu_filter) <- "RNA_snn_res.1"
markers <- c("CD3E", "IL7R", "CSF1R", "ITGAX", "CDH5", "VWF", "DCN", "LUM", "ACTA2", "TAGLN", "CSPG4", "RGS5", "CADM2", "NCAM1", "KRT5", "TP63", "PCA3", "FOLH1", "LTF", "PIGR")
DotPlot(seu_filter, features = markers,dot.scale = 4) + coord_flip()+ RotatedAxis()+
  theme(axis.text.y = element_text(size = 6))
scCustomize::Clustered_DotPlot(seu_filter, features = markers)
seu_filter <- RenameIdents(seu_filter,
                      "0"="Luminal Cell",     
                      "1"="Luminal Cell",
                      "2"="Luminal Cell",
                      "3"="Luminal Cell",
                      "4"="Luminal Cell",
                      "5_0"="Luminal Cell",
                      "5_1"="Club Cell",
                      "5_2"="Club Cell",
                      "6"="Endothelial Cell",
                      "7"="Luminal Cell", 
                      "8"="Fibroblast", 
                      "9"="Luminal Cell",
                      "10"="Luminal Cell", 
                      "11"="Luminal Cell", 
                      "12"="Luminal Cell", 
                      "13"="Basal Cell", 
                      "14"="Myeloid Cell", 
                      "15"="Luminal Cell", 
                      "16"="Smooth Muscle Cell", 
                      "17"="Luminal Cell", 
                      "18_0"= "Basal Cell",
                      "18_1"= "Club Cell",
                      "18_2"= "Luminal Cell",
                      "18_3"= "Club Cell",
                      "18_4"= "Basal Cell",
                      "19"= "Luminal Cell", 
                      "20_0"= "Luminal Cell",
                      "20_1"= "Luminal Cell",
                      "20_2"= "Club Cell",
                      "21"= "Luminal Cell", 
                      "22"= "T Cell", 
                      "23"= "Pericyte", 
                      "24"= "Schwann Cell", 
                      "25"= "Luminal Cell", 
                      "26"= "Luminal Cell", 
                      "27"= "Luminal Cell"
)

DimPlot(seu_filter,label = T,reduction = "umap")
metadata <- seu_filter@meta.data
seu_filter@meta.data$celltype = Idents(seu_filter)
DimPlot(seu_filter,label = T,reduction = "umap",group.by = "celltype")+ggsci::scale_color_d3("category20", na.value = "grey80")

qs::qsave(seu_filter,"E:/task/PNI/文章1/output/tmp.qs")

## 上皮细胞
seu_ep = subset(seu_filter,subset = celltype == "Luminal Cell")
seu_ep <- NormalizeData(seu_ep)
seu_ep <- FindVariableFeatures(seu_ep, nfeatures = 2000)
seu_ep <- ScaleData(seu_ep)
seu_ep <- RunPCA(seu_ep)

## 矫正细胞周期的影响
seu_ep<- ScaleData(seu_ep, vars.to.regress = c("S.Score", "G2M.Score"))
seu_ep <- RunPCA(seu_ep)

## 检查批次效应
seu_ep <- RunUMAP(seu_ep, reduction = "pca", dims = 1:30)

## 去除批次效应
library(harmony)
seu_ep <- RunHarmony(seu_ep, group.by.vars = c("orig.ident"), reduction.use = "pca", dims.use = 1:50)

seu_ep <- RunUMAP(seu_ep, reduction = "harmony", dims = 1:30)
Idents(seu_ep) = "RNA_snn_res.1"
DimPlot(seu_ep, reduction = "umap", group.by = "orig.ident",raster = F) + ggsci::scale_color_d3()
DimPlot(seu_ep, reduction = "umap", raster = F,label = T) + ggsci::scale_color_d3("category20")
DimPlot(seu_ep, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident",raster = F) + ggsci::scale_color_d3("category20")

## 聚类
resolutions <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1)
seu_ep <- FindNeighbors(seu_ep, reduction = "pca", dims = 1:30, k.param = 20)
seu_ep <- FindClusters(seu_ep, resolution = resolutions)

DimPlot(seu_ep, reduction = "umap", group.by = paste0("RNA_snn_res.", resolutions), ncol = 3, label = T,raster = F) & NoLegend()
Idents(seu_ep) <- "RNA_snn_res.1"
DimPlot(seu_ep,label = T,reduction = "umap",raster = F)
DimPlot(seu_ep,reduction ="umap",cells.highlight = colnames(subset(seu_ep,subset= RNA_snn_res.1%in% "19")) )

seu_ep = subset(seu_ep,subset = RNA_snn_res.1 != "20")

## 细胞群的占比
library(dplyr)
library(ggplot2)
library(RColorBrewer)

metadata <- seu_ep@meta.data

df <- metadata %>%
  dplyr::count(RNA_snn_res.1, orig.ident) %>%
  group_by(RNA_snn_res.1) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

idents <- sort(unique(df$orig.ident))
pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(idents))

ggplot(df, aes(x = RNA_snn_res.1, y = prop, fill = orig.ident)) +
  geom_col(color = "grey20", width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = pal, name = "orig.ident") +
  labs(x = "RNA_snn_res.1", y = "Sample proportion",
       title = "Sample composition per cluster") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

seu_ep <- RenameIdents(seu_ep,
                           "0"="A1842 Specific Luminal Cell",     
                           "1"="A1842 Specific Luminal Cell",
                           "2"="A2516 Specific Luminal Cell",
                           "3"="A2518 Specific Luminal Cell",
                           "4"="A2516 Specific Luminal Cell",
                           "5"="A1842 Specific Luminal Cell",
                           "6"="A3289 Specific Luminal Cell",
                           "7"="A3290 Specific Luminal Cell", 
                           "8"="A3287 Specific Luminal Cell", 
                           "9"="A3289 Specific Luminal Cell",
                           "10"="A3291 Specific Luminal Cell", 
                           "11"="A3291 Specific Luminal Cell", 
                           "12"="A1842 Specific Luminal Cell", 
                           "13"="A3291 Specific Luminal Cell", 
                           "14"="A3289 Specific Luminal Cell", 
                           "15"="A2516 Specific Luminal Cell", 
                           "16"="A3291 Specific Luminal Cell", 
                           "17"="A3290 Specific Luminal Cell", 
                           "18"= "A3290 Specific Luminal Cell",
                           "19"= "Shared Luminal Cell"
)

DimPlot(seu_ep,label = T,reduction = "umap")
metadata <- seu_ep@meta.data
seu_ep@meta.data$celltype_lum = Idents(seu_ep)
DimPlot(
  seu_ep,
  reduction = "umap",
  group.by = "celltype_lum",
  label = TRUE,
  repel = TRUE,          # 关键：标签避让
  label.size = 3.5       # 可按需调整字体大小
) + ggsci::scale_color_d3("category20", na.value = "grey80")

## 差异分析
metadata=seu_ep@meta.data
metadata$state=ifelse(metadata$celltype_lum=="Shared Luminal Cell","1","0")
seu_ep@meta.data=metadata

Idents(seu_ep) <- "state"
diffgenes <- FindMarkers(seu_ep, 
                         ident.1 = 1, ## 结果正的在这里升高，负的下降
                         ident.2 = 0, 
                         logfc.threshold = 0)
markers <- diffgenes %>%
  filter(p_val_adj < 0.05)

## 火山图
library(ggrepel)
data <- markers
colnames(data) <- c("P.Value","logFC","pct.1","pct.2","adj.P.Val")
data$gene_id <- rownames(data)

##设置竖线和标注颜色
logFCfilter = 1
logFCcolor =1.5

### 标记上下调
index = data$adj.P.Val <0.05 & abs(data$logFC) > logFCfilter
data$group <- 0
data$group[index & data$logFC>0] = 1
data$group[index & data$logFC<0] = -1
data$group <- factor(data$group,levels = c(1,0,-1),labels =c("Up","NS","Down") )

ggplot(data=data, aes(x=logFC, y =-log10(adj.P.Val),color=group)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values = c("red", "grey50", "blue4"))+
  labs(x="log2 (fold change)",y="-log10 (adj.P.Val)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(-logFCfilter,logFCfilter),lty=4,lwd=0.6,alpha=0.8)+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.position="top")+
  geom_point(data=subset(data, abs(logFC) >= logFCcolor & adj.P.Val <0.05),alpha=0.8, size=3,col="green4")+
  geom_text_repel(data=subset(data, abs(logFC) >= logFCcolor & adj.P.Val <0.05),
                  aes(label=gene_id),col="black",alpha = 0.8)

## GSEA
library(clusterProfiler)
diffgenes$gene <- rownames(diffgenes)
gid <- bitr(unique(diffgenes$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
colnames(gid)[1] <- "gene"
diffgenes <- merge(diffgenes, gid, by='gene')

gene_df <- diffgenes
geneList <- gene_df$avg_log2FC
names(geneList) =  gene_df$gene
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
genesets <- read.gmt("resource/h.all.v2022.1.Hs.symbols.gmt")

y <- GSEA(geneList,TERM2GENE = genesets)
yd <- as.data.frame(y)

##炫酷一点
ggplot(y, showCategory = 30, aes(NES, forcats::fct_reorder(Description, NES))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  ##scale_color_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL)+
  ggtitle("Shared Luminal Cell vs other Luminal Cell")+
  theme(plot.title = element_text(hjust = 0.5))

library(GseaVis)

setid <- c("HALLMARK_ANDROGEN_RESPONSE")
gseaNb(object = y,
       geneSetID = setid,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.8,
       pvalY = 0.5,
       pFill = "white",
       addGene = T,
       markTopgene = T,
       geneCol = "#009933",
       topGeneN = 10)

qs::qsave(seu_filter,"E:/task/PNI/文章1/output/filter.qs")
qs::qsave(seu_ep,"E:/task/PNI/文章1/output/ep.qs")

## 细胞通讯
library(CellChat)
library(patchwork)
library(Seurat)

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

##根据实际情况选择互作模式
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

##构建cellchat对象
seu_SC = subset(seu_filter,subset = celltype == "Schwann Cell")
seu_combined <- merge(seu_ep, y = seu_SC)
seu_combined$celltype_lum[seu_combined$celltype_lum=="NA"] = "Schwann Cell"
  
cellchat <- createCellChat(object = seu_combined, 
                           group.by = "celltype_lum", assay = "RNA")
cellchat <- setIdent(cellchat, ident.use = "celltype_lum") 
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

##插入互作信息
cellchat@DB <- CellChatDB.use

##提取配受体库基因的信息
cellchat@data.signaling
cellchat <- subsetData(cellchat)
cellchat@data.signaling

##计算差异基因并构建交互信息构建
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
OEI_df <- identifyOverExpressedInteractions(cellchat,return.object = F)

##算法核心
## cellchat@idents <- droplevels(cellchat@idents)(subset细胞时候用)
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = T,population.size = TRUE)

##过滤细胞数目少的cluster
cellchat <- filterCommunication(cellchat, min.cells = 10)

##提取cellchat结果
df.net <- subsetCommunication(cellchat) 

##计数通路内的配受体对
cellchat <- computeCommunProbPathway(cellchat)

##计数细胞间交互作用
cellchat <- aggregateNet(cellchat)

##可视化
##  所有细胞交互作用的可视化
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,  vertex.label.cex = 0.8,weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,  vertex.label.cex = 0.8,weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

##  单个细胞发送的信号可视化
mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in seq(1,9)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

##整体信号通路
# Babble plot
netVisual_bubble(cellchat, sources.use = 1:13, targets.use = 1:13, remove.isolate = FALSE,font.size = 8)

## Shared lum cell 基因选择 logFC>1.5
diffgenes = arrange(diffgenes,desc(avg_log2FC))
gene = diffgenes$gene[1:17]

## 验证
## 生存
library(data.table)
clinical <- fread(file = "C:/Users/Administrator/Desktop/R_Date/TCGA/survival_PRAD_survival.txt",data.table = F)
clinical$num <- ifelse(substr(clinical$sample,14,15)=="01","Tum","Nor")

load(file = "C:/Users/Administrator/Desktop/R_Date/TCGA/TCGA_PRAD_TPM_CancerOnly.Rdata")
library(dplyr)
library(tidyr)
tmp = select(data,gene)
tmp$sample = substr(rownames(tmp),1,15)
clinical <- merge(tmp,clinical,by="sample")
clin <- dplyr::filter(clinical,num == "Tum")
library(survival)
library(survminer)
df_clean <- clin %>%
  filter(!is.na(PFI))
group <- ifelse(df_clean$ANLN>median(df_clean$ANLN),'high','low')
sfit <- survfit(Surv(PFI.time, PFI)~group, data=df_clean)
ggsurvplot(sfit, conf.int=F, pval=TRUE,title="ANLN_PFI")

PNI = fread("E:/task/PNI/data_mining/resource/PNI分级最新版.csv",data.table = F)
PNI = PNI[-1,-1]
colnames(PNI) = c("PNI","ID")
clin$ID = substr(clin$sample,1,12)
clin_PNI = merge(clin,PNI,by = "ID")
df_clean <- clin_PNI %>%
  filter(!is.na(PFI))
group <- df_clean$PNI
sfit <- survfit(Surv(PFI.time, PFI)~group, data=df_clean)
ggsurvplot(sfit, conf.int=F, pval=TRUE,title="PNI")

## 表达
data$ID = substr(rownames(data),1,12)
data_PNI=merge(PNI,data,by = "ID")
tmp =data_PNI[,c("PNI",gene)]

tmp <- tmp %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")

library(ggpubr)
my_comparisons <- list(c("PNI","noPNI"))
ggplot(data = tmp,aes(x=PNI,y=expression,fill=PNI))+
  geom_boxplot()+
  geom_jitter(size = 0.5, width = 0.1)+
  theme_bw()+
  facet_grid(.~gene)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")+
  theme(axis.text =element_text(size = 12),
        legend.title = element_text(size = 0.5))+ggsci::scale_fill_d3("category20")   ## x轴字体大小

##相关性
tmp =data_PNI[,c("PNI",gene)]
tmp <- tmp %>%
  group_by(PNI) %>% 
  mutate(PNI_id = paste0(PNI, row_number())) %>% 
  ungroup()
tmp = as.data.frame(tmp)
rownames(tmp) = tmp$PNI_id
tmp = tmp[,-c(1,17)]
tmp = as.data.frame(t(tmp))
heatdata <- tmp
heatdata = select(heatdata,c(paste0("PNI",seq(1,166)),paste0("noPNI",seq(1,48))))

##列聚类
group = colnames(heatdata)
group = gsub("[0-9]+", "", group)
annotation_col <- data.frame(group)
rownames(annotation_col) <-colnames(heatdata)

library(pheatmap)
library(viridisLite)
custom_colors <- colorRampPalette(c("navy", "white", "red"))(100)
pheatmap(heatdata, # 热图的数据
         cluster_rows = TRUE, # 行聚类
         cluster_cols = F, # 列聚类
         annotation_col = annotation_col, # 标注样本分类
         annotation_legend = TRUE, # 显示注释
         show_rownames = TRUE, # 显示行名
         show_colnames = FALSE, # 不显示列名
         scale = "row", # 以行来标准化
         color = custom_colors, # 使用自定义颜色
         cellwidth = 5, # 格子宽度
         cellheight = 20, # 格子高度
         fontsize = 10,# 字体大小
         main = "Shared Luminal Cell DEG"
)

### 选取数据相关性分析
M<-cor(as.data.frame(t(heatdata)))
### 作图
library(corrplot)
### method是展示形式，order是是否聚类，type限定上下三角
corrplot(M, method="pie", order="hclust",tl.col = "black",col = rev(COL2("RdBu",200)))

## 差异分析
load(file = "C:/Users/Administrator/Desktop/R_Date/TCGA/TCGA_PRAD_count.Rdata")
expr_df <- count_matrix
expr_df <- as.data.frame(expr_df)
expr_df <- rownames_to_column(expr_df,"gene_id")

####去掉点
expr_df_nopoint <- expr_df %>%
  tidyr::separate(gene_id,into = c("gene_id"),sep="\\.")

# gtf1 <- rtracklayer::import('resource/Homo_sapiens.GRCh38.104.chr.gtf')
# gtf_df <- as.data.frame(gtf1)
# save(gtf_df,file = "resource/gtf_df.Rda")

####ID转换
load(file = "C:/Users/Administrator/Desktop/R_Date/resource/gtf_df.Rda")
mRNA_exprSet <- gtf_df %>%
  dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>% #筛选gene,和编码指标
  dplyr::select(gene_name,gene_id) %>%
  dplyr::inner_join(expr_df_nopoint,by ="gene_id") %>% 
  select(-gene_id)

###重复基因名
library(dplyr)
library(tidyr)
mRNA_exprSet <- mRNA_exprSet %>% mutate(AEV=rowMeans(.[,-1])) %>% 
  arrange(desc(AEV)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  select(-"AEV") 
mRNA_exprSet <- mRNA_exprSet[rowSums(is.na(mRNA_exprSet)) == 0, ]
mycounts = mRNA_exprSet
rownames(mycounts) = mycounts$gene_name
mycounts = mycounts[,-1]
mycounts = as.data.frame(t(mycounts))
mycounts$state = ifelse(substr(rownames(mycounts),14,15)=="01","Tum","Nor")
mycounts = filter(mycounts,state == "Tum")
mycounts$ID = substr(rownames(mycounts),1,12)
mycounts = merge(mycounts,PNI,by = "ID")

##metadata创建
metadata <- select(mycounts,c("PNI","ID"))

mycounts = select(mycounts,-c("state","PNI"))
rownames(mycounts) <- make.unique(mycounts$ID, sep = "_")
mycounts = mycounts[,-1]
mycounts = as.data.frame(t(mycounts))
mycounts = rownames_to_column(mycounts,"gene_name")

##差异分析
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=mycounts,
                             colData=metadata,
                             design=~PNI,
                             tidy=TRUE)
dds <- DESeq(dds)
nrow(dds)
# save(dds,file="output/mRNA_exprSet_DEseq2_PC.Rdata")

dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd,"PNI")
exprSet_vst <- as.data.frame(assay(vsd))
dds <- DESeq(dds)
### 依次是，1.分组信息(metadata中的列) 2.处理组，3.对照组
contrast=c("PNI", "PNI", "noPNI")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") 
##res <- separate(res,gene_id,into = c("gene_id","others"),sep = " | ")
colnames(res) <- c("gene_id","AveExpr","logFC","lfcSE","P.Value","adj.P.Val")

## IHC基因展示
tmp =data_PNI[,c("PNI",gene)]
library(ggstatsplot)
ggbetweenstats(
  data = tmp,
  x = PNI,
  y = TOP2A)+ggsci::scale_fill_d3("category20") +ggtitle("TOP2A")+
  theme(plot.title = element_text(hjust = 0.5) )

##施万细胞评分
SC = fread("E:/task/PNI/data_mining/resource/TCGA_PNI_total_score.csv",data.table = F)
SC = SC[,-1]
SC = merge(SC,PNI,by = "ID")
colnames(SC)[5] = "Schwann_Cell_Score"

ggbetweenstats(
  data = SC,
  x = PNI,
  y = Schwann_Cell_Score)+ggsci::scale_fill_d3("category20") +ggtitle("Schwann_Cell_Score")+
  theme(plot.title = element_text(hjust = 0.5) )

SC = fread("E:/task/PNI/data_mining/resource/TCGA_PNI_total_score.csv",data.table = F)
SC = SC[,-1]
clin = merge(clin,SC,by = "ID")
df_clean <- clin %>%
  filter(!is.na(PFI))
group <- ifelse(df_clean$Schwann_Cell_total>median(df_clean$Schwann_Cell_total),'high','low')
sfit <- survfit(Surv(PFI.time, PFI)~group, data=df_clean)
ggsurvplot(sfit, conf.int=F, pval=TRUE,title="Schwann_Cell_Score")