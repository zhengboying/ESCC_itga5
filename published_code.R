#big_theme####
#上皮亚群rds---
scObject <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
theme_set(theme_classic(base_size = 22))
big_theme <- theme(
  plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
  axis.title = element_text(size = 22, face = "bold"),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 22, face = "bold"),
  legend.text = element_text(size = 22)
)
source("./项目二：三代单细胞食管癌新辅机制研究/project_proNCI/03.Function/PlotTheme.R")
colors = PlotTheme$Color

#Fig1C####
library(tidyverse)
library(ggplot2)
library(cowplot)  # 用于组合图形
library(ggrepel)  # 用于标签防重叠
library("RColorBrewer")
stats <- read.csv("/Users/user/云平台/项目二：三代单细胞食管癌新辅机制研究/project_new/统计/barcode_stats/stats.csv")
p3 <- ggplot(stats, aes(x = G2_Total, y = G3_Gene_Total)) +
  geom_point(aes(color = Overlap_Ratio), size = 5, alpha = 0.8) +
  xlim(0,16000)+ylim(0,15000)+
  geom_smooth(method = "lm", se = FALSE, color = "darkgray", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dotted") +
  scale_color_gradient(low = "blue", high = "red",
                       name = "Overlap Ratio",
                       limits = c(0, 1)) +
  geom_text_repel(
    aes(label = Sample), 
    size = 5,
    max.overlaps = Inf,       # 最关键修复：取消默认10个重叠标签的限制，允许全量展示
    force = 5,                # 放大标签间排斥力（默认1），强制重叠标签互相分离
    force_pull = 1,           # 平衡标签向数据点的拉力，避免标签离点过远
    box.padding = 1,          # 标签间安全间距翻倍（原0.5），避免标签边缘贴合
    point.padding = 0.5,      # 新增标签与散点的间距，避免标签覆盖数据点
    max.iter = 10000,         # 位置优化迭代次数从2000提升至10000，保证布局最优
    seed = 123,               # 固定随机种子，每次运行标签位置完全一致可复现
    segment.color = "gray50", # 优化标签连线颜色，避免干扰主图
    segment.size = 0.3
  )+
  labs(title = NULL,
       x = "2nd Cell Count", y = "3nd Cell Count") +
  theme_minimal() +
  theme(plot.margin = margin(10,10,10,10),legend.background = element_blank()) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line= element_line(colour = "black",linewidth  = 1)) + 
  theme(axis.text.x = element_text(size = 18,face = "bold"),
        axis.text.y = element_text(face = "bold",size = 18),
        axis.title.y = element_text(size = 20,face = "bold"),
        axis.title.x = element_text(size = 20,face = "bold"),
        legend.text = element_text(size = 15,face = "bold"),
        legend.title = element_blank())
p3

#Fig1D####
percent_data <- stats %>%
  mutate(Total = G2_Total + G3_Gene_Total - Common,
         Common_Pct = Common,
         G2_Only_Pct = G2_Only,
         G3_Only_Pct = G3_Only) %>%
  pivot_longer(cols = ends_with("Pct"),
               names_to = "Type",
               values_to = "Percentage") %>%
  mutate(Type = str_remove(Type, "_Pct"))
my_colors <- brewer.pal(3, "Set2")
names(my_colors) <- c("Common","G2_Only", "G3_Only")
p4 <- ggplot(percent_data, aes(x = Sample, y = Percentage, fill = Type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = my_colors,
                    name = "Cell Type",
                    labels = c("Common", "G2_Only", "G3_Only")) +
  labs(title = "Barcode Overlap Between 2nd and 3nd",
       x = NULL,
       y = "Cells Count") +
  geom_text(aes(label = Percentage),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "black") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18)) + 
  theme(plot.margin = margin(10,10,10,10),legend.background = element_blank()) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line= element_line(colour = "black",linewidth  = 1)) + 
  theme(axis.text.x = element_text(size = 15,face = "bold",angle = 315,hjust = 0.05),
        axis.text.y = element_text(face = "bold",size = 15),
        axis.title.y = element_text(size = 18,face = "bold"),
        axis.title.x = element_text(size = 18,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),legend.title = element_blank())
p4

#Fig1F:总群marker####
CanonicalMarker <- c("PTPRC", "CD3D", "CD3E", # T/NK
                     "CD79A", "MS4A1", # B/Plasma
                     "CD14", "LYZ", # Myeloid
                     "CPA3", "TPSAB1", "KIT", #Mast
                     "PECAM1", "CLDN5", "VWF", # Endothelial
                     "EPCAM", "KRT19", "CDH1",# Epithelial
                     "COL1A1","DCN")# Fibroblast SMC

scObject@meta.data$cell_type <- factor(
  scObject@meta.data$cell_type,
  levels = c(
    "B cells",
    "Endothelial",
    "Epithelial",
    "Fibroblast",
    "Mast",
    "Myeloid",
    "T cells"
  )
)

marker_group <- data.frame(
  gene = CanonicalMarker,
  group = rep(
    c("T cells", "B cells", "Myeloid", "Mast", "Endothelial", "Epithelial", "Fibroblast"),
    times = c(3, 2, 2, 3, 3, 3, 2)
  ),
  stringsAsFactors = FALSE
)

marker_group$gene <- factor(marker_group$gene, levels = CanonicalMarker)

dot_data <- DotPlot(scObject, 
                    features = CanonicalMarker,
                    group.by = "cell_type")$data

dot_data <- merge(dot_data, marker_group, by.x = "features.plot", by.y = "gene")
dot_data$features.plot <- factor(dot_data$features.plot, levels = CanonicalMarker)

cell_order <- rev(c("B cells","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","T cells"))
dot_data$id <- factor(dot_data$id, levels = cell_order) # 关键：固定顺序

group_n <- table(marker_group$group)[unique(marker_group$group)]
split_pos <- cumsum(group_n) + 0.5
split_pos <- c(cumsum(group_n)[-length(cumsum(group_n))] + 0.5, tail(cumsum(group_n),1) + 0.7)
library(ggplot2)
library(ggnewscale)

p <- ggplot(dot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled), 
             shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(0, 9)) +
  geom_vline(xintercept = split_pos, color = "black", linewidth = 0.5) +
  
  
  labs(x = NULL, y = NULL) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 17, angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin = margin(10,10,10,30),
    legend.position = "none"
  ) +
  scale_y_discrete(position = "right", expand = expansion(add = 0.5))

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig1D.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)
#Fig1I####
#细胞数目统计柱形图
celltype_colors <- c(
  "T cells" = "#377EB8",       # 蓝色
  "B cells" = "#4DAF4A",       # 绿色
  "Myeloid" = "#984EA3",       # 紫色
  "Epithelial" = "#bcbd22",       # 橙色
  "Endothelial" = "#FFFF33",   # 黄色
  "Fibroblast" = "#A65628",   # 棕色
  "Mast" = "#17becf"
)
cluster_table <- table(scObject$group_1, scObject$cell_type)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("N", "T")
cluster_order <- c("B cells","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","T cells")
# cluster_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10","11","12")
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
colors = PlotTheme$Color
colors = celltype_colors
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  coord_flip()+
  labs(x = NULL, y = "Proportion", fill = "Cell type") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig1G.pdf", 
       p, bg = "white",width = 8, height = 4, dpi = 300)

#Fig1H####
#细胞数目统计柱形图
celltype_colors <- c(
  "T cells" = "#377EB8",       # 蓝色
  "B cells" = "#4DAF4A",       # 绿色
  "Myeloid" = "#984EA3",       # 紫色
  "Epithelial" = "#bcbd22",       # 橙色
  "Endothelial" = "#FFFF33",   # 黄色
  "Fibroblast" = "#A65628",   # 棕色
  "Mast" = "#17becf"
)
cluster_table <- table(scObject$Sample, scObject$cell_type)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("10820783-N", "10820783-T", "10895817-N", "10895817-T", 
                    "10979579-N", "10979579-T", "11011381-N","11011381-T")
cluster_order <- c("B cells","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","T cells")
# cluster_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10","11","12")
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
# colors = PlotTheme$Color
colors = celltype_colors
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  labs(x = NULL, y = "Proportion", fill = "Cell type") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20,angle = 45,hjust = 1),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.margin = margin(10,10,10,20))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig1F.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)

#Fig1J####
library(Seurat)
library(ggplot2)

# 1. 提取每个细胞的 nFeature（Transcript 矩阵的基因数）
# 注意：这里的nFeature_Transcript 是每个细胞在Transcript assay中检测到的基因数
scObject$nFeature_Transcript <- Matrix::colSums(GetAssayData(scObject, assay = "Transcript", layer = "counts") > 0)

# 2. 提取绘图数据
plot_data <- data.frame(
  cell_type = scObject$cell_type,
  nFeature = scObject$nFeature_Transcript
)

# 3. 按你示例图的顺序定义细胞类型（和图例顺序对应）
cell_order <- c("B cells","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","T cells")
plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_order)

# 4. 示例图同款配色
celltype_colors <- c(
  "T cells" = "#377EB8",       # 蓝色
  "B cells" = "#4DAF4A",       # 绿色
  "Myeloid" = "#984EA3",       # 紫色
  "Epithelial" = "#bcbd22",       # 橙色
  "Endothelial" = "#FFFF33",   # 黄色
  "Fibroblast" = "#A65628",   # 棕色
  "Mast" = "#17becf"
)

# 5. 绘制和示例完全一致的小提琴图
p <- ggplot(plot_data, aes(x = cell_type, y = nFeature, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5) + # 小提琴轮廓线
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "nFeature_RNA") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20,angle = 45,hjust = 1),
    legend.text = element_text(size = 18),
    plot.margin = margin(10,10,10,10),
    legend.position = "right",
    legend.title = element_blank()
  )

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig1H.pdf", 
       p, bg = "white",width = 7, height = 6, dpi = 300)
#Fig1K####
library(Seurat)
library(ggplot2)

# 1. 提取Transcript assay的基因数
scObject$nFeature_Transcript <- Matrix::colSums(GetAssayData(scObject, assay = "Transcript", layer = "counts") > 0)

# 2. 筛选Epithelial细胞，并提取绘图数据
plot_data <- scObject@meta.data %>%
  filter(cell_type == "Epithelial") %>%
  select(cell_type, group_1, nFeature_Transcript)

# 3. 确保分组顺序和示例图一致（N在上，T在下）
plot_data$group_1 <- factor(plot_data$group_1, levels = c("N", "T"))

# 4. 示例图同款配色
group_colors <- c("N" = "#00A087", "T" = "#E64B35")

# 5. 绘制小提琴图（和示例图格式完全匹配）
p <- ggplot(plot_data, aes(x = cell_type, y = nFeature_Transcript, fill = group_1)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.5, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "nFeature_RNA", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 20, face = "bold")
  )

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig1I.pdf", 
       p, bg = "white",width = 4, height = 5, dpi = 300)

#FigS1####
library(Seurat)
library(ggplot2)
library(dplyr)

# 提取元数据---
plot_df <- scObject@meta.data %>%
  select(Sample, nFeature_RNA, nCount_RNA, percent.mt, percent.rib, GenesPerUMI)

# 固定Sample顺序（按你的分组顺序）
plot_df$Sample <- factor(plot_df$Sample, levels = unique(plot_df$Sample))

# 你指定的颜色
my_color <- PlotTheme$Color

my_plot_theme <- theme_classic() +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA), # 完整外边框
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = "none"
  )

# ===================== 1. nFeature_RNA 小提琴图 =====================
p1 <- ggplot(plot_df, aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
  geom_violin(scale = "width", trim = T, linewidth = 0.5) +
  scale_fill_manual(values = my_color) +
  labs(x = NULL, title = "nFeature_RNA") +
  my_plot_theme
p1
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1C.pdf", 
       p1, bg = "white",width = 6, height = 5, dpi = 300)
# ===================== 2. nCount_RNA 小提琴图 =====================
p2 <- ggplot(plot_df, aes(x = Sample, y = nCount_RNA, fill = Sample)) +
  geom_violin(scale = "width", trim = T, linewidth = 0.5) +
  scale_fill_manual(values = my_color) +
  labs(x = NULL, title = "nCount_RNA") +
  my_plot_theme
p2
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1D.pdf", 
       p2, bg = "white",width = 6, height = 5, dpi = 300)
# ===================== 3. percent.mt 箱线图（无点） =====================
p3 <- ggplot(plot_df, aes(x = Sample, y = percent.mt, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = my_color) +
  labs(x = NULL, title = "percent.mt") +
  my_plot_theme+
  theme(plot.margin = margin(10,10,10,30))
p3
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1E.pdf", 
       p3, bg = "white",width = 6, height = 5, dpi = 300)
# ===================== 4. percent.rib 箱线图（无点） =====================
p4 <- ggplot(plot_df, aes(x = Sample, y = percent.rib, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = my_color) +
  labs(x = NULL, title = "percent.hgb") +
  my_plot_theme+
  theme(plot.margin = margin(10,10,10,30))
p4
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1F.pdf", 
       p4, bg = "white",width = 6, height = 5, dpi = 300)
# ===================== 5. GenesPerUMI 箱线图（无点） =====================
p5 <- ggplot(plot_df, aes(x = Sample, y = GenesPerUMI, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = my_color) +
  labs(x = NULL, title = "GenesPerUMI") +
  my_plot_theme+
  theme(plot.margin = margin(10,10,10,30))
p5
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1G.pdf", 
       p5, bg = "white",width = 6, height = 5, dpi = 300)

# ===================== 6. Filter_nCount_RNA-nFeature_RNA_Correlation====================
library(Seurat)
library(ggplot2)
library(dplyr)
library(corrplot)

# 1. 提取绘图数据
plot_df <- scObject@meta.data %>%
  select(Sample, nCount_RNA, nFeature_RNA)

# 2. 计算相关系数（Pearson，和示例图顶部的数字一致）
cor_val <- cor(plot_df$nCount_RNA, plot_df$nFeature_RNA, method = "pearson")
cor_label <- round(cor_val, 2)

# 3. 绘图（完全参照示例图）
p <- ggplot(plot_df, aes(x = nCount_RNA, y = nFeature_RNA, color = Sample)) +
  geom_point(alpha = 0.5, size = 0.3) +  # 点透明度+大小和示例一致
  scale_color_manual(values = PlotTheme$Color) +  # 用你之前的颜色方案
  labs(x = "nCount_RNA", y = "nFeature_RNA", title = cor_label) +  # 标题为相关系数
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = "right",  # 图例在右侧，和示例一致
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA) # 完整外边框
  )+
  guides(color = guide_legend(override.aes = list(size = 4)))

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1H.pdf", 
       p, bg = "white",width = 7, height = 5, dpi = 300)

#FigS1I####
# 统计每个class_code的数量，并按频数降序排序（可选）
library(ggplot2)
library(dplyr)
data1 <- read.table("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/all_samples_gffcompare_strict-match.combined_transcripts.txt",
                    header = T,      # 看你文件应该没有列名
                    fill = TRUE,         # 自动补齐缺失的列（关键修复）
                    comment.char = "",   # 不把 # 当注释（避免误读）
                    sep = "\t"  )
class_code_counts <- data1 %>%
  # 2. 按class_code分组统计频数
  group_by(class_code) %>%
  summarise(count = n(), .groups = "drop") %>%
  # 3. 按频数降序排序
  arrange(desc(count))

# 查看统计结果（仅包含annotion=TRUE的class_code分布）
print(class_code_counts)

# 绘制基础柱形图（ggplot2
p <- ggplot(class_code_counts, 
            aes(x = factor(class_code, levels = class_code),  # 强制按统计顺序排序
                y = count)) +
  geom_col(fill = "#377EB8", width = 0.7) +
  geom_text(aes(label = count), vjust = -0.3, size = 5, fontface = "bold") +
  labs(
    x = "Class Code Type",  
    y = "Number",           
    title = "Distribution of Class Code"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),)


# 显示图表
print(p)
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1I.pdf",
       p,bg = "white",width = 8,height = 6)

#FigS1K####
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 读入你那个class_code统计表格（第一行是表头）
count_df <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/统计/转录本类型统计/isoquant_transcript_class_statistics_wide.csv", header = TRUE)  # 换成你自己的文件路径

# 2. 计算每个样本的总转录本数目（每行加和，排除sample列）
count_df$total_transcripts <- rowSums(count_df[, -1])  # -1是去掉第一列sample

# 3. 从scObject中统计每个样本的细胞数
cell_count_df <- scObject@meta.data %>%
  count(Sample, name = "cell_count")

# 4. 合并转录本总数和细胞数，计算平均值
plot_df <- count_df %>%
  select(sample, total_transcripts) %>%
  left_join(cell_count_df, by = c("sample" = "Sample")) %>%
  mutate(mean_transcripts_per_cell = total_transcripts / cell_count)

# 5. 拆分样本名，分组（按N/T着色）
plot_df <- plot_df %>%
  separate(sample, into = c("patient", "group"), sep = "-", remove = FALSE) %>%
  mutate(group = factor(group, levels = c("N", "T")))  # 控制图例顺序

# 6. 绘制和示例一样的分组柱状图
p <- ggplot(plot_df, aes(x = sample, y = mean_transcripts_per_cell, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("N" = "#377EB8", "T" = "#E41A1C")) +  # 蓝/红配色，和示例一致
  labs(
    x = NULL,
    title = "Mean transcripts per cell"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 20),
    plot.margin = margin(10,10,10,20)
  )

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS1K.pdf",
       p,bg = "white",width = 6,height = 6,dpi = 600)
#FigS2A####
colors = PlotTheme$Color
DefaultAssay(scObject)<-"RNA"
p<-DimPlot(scObject, reduction = "UMAP", group.by = "Sample",
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
           pt.size = 0.5, label.size = 5, seed = 16) +
  scale_color_manual(values = colors)+
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 18)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2A.pdf",
       p,bg = "white",width = 8,height = 6,dpi = 300)
#FigS2B####
colors = PlotTheme$Color
DefaultAssay(scObject)<-"RNA"
p<-DimPlot(scObject, reduction = "UMAP", group.by = "Phase",
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
           pt.size = 0.5, label.size = 5, seed = 16) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 18)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2B.pdf",
       p,bg = "white",width = 6,height = 5,dpi = 300)

#FigS2C####
CanonicalMarker <- c("PTPRC", "CD3D", "CD3E", "CD4", "CD8A", "NKG7", "KLRB1", # T/NK
                     "CD79A", "MS4A1", "IGHD", "JCHAIN", "XBP1", # B/Plasma
                     "CD14", "FCGR3A", "LYZ", "CSF3R", # Myeloid
                     "CPA3", "TPSAB1", "KIT", #Mast
                     "PECAM1", "CLDN5", "VWF", # Endothelial
                     "EPCAM", "KRT19", "CDH1",# Epithelial
                     "COL1A1", "ACTA2", "TAGLN", "DCN", "RGS5") # Fibroblast 
cell_order_plot <- c("T cells", "B cells", "Myeloid", "Mast", "Endothelial", "Epithelial", "Fibroblast")
# 提前把 cell_type 转成因子并固定顺序
scObject$cell_type <- factor(scObject$cell_type, levels = cell_order_plot)
p<-DotPlot(scObject,
           features = CanonicalMarker,
           group.by = "cell_type",
           cols = c("blue", "red"),
           dot.scale = 4) +
  ggtitle("scRNA 2nd genes")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 20, face = "bold"),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        plot.margin = margin(10,10,10,10))
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2C.pdf",
       p,bg = "white",width = 12,height = 6,dpi = 300)

#FigS2D####
scObject_1 <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/5.scRNA_3nd_gene_umap_tsne_nFeature2000_res0.4_cellannotion.Rds"))
p<-DimPlot(scObject_1, reduction = "UMAP", group.by = "Sample",
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
           pt.size = 0.5, label.size = 5, seed = 16) +
  scale_color_manual(values = colors)+
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 18)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2D.pdf",
       p,bg = "white",width = 8,height = 6,dpi = 300)
#FigS2E####
scObject_1$Phase <- scObject$Phase
p<-DimPlot(scObject_1, reduction = "UMAP", group.by = "Phase",
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
           pt.size = 0.5, label.size = 5, seed = 16) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 18)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2E.pdf",
       p,bg = "white",width = 6,height = 5,dpi = 300)

#FigS2F####
scObject_1 <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/5.scRNA_3nd_isoform_umap_tsne_nFeature2000_res0.4_cellannotion.Rds"))
p<-DimPlot(scObject_1, reduction = "UMAP", group.by = "Sample",
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
           pt.size = 0.5, label.size = 5, seed = 16) +
  scale_color_manual(values = colors)+
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 18)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2G.pdf",
       p,bg = "white",width = 8,height = 6,dpi = 300)

#FigS2G####
scObject_1$Phase <- scObject$Phase
p<-DimPlot(scObject_1, reduction = "UMAP", group.by = "Phase",
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
           pt.size = 0.5, label.size = 5, seed = 16) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.title = element_blank(),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 18)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2H.pdf",
       p,bg = "white",width = 6,height = 5,dpi = 300)

#FigS2H####
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 你的颜色
my_color <- PlotTheme$Color

# 2. 生成表格（你原来的代码）
cluster_table <- table(scObject$Sample, scObject$cell_type)

# 3. 转为长数据（ggplot 必须用长数据）
df <- as.data.frame(cluster_table)
colnames(df) <- c("Sample", "cell_type", "count")

# 4. 绘图：横轴=细胞类型，纵轴=数量，堆叠=Sample
p <- ggplot(df, aes(x = cell_type, y = count, fill = Sample)) +
  geom_bar(stat = "identity") +
  # 使用你指定的颜色
  scale_fill_manual(values = my_color) +
  # 标签
  labs(x = "Cell Type", y = "Cell Number", fill = "Sample") +
  # 主题
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# 出图
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2J.pdf",
       p,bg = "white",width = 6,height = 5,dpi = 300)

#FigS2I####
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 你的颜色
my_color <- PlotTheme$Color

# 2. 生成表格（你原来的代码）
cluster_table <- table(scObject$Sample, scObject$cell_type)

# 3. 转为长数据（ggplot 必须用长数据）
df <- as.data.frame(cluster_table)
colnames(df) <- c("Sample", "cell_type", "count")

# 4. 绘图：横轴=细胞类型，纵轴=数量，堆叠=Sample
p <- ggplot(df, aes(x = cell_type, y = count, fill = Sample)) +
  geom_bar(stat = "identity",position = "fill") +
  # 使用你指定的颜色
  scale_fill_manual(values = my_color) +
  # 标签
  labs(x = "Cell Type", y = "Cell Number", fill = "Sample") +
  # 主题
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title = element_blank(),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# 出图
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/figS2K.pdf",
       p,bg = "white",width = 6,height = 5,dpi = 300)




#Fig2D####
scObject<-readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
p<-DimPlot(object = scObject,group.by = "RNA_snn_res.0.2Harmony",
           split.by = "group_1",ncol = 1,pt.size = 0.05)+
  scale_color_manual(values = colors)+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig2B.pdf",
       p,bg = "white",width = 12,height = 20,dpi = 300)
#Fig2B####
cluster_table <- table(scObject$group_1, scObject$RNA_snn_res.0.2Harmony)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("N", "T")
# cluster_order <- c("B cells","Endothelial","Epithelial","Fibroblast","Mast","Myeloid","T cells")
cluster_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10","11","12")
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
colors = PlotTheme$Color
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Proportion", fill = "Cluster") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x  = element_text(size = 22,face = "bold"),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/fig2C.pdf", 
       p, bg = "white",width = 4, height = 4, dpi = 300)


#Fig2C####
cell_counts <- scObject@meta.data %>%
  dplyr::select(RNA_snn_res.0.2Harmony, Sample, group_1) %>%
  group_by(CellType = RNA_snn_res.0.2Harmony, Group = Sample, Group1 = group_1) %>%
  summarise(Count = n(), .groups = 'drop')

# 计算每个CellType中Group的百分比（按CellType分组）
cell_counts_percent <- cell_counts %>%
  group_by(CellType) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# 设置堆叠顺序：N在下，T在上
cell_counts_percent <- cell_counts_percent %>%
  mutate(Group1 = factor(Group1, levels = c("T", "N"))) %>%
  arrange(CellType, Group1)  # 先按CellType，再按Group1排序

# 设置Group因子的水平，确保N样本在前，T样本在后
n_samples <- cell_counts_percent %>% filter(Group1 == "N") %>% pull(Group) %>% unique()
t_samples <- cell_counts_percent %>% filter(Group1 == "T") %>% pull(Group) %>% unique()
sample_levels <- c(t_samples, n_samples)

cell_counts_percent$Group <- factor(cell_counts_percent$Group, levels = sample_levels)

# 设置颜色方案（N组用蓝色系，T组用红色系）
n_colors <- c("#1f77b4", "#4292c6", "#6baed6", "#9ecae1")[1:length(n_samples)]
t_colors <- c("#e31a1c", "#ef3b2c", "#fb6a4a", "#fc9272")[1:length(t_samples)]
group_colors <- setNames(c(t_colors, n_colors), sample_levels)

# 绘制柱形图（N在下，T在上）
p <- ggplot(cell_counts_percent, aes(x = Percentage, y = CellType, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.7, width = 0.8) +
  scale_fill_manual(values = group_colors) +
  labs(y = "Cell Cluster", x = "Percentage (%)", fill = "Sample",
       title = "Sample Composition by Cell Cluster") +
  theme_classic(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    plot.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.line = element_line(linewidth = 1, color = "black"),
    # 去除与坐标轴的缝隙
    plot.margin = margin(5, 5, 5, 5)
  ) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +  # 去除x轴缝隙
  scale_y_discrete(expand = c(0, 0))  # 去除y轴缝隙

print(p)
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/fig2D.pdf",
       p,bg = "white",width = 8,height = 5,dpi = 300)
#Fig2E####
library(dplyr)
library(ggplot2)
library(ggpubr)
cell_counts <- scObject@meta.data %>%
  dplyr::select(RNA_snn_res.0.2Harmony, Sample, group_1) %>%
  group_by(CellType = RNA_snn_res.0.2Harmony, Group = Sample, Group1 = group_1) %>%
  summarise(Count = n(), .groups = 'drop')
#1. 按样本分别计算【CellType %in% c(3,5,8) 总count / 全样本总count】
ratio_df <- cell_counts %>%
  #先按样本Group汇总全样本总数
  group_by(Group) %>%
  mutate(total_all = sum(Count)) %>%
  ungroup() %>%
  #筛选目标亚群3/5/8，按样本求和
  filter(CellType %in% c(3,5,8)) %>%
  group_by(Group, Group1, total_all) %>%
  summarise(sub_total = sum(Count), .groups = "drop") %>%
  #计算比例
  mutate(ratio = sub_total / total_all)

#2. 绘制箱线图，Group1分T/N两组 + 散点 + 组间显著性
p<-ggplot(ratio_df, aes(x = Group1, y = ratio, fill = Group1)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, staplewidth = 0.5,
               colour = "black") +
  ylim(0,0.7)+
  geom_jitter(width = 0.1,size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("T"="#e74c3c","N"="#3498db")) +
  # 显著性变大
  stat_compare_means(comparisons = list(c("T","N")), method = "wilcox.test", label = "p.signif", size = 6) +
  labs(y = "Proportion of C358 epithelial", x = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 完整外边框
    axis.line = element_blank()  # 配合完整边框使用
  ) +
  scale_x_discrete(expand = c(0.3,0.3)) +  # 缩短两组间距
  big_theme
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/fig2E.pdf",
       p,bg = "white",width = 3,height = 5,dpi = 300)


#Fig2I####
library(ggplot2)
library(dplyr)
df1 <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.Epi_sub358_assays_diff/RNA/Epi_RBP1+_T_vs_N_RNA_full.csv")
df2 <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.Epi_sub358_assays_diff/RBP1+_diff_RNA.csv")

merge_df <- inner_join(
  df1 %>% dplyr::select(gene, avg_log2FC, p_val_adj) %>% dplyr::rename(xFC = avg_log2FC, x_padj = p_val_adj),
  df2 %>% dplyr::select(gene, avg_log2FC, p_val_adj) %>% dplyr::rename(yFC = avg_log2FC, y_padj = p_val_adj),
  by = "gene"
)

cutoff = 1.5
merge_df <- merge_df %>%
  mutate(
    group = case_when(
      xFC >= cutoff & yFC >= cutoff & x_padj < 0.05 & y_padj < 0.05 ~ "Up",
      xFC <= -cutoff & yFC <= -cutoff & x_padj < 0.05 & y_padj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

merge_df$group <- factor(merge_df$group, levels = c("Up", "Down", "NS"))

stat_tab <- merge_df %>% count(group)
nUp    = ifelse("Up" %in% stat_tab$group, stat_tab$n[stat_tab$group=="Up"], 0)
nDown  = ifelse("Down" %in% stat_tab$group, stat_tab$n[stat_tab$group=="Down"], 0)
nNS    = ifelse("NS" %in% stat_tab$group, stat_tab$n[stat_tab$group=="NS"], 0)

cor_res <- cor.test(merge_df$xFC, merge_df$yFC, method = "pearson")
r_val = round(cor_res$estimate, 2)
p_val <- cor_res$p.value
if (p_val < 2.225074e-300) {
  p_val <- "< 2.2e-300"
} else {
  p_val <- format(p_val, digits = 2, scientific = TRUE)
}

# 提取 ITGA5 坐标
itga5 <- merge_df %>% filter(gene == "ITGA5")

# 开始绘图
p <- ggplot(merge_df, aes(x = xFC, y = yFC, color = group)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(xintercept = c(-cutoff, cutoff), lty = 2, col = "black", linewidth = 0.8) +
  geom_hline(yintercept = c(-cutoff, cutoff), lty = 2, col = "black", linewidth = 0.8) +
  
  scale_color_manual(
    values = c("Up"="#e72e2e","Down"="#3498db","NS"="gray70"),
    labels = c(paste0("Up (n=",nUp,")"),paste0("Down (n=",nDown,")"),paste0("NS (n=",nNS,")")),
    limits = c("Up","Down","NS")
  ) +
  xlim(-10,12) + ylim(-5,5) +
  labs(x = "TvsN in C358", y = "C358vsOthers") +
  annotate("text",x=5,y=-2.5,label=paste0("r = ",r_val,"\np = ",p_val),hjust=0,size=6) +
  
  
  geom_segment(
    data = itga5,
    aes(x = xFC, y = yFC, xend = xFC + 1.5, yend = yFC + 1.5),
    color = "black", linewidth = 0.8, linetype = 1
  ) +
  geom_text(
    data = itga5,
    aes(x = xFC + 1.6, y = yFC + 1.6, label = gene),
    color = "black", size = 6, fontface = "bold"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.border = element_rect(linewidth = 1),
    plot.title = element_text(size=28, face="bold", hjust=0.5),
    axis.title = element_text(size=22, face="bold"),
    axis.text = element_text(size=20),
    legend.title = element_blank(),
    legend.text = element_text(size=18)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/fig2J.pdf",
       p,bg = "white",width = 7,height = 6,dpi = 300)

#Fig2J####
hallmarks_file<-"/Users/zhengboying/Documents/project/PJ12233-ESCC-scRNA-zhengby/2rd+3rd_N_T/published/h.all.v2024.1.Hs.symbols.gmt"
c358up_genelist_file<-"/Users/zhengboying/Documents/project/PJ12233-ESCC-scRNA-zhengby/2rd+3rd_N_T/published/c358_up_genelist.tsv"
# 
c358up_genelist<-read.csv(c358up_genelist_file, sep = "\t")
c358up_genelist <- c358up_genelist$Gene 


hallmarks <- read.csv(hallmarks_file, sep = "\t", header = FALSE)
term2gene <- do.call(rbind, lapply(1:nrow(hallmarks), function(i) data.frame(Term = hallmarks[i,1], Gene = trimws(as.character(hallmarks[i,3:ncol(hallmarks)])), stringsAsFactors = FALSE)))

# 清理空值
term2gene <- unique(term2gene[!is.na(term2gene$Gene) & term2gene$Gene != "", ])

# 查看结果
head(term2gene)



# 准备背景基因（所有人类Gene Symbol）
library(org.Hs.eg.db)
all_human_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")



result <- enricher(
  gene = c358up_genelist,
  TERM2GENE = term2gene,
  universe = all_human_symbols,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)


dotplot(result, 
        showCategory = 15,           # 一共就10个
        title = "Hallmark Gene Sets Enrichment")

ggsave("Fig2J_hallmark_enrichment_dotplot.pdf", width = 8, height = 6)
result_df <- as.data.frame(result)

# 查看前几行
head(result_df)
write.csv(result_df, "Fig2J_hallmark_enrichment_results.csv", row.names = FALSE)
#FigS3I####
library(ggplot2)
library(dplyr)
library(RColorBrewer)
gsea_res <- read.csv("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.Epi_sub358_assays_diff/enrichment_plots/RBP1+_diff_RNA_GSEA_full.csv",
                     stringsAsFactors = FALSE)

gsea_res <- gsea_res %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%  # 去掉前缀
  filter(padj < 0.05) %>%                              # 显著通路
  mutate(
    group = case_when(
      NES > 0 ~ "activated",
      NES < 0 ~ "suppressed"
    ),
    group = factor(group, levels = c("suppressed", "activated"))  # 固定左右顺序
  ) %>%
  arrange(desc(abs(NES)))  # 按NES排序

p <- ggplot(gsea_res) +
  # 棒棒糖竖线（从0到气泡）
  geom_segment(aes(x = 0, xend = NES, y = reorder(pathway, NES), yend = pathway),
               color = "gray50", linewidth = 1) +
  # 气泡
  geom_point(aes(x = NES, y = pathway, size = size, color = padj),
             alpha = 0.92) +
  
  facet_wrap(~group, nrow = 1, scales = "free_x") +
  
  # 颜色梯度（红=更显著，蓝=低显著）
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "Spectral")),
    trans = "log10",
    name = "P-value"
  ) +
  
  # 气泡大小
  scale_size_continuous(range = c(3,6), name = "Gene Count") +
  
  # 标签
  labs(
    x = "NES",
    y = "Hallmark Pathway",
    title = "GSEA Hallmark Enrichment: C358 vs Others"
  ) +
  theme_bw(base_size = 18) +  # 全局字体变大
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),  # 分面标题字体
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.grid = element_line(color = "gray90"),
    legend.position = "right"
  )
p

ggsave("./项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3J.pdf", 
       p, width = 12, height = 8, dpi=300)


#FigS3J####
library(ggplot2)
library(dplyr)
library(RColorBrewer)
gsea_res <- read.csv("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.Epi_sub358_assays_diff/RNA/enrichment_plots/Epi_sub358_T_vs_N_RNA_GSEA_full.csv",
                     stringsAsFactors = FALSE)

gsea_res <- gsea_res %>%
  mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%  # 去掉前缀
  filter(padj < 0.05) %>%                              # 显著通路
  mutate(
    group = case_when(
      NES > 0 ~ "activated",
      NES < 0 ~ "suppressed"
    ),
    group = factor(group, levels = c("suppressed", "activated"))  # 固定左右顺序
  ) %>%
  arrange(desc(abs(NES)))  # 按NES排序

p <- ggplot(gsea_res) +
  # 棒棒糖竖线（从0到气泡）
  geom_segment(aes(x = 0, xend = NES, y = reorder(pathway, NES), yend = pathway),
               color = "gray50", linewidth = 1) +
  # 气泡
  geom_point(aes(x = NES, y = pathway, size = size, color = padj),
             alpha = 0.92) +
  
  facet_wrap(~group, nrow = 1, scales = "free_x") +
  
  # 颜色梯度（红=更显著，蓝=低显著）
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "Spectral")),
    trans = "log10",
    name = "P-value"
  ) +
  
  # 气泡大小
  scale_size_continuous(range = c(3,6), name = "Gene Count") +
  
  # 标签
  labs(
    x = "NES",
    y = "Hallmark Pathway",
    title = "GSEA Hallmark Enrichment: TvsN in C358"
  ) +
  theme_bw(base_size = 18) +  # 全局字体变大
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),  # 分面标题字体
    strip.background = element_rect(fill = "#f0f0f0"),
    panel.grid = element_line(color = "gray90"),
    legend.position = "right"
  )
p

ggsave("./项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3K.pdf", 
       p, width = 12, height = 8, dpi=300)

#FigS3A####
p<-DimPlot(object = scObject,group.by = "Sample",ncol = 1,pt.size = 0.6,
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE)+
  scale_color_manual(values = colors)+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3A.pdf",
       p,bg = "white",width = 8,height = 6,dpi = 300)
#FigS3B####
p<-DimPlot(object = scObject,group.by = "Phase",ncol = 1,pt.size = 0.6,
           label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE)+
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22)
  )
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3B.pdf",
       p,bg = "white",width = 6,height = 5,dpi = 300)
#FigS3C####
cluster_table <- table(scObject$Sample, scObject$RNA_snn_res.0.2Harmony)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("10820783-N", "10820783-T", "10895817-N", "10895817-T", 
                    "10979579-N", "10979579-T", "11011381-N","11011381-T")
cluster_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9","10","11","12")
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
colors = PlotTheme$Color
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Proportion", fill = "Cluster") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x  = element_text(size = 22,face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    plot.margin = margin(10,10,10,20))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3C.pdf", 
       p, bg = "white",width = 7, height = 6, dpi = 300)
#FigS3D####
all_markers <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.epithelial_subcluster/Epi/HarmonyIntegration2000/diff/res0.2/RNA/cluster_markers_RNA_all.csv")

# 提取每个cluster的top 5 marker基因
top_markers <- all_markers %>%
  filter(
    !grepl("^ENSG", gene, ignore.case = TRUE),    # 排除ENSG开头的
    !grepl("^LINC", gene, ignore.case = TRUE),    # 排除LINC开头的
    !grepl("^IG[HKL]?", gene, ignore.case = TRUE), # 排除IG开头的（包括IGH、IGK、IGL）
    !grepl("^AC\\d", gene, ignore.case = TRUE),   # 排除AC后接数字的
    !grepl("^AP\\d", gene, ignore.case = TRUE),   # 排除AP后接数字的
    !grepl("^AL\\d", gene, ignore.case = TRUE),
    !grepl("-DT$", gene, ignore.case = TRUE),
    !grepl("-AS1$", gene, ignore.case = TRUE)
  ) %>%
  filter(pct.1 > 0.5&pct.2<0.25) %>%  # 先筛选pct.2 < 0.1的基因
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
unique_genes <- unique(top_markers$gene)

p <- DotPlot(scObject, 
             # features = CanonicalMarker,
             features = unique_genes,
             cols = c("blue", "red"),
             dot.scale = 4,
             group.by = "RNA_snn_res.0.2Harmony") +
  coord_flip()+
  # labs(title = paste("Cluster Top5 markers"), x = NULL, y = "Clusters") +
  labs(title = paste("Top5 Markers"), x = NULL, y = "clusters") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank()) +
  big_theme+
  theme(plot.margin = margin(l = 10, r = 10, t = 10, b = 10))  # 增加左侧边距
p  
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3D.pdf", 
       p, bg = "white",width = 6, height = 8, dpi = 300)
#FigS3E####
p<-FeaturePlot(scObject, features = c("RBP1","TNFSF10","VEGFA","TNFAIP2"), coord.fixed = T, 
               cols = c("grey", "yellow", "red"), order = T,
               raster = FALSE,pt.size=0.3)+
  theme(strip.text = element_text(size = 28, face = "bold")) &
  big_theme
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3E.pdf", 
       p, bg = "white",width = 10, height = 8, dpi = 300)

#FigS3F####
cell_scores_CNV <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/6.infercnv_epi_mast_T_ref/cell_scores_CNV.csv", row.names = 1)
common_cells <- intersect(rownames(cell_scores_CNV), colnames(scObject_all))
cell_scores_CNV <- cell_scores_CNV[common_cells, , drop = FALSE]

# # 创建专门的分析列
scObject_all@meta.data$CNV_score_epi_vs_mast_T <- NA
scObject_all@meta.data[rownames(cell_scores_CNV), "CNV_score_epi_vs_mast_T"] <- cell_scores_CNV$cnv_score


# 1. 筛选和分组数据（保持不变）
meta_all <- scObject_all@meta.data

meta_filtered <- meta_all %>%
  filter(cell_type %in% c("T cells", "Mast", "Epithelial"))

meta_filtered <- meta_filtered %>%
  mutate(
    CNV_group = case_when(
      cell_type %in% c("T cells", "Mast") ~ "reference",
      cell_type == "Epithelial" ~ "Epithelial",
      TRUE ~ "Other"
    )
  )

# 2. 生成 fill_color，并转为有序因子（关键步骤！）
meta_plot <- meta_filtered %>%
  mutate(
    fill_color = case_when(
      CNV_group == "reference" ~ "reference",
      CNV_group == "Epithelial" & group_1 == "N" ~ "N",
      CNV_group == "Epithelial" & group_1 == "T" ~ "T",
      TRUE ~ "Other"
    ),
    # 强制绘图顺序：reference → N → T
    fill_color = factor(fill_color, 
                        levels = c("reference", "N", "T"), 
                        ordered = TRUE)
  )

# 3. 定义颜色映射
color_map <- c(
  "reference" = "gray50",
  "N" = "blue",
  "T" = "red"
)

# 4. 绘制箱线图
p <- ggplot(meta_plot, aes(x = Sample, y = CNV_score_epi_vs_mast_T)) +
  geom_boxplot(
    aes(fill = fill_color),
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.3
  ) +
  scale_fill_manual(
    values = color_map,
    breaks = c("reference", "N", "T"),
    labels = c("reference", "N", "T")
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "CNV score",
    fill = NULL
  )+
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.y  = element_text(size = 20),
    axis.text.x  = element_text(size = 20,angle = 45,hjust = 1),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 18)
  )

print(p)
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3F.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)


#FigS3G####
cell_scores_CNV <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/6.infercnv_epi_mast_T_ref/cell_scores_CNV.csv", row.names = 1)
common_cells <- intersect(rownames(cell_scores_CNV), colnames(scObject_all))
cell_scores_CNV <- cell_scores_CNV[common_cells, , drop = FALSE]

# # 创建专门的分析列
scObject_all <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/5.Final_umap_nFeature2000_res0.4_cellannotion.Rds"))
scObject_all@meta.data$CNV_score_epi_vs_mast_T <- NA
scObject_all@meta.data[rownames(cell_scores_CNV), "CNV_score_epi_vs_mast_T"] <- cell_scores_CNV$cnv_score

scObject <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
target_barcode <- scObject@meta.data %>%
  filter(RNA_snn_res.0.2Harmony %in% c(3, 5, 8)) %>%
  rownames()
meta_all <- scObject_all@meta.data
meta_filtered <- meta_all %>%
  filter(
    # 满足二选一：细胞类型是T cells/Mast  OR  细胞名在目标barcode里
    cell_type %in% c("T cells", "Mast") | rownames(.) %in% target_barcode
  )

meta_filtered <- meta_filtered %>%
  mutate(
    CNV_group = case_when(
      cell_type %in% c("T cells", "Mast") ~ "reference",
      cell_type == "Epithelial" ~ "Epithelial",
      TRUE ~ "Other"
    )
  )

# 2. 生成 fill_color，并转为有序因子（关键步骤！）
meta_plot <- meta_filtered %>%
  mutate(
    fill_color = case_when(
      CNV_group == "reference" ~ "reference",
      CNV_group == "Epithelial" & group_1 == "N" ~ "N",
      CNV_group == "Epithelial" & group_1 == "T" ~ "T",
      TRUE ~ "Other"
    ),
    # 强制绘图顺序：reference → N → T
    fill_color = factor(fill_color, 
                        levels = c("reference", "N", "T"), 
                        ordered = TRUE)
  )

# 3. 定义颜色映射
color_map <- c(
  "reference" = "gray50",
  "N" = "blue",
  "T" = "red"
)

# 4. 绘制箱线图
p <- ggplot(meta_plot, aes(x = Sample, y = CNV_score_epi_vs_mast_T)) +
  geom_boxplot(
    aes(fill = fill_color),
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.size = 0.5,
    outlier.alpha = 0.3
  ) +
  scale_fill_manual(
    values = color_map,
    breaks = c("reference", "N", "T"),
    labels = c("reference", "N", "T")
  ) +
  theme_bw() +
  labs(
    title = "C358 Epithelial CNV score",
    x = NULL,
    y = "CNV score",
    fill = NULL
  )+
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text.y  = element_text(size = 20),
    axis.text.x  = element_text(size = 20,angle = 45,hjust = 1),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 18)
  )

print(p)
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3G.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)

#FigS3H####
library(pheatmap)
corr_matrix <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.epithelial_subcluster/Epi/HarmonyIntegration2000/aucell/spearman_correlation_all_genesets.csv", 
                        header = TRUE, row.names = 1, check.names = FALSE)

# 转成矩阵（防止数据框格式导致问题）
corr_matrix <- as.matrix(corr_matrix)

color_palette <- colorRampPalette(c(
  "#4A7BB7",  # 深蓝色（负相关端）
  "#A6CBE3",  # 浅蓝
  "#FFFFFF",  # 白色（中间0）
  "#FFF6CC",  # 浅黄
  "#FFB366",  # 橘黄
  "#E64B35"   # 橘红（正相关端）
))(100)

# 3. 绘制热图
p<-pheatmap(
  corr_matrix,
  color = color_palette,
  breaks = seq(-1, 1, length.out = 101),  # 固定色阶范围-1~1
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 30,
  treeheight_col = 30,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 270,
  main = "Spearman corplot",
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  legend = TRUE,
  legend_breaks = seq(-1, 1, 0.2),
  legend_labels = seq(-1, 1, 0.2),
  legend_title = ""
)

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2-6/figS3H.pdf", 
       p, bg = "white",width = 10, height = 8, dpi = 300)


#Fig3A/B####

data_fusion <- read.csv("fusion_data.csv", sep = ",", header = TRUE,row.names = 1)
##
##
#all_values <- unlist(data_fusion[, 4:ncol(data_fusion)])
#> table(all_values)
#all_values
#0       1       2       3       4       5       6       7       8       9      10      11 
#8568956    8620    1520     468     183      93      75      24      14       3       2       1 
#12      13      14      15      16      17      20      21      22      25      28      29 
#3       1       4       1       1       1       1       1       1       1       1       3 
#31      36      38      39      44      48      53     166     227     449     481    1961 
#1       1       1       1       1       1       1       1       1       1       1       1 
#1996 
#1


####data_fusion所有值为1,2替换成0
data_fusion_clean<-data_fusion
data_fusion_clean[data_fusion_clean == 1 | data_fusion_clean == 2] <- 0

###计算每行 4:ncol(data_fusion_clean)不为零的值的数量，添加新列fusion_count
data_fusion_clean$fusion_count <- apply(data_fusion_clean[, 4:ncol(data_fusion_clean)], 1, function(x) sum(x != 0))
##提取Sample列,cell_barcode,scRNA_2nd_cell_type,fusion_count列，保存为data_fusion_df
data_fusion_df <- data_fusion_clean[, c("Sample", "cell_barcode", "scRNA_2nd_cell_type", "fusion_count")]

##提取data_fusion_df的fusion_count不为0的行
data_fusion_nonzero <- data_fusion_df[data_fusion_df$fusion_count != 0, ]
#> table(data_fusion_nonzero$Sample)

#10820783-N 10820783-T 10895817-N 10895817-T 10979579-N 10979579-T 11011381-N 11011381-T 
#156         35         51        119         40         56         48        288 

##以Sample为横坐标，scRNA_2nd_cell_type做堆叠柱状图
celltype_colors <- c(
  "T cells" = "#377EB8",  
  "B cells" = "#4DAF4A",       
  "Myeloid" = "#984EA3",       
  "Epithelial" = "#bcbd22",       
  "Endothelial" = "#FFFF33",   
  "Fibroblast" = "#A65628",   
  "Mast" = "#17becf"
)
##count数打在图上
ggplot(data_fusion_nonzero, aes(x = Sample, fill = scRNA_2nd_cell_type)) +
  geom_bar(position = "stack") +
  geom_text(stat = "count",
            aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "white") +  # 白色文字更清晰
  labs(x = "Sample", y = "Count", fill = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = celltype_colors)

##
ggsave("Fig3_fusion_stacked_barplot.pdf", width = 7, height = 4)

##增加group列，Sample的后缀为N 定为N组，否则为T组
data_fusion_nonzero$Group <- ifelse(grepl("-N$", data_fusion_nonzero$Sample), "N", "T")

##以Group为横坐标，scRNA_2nd_cell_type做堆叠柱状图
ggplot(data_fusion_nonzero, aes(x = Group, fill = scRNA_2nd_cell_type)) +
  geom_bar(position = "stack") +
  geom_text(stat = "count",
            aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "white") +  # 白色文字更清晰
  labs(x = "Group", y = "Count", fill = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text()) +
  scale_fill_manual(values = celltype_colors)
ggsave("Fig3_fusion_stacked_barplot_group.pdf", width = 3, height = 4)


#Fig3D####
scObject <- readRDS("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/8.scRNA_Fusion_gene_umap_tsne_nFeature2000_res0.4_cellannotion.Rds")

volcano_data <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/10.fusion_gene/2nd_gene/diff_Epithelial_T_vs_N_full.csv")
top_genes <- volcano_data %>%
  filter(significant == "Significant") %>%
  top_n(n = 25, wt = avg_log2FC)
significant_genes <- unique(top_genes$gene)  # 确保唯一性
# 检查融合基因是否存在于Seurat对象中
available_genes <- significant_genes[significant_genes %in% rownames(scObject)]

# 绘制点图
p1 <- DotPlot(scObject,
              features = available_genes,
              group.by = "Sample",dot.scale = 10) +
  ylab("Sample Group")+
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x  = element_text(angle = 45,size = 20,hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18)
  ) +
  coord_flip()+
  ggtitle("Fusion Genes in Epithelial cells")
p1
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig3D.pdf",
       p1, bg = "white", width = 10, height = 8, dpi = 300)

#Fig4A####
library(ggplot2)
library(dplyr)
df1 <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.Epi_sub358_assays_diff/Transcript/Epi_RBP1+_T_vs_N_Transcript_full.csv")
df2 <- read.csv("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/5.Epi_sub358_assays_diff/RBP1+_diff_Transcript.csv")

merge_df <- inner_join(
  df1 %>% dplyr::select(gene, avg_log2FC, p_val_adj) %>% dplyr::rename(xFC = avg_log2FC, x_padj = p_val_adj),
  df2 %>% dplyr::select(gene, avg_log2FC, p_val_adj) %>% dplyr::rename(yFC = avg_log2FC, y_padj = p_val_adj),
  by = "gene"
)

cutoff = 1.5
merge_df <- merge_df %>%
  mutate(
    group = case_when(
      xFC >= cutoff & yFC >= cutoff & x_padj < 0.05 & y_padj < 0.05 ~ "Up",
      xFC <= -cutoff & yFC <= -cutoff & x_padj < 0.05 & y_padj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

merge_df$group <- factor(merge_df$group, levels = c("Up", "Down", "NS"))

stat_tab <- merge_df %>% count(group)
nUp    = ifelse("Up" %in% stat_tab$group, stat_tab$n[stat_tab$group=="Up"], 0)
nDown  = ifelse("Down" %in% stat_tab$group, stat_tab$n[stat_tab$group=="Down"], 0)
nNS    = ifelse("NS" %in% stat_tab$group, stat_tab$n[stat_tab$group=="NS"], 0)

cor_res <- cor.test(merge_df$xFC, merge_df$yFC, method = "pearson")
r_val = round(cor_res$estimate, 2)
p_val <- cor_res$p.value
if (p_val < 2.225074e-300) {
  p_val <- "< 2.2e-300"
} else {
  p_val <- format(p_val, digits = 2, scientific = TRUE)
}

# 提取 ITGA5 坐标
itga5 <- merge_df %>% filter(gene == "TCONS-00035775")

# 开始绘图
p <- ggplot(merge_df, aes(x = xFC, y = yFC, color = group)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(xintercept = c(-cutoff, cutoff), lty = 2, col = "black", linewidth = 0.8) +
  geom_hline(yintercept = c(-cutoff, cutoff), lty = 2, col = "black", linewidth = 0.8) +
  
  scale_color_manual(
    values = c("Up"="#e72e2e","Down"="#3498db","NS"="gray70"),
    labels = c(paste0("Up (n=",nUp,")"),paste0("Down (n=",nDown,")"),paste0("NS (n=",nNS,")")),
    limits = c("Up","Down","NS")
  ) +
  xlim(-10,12) + ylim(-5,5) +
  labs(x = "TvsN in C358", y = "C358vsOthers") +
  annotate("text",x=5,y=-2.5,label=paste0("r = ",r_val,"\np = ",p_val),hjust=0,size=6) +
  
  
  geom_segment(
    data = itga5,
    aes(x = xFC, y = yFC, xend = xFC + 1.5, yend = yFC + 1.5),
    color = "black", linewidth = 0.8, linetype = 1
  ) +
  geom_text(
    data = itga5,
    aes(x = xFC, y = yFC + 1.6, label = gene),
    color = "black", size = 6, fontface = "bold"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.border = element_rect(linewidth = 1),
    plot.title = element_text(size=28, face="bold", hjust=0.5),
    axis.title = element_text(size=22, face="bold"),
    axis.text = element_text(size=20),
    legend.title = element_blank(),
    legend.text = element_text(size=18)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig4A.pdf",
       p,bg = "white",width = 7,height = 6,dpi = 300)

#Fig3F####
DefaultAssay(scObject_Fusion)<-"Fusion"
p <- VlnPlot(scObject_Fusion, features = "SF3A3-MMACHC", 
             group.by = "Sample",pt.size = 0.1, ncol = 1,cols = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 18),
        axis.title.x = element_blank())
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig3F.pdf",
       p,bg = "white",width = 7,height = 4,dpi = 300)

#Fig3G####
celltype_colors <- c(
  "B cells" = "#4DAF4A",       # 绿色
  "Endothelial" = "#FFFF33",   # 黄色
  "Epithelial" = "#bcbd22",       # 橙色  
  "Fibroblast" = "#A65628",   # 棕色
  "Mast" = "#17becf",
  "Myeloid" = "#984EA3",       # 紫色
  "T cells" = "#377EB8"      # 蓝色
)
# scObject <- readRDS("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/8.scRNA_Fusion_gene_umap_tsne_nFeature2000_res0.4_cellannotion.Rds")
DefaultAssay(scObject)<-"Fusion"
scObject$cell_type <- factor(scObject$cell_type, levels = names(celltype_colors))

p <- VlnPlot(
  scObject, 
  features = "SF3A3-MMACHC", 
  group.by = "cell_type",
  pt.size = 0.1, 
  ncol = 1,
  cols = celltype_colors
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x = element_blank()
  )
print(p)
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig3K.pdf",
       p,bg = "white",width = 7,height = 4,dpi = 300)


#Fig4C####

# 读取数据
data <- read.csv("Table_S6_diff_isoform.tsv", sep = "\t", header = TRUE)

# 绘制柱状图（横放），data$Chr数量，排序
ggplot(data, aes(x = fct_rev(fct_infreq(Chr)))) +
  geom_bar(fill = "steelblue") +
  coord_flip() +
  labs(x = "Chromosome", y = "Number of Isoforms") +
  theme_classic()
ggsave("Fig4_diff_isoform_barplot.pdf", width = 6, height = 12)

#Fig4D####
#class_code的Ensembl替换成IsoQuant
data$Source <- gsub("ENSEMBL", "IsoQuant", data$Source)
##堆叠柱状图横轴Source 纵轴堆叠分组class_code
# 使用 hue_pal 增加饱和度
library(scales)
ggplot(data, aes(x = Source, fill = class_code)) +
  geom_bar(position = "stack") +
  labs(x = "Source", y = "Count", fill = "Class Code") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_hue(saturation = 0.8, l = 55)  # 提高饱和度

ggsave("Fig4_diff_isoform_stacked_barplot.pdf", width = 2.5, height = 5)

#Fig4E####
##Protein_coding 和Source 饼图
# 计算每个Protein_coding的class_code，
library(dplyr)
pie_data <- data %>%
  group_by(class_code,Protein_coding) %>%
  summarise(count = n()) %>%
  ungroup()
# 绘制饼图
ggplot(pie_data, aes(x = "", y = count, fill = class_code)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ Protein_coding) +
  labs(x = NULL, y = NULL, fill = "Protein Coding") +
  theme_classic() +
  scale_fill_brewer(palette = "Paired")
theme(legend.position = "default")

ggsave("Fig4_diff_isoform_pie_chart.pdf", width = 6, height = 3)

#Fig4G####
# 读取数据
#data的gene_name和Isoform以“_”串联生成新列ID
data$ID <- paste(data$gene_name, data$Isoform, sep = "_")

##MS质谱有数据检出的isoform(非特异性)
gene_list <- c(
  "ITGA5_TCONS-00035775",
  "KRTDAP_TCONS-00068213",
  "CEBPB_TCONS-00087975",
  "RBP1_TCONS-00096269",
  "TPM4_TCONS-00069311",
  "SDC1_TCONS-00081659",
  "HSPA8_TCONS-00026698",
  "HSP90B1_TCONS-00031184",
  "KRT6B_TCONS-00034309",
  "GSTA1_TCONS-00116840",
  "TUBA4A_TCONS-00084661",
  "ARPC2_TCONS-00084160",
  "SERPINB4_TCONS-00064062",
  "GJB2_TCONS-00038241",
  "CDH3_TCONS-00052884",
  "HLA-E_TCONS-00119443",
  "HLA-B_TCONS-00120382",
  "NFKBIA_TCONS-00039779",
  "CTSB_TCONS-00131009"
)


# 选择ID是gene_list中的数据
top19_data <- data %>% filter(ID %in% gene_list)
# 绘制棒棒糖图,gene_name是y轴
##左边棒棒糖长度avg_log2FC_C358vsAllCluster	圆点大小pct.1_C358vsAllCluster 颜色p_val_adj_C358vsAllCluster
##右边棒棒糖长度avg_log2FC_NvsT	圆点大小pct.1_NvsT 颜色p_val_adj_NvsT
ggplot(top19_data) +
  geom_segment(aes(x = -avg_log2FC_C358vsAllCluster, xend = 0, y = ID), color = "grey") +
  geom_point(aes(x = -avg_log2FC_C358vsAllCluster, y = ID, 
                 size = pct.1_C358vsAllCluster, 
                 color = p_val_adj_C358vsAllCluster)) +
  geom_segment(aes(x = 0, xend = avg_log2FC_NvsT, y = ID), color = "grey") +
  geom_point(aes(x = avg_log2FC_NvsT, y = ID, 
                 size = pct.1_NvsT, 
                 color = p_val_adj_NvsT)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient(low = "blue", high = "red", trans = "log10") +
  labs(x = "Log2 Fold Change", y = "Isoforms", 
       color = "Adjusted P-value", 
       size = "pct.1") +  # 直接使用 pct.1，不翻译
  theme_classic()
##
ggsave("Fig4_diff_isoform_lollipop_plot.pdf", width = 8, height = 6)


#Fig4B####
##overlap的isoform>genes,功能富集

hallmarks_file<-"/Users/zhengboying/Documents/project/PJ12233-ESCC-scRNA-zhengby/2rd+3rd_N_T/published/h.all.v2024.1.Hs.symbols.gmt"
c358_up_isoform2genelist_file<-"/Users/zhengboying/Documents/project/PJ12233-ESCC-scRNA-zhengby/2rd+3rd_N_T/published/c358_up_isoform2genelist.tsv"
# 
c358_up_isoform2genelist<-read.csv(c358_up_isoform2genelist_file, sep = "\t")
c358_up_isoform2genelist <- c358_up_isoform2genelist$Gene 


hallmarks <- read.csv(hallmarks_file, sep = "\t", header = FALSE)
term2gene <- do.call(rbind, lapply(1:nrow(hallmarks), function(i) data.frame(Term = hallmarks[i,1], Gene = trimws(as.character(hallmarks[i,3:ncol(hallmarks)])), stringsAsFactors = FALSE)))

# 清理空值
term2gene <- unique(term2gene[!is.na(term2gene$Gene) & term2gene$Gene != "", ])

# 查看结果
head(term2gene)

# 准备背景基因（所有人类Gene Symbol）
library(org.Hs.eg.db)
all_human_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

result_iso <- enricher(
  gene = c358_up_isoform2genelist,
  TERM2GENE = term2gene,
  universe = all_human_symbols,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)


dotplot(result_iso, 
        showCategory = 15,           # 一共就10个
        title = "Hallmark Gene Sets Enrichment")

ggsave("Fig4B_hallmark_iso2gene_enrichment_dotplot.pdf", width = 7, height = 6)
result_iso_df <- as.data.frame(result_iso)

# 查看前几行
head(result_iso_df)
write.csv(result_iso_df, "Fig4B_hallmark_iso2gene_enrichment_results.csv", row.names = FALSE)


#Fig4H####
DefaultAssay(scObject_Fusion)<-"Transcript"
p <- VlnPlot(scObject_Fusion, features = "TCONS-00035775", 
             group.by = "Sample",pt.size = 0.1, ncol = 1,cols = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 18),
        axis.title.x = element_blank())
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig4H.pdf",
       p,bg = "white",width = 7,height = 4,dpi = 300)


#Fig4I####
# scObject<-readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
p<-FeaturePlot(scObject_Fusion, features = c("TCONS-00035775"),reduction = "UMAP",
               split.by = "group_1",
               coord.fixed = T, cols = c("grey", "yellow", "red"), order = T,
               ncol = 1,raster = FALSE,pt.size=1.2)
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/fig4I.pdf",
       p,bg = "white",width = 8,height = 4,dpi = 300)

#FigS4B####
# scObject<-readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
DefaultAssay(scObject_Fusion)<-"RNA"
p<-FeaturePlot(scObject_Fusion, features = c("ITGA5"),reduction = "UMAP",
               split.by = "group_1",
               coord.fixed = T, cols = c("grey", "yellow", "red"), order = T,
               ncol = 1,raster = FALSE,pt.size=1.2)
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/figS5B.pdf",
       p,bg = "white",width = 8,height = 4,dpi = 300)

#FigS4C####
# scObject<-readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
p<-FeaturePlot(scObject, features = c("TCONS-00035775"),reduction = "UMAP",
               split.by = "group_1",
               coord.fixed = T, cols = c("grey", "yellow", "red"), order = T,
               ncol = 1,raster = FALSE,pt.size=0.8)
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/figS5C.pdf",
       p,bg = "white",width = 8,height = 4,dpi = 300)

#FigS4D####
# scObject<-readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/Rds/6.umap_tsne_Epi_nFeature2000_with_score.Rds"))
p<-FeaturePlot(scObject, features = c("ITGA5"),reduction = "UMAP",
               split.by = "group_1",
               coord.fixed = T, cols = c("grey", "yellow", "red"), order = T,
               ncol = 1,raster = FALSE,pt.size=0.8)
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig2_S3/figS5D.pdf",
       p,bg = "white",width = 8,height = 4,dpi = 300)


#Fig5B####
library(pheatmap)
##the person——cor is from GEPIA2 WEBSET
TF4_cor <- read.csv("/Users/zhengboying/Documents/project/PJ12233-ESCC-scRNA-zhengby/2rd+3rd_N_T/published/TF4_cor_with_ITGA5_TCGA.tsv",sep = "\t", row.names = 1)

TF4_cor[] <- lapply(TF4_cor, as.numeric)
# 提取4个基因的相关系数和P值
#genes <- c("RUNX1", "SREBF1", "TCF4", "PRDM1")
genes <- c("RUNX1", "PRDM1")
corr <- TF4_cor[, paste0(genes, "_person")]
pval <- TF4_cor[, paste0(genes, "_P.value")]

my_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
sig <- matrix("", nrow = nrow(pval), ncol = ncol(pval))
colnames(sig) <- colnames(pval)
rownames(sig) <- rownames(pval)

# 添加星号标记（不显示数值）
sig[pval < 0.05] <- "*"
sig[pval < 0.01] <- "**"
sig[pval < 0.001] <- "***"

# 定义颜色断点（确保0对应白色）
max_abs <- max(abs(corr), na.rm = TRUE)  # 获取最大绝对值
breaks <- seq(-max_abs, max_abs, length.out = 101)  # 对称断点，0在中间
# 直接转置 corr 矩阵
pheatmap(t(corr),
         main = "Correlation Heatmap of Transcription Factors",
         display_numbers = t(sig),  # 注意：sig 也要转置
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = my_colors,
         fontsize_number = 8,
         angle_col = 45,
         breaks = breaks,  # 关键：设置对称断点
         filename = "Fig5B_TF4_ITGA5_correlation_heatmap.pdf",  # 直接保存
         width = 10,
         height = 1.25)


#Fig5C####
epi<-readRDS("/Users/zhengboying/Documents/project/PJ12233-ESCC-scRNA-zhengby/2rd+3rd_N_T/epi/6.umap_tsne_Epi_nFeature2000.Rds")

##C1-8
DimPlot(epi,group.by = "RNA_snn_res.0.2Harmony",label=T)
##新增一列meta.data ,RNA_snn_res.0.2Harmony是3,5,8的时候，标记为ITGA5+,否则标记为ITGA5-
epi$ITGA5_status <- ifelse(epi$RNA_snn_res.0.2Harmony %in% c(3, 5, 8), "ITGA5+", "ITGA5-")


FeaturePlot(epi,features = c("RUNX1"))
ggsave("Fig5C_RUNX1.pdf",width = 7,height = 6)

#Fig5D####
VlnPlot(epi, features = c("RUNX1"),
        group.by = "ITGA5_status", 
        split.by = "group_1",pt.size = 0) 
ggsave("Fig5D_RUNX1_violin_plot.pdf", width = 4, height = 6)


#FigS5B####
FeaturePlot(epi,features = c("PRDM1"))
ggsave("FigS5_PRDM1.pdf",width = 7,height = 6)

VlnPlot(epi, features = c("PRDM1"),
        group.by = "ITGA5_status", 
        split.by = "group_1",pt.size = 0) 
ggsave("FigS5B_PRDM1_violin_plot.pdf", width = 4, height = 6)


#Fig5I(髓系)####
library(ggplot2)
library(dplyr)
library(ggpubr)
# 读取数据
scObject <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Myeloid/Myeloid_feature2000.rds"))

# 提取meta.data
meta <- scObject@meta.data

# 筛选需要的分组
meta_filtered <- meta[meta$group_1 %in% c("N", "T"), ]

# 计算每个样本中每个亚群占该样本所有T细胞总数的比例
plot_data <- meta_filtered %>%
  group_by(Sample, group_1, cluster_annotation) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(total_cells = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total_cells * 100)
# 确保分组顺序
plot_data$group_1 <- factor(plot_data$group_1, 
                            levels = c("N", "T"))

sci_colors <- c("N" = "#0072B2",   
                "T" = "#D55E00")
# 定义比较组
comparisons <- list(
  c("N", "T")
)

subclusters <- unique(plot_data$cluster_annotation)

# 获取所有亚群名称
for (sub in subclusters) {
  # 筛选当前亚群数据
  sub_data <- plot_data[plot_data$cluster_annotation == sub, ]
  
  # 检测每个比较组是否有足够的观测值（每组至少2个样本）
  valid_comparisons <- list()
  for (comp in comparisons) {
    group1_data <- sub_data$proportion[sub_data$group_1 == comp[1]]
    group2_data <- sub_data$proportion[sub_data$group_1 == comp[2]]
    
    # 检查每组是否有至少2个观测值
    if (length(group1_data) >= 2 & length(group2_data) >= 2) {
      valid_comparisons <- c(valid_comparisons, list(comp))
    }
  }
  
  # 打印有效比较组信息
  cat("\n========== ", sub, " ==========\n")
  cat("有效比较组数量:", length(valid_comparisons), "\n")
  for (comp in valid_comparisons) {
    cat(paste(comp[1], "vs", comp[2], "- 有足够观测值\n"))
  }
  
  # 绘制箱线图
  p <- ggplot(sub_data, aes(x = group_1, y = proportion, fill = group_1)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
    scale_fill_manual(values = sci_colors) +
    theme_classic() +
    labs(title = sub,
         x = "", y = "Proportion (%)") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 18, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 18)
    )
  
  # 只有当有有效比较组时才添加p值
  if (length(valid_comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = valid_comparisons, 
                                method = "t.test",
                                label = "p.format", 
                                label.y.npc = "top", 
                                tip.length = 0.03, 
                                size = 5)
  } else {
    cat("警告: 没有足够的观测值进行任何比较\n")
  }
  
  # 保存图片
  ggsave(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Myeloid/boxplot/", sub, "_响应_boxplot.png"), 
         p, width = 6, height = 6, dpi = 300)
  ggsave(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Myeloid/boxplot/", sub, "_响应_boxplot.pdf"), 
         p, width = 6, height = 6, dpi = 300)
  
  # 打印图片
  print(p)
}

#Fig5I(成纤维)####
library(ggplot2)
library(dplyr)
library(ggpubr)
# 读取数据
scObject <- readRDS(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Myeloid/Myeloid_feature2000.rds"))

# 提取meta.data
meta <- scObject@meta.data

# 筛选需要的分组
meta_filtered <- meta[meta$group_1 %in% c("N", "T"), ]

# 计算每个样本中每个亚群占该样本所有T细胞总数的比例
plot_data <- meta_filtered %>%
  group_by(Sample, group_1, cluster_annotation) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(total_cells = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total_cells * 100)
# 确保分组顺序
plot_data$group_1 <- factor(plot_data$group_1, 
                            levels = c("N", "T"))

sci_colors <- c("N" = "#0072B2",   
                "T" = "#D55E00")
# 定义比较组
comparisons <- list(
  c("N", "T")
)

subclusters <- unique(plot_data$cluster_annotation)

# 获取所有亚群名称
for (sub in subclusters) {
  # 筛选当前亚群数据
  sub_data <- plot_data[plot_data$cluster_annotation == sub, ]
  
  # 检测每个比较组是否有足够的观测值（每组至少2个样本）
  valid_comparisons <- list()
  for (comp in comparisons) {
    group1_data <- sub_data$proportion[sub_data$group_1 == comp[1]]
    group2_data <- sub_data$proportion[sub_data$group_1 == comp[2]]
    
    # 检查每组是否有至少2个观测值
    if (length(group1_data) >= 2 & length(group2_data) >= 2) {
      valid_comparisons <- c(valid_comparisons, list(comp))
    }
  }
  
  # 打印有效比较组信息
  cat("\n========== ", sub, " ==========\n")
  cat("有效比较组数量:", length(valid_comparisons), "\n")
  for (comp in valid_comparisons) {
    cat(paste(comp[1], "vs", comp[2], "- 有足够观测值\n"))
  }
  
  # 绘制箱线图
  p <- ggplot(sub_data, aes(x = group_1, y = proportion, fill = group_1)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
    scale_fill_manual(values = sci_colors) +
    theme_classic() +
    labs(title = sub,
         x = "", y = "Proportion (%)") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 18, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 18)
    )
  
  # 只有当有有效比较组时才添加p值
  if (length(valid_comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = valid_comparisons, 
                                method = "t.test",
                                label = "p.format", 
                                label.y.npc = "top", 
                                tip.length = 0.03, 
                                size = 5)
  } else {
    cat("警告: 没有足够的观测值进行任何比较\n")
  }
  
  # 保存图片
  ggsave(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Fibroblast/boxplot/", sub, "_响应_boxplot.png"), 
         p, width = 6, height = 6, dpi = 300)
  ggsave(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Fibroblast/boxplot/", sub, "_响应_boxplot.pdf"), 
         p, width = 6, height = 6, dpi = 300)
  
  # 打印图片
  print(p)
}


#FigS7C####
cluster_table <- table(scObject$Sample, scObject$cluster_annotation)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("10820783-N", "10820783-T", "10895817-N", "10895817-T", 
                    "10979579-N", "10979579-T", "11011381-N","11011381-T")
cluster_order <- c("APOE_Macrophage",
                   "SPP1_Macrophage",
                   "MKI67_Mcyc",
                   "FCER1A_DC",
                   "LAMP3_DC",
                   "LILRA4_DC",
                   "CLEC9A_DC")
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
colors = PlotTheme$Color
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Proportion", fill = "Cluster") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x  = element_text(size = 22,face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    plot.margin = margin(10,10,10,20))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig5_S6_S7/figS6A.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)

#Fig7E####
scObject<-readRDS("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/Fibroblast/Fibroblast_feature2000.rds")
cluster_table <- table(scObject$Sample, scObject$cluster_annotation)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("10820783-N", "10820783-T", "10895817-N", "10895817-T", 
                    "10979579-N", "10979579-T", "11011381-N","11011381-T")
cluster_order <- c("CFD_iCAF",
                   "MMP11_myCAF",
                   "RGS5_myCAF",
                   "STMN1_CAF",
                   "MUSTN1_myCAF",
                   "DES_myCAF",
                   "SFRP4_myCAF",
                   "CD74_apCAF")
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
colors = PlotTheme$Color
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Proportion", fill = "Cluster") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x  = element_text(size = 22,face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    plot.margin = margin(10,10,10,20))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig5_S6_S7/figS6D.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)

#FigS8B####
scObject<-readRDS("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/T cells/T cells_feature2000.rds")
cluster_table <- table(scObject$Sample, scObject$cluster_annotation)

# 转换为数据框并绘图
library(ggplot2)
library(tidyr)
library(RColorBrewer)
# 转换为适合ggplot的数据格式
plot_data <- as.data.frame.matrix(cluster_table)
plot_data$Sample <- rownames(plot_data)

# 转换为长格式
# 定义Sample和Cluster的显示顺序
sample_order_R <- c("10820783-N", "10820783-T", "10895817-N", "10895817-T", 
                    "10979579-N", "10979579-T", "11011381-N","11011381-T")
cluster_order <- c("CD4_Treg_FOXP3",        
                   "CD4_Tn_IL7R",      
                   "CD8_Tcyc_MKI67",     
                   "CD8_Tex_CXCL13" )
data_long <- pivot_longer(plot_data, 
                          cols = -Sample, 
                          names_to = "Cluster", 
                          values_to = "Count")
# 将data_long中的因子按照指定顺序排序
data_long <- data_long %>%
  mutate(
    Sample = factor(Sample, levels = sample_order_R),
    Cluster = factor(Cluster, levels = cluster_order)
  ) %>%
  arrange(Sample, Cluster)

# 绘制堆叠柱形图（比例显示）
colors = PlotTheme$Color
p <- ggplot(data_long, aes(x = Sample, y = Count, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Proportion", fill = "Cluster") +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.x  = element_text(size = 22,face = "bold",angle = 45,hjust = 1),
    axis.text.y  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    plot.margin = margin(10,10,10,20))

p

ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig5_S6_S7/FigS7B.pdf", 
       p, bg = "white",width = 8, height = 6, dpi = 300)

#FigS8E####
# 先定义配色 + 换行拆分后的图例标签
colors <- c(
  "Hyperexpanded (0.1 < X <= 1)" = "#377eb8",
  "Large (0.01 < X <= 0.1)"       = "#ff9900",
  "Medium (0.001 < X <= 0.01)"    = "#4daf4a",
  "NA"                            = "#999999"
)

# 换行：\n 实现两行显示
new_labels <- c(
  "Hyperexpanded\n(0.1 < X <= 1)",
  "Large\n(0.01 < X <= 0.1)",
  "Medium\n(0.001 < X <= 0.01)",
  "NA"
)

p <- DimPlot(scObject,pt.size = 1.2,reduction = "umap",group.by = "cloneSize")+
  scale_color_manual(
    values = colors,
    labels = new_labels,   # 关键：替换成两行标签
    name = "cloneSize"
  )+
  big_theme
p
ggsave("项目二：三代单细胞食管癌新辅机制研究/project_new/文章主图/fig5_S6_S7/FigS7E.pdf",
       p,width = 8, height = 6, dpi = 300)

#FigS8D####
library(ggplot2)
library(dplyr)
library(ggpubr)
# 读取数据
scObject<-readRDS("项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/T cells/T cells_feature2000.rds")

# 提取meta.data
meta <- scObject@meta.data

# 筛选需要的分组
meta_filtered <- meta[meta$group_1 %in% c("N", "T"), ]

# 计算每个样本中每个亚群占该样本所有T细胞总数的比例
plot_data <- meta_filtered %>%
  group_by(Sample, group_1, cluster_annotation) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(total_cells = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total_cells * 100)
# 确保分组顺序
plot_data$group_1 <- factor(plot_data$group_1, 
                            levels = c("N", "T"))

sci_colors <- c("N" = "#0072B2",   
                "T" = "#D55E00")
# 定义比较组
comparisons <- list(
  c("N", "T")
)

subclusters <- unique(plot_data$cluster_annotation)

# 获取所有亚群名称
for (sub in subclusters) {
  # 筛选当前亚群数据
  sub_data <- plot_data[plot_data$cluster_annotation == sub, ]
  
  # 检测每个比较组是否有足够的观测值（每组至少2个样本）
  valid_comparisons <- list()
  for (comp in comparisons) {
    group1_data <- sub_data$proportion[sub_data$group_1 == comp[1]]
    group2_data <- sub_data$proportion[sub_data$group_1 == comp[2]]
    
    # 检查每组是否有至少2个观测值
    if (length(group1_data) >= 2 & length(group2_data) >= 2) {
      valid_comparisons <- c(valid_comparisons, list(comp))
    }
  }
  
  # 打印有效比较组信息
  cat("\n========== ", sub, " ==========\n")
  cat("有效比较组数量:", length(valid_comparisons), "\n")
  for (comp in valid_comparisons) {
    cat(paste(comp[1], "vs", comp[2], "- 有足够观测值\n"))
  }
  
  # 绘制箱线图
  p <- ggplot(sub_data, aes(x = group_1, y = proportion, fill = group_1)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, color = "black") +
    scale_fill_manual(values = sci_colors) +
    theme_classic() +
    labs(title = sub,
         x = "", y = "Proportion (%)") +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 18, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 18)
    )
  
  # 只有当有有效比较组时才添加p值
  if (length(valid_comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = valid_comparisons, 
                                method = "t.test",
                                label = "p.format", 
                                label.y.npc = "top", 
                                tip.length = 0.03, 
                                size = 5)
  } else {
    cat("警告: 没有足够的观测值进行任何比较\n")
  }
  
  # 保存图片
  ggsave(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/T cells//boxplot/", sub, "_boxplot.png"), 
         p, width = 6, height = 6, dpi = 300)
  ggsave(paste0("./项目二：三代单细胞食管癌新辅机制研究/project_new/02.Analysis—all/11.cell_type_subcluster_new/T cells//boxplot/", sub, "_boxplot.pdf"), 
         p, width = 6, height = 6, dpi = 300)
  
  # 打印图片
  print(p)
}

