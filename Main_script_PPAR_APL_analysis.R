library(readxl)
library(GSVA)
library(GSEABase)
library(limma)
geneSetlist <- read_xlsx("c2.cp.kegg.v7.2.symbols.xlsx")
geneSetlist[is.na(geneSetlist)] <- ""
geneSets <- list()
for (i in 1:186) {
  genes <- data.frame(geneSetlist[geneSetlist[,i] != "", i])
  colnames(genes) <- "v1"
  geneSets[[i]] <- genes$v1
}
names(geneSets) <- colnames(geneSetlist)

# Run GSVA
exp <- read.table("low20high20_CodingGenes_TPM.txt",header=T, row.names = 1)
exp <- as.matrix(exp)
res_es <- gsva(exp, geneSets, min.sz=10, max.sz=500, 
               verbose=FALSE,
               parallel.sz=0)

# Save the results
write.csv(res_es, "186_KEGG_GSVA_score.csv")

##########################
### Limma DE analysis  ###
##########################
## Differential pathway activity analysis of the specific pathway between APL patients with high or low triglyceride levels

library(limma)
mat_DE <- read.csv("186_KEGG_GSVA_score.csv", header = T, row.names = 1)
mat <- mat_DE
design <- c(rep("Normal",20), rep("High",20))
design <- factor(design, levels = c("High", "Normal"))
design.mat <- model.matrix(~0+design)
colnames(design.mat) <- levels(design)
fit <- lmFit(mat,design.mat)
cont.matrix <- makeContrasts(High-Normal, levels = design.mat)
fit1 <- contrasts.fit(fit, cont.matrix)
fit1 <- eBayes(fit1)
DEout <- topTable(fit1,coef=1,adjust.method="BH",number = Inf, p.value = 1, lfc = 0)
write.csv(DEout, "186_pathway_limma.out.csv")

##########################
###  pathway valcano   ###
##########################

library(ggplot2)
library(RColorBrewer)
library(ggmap)
library(ggrepel)
library(dplyr)
source("plotVolcano3.R")
mat_DE <- read.csv("186_pathway_limma.out.csv", header = T)
sigs <- read.csv("selected_pathways.csv", header = T)
sigs <- as.vector(sigs$Pathway)
gene_selecte <- sigs
p <- plotVolcano(mat=mat_DE, gene.col="Pathway", x.col="logFC", 
                 y.col="adj.P.Val", labx = bquote(~Log[2]~"(fold change)"), 
                 laby = bquote(~-Log[10]~italic("Q-value")), 
                 x_cut1 = 0, x_cut2 = 0, 
                 y_cut1=0.3, y_cut2 = 0.4537, 
                 x.lim=0.3, y.lim=3.5,label=TRUE, 
                 title = "", 
                 selected_genes = gene_selecte)

ggsave(p, filename = "01_DE_pathways.pdf", width = 5, height = 5)


##################
###  heatmap   ###
##################

library(pheatmap)
t1 <- read.csv("Top20_GSVA_low20high20_pathway_matrix.csv", header = T, row.names = 1)
bk = unique(c(seq(-3,3, length=100)))

annotation_col = data.frame(
  group = factor(c(rep("LowTG",20), rep("HighTG", 20))))
rownames(annotation_col) = colnames(t1)
ann_colors = list(
  group = c(LowTG = "#7570B3", HighTG = "#1B9E77"))

pdf("02_GSVA_top20_heatmap.pdf",10,4)
pheatmap(t1,breaks = bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         scale = "row", 
         angle_col = 45,
         annotation_colors = ann_colors,
         gaps_col = 20,
         gaps_row = 10,
         annotation_col = annotation_col,
         show_colnames = F,
         cluster_cols = F, cluster_rows = F
)
dev.off()


### DE genes analysis ###
options(warn=-1) # Ignore the warning messages
library(DESeq2)
setwd("D:\\Documents\\Desktop\\pts\\shi\\2021.4.21\\plot")
tpm <- read.csv("low20high20_CodingGenes_TPM.csv", header=T)
mat <- read.csv("low20high20_count.csv",header=T, row.names = 1)
mat <- mat[rowMeans(mat[,2:ncol(mat)]) >= 2, ]
mat <- mat[,-1]
source("C:\\JP_R.packages\\JPComUse\\R\\RunDESeq2.R")

RunDESeq2(count_mat = mat, n.cont = 20, n.treat = 20,
          prefix = "DESeq2_out", sort.p = FALSE, merge.normalized = TRUE,
          normalized_mat = tpm)

##################
### GSEA plot ###
##################
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
setwd("D:\\02.files\\pts\\shi\\2021.4.21\\plot")
DEgene_mat <- read.table("DESeq2_out_with_normalized_mat.txt", header = T, sep = "\t")
DE_list <- DEgene_mat[,c(1,3)]
DE_list <- dplyr::distinct(DE_list, ID, .keep_all=TRUE)
DE_list <- na.omit(DE_list)
geneList <- DE_list$log2FoldChange
names(geneList)=DE_list$ID
geneList <- sort(geneList,decreasing = T)
geneSet <- read.gmt("PR_SE_206.gmt")
#geneSet <- read.gmt("c2.cp.kegg.v7.2.symbols.gmt")
source("D:\\03.R_functions\\JPComUse\\R\\RunGSEA.R")
p <- RunGSEA(geneList, geneSet, color="#6BAED6", save=FALSE, 
             scoreTypePos = TRUE)

p
gsea.enrich <- GSEA(geneList,TERM2GENE = geneSet, pvalueCutoff = 1,
                    scoreType = "pos", pAdjustMethod = "BH", seed = TRUE)
ggsave(p, filename = "03_GSEA_PR_PPARa_SE.pdf", width = 4, height = 3.5)

### heatmap
library(pheatmap)
library(grid)
source("C:\\JP_R.packages\\JPComUse\\R\\pheatmap_add_flag.R")

core_mat <- read.table("PR_PPARa_SE_normalized_mat.txt", header = T, row.names = 1)
core_mat <- log2(core_mat +1)
core_mat <- as.matrix(core_mat)
# core_mat <- core_mat - apply(core_mat, 1, mean)
# core_mat[which(core_mat < -5)] <- -5
# core_mat[which(core_mat > 5)] <- 5
bk = unique(c(seq(-2,2, length=100)))
annotation_col = data.frame(
  Samples = factor(c(rep("Normal_TG",20),rep("High_TG",20))))
rownames(annotation_col) = colnames(core_mat)
ann_colors = list(
  Samples = c(Normal_TG = "#7FC97F", High_TG = "#BEAED4")
)

p1 <- pheatmap(core_mat,
               breaks = bk,
               color = colorRampPalette(c("#4292C6", "#F0F0F0", "#A50F15"))(100),
               scale = "row",
               border_color = "black",
               show_rownames = TRUE,
               show_colnames = FALSE,
               #clustering_method = "ward.D2",
               annotation_colors = ann_colors,
               annotation_col = annotation_col,
               cluster_cols = FALSE, cluster_rows = T
)

####################################
#### Dot plot for GO enrichment ####
####################################

library(readxl)
library(ggplot2)
library(ggthemr)
ggthemr("fresh")
library(paletteer)
paletteer_d("RColorBrewer::Dark2")
table_go <- read.table("GO_plot.txt", header = T)
table_go$Description <- factor(table_go$Description, levels = rev(table_go$Description))
p <- ggplot(table_go, aes(x=NA, y=Description, color=negLopP, size=Z.score)) + 
  geom_point(stroke = 1) + labs(x="", y="") +
  theme_classic() +
  scale_color_continuous(low=paletteer_d("RColorBrewer::OrRd")[3], high=paletteer_d("RColorBrewer::OrRd")[8]) + 
  scale_size("Z-score",limits = c(3,11), range=c(3,5))
p
ggsave(p, filename = "Overlap_gene_GO.pdf", width = 6, height = 5)

### Correlation analysis

# library(ggplot2)
# library(ggExtra)
# 
# t1 <- read.table("4Gene_40samples.txt", header = T)
# col <- brewer.pal(3, "Dark2")
# t1$Group <- factor(t1$Group, levels = c("LowTG","HighTG"))
# p1 <- ggplot(data = t1, aes(PPARA, FLT3, color = Group)) + 
#   geom_point(size = 3, alpha = .65, shape = 16, stroke = 1) + 
#   labs(x ="PPARA expression", y = "FLT3 expression") +
#   scale_color_manual(values = col) + theme_bw() + 
#   geom_smooth(method = "lm", se = TRUE, color = col[3])
# stat <- cor.test(t1$PPARA, t1$FLT3, method = "pearson")
# anno <- c(stat$estimate, stat$p.value)
# names(anno) <- c("R", "P")
# lab <- paste(names(anno), "=",  round(anno, 3), collapse="\n")
# p2 <- p1 + annotate("text", 5, 1.5, label = lab)
# p3 <- ggMarginal(p2, type = "density", groupColour = TRUE, groupFill = TRUE)
# 
# ggsave(p3, filename = "04_PPARA_FLT3.pdf", width = 5, height = 4)
# 
# 
# ## Violinplot
# library(patchwork)
# library(ggpubr)
# 
# t1 <- read.table("4Gene_40samples.txt", header = T)
# t1$Group <- factor(t1$Group, levels = c("HighTG", "LowTG"))
# p1 <- ggviolin(t1, x="Group", y="PPARA", fill = "Group", 
#                palette = "npg", 
#                add = "boxplot", add.params = list(fill="white")) +
#   stat_compare_means() + theme(legend.position = "")
# p2 <- ggviolin(t1, x="Group", y="PPARD", fill = "Group", 
#                palette = "npg", 
#                add = "boxplot", add.params = list(fill="white")) +
#   stat_compare_means() + theme(legend.position = "") 
# p3 <- ggviolin(t1, x="Group", y="PPARG", fill = "Group", 
#                palette = "npg", 
#                add = "boxplot", add.params = list(fill="white")) +
#   stat_compare_means() + theme(legend.position = "")
# p4 <- ggviolin(t1, x="Group", y="FLT3", fill = "Group", 
#                palette = "npg", 
#                add = "boxplot", add.params = list(fill="white")) +
#   stat_compare_means() + theme(legend.position = "")
# p4
# 
# p <- p1 + p2 + p3 + p4
# ggsave(p, filename = "05_ViolinPlot.pdf", width = 5, height = 5)
# 
# library(clusterProfiler)
# t1 <- read.table("PR_PPARa_SE_gene.txt", header = T, row.names = 1)
# x <- as.character(t1$x)
# eg <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
# genelist <- eg$ENTREZID
# genelist[duplicated(genelist)]
# 
# kegg <- enrichKEGG(genelist, organism = "hsa",
#                    keyType = "kegg",
#                    pvalueCutoff = 0.05,
#                    pAdjustMethod = "BH", 
#                    minGSSize = 10,maxGSSize = 500,
#                    qvalueCutoff = 0.2,use_internal_data = FALSE)
# 
# dotplot(kegg, showCategory=30)
# 
# table_go <- read.table("dotplot.txt", header = T)
# table_go$Description <- factor(table_go$Description, levels = rev(table_go$Description))
# p <- ggplot(table_go, aes(x=NA, y=Description, color=negLog2FDR, size=Count)) + 
#   geom_point(stroke = 1) + labs(x="", y="") +
#   theme_classic() +
#   scale_color_continuous(low=paletteer_d("RColorBrewer::OrRd")[3], high=paletteer_d("RColorBrewer::OrRd")[8]) + 
#   scale_size("Count",limits = c(4,9), range=c(3,5))
# p
# ggsave(p, filename = "kegg.pdf", width = 4.5, height = 3.2)


# CIBERsortx

#####################
### heatmap plot  ###
#####################

core_mat <- read.csv("CIBERSORTx_Low20_High20.csv", header = T, row.names = 1)
library(pheatmap)
core_mat <- as.matrix(core_mat)
core_mat <- log2(core_mat + 1)
# core_mat <- core_mat - apply(core_mat, 1, mean)
# core_mat[which(core_mat < -0.15)] <- -0.15
# core_mat[which(core_mat > 0.3)] <- 0.3
bk = unique(c(seq(-2,2, length=100)))
annotation_col = data.frame(
  group = factor(c(rep("Low_TG",20),rep("High_TG",20))))
rownames(annotation_col) = colnames(core_mat)
ann_colors = list(
  group = c(Low_TG = "#D95F02", High_TG = "grey")
)
pdf("Cibersort_Heatmap.pdf",11, 5)
out <- pheatmap(core_mat,
                breaks = bk,
                cutree_cols = 2,
                color = colorRampPalette(c("#4292C6", "#F0F0F0", "#A50F15"))(100),
                scale = "row", 
                border_color = "black",
                angle_col = 90,
                show_rownames = TRUE,
                show_colnames = FALSE,
                clustering_method = "ward.D2",
                annotation_colors = ann_colors,
                annotation_col = annotation_col,
                cluster_cols = T, cluster_rows = T
)
dev.off()

# Boxplot
library(ggplot2)
library(reshape2)
library(RColorBrewer)
mat <- read.table("CIBERSORTx_Low20_High20.txt", header = T, row.names = 1)
mat1 <- melt(mat)
p <- ggplot(mat1, aes(x=variable, y= value)) +
  geom_boxplot(aes(fill = Group), position = "dodge2", outlier.alpha = 0) +
  scale_fill_manual(values = c("#FC8D62","#66C2A5")) +
  theme_classic2() + 
  labs(x = "", y = "Relative immune cell composition") +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(aes(group = Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = TRUE)
p
ggsave(p, filename = "cibersortx-boxplot.pdf", width = 10, height = 6)
