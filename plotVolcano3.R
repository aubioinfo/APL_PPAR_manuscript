#' Function to visualise a valcanoplot
#'
#' @param mat The input file for visualization, which is commonly derived from the differentially expressed analysis using DESeq2, edgeR or limma packages. The mat reqiures a column "Gene", "logFC", "P.Value"
#' @param x_cut The cut-off value of x-tile, e.g., Log2fold-change obtained from DE analysis
#' @param y_cut The cut-off value of y-tile, e.g., P.values
#' @param x.col The x column used for plot.
#' @param y.col The y column used for plot.
#' @param gene.col The column of the genes.
#' @param label Whether to show labels of genes of interest. The default value is TRUE.
#' @param selected_genes A string indicates selected genes for plot.
#' @param x.lab The x lab content. Used when the label is TRUE.
#' @param y.lab The y lab content.
#' @param x.lim The max/min value of x.
#' @return a ggplot object (valcanoplot)
#' @export
#'
#' @examples plotVolcano(mat=m1, gene.col=Symbol, x.col="log2FoldChange", y.col="p.adj",
#' x_cut = 1, y_cut=0.05, x.lab="log2FC", y.lab="FDR", x.lim=5,
#' label=TRUE, selected_genes=genes)
#'
library(ggrastr)
library(ggrepel)
plotVolcano <- function(mat, gene.col, x.col, y.col, x_cut1, x_cut2,
                        y_cut1, y_cut2, labx, laby,
                        x.lim, y.lim, title="", label=FALSE, selected_genes=NULL) {
  colnames(mat)[colnames(mat) %in% x.col] = 'x.col'
  colnames(mat)[colnames(mat) %in% y.col] = 'y.col'
  colnames(mat)[colnames(mat) %in% gene.col] = 'Gene'

  mat$label <- ""
  genes_label <- intersect(mat$Gene, selected_genes)
  mat[which(mat$Gene %in% selected_genes), "label"] <- as.character(mat[which(mat$Gene %in% genes_label), "Gene"])

  xmin <- (range(mat$x.col)[1]- (range(mat$x.col)[1] + x.lim))
  xmax <- (range(mat$x.col)[1]+ (x.lim - range(mat$x.col)[1]))
  ymin <- 0
  ymax <- y.lim

  mycol <- c("#223D6C","#D20A13","#088247","darkgreen","chocolate4","blueviolet","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

  n1 <- length(mat[, 1])
  cols <- rep("grey", n1)
  names(cols)<- rownames(mat)
  selected_genes <- selected_genes
  #

  #set colors
  cols[mat$y.col < y_cut1 & mat$x.col > x_cut1]<- "#9970AB"
  cols[mat$y.col < y_cut2 & mat$x.col > x_cut2]<- "#762A83"
  cols[mat$y.col < y_cut1 & mat$x.col < -x_cut1]<- "#5AAE61"
  cols[mat$y.col < y_cut2 & mat$x.col < -x_cut2]<- "#1B7837"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  mat$color_transparent <- color_transparent

  # set sizes
  n1 <- length(mat[, 1])
  size <- rep(1, n1)

  size[mat$y.col < y_cut1 & mat$x.col > x_cut1]<- 2
  size[mat$y.col < y_cut2 & mat$x.col > x_cut2]<- 4
  size[mat$y.col < y_cut1 & mat$x.col < -x_cut1]<- 2
  size[mat$y.col < y_cut2 & mat$x.col < -x_cut2]<- 4

  p <- ggplot(data=mat, aes(x.col, -log10(y.col))) +
    geom_point_rast(alpha = 0.6, size = size, colour = mat$color_transparent) +
    labs(x=labx, y=laby, title=title) +
    ylim(c(ymin,ymax)) +
    xlim(c(xmin,xmax)) +
    geom_vline(xintercept = c(-x_cut1, x_cut1), color="grey40",
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(y_cut1), color="grey40",
               linetype="longdash", lwd = 0.5) +
    geom_vline(xintercept = c(-x_cut2, x_cut2), color="grey40",
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(y_cut2), color="grey40",
               linetype="longdash", lwd = 0.5)
    theme_classic(base_size = 12)

  p1 <- p + geom_label_repel(aes(label = label),
                     max.overlaps = 100000, size = 2.5,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"),
                     segment.color = "black", segment.size = 0.3)
    # geom_point_rast(data = selected_genes, alpha = 1, size = 4.6, shape = 1,
    #            stroke = 1,
    #            color = "black") +
    # scale_color_manual(values = mycol) +
    # geom_text_repel(data = selected_genes,
    #                 show.legend = FALSE,
    #                 size = 5, box.padding = unit(0.35, "lines"),
    #                 point.padding = unit(0.3, "lines")) +
    # guides(color=guide_legend(title = NULL))

  if(label){
      return(p1)
      }else{
        return(p)
      }
    }

