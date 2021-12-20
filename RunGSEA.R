#' Function to perfom GSEA analysis and plot
#'
#' @param gene.List A pre-ranked gene list
#' @param gene.Set The gene sets used for enrichment analysis
#' @param color The color of the line, the default is "#6BAED6"
#' @param save A logic value to indicate whether to save the plot, the default is TRUE
#' @param prefix A string value to indicate the prefix of output file. Used only the save=TRUE
#' @param width The width of the plot
#' @param height The height of the plot
#'
#' @return A .pdf file for enrichment plot
#' @export
#'
#' @examples RunGSEA(geneList, geneSet, color="#6BAED6", save=TRUE, prefix="ISR",width=5, height=4)
#'
#'

RunGSEA <- function(gene.List, gene.Set, color="#6BAED6", save=TRUE,
                    prefix, width=5, height=4, scoreTypePos = TRUE){
  set.seed(123456)
  if (scoreTypePos) {
  gsea.enrich <- GSEA(gene.List,TERM2GENE = gene.Set, pvalueCutoff = 1,
                      scoreType = "pos", pAdjustMethod = "BH", seed = TRUE)
  anno <- gsea.enrich[1, c("NES", "pvalue", "p.adjust")]
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  }else{
    gsea.enrich <- GSEA(gene.List,TERM2GENE = gene.Set, pvalueCutoff = 1,
                        pAdjustMethod = "BH", seed = TRUE)
    anno <- gsea.enrich[1, c("NES", "pvalue", "p.adjust")]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  }
  p <- gseaplot2(gsea.enrich, 1, color=color, title=gsea.enrich$Description[1],
                 ES_geom="line", rel_heights = c(1.5, 0.5, 1),
                 subplots = 1:3, pvalue_table = F) +
    annotate("text", 0.8, 0.85, label = lab)
  if (save) {
    ggsave(p, filename = paste0(prefix, "_GSEA.pdf"), width = width, height = height)
  }else{

    return(p)

    }
}
