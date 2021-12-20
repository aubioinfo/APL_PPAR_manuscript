#' Function to visualise a K-M survival curve
#'
#' @param count_mat A matrix of RNA-Seq raw count data with rows for genes and columns for samples.
#' @param n.cont The number of samples in the control group.
#' @param n.treat The number of samples in the treatment group.
#' @param prefix A string value to indicate the prefix of output file.
#' @param sort.p A logic value to indicate whether to sort adjusted p value for output table; TRUE by default.
#' @param merge.normalized A logic value to indicate whether to merge normalized data (tpm, fpkm) for output table; TRUE by default.
#' @param normalized_mat If TRUE for the merge.normalized, provide the normalized matrix.
#'
#' @return A .txt file storing the results of differential expression analysis by DESeq2
#' @export
#'
#' @examples RunDESeq2(count_mat, n.cont, n.treat, prefix, sort.p=TRUE, merge.normalized=TRUE, normalized_mat)
#'

RunDESeq2 <- function(count_mat, n.cont, n.treat, prefix, sort.p=TRUE, merge.normalized=TRUE, normalized_mat){
  count_mat <- round(count_mat)
  condition <- factor(c(rep("control", n.cont),rep("treatment", n.treat)))
  coldata <- data.frame(row.names=colnames(count_mat),condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat,
                                        colData = coldata,
                                        design = ~condition)
  dds$condition <- relevel(dds$condition,ref = "control")
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("condition","treatment","control"))
  if(sort.p) {
    resData <- as.data.frame(res[order(res$padj),])
  } else {
    resData <- as.data.frame(res)
  }
  resData$ID <- rownames(resData)
  resData <- resData[,c("ID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
  resData$FoldChange <- 2^resData$log2FoldChange
  if(merge.normalized) {
    outfile <- paste0(prefix, "_with_normalized_mat.txt")
    tpm <- normalized_mat
    colnames(tpm)[1] <- "ID"
    resData <- merge(resData, tpm, by="ID")
  }else {
    outfile <- paste0(prefix, "_statistic_only.txt")
  }

  write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
}





