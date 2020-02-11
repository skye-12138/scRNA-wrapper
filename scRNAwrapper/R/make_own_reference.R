
#'  Transform matrix to a standared reference format for next analysis
#'
#' @param mat a matrix of gene expression value , colnames are classified variables , rownames are genes
#' @param ref.count logic value, if the expression values in mat are count data, FALSE means the mat has been lognormalized.
#'
#' @return a standard reference data for next analysis
#' @export
#'
#' @examples
#'
ref.make<-function(mat,ref.count=TRUE){
  if(ref.count == TRUE){
    ref<-SummarizedExperiment(assays = list(counts=mat))
    ref<-logNormCounts(ref)
    ref<-ref@assays@data$logcounts
    ref<-SummarizedExperiment(assays = list(logcounts=ref))
    ref$label<-ref@colData@rownames
  }
  else if(ref.count == FALSE){
    ref<-SummarizedExperiment(assays = list(logcounts=ref))
    ref$label<-ref@colData@rownames
  }
  return(ref)
}
