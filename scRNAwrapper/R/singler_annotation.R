#' combine seurat with singler , make auto annotation for seurat cluster
#'
#' @param seuratObject a seurat object
#' @param ref matrix of gene exression or the output of ref.make(), or you can add one of database built in singleR. default is DatabaseImmuneCellExpressionData().
#' @param ref.count logic value, if the expression values in ref are count data, FALSE means values in ref has been lognormalized.
#' @param method "cluster" or "single", used for singler annotation strategy.
#'
#' @return dataframe of seurat cluster annotaion
#' @export
#'
#' @examples
#' library(scRNAwrapper)
#' sample<-readRDS("~/work/2020.1.7/data/srt.final.RDS")
#' sample<-subset(sample,ident=0:3)
#' anno<-seu.singler(sample)
#' sample$singler<-Idents(sample)
#' levels(sample$singler)<-anno[,2]
#' DimPlot(sample,group.by = "singler",label = TRUE)
#'
seu.singler<-function(seuratObject,ref=NULL,ref.count=TRUE,method="cluster"){
  if(is.null(ref)){
    ##HumanPrimaryCellAtlasData() is not avialible at 2020/2/11 1.1.9
    ref <- DatabaseImmuneCellExpressionData()
    ref$label=ref$label.main
  }
  else if(is.matrix(ref) == TRUE){
    ref<-ref.make(ref,ref.count = ref.count)
  }
  data<-as.SingleCellExperiment(seuratObject)
  data<-logNormCounts(data)
  if(method == "cluster" || "single"){
    if(method == "cluster" ){
      cluster<-seuratObject@meta.data$seurat_clusters
    }
    else{
      cluster<-""
    }
    singler = SingleR(test= data , ref = ref ,labels = ref$label , method = method,
                      clusters = cluster, genes = "de", quantile = 0.8, fine.tune = F,
                      tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                      assay.type.test = "logcounts", assay.type.ref = "logcounts",
                      check.missing = TRUE,
                      BPPARAM = SerialParam())
  }
  ###make annotation file in the same order with seurat object
  anno_mat<-cbind(rownames(singler),singler$labels)
  anno_mat<-as.data.frame(anno_mat)
  anno_mat[,1]<-as.numeric(as.character(anno_mat[,1]))
  anno_mat<-anno_mat[order(anno_mat[,1]),]
  colnames(anno_mat)<-c("seurat.cluster","singler.annotation")
  return(anno_mat)
}
