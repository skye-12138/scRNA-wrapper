
#' make quick heatmap plot for a seurat object
#'
#' @param seuratObject a Seurat object, need to heve done the cluster step.
#' @param anno.col data.frame, annotation for clusters, rows are clusters same as seuratObject, columns are classification variables
#' @param gene.mark a series of genes, which will significantly marked in the final heatmap.
#' @param marker a series of genes, shows the row value of heatmap.
#' @param duplicate logic value. can the marker in all groups be duplicated, accept TRUE/FALSE. if FALSE, the program will maintain the only marker with maixmum avg_logFC.
#' @param top filter paragram for default marker generate.
#'
#' @return a heatmap of interesting markers in all clusters
#' @export plot.heatmap
#'
#' @examples
#'
plot.heatmap<-function(seuratObject,anno.col=NULL,gene.mark=NULL,marker=NULL,duplicate=FALSE,top=20){
  data<-as.SingleCellExperiment(seuratObject)
  data<-logNormCounts(data)
  if(is.null(marker)){
    marker <- FindAllMarkers(seuratObject,
                             only.pos = TRUE)
    marker <- marker %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    marker$order<-1:length(marker$gene)
    if(duplicate == FALSE){
      marker.dup<-names(which(table(marker$gene) >1))
      for(i in marker.dup){
        dup<-marker[which(marker$gene == i),]
        del.col<-dup[-which.max(dup$avg_logFC),"order"]
        marker<-marker[-which(marker$order %in% unlist(del.col)),]
      }
    }
  }
  else{
    marker<-as.data.frame(marker)
    colnames(marker)<-gene
  }
  ##extract marker gene lognormcounts
  marker_order<-c()
  for (a in 1:length(marker$gene)) {
    marker_order[a]<-which(rownames(data) == marker$gene[a])
  }
  marker.lognorcount<-data[marker_order,]
  ##generate average expression in each cluster for each gene
  table<-matrix(ncol = length(unique(marker.lognorcount$ident)),nrow = length(marker$gene))
  num=1
  for(i in unique(marker.lognorcount$ident)){
    sub<-marker.lognorcount[,which(marker.lognorcount$ident == i)]
    sub<-sub@assays@data$logcounts
    a<-rowMeans(as.data.frame(sub))
    table[,num]<-a
    num<-num+1
  }
  colnames(table)<-unique(marker.lognorcount$ident)
  rownames(table)<-rownames(marker.lognorcount)
  ##make scaled data for plot
  mat<-apply(table,1,scale)
  rownames(mat)<-colnames(table)
  mat<-mat[order(as.numeric(rownames(mat))),]
  mat<-t(mat)
  rownames(mat)<-rownames(table)
  ##make annotation
  if(is.null(gene.mark) == FALSE){
    gene_pos<-which(rownames(mat) %in% gene.mark)
    gene.mark<-rownames(mat)[gene_pos]
    gene.mark<-rowAnnotation(mark_gene=anno_mark(at= gene_pos,labels=gene.mark))
  }
  if(is.null(anno.col) == FALSE){
    anno.col <- HeatmapAnnotation(
      df=anno.col
    )
  }
  Heatmap(mat,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = TRUE,show_row_names = FALSE,top_annotation=anno.col,right_annotation=gene.mark,name = "level")
}
