#' Perform integrative unsupervised clustering for making PRPCs without known biological labels
#' @description given a seurat object with multi-modal data (each as an assay), do integrative unsupervised clustering
#'
#' @param seurat_obj A seurat object with multi-modal data (each as an assay)
#' @param assays A character indicating the names of the assays to use.
#' @param uv_variables A character indicating the names of the unwanted variables you do not want to be affected by in your clustering
#' @param npcs A numeric vector indicating how many PCs to use for clustering for each modality. Must match length and order of assays.
#' @param graph.name Name of the integrated wnn graph to save as.
#' @param normalization_methods Names of the quick and fast normalizations to run before clustering. Must match length and order of assays.
#' @return A Seurat object with the integrated shared neighbour graph. Clustering can be done on this graph separately to tune for resolution.
#' @export

FindCorrectedMultimodalNeighbours <- function(
    seurat_obj,
    assays,
    uv_variables,
    npcs,
    graph.name = 'harmony_wsnn',
    normalization_methods = c('LogNormalize', 'CLR')){

  if(length(assays) != length(normalization_methods))
  {stop('The number of assays and normalization methods must be the same (and in the same order)')}

  for(i in 1:length(assays)){
    Seurat::DefaultAssay(seurat_obj) <- assays[[i]]
    seurat_obj <- Seurat::NormalizeData(
      seurat_obj,
      normalization.method = normalization_methods[[i]],
      margin = 2) %>%
      Seurat::FindVariableFeatures() %>%
      Seurat::ScaleData() %>%
      Seurat::RunPCA(
        reduction.name = paste0('pca_', assays[[i]]),
        npcs = npcs[[i]])
    message('Normalization & PCA completed, removing batch effect with harmony')
    seurat_obj <- harmony::RunHarmony(
      seurat_obj, group.by.vars = uv_variables,
      reduction.use = paste0('pca_', assays[[i]]),
      reduction.save = paste0('harmony_',assays[[i]]))
  }

  reductions <- paste0('harmony_', assays)
  dims <- list(1:ncol(seurat_obj[[reductions[[1]]]]),
               1:ncol(seurat_obj[[reductions[[2]]]]))

  message('Computing multi-modal neighbours with')
  seurat_obj <- Seurat::FindMultiModalNeighbors(
    seurat_obj,
    reduction.list = reductions,
    dims.list = dims,
    snn.graph.name = graph.name)

  return(seurat_obj)
}
