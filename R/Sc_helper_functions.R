


#' Run RUVIII for benchmarkmetrics object
#' @description Runs RUVIII, the PCA for its adjusted data, and records runtime.
#'  All stored back into the BenchmarkMetrics object.
#'
#' @param obj A BenchmarkMetrics object with populated raw_data and metadata
#' @param prpc A matrix of your pseudo replicates
#' @param ncgs A logical indicating which genes are control genes
#' @param celltype A vector indicating the biological groups
#' @param M Replicate matrix
#' @param pcs Number of PCs to calculate for the adjusted data
#' @param pseudo_count The number to add to the raw counts before log transforming
#' @return A BenchmarkMetrics object with the adjusted data, their PCs and runtimes added to the corresponding slots.
#' @export
RunRUVIII <- function(obj, prpc, ncgs, k, name = 'RUVIII', pcs = 30){
  start <- Sys.time()
  results <- Sparse_RUV_III(
    log2_sparse(t(obj@Raw_data)),
    Yrep = t(prpc),
    M = ruv::replicate.matrix(colnames(prpc)),
    ctl = ncgs,
    k = k)
  obj@RunningTime[[name]] <- difftime(Sys.time(),start, units = 'mins')
  obj@Adj_data[[name]] <- t(results)
  message('\n Computing PCA')
  obj@PCs[[name]] <- run_PCA(obj@Adj_data[[name]], pcs)$u
  message('\n Completed! \U1F483')
  return(obj)
}




