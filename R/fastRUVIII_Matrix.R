#' Run RUVIII for a data matrix
#' @description Runs RUVIII and outputs an adjusted matrix
#'
#' @param matrix A matrix with rows as features/genes and columns as cells.
#' @param prpc A matrix of your pseudo replicates, obtained using the CreatePRPC function.
#' @param ncgs A logical indicating which genes are control genes
#' @param k How much unwanted variation to remove. Lower k = remove less, high k = remove more.
#' @return An adjusted (log transformed) matrix with cells as columns and features as rows.
#' @export
fastRUVIII_Matrix <- function(matrix, prpc, ncgs, k){
  matrix <- as(matrix, 'CSparseMatrix')
  results <- Sparse_RUV_III(
    log2_sparse(t(matrix)),
    Yrep = t(prpc),
    M = ruv::replicate.matrix(colnames(prpc)),
    ctl = ncgs,
    k = k,
    return.info = TRUE)

  return(results)
}
