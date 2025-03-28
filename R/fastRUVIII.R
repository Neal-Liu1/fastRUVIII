#' Run fastRUVIII to remove unwanted variation
#'
#' @description
#' This function runs a fast implementation of RUVIII (Remove Unwanted Variation) on
#' various input types including matrices, Seurat objects, and SingleCellExperiment objects.
#'
#' @param object A matrix, Seurat object, or SingleCellExperiment object containing gene expression data
#' @param k Integer indicating how much unwanted variation to remove. Lower k = remove less, high k = remove more
#' @param replicates A matrix of pseudo replicates, obtained using the CreatePRPC function
#' @param control_genes A logical vector indicating which genes are control genes
#' @param assay String specifying which assay to use (for Seurat and SingleCellExperiment objects)
#' @param apply_log Logical, whether to apply log transformation to the data
#'
#' @return An adjusted expression matrix or an updated Seurat/SingleCellExperiment object with the
#'         corrected expression data in a new assay called "RUVIII"
#'
#' @importFrom methods setGeneric setMethod as
#' @importFrom Matrix t
#' @importFrom ruv replicate.matrix
#'
#' @export
#'
setGeneric(
  name = 'fastRUVIII',
  function(
    object,
    k,
    replicates = NULL,
    control_genes = NULL,
    assay = NULL,
    apply_log = TRUE)
  {standardGeneric('fastRUVIII')})

#' @rdname fastRUVIII
#' @method fastRUVIII Matrix
#' @export
setMethod(
  'fastRUVIII',
  signature = c(object = 'Matrix'),
  function(object,
           k,
           replicates,
           control_genes,
           assay,
           apply_log) {

    matrix <- as(object, 'dgCMatrix')

    # Apply log transformation if requested
    if (apply_log) {
      matrix_data <- log2_sparse(t(matrix))
    } else {
      matrix_data <- t(matrix)
    }

    results <- Sparse_RUV_III(
      Y = matrix_data,
      Yrep = t(replicates),
      M = ruv::replicate.matrix(colnames(replicates)),
      ctl = control_genes,
      k = k,
      return.info = TRUE)

    return(results)
  }
)

#' @rdname fastRUVIII
#' @method fastRUVIII Seurat
#' @importFrom Seurat GetAssayData DefaultAssay CreateAssayObject
#' @export
setMethod(
  'fastRUVIII',
  signature = c(object = 'Seurat'),
  function(object,
           k,
           replicates = NULL,
           control_genes = NULL,
           assay = NULL,
           apply_log = TRUE) {

    # Use default assay if not specified
    if (is.null(assay)) {
      assay <- DefaultAssay(object)
    }

    # Extract expression matrix
    matrix <- GetAssayData(object, assay = assay, slot = "counts")

    # Apply log transformation if requested
    if (apply_log) {
      matrix_data <- log2_sparse(t(matrix))
    } else {
      matrix_data <- t(matrix)
    }

    # Run RUVIII
    results <- Sparse_RUV_III(
      Y = matrix_data,
      Yrep = t(replicates),
      M = ruv::replicate.matrix(colnames(replicates)),
      ctl = control_genes,
      k = k,
      return.info = TRUE)

    # Create a new assay with the corrected data
    corrected_assay <- CreateAssayObject(counts = t(results$newY))

    # Add the corrected assay to the Seurat object
    object[["RUVIII"]] <- corrected_assay

    return(object)
  }
)

#' @rdname fastRUVIII
#' @method fastRUVIII SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment counts
#' @importFrom SummarizedExperiment assays assay
#' @export
setMethod(
  'fastRUVIII',
  signature = c(object = 'SingleCellExperiment'),
  function(object,
           k,
           replicates = NULL,
           control_genes = NULL,
           assay = NULL,
           apply_log = TRUE) {

    # Use default assay if not specified
    if (is.null(assay)) {
      assay <- "counts"
    }

    # Extract expression matrix
    matrix <- SummarizedExperiment::assay(object, assay)

    # Apply log transformation if requested
    if (apply_log) {
      matrix_data <- log2_sparse(t(matrix))
    } else {
      matrix_data <- t(matrix)
    }

    # Run RUVIII
    results <- Sparse_RUV_III(
      Y = matrix_data,
      Yrep = t(replicates),
      M = ruv::replicate.matrix(colnames(replicates)),
      ctl = control_genes,
      k = k,
      return.info = TRUE)

    # Add corrected assay to the SingleCellExperiment object
    assays(object)[["RUVIII"]] <- t(results$newY)

    return(object)
  }
)

Sparse_RUV_III <- function (
    Y,
    Yrep,
    M,
    ctl,
    k = NULL,
    eta = NULL,
    Ynord = NULL,
    eigVec = NULL,
    include.intercept = TRUE,
    average = FALSE,
    return.info = FALSE,
    inputcheck = FALSE)
{
  require(ruv)
  require(Matrix)
  m <- nrow(Y)
  n <- ncol(Y)
  ctl <- tological(ctl, n)
  message('check the inputs finished')
  ############# RUV-I
  time_eval({
    Y <- ruv::RUV1(Y, eta, ctl, include.intercept = include.intercept)
    Yrep <- ruv::RUV1(Yrep, eta, ctl, include.intercept = include.intercept)
  }, "RUV1 on Y and Yrep separately")
  if (ncol(M) >= m)
    newY = Y
  else if (is.null(k)) {
    ycyctinv = solve(Y[, ctl] %*% t(Y[, ctl]))
    newY = M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv) %*% Y
    fullalpha = NULL
  }
  else if (k == 0) {
    newY = Y
    fullalpha = NULL
  }
  else {
    if (is.null(Ynord) & is.null(eigVec)) {
      ############# residual operators
      message('Y0 and eigVec are not provided')

      time_eval({
        Y0 <- residop(Yrep, M)
      }, "residual operator on Yrep")

      ############# eigen vectors
      time_eval({
        eigenvector <- BiocSingular::runSVD(
          x = Y0,
          k = k,
          BSPARAM = BiocSingular::FastAutoParam(),
          center = FALSE,
          scale = FALSE
        )$u
        if (!return.info){
          rm(Y0)
          gc()
        }
      }, "eigendecomposition on Y0")
      ############# fullalpha
      time_eval({
        fullalpha <- t(eigenvector[, 1:k, drop = FALSE]) %*% Yrep
      }, "calculation of fullalpha ram")
    }
    if (!is.null(Ynord) & is.null(eigVec)) {
      ############# eigen vectors
      message('Ynord is provided')
      time_eval({
        eigenvector <- BiocSingular::runSVD(
          x = Ynord %*% t(Ynord),
          k = k,
          BSPARAM = BiocSingular::bsparam(),
          center = TRUE,
          scale = FALSE
        )$u
      }, "eigendecomposition on Ynord")
      ############# fullalpha
      time_eval({
        fullalpha <- t(eigenvector[, 1:k, drop = FALSE]) %*% Yrep
      }, "obtaining fullalpha")
    }
    if (is.null(Ynord) & !is.null(eigVec)) {
      message('eigVec is provided')
      time_eval({
        fullalpha <- t(eigVec[, 1:k, drop = FALSE]) %*% Yrep
      }, "obtaining fullalpha")
    }
    ############# alpha
    time_eval({
      alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
      ac <- alpha[, ctl, drop = FALSE]
    }, "obtaining alpha")
    ############# Standardize

    # THIS MIGHT BE A PROBLEM
    time_eval({
      Y_stand <- scale(Y[ , ctl], center=TRUE, scale=FALSE)
    }, 'standardization of the data')


    ############# W
    time_eval({
      W <- Y_stand %*% t(ac) %*% solve(ac %*% t(ac))
      rm(Y_stand)
      gc()
    }, 'calculationg of W')
    #####################################  data adjustment
    time_eval({
      newY <- Y - (W %*% alpha)
    }, "adjustment of Y")
  }
  if (average)
    newY = ((1 / apply(M, 2, sum)) * t(M)) %*% newY
  if (!return.info)
    return(newY)
  if (is.null(Ynord) & is.null(eigVec))
    return(list(
      newY = newY,
      eigenvector = eigenvector,
      W = W,
      Ynord = Y0
    ))
  if (is.null(eigVec))
    return(list(
      newY = newY,
      eigenvector = eigenvector,
      W = W
    ))
  else
    return(list(
      newY = newY,
      W = W))
}
