# Fast log2 transform for sparse matrices
log2_sparse <- function(matrix, pseudocount = 1){
  require(Matrix)
  if(!is(matrix, "sparseMatrix")){stop('Your matrix is not a sparse matrix.')}
  matrix@x <- log2(matrix@x + pseudocount)
  matrix <- Matrix::drop0(matrix)
  return(matrix)
}

# Fast averaging for sparse matrices
average_celltypes_sparse <- function(expr_sparse, batch_info, celltype_info, min_cells = 3) {
  require(Matrix)
  # Create an index to identify unique (celltype, batch) pairs
  group_indices <- paste(celltype_info, batch_info, sep = ".")
  # Initialize a sparse matrix to store the sums
  gene_names <- rownames(expr_sparse)
  unique_groups <- unique(group_indices)

  sum_matrix <- Matrix(0, nrow = nrow(expr_sparse), ncol = length(unique_groups), sparse = TRUE)
  colnames(sum_matrix) <- unique_groups
  rownames(sum_matrix) <- gene_names
  # Loop over each group to sum the expression values
  for (i in seq_along(unique_groups)) {
    group <- unique_groups[i]
    group_cells <- which(group_indices == group)
    if (length(group_cells) < min_cells){next}
    sum_matrix[, i] <- Matrix::rowSums(expr_sparse[, group_cells, drop = FALSE])
  }

  # Create a vector of cell counts per group
  cell_counts <- table(group_indices)
  cell_counts <- cell_counts[colnames(sum_matrix)]

  # Divide the sums by the cell counts to get the averages
  avg_matrix <- sum_matrix %*% Diagonal(x = 1 / as.numeric(cell_counts))

  #celltype_names <- sapply(colnames(sum_matrix), function(x){unlist(strsplit(x, '[.]'))[1]})
  #colnames(avg_matrix) <- celltype_names
  avg_matrix <- as.matrix(avg_matrix)
  colnames(avg_matrix) <- colnames(sum_matrix)
  rownames(avg_matrix) <- rownames(sum_matrix)

  if(length(which(colSums(avg_matrix) == 0)) == 0)
  {return(avg_matrix)}

  avg_matrix <- avg_matrix[,-which(colSums(avg_matrix) == 0)]

  return(avg_matrix)
}


residop = function(A, B){
  decomp = qr(B)
  qr.resid(decomp, A)
}

#' Runs an expression, times it, and reports the time taken
#' @param expr An expression to run
#' @param description A character scalar, describing the operation as a noun
time_eval <- function(
    expr,
    description
){
  message(paste0('Performing ', description))
  elapsed = system.time(expr)[["elapsed"]]
  message(paste0(
    'Finished performing ',
    description,
    ' in ',
    round(elapsed, digits = 2),
    ' seconds.'
  ))
}

tological <-function (ctl, n)
{
  ctl2 = rep(FALSE, n)
  ctl2[ctl] = TRUE
  return(ctl2)
}

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
    Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
    Yrep <- RUV1(Yrep, eta, ctl, include.intercept = include.intercept)
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


setGeneric(name = 'find_corr',
           function(data, variable_vector){standardGeneric('find_corr')})

setMethod('find_corr',
          signature = c(variable_vector = 'numeric'),
          function(data, variable_vector){
            results <- Rfast::correls(y=variable_vector, x=t(data), type = 'spearman') %>% as.data.frame()
            return(as.vector(results$correlation))
          })

setMethod('find_corr',
          signature = c(variable_vector = 'character'),
          function(data, variable_vector){
            #results <- ftest(data, variable_vector)
            results <- parallel::mclapply(1:nrow(data),
                                          function(x){summary(aov(as.vector(data[x,])~variable_vector))[[1]][4][1,1]})
            return((unlist(results)))
          })



