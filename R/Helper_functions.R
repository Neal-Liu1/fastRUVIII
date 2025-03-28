
# Fast log2 transform for sparse matrices
log2_sparse <- function(matrix, pseudocount = 1){
  require(Matrix)
  if(!is(matrix, "dgCMatrix")){matrix <- as(matrix, 'dgCMatrix')}
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

  sum_matrix <- Matrix::Matrix(0, nrow = nrow(expr_sparse), ncol = length(unique_groups), sparse = TRUE)
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
  avg_matrix <- sum_matrix %*% Matrix::Diagonal(x = 1 / as.numeric(cell_counts))

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



