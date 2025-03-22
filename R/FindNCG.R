
# NEED TO FIX ASSAY SPECIFICATION

setGeneric(name = 'FindNCG',
           function(object,
                    unwanted_variables,
                    bio_variables,
                    assay = NULL,
                    metadata = NULL,
                    no.ncg = 1000,
                    apply_log = T,
                    sample_fraction = 0.1)
           {standardGeneric('FindNCG')})



setMethod(
  'FindNCG',
  signature = c(
    object = 'Matrix',
    unwanted_variables = 'character',
    bio_variables = 'character'),
  function(
    object,
    unwanted_variables,
    bio_variables,
    assay,
    metadata,
    no.ncg,
    apply_log,
    sample_fraction)
  {
    data <- as(object, 'CSparseMatrix')
    if(apply_log)
    {data <- log2_sparse(
      data,
      pseudocount = 1)}

    final_gene_ranks <- findNCG_default(
      matrix = data,
      unwanted_variables = unwanted_variables,
      bio_variables = bio_variables,
      metadata = metadata,
      sample_fraction = sample_fraction)

    return(names(final_gene_ranks)[1:no.ncg])

  })





"
#' @importClassesFrom Seurat Seurat
setMethod('FindNCG',
          signature = c(object = 'Seurat',
                        unwanted_variables = 'character',
                        bio_variables = 'character'),
          function(object,
                   unwanted_variables,
                   bio_variables,
                   assay,
                   metadata,
                   no.ncg,
                   apply_log,
                   sample_fraction){

            sample_num <- sample_fraction * ncol(object)
            sample_ <- sample(ncol(object), sample_num)
            message(paste0('Sampling ', sample_num,' cells from your data'))

            # Hardcoding main data matrix to assays > RNA > layers > counts
            if(apply_log){data <- log2_sparse(object@assays$RNA@layers$counts, pseudocount = 1)}
            else{data <- object@assays$RNA@layers$counts}

            data <- data[,sample_] %>% as.matrix()
            metadata <- object@meta.data[sample_,]

            message('Calculating Spearman correlation for continuous variables & F score for categorical variables')
            corr_data <- lapply(c(unwanted_variables, bio_variables),
                                function(x){find_corr(data, metadata[[x]])})

            names(corr_data) <- c(unwanted_variables, bio_variables)

            # Flip bio variables
            message('Finding genes highly affected by unwanted variables but not affected by biological variables')
            for(name in names(corr_data)){
              if(name %in% bio_variables){
                corr_data[[name]] <- corr_data[[name]] * (-1)
              }
            }
            corr_data <- do.call(cbind, corr_data) %>% as.data.frame()
            rownames(corr_data) <- rownames(object[['RNA']])
            ranks <- apply(corr_data, 2, rank) %>% as.data.frame()
            ranks$avg_expr <- Rfast::rowmeans(data)
            prod_ranks <- apply(ranks, 1, prod)
            final_gene_ranks <- prod_ranks[base::order(-prod_ranks)]
            message('FindNCG Completed!')

            return(names(final_gene_ranks)[1:no.ncg])

          })
"

#' @importClassesFrom Seurat Seurat
setMethod('FindNCG',
          signature = c(object = 'Seurat',
                        unwanted_variables = 'character',
                        bio_variables = 'character'),
          function(object,
                   unwanted_variables,
                   bio_variables,
                   assay,
                   metadata,
                   no.ncg,
                   apply_log,
                   sample_fraction){

            # Hardcoding main data matrix to assays > RNA > layers > counts
            if(apply_log){data <- log2_sparse(Seurat::GetAssay(object, assay = assay)$counts, pseudocount = 1)}
            else{data <- Seurat::GetAssay(object, assay = assay)$counts}

            final_gene_ranks <- findNCG_default(
              matrix = data,
              unwanted_variables = unwanted_variables,
              bio_variables = bio_variables, metadata = object@meta.data,
              sample_fraction = sample_fraction)

            return(names(final_gene_ranks)[1:no.ncg])

          })





#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod('FindNCG',
          signature = c(object = 'SingleCellExperiment',
                        unwanted_variables = 'character',
                        bio_variables = 'character'),
          function(object,
                   unwanted_variables,
                   bio_variables,
                   assay,
                   metadata,
                   no.ncg,
                   apply_log,
                   sample_fraction){
            if(is.null(assay))
              {stop("You're using a SingleCellExperiment but didn't specify an assay.")}
            if(apply_log)
            {data <- log2_sparse(
                SummarizedExperiment::assay(object, assay) %>%
                  as('CSparseMatrix'), pseudocount = 1)}
            else{data <- SummarizedExperiment::assay(object, assay)}

            final_gene_ranks <- findNCG_default(
              matrix = data,
              unwanted_variables = unwanted_variables,
              bio_variables = bio_variables,
              metadata = SummarizedExperiment::colData(object),
              sample_fraction = sample_fraction)

            return(names(final_gene_ranks)[1:no.ncg])

          })



findNCG_default <- function(
    matrix,
    unwanted_variables,
    bio_variables,
    metadata,
    sample_fraction){

  sample_num <- sample_fraction * ncol(matrix)
  sample_ <- sample(ncol(matrix), sample_num)
  message(paste0('Sampling ', sample_num,' cells from your data'))

  data <- matrix[,sample_] %>% as.matrix()
  metadata <- metadata[sample_,]

  message('Calculating Spearman correlation for continuous variables & F score for categorical variables')
  corr_data <- lapply(c(unwanted_variables, bio_variables),
                      function(x){find_corr(data, metadata[[x]])})
  names(corr_data) <- c(unwanted_variables, bio_variables)

  # Flip bio variables
  message('Finding genes highly affected by unwanted variables but not affected by biological variables')
  for(name in names(corr_data)){
    if(name %in% bio_variables){
      corr_data[[name]] <- corr_data[[name]] * (-1)
    }
  }
  corr_data <- do.call(cbind, corr_data) %>% as.data.frame()
  rownames(corr_data) <- rownames(object[['RNA']])
  ranks <- apply(corr_data, 2, rank) %>% as.data.frame()
  ranks$avg_expr <- Rfast::rowmeans(data)
  prod_ranks <- apply(ranks, 1, prod)
  final_gene_ranks <- prod_ranks[base::order(-prod_ranks)]
  return(final_gene_ranks)
}





