#' Filter Lowly Expressed Genes Across Multiple Batches
#'
#' @description This function identifies highly variable genes across multiple batches of cells
#' and filters the data to retain only these genes. It uses Seurat's variance stabilizing
#' transformation (VST) to identify variable genes separately in each batch, then finds their
#' intersection.
#'
#' @param obj An object of class 'BenchmarkMetrics', 'Matrix', 'Seurat', or 'SingleCellExperiment'
#'        containing gene expression data
#' @param batch_variable For 'BenchmarkMetrics' and 'Seurat'/'SingleCellExperiment': A character
#'        string specifying the column name in metadata that identifies batches.
#'        For 'Matrix': A vector of batch identifiers with length equal to the number of columns
#'        in the matrix
#' @param ngenes_per_batch Integer specifying the number of top variable genes to select from
#'        each batch (default: 10000)
#' @param n_not_detected_batch_to_permit Integer specifying the number of batches a gene can be
#'        missing from and still be included (default: 1)
#' @param assay Assay to use if using Seurat or SingleCellExperiment
#'
#' @return A filtered version of the input object, containing only the highly variable genes
#'
#' @details This function implements a batch-aware approach to gene filtering that helps retain
#' genes that are consistently variable across experimental batches. It works by:
#' 1. Splitting the data by batch
#' 2. Identifying the top variable genes in each batch using Seurat's VST method
#' 3. Finding genes that are variable in at least (n_batches - n_not_detected_batch_to_permit) batches
#' 4. Filtering the input object to retain only these genes
#'
#' @examples
#' \dontrun{
#' # For BenchmarkMetrics object
#' filtered_obj <- FilterLowlyExpressedGenes(
#'   benchmark_obj,
#'   batch_variable = "batch",
#'   ngenes_per_batch = 5000
#' )
#'
#' # For Matrix
#' batch_ids <- rep(c("batch1", "batch2"), each = 500)
#' filtered_matrix <- FilterLowlyExpressedGenes(
#'   expression_matrix,
#'   batch_variable = batch_ids,
#'   ngenes_per_batch = 3000
#' )
#'
#' # For Seurat
#' filtered_seurat <- FilterLowlyExpressedGenes(
#'   seurat_obj,
#'   batch_variable = "orig.ident",
#'   ngenes_per_batch = 2000
#' )
#'
#' # For SingleCellExperiment
#' filtered_sce <- FilterLowlyExpressedGenes(
#'   sce_obj,
#'   batch_variable = "batch",
#'   ngenes_per_batch = 2000,
#'   assay_name = "counts"
#' )
#' }
#'
#' @export
setGeneric('FilterLowlyExpressedGenes',
           function(obj,
                    batch_variable,
                    assay = NULL,
                    ngenes_per_batch = 10000,
                    n_not_detected_batch_to_permit = 1){
             standardGeneric('FilterLowlyExpressedGenes')})

#' @rdname FilterLowlyExpressedGenes
#' @export
setMethod(
  'FilterLowlyExpressedGenes',
  signature = c(obj = 'BenchmarkMetrics'),
  function(obj,
           batch_variable,
           assay,
           ngenes_per_batch,
           n_not_detected_batch_to_permit){

    # Input validation
    if(missing(batch_variable)) {
      stop("batch_variable must be provided")
    }

    if(!batch_variable %in% colnames(obj@Metadata)) {
      stop(paste0("batch_variable '", batch_variable, "' not found in BenchmarkMetrics object metadata"))
    }

    if(!is.numeric(ngenes_per_batch) || ngenes_per_batch <= 0) {
      stop("ngenes_per_batch must be a positive integer")
    }

    if(!is.numeric(n_not_detected_batch_to_permit) || n_not_detected_batch_to_permit < 0) {
      stop("n_not_detected_batch_to_permit must be a non-negative integer")
    }

    # Check if there's enough data
    if(nrow(obj@Raw_data) == 0 || ncol(obj@Raw_data) == 0) {
      stop("The data matrix in the BenchmarkMetrics object is empty")
    }

    message('Splitting data by batch \U0001F92F')
    batch_data_list <- lapply(
      unique(obj@Metadata[[batch_variable]]),
      function(batch){
        obj@Raw_data[, obj@Metadata[[batch_variable]] == batch]})

    message(paste0('Running VST and selecting top ',ngenes_per_batch,' genes per batch \U0001F92F'))
    batch_hvgs <- lapply(
      batch_data_list,
      function(x){
        Seurat::FindVariableFeatures(
          as.matrix(x),
          selection.method = 'vst',
          nfeatures = ngenes_per_batch)[['variable']]
      }) %>% as.data.frame()

    message(paste0('Finding intersection between the ',ngenes_per_batch,' genes per batch \U0001F92F'))
    common_hvgs <- rowSums(batch_hvgs) >= length(unique(obj@Metadata[[batch_variable]]))- n_not_detected_batch_to_permit
    message(paste0('Found ', sum(common_hvgs), ' HVGs across all batches \U0001F92F'))

    if(sum(common_hvgs) == 0) {
      warning("No common highly variable genes found across batches. Consider increasing ngenes_per_batch or n_not_detected_batch_to_permit")
      return(obj)
    }

    obj@Raw_data <- obj@Raw_data[common_hvgs,]
    return(obj)
  })

#' @rdname FilterLowlyExpressedGenes
#' @param assay Character string specifying the assay to use (for Seurat objects only, default: "RNA")
#' @export
setMethod(
  'FilterLowlyExpressedGenes',
  signature = c(obj = 'Matrix'),
  function(obj,
           batch_variable,
           assay,
           ngenes_per_batch,
           n_not_detected_batch_to_permit){

    # Input validation
    if(missing(batch_variable)) {
      stop("batch_variable must be provided")
    }

    if(is.null(colnames(obj))){
      stop("Matrix must have column names")
    }

    if(length(batch_variable) != ncol(obj)) {
      stop("Length of batch_variable must equal the number of columns in the matrix")
    }

    if(!is.numeric(ngenes_per_batch) || ngenes_per_batch <= 0) {
      stop("ngenes_per_batch must be a positive integer")
    }

    if(!is.numeric(n_not_detected_batch_to_permit) || n_not_detected_batch_to_permit < 0) {
      stop("n_not_detected_batch_to_permit must be a non-negative integer")
    }

    # Check if there's enough data
    if(nrow(obj) == 0 || ncol(obj) == 0) {
      stop("The input matrix is empty")
    }

    message('Splitting data by batch \U0001F92F')
    unique_batches <- unique(batch_variable)

    if(length(unique_batches) < 2) {
      warning("Only one batch detected. Consider using Seurat::FindVariableFeatures directly.")
    }

    batch_data_list <- lapply(
      unique_batches,
      function(batch){
        obj[, batch_variable == batch, drop = FALSE]})

    message(paste0('Running VST and selecting top ', ngenes_per_batch, ' genes per batch \U0001F92F'))
    batch_hvgs <- lapply(
      batch_data_list,
      function(x){
        if(ncol(x) < 3) {
          warning(paste0("Batch has fewer than 3 cells. Skipping VST for this batch."))
          return(character(0))
        }
        Seurat::FindVariableFeatures(
          as.matrix(x),
          selection.method = 'vst',
          nfeatures = min(ngenes_per_batch, nrow(x)))[['variable']]
      }) %>% as.data.frame()

    message(paste0('Finding intersection between the ', ngenes_per_batch, ' genes per batch \U0001F92F'))
    common_hvgs <- rowSums(batch_hvgs) >= length(unique_batches) - n_not_detected_batch_to_permit
    message(paste0('Found ', sum(common_hvgs), ' HVGs across all batches \U0001F92F'))

    if(sum(common_hvgs) == 0) {
      warning("No common highly variable genes found across batches. Consider increasing ngenes_per_batch or n_not_detected_batch_to_permit")
      return(obj)
    }

    # Return filtered matrix
    filtered_matrix <- obj[rownames(obj)[common_hvgs], , drop = FALSE]
    return(filtered_matrix)
  })

#' @rdname FilterLowlyExpressedGenes
#' @param assay Character string specifying the assay to use (for Seurat objects only, default: "RNA")
#' @export
setMethod(
  'FilterLowlyExpressedGenes',
  signature = c(obj = 'Seurat'),
  function(obj,
           batch_variable,
           assay = 'RNA',
           ngenes_per_batch,
           n_not_detected_batch_to_permit){

    # Input validation
    if(missing(batch_variable)) {
      stop("batch_variable must be provided")
    }

    if(!batch_variable %in% colnames(obj@meta.data)){
      stop(paste0("batch_variable '", batch_variable, "' not found in Seurat object metadata"))
    }

    if(!assay %in% Seurat::Assays(obj)) {
      stop(paste0("Assay '", assay, "' not found in Seurat object"))
    }

    if(!is.numeric(ngenes_per_batch) || ngenes_per_batch <= 0) {
      stop("ngenes_per_batch must be a positive integer")
    }

    if(!is.numeric(n_not_detected_batch_to_permit) || n_not_detected_batch_to_permit < 0) {
      stop("n_not_detected_batch_to_permit must be a non-negative integer")
    }

    # Get the expression matrix from the Seurat object
    count_matrix <- Seurat::GetAssayData(obj, slot = "counts", assay = assay)

    # Check if there's enough data
    if(nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
      stop("The count matrix in the Seurat object is empty")
    }

    # Get batch information from metadata
    batch_info <- obj@meta.data[[batch_variable]]

    message('Splitting data by batch \U0001F92F')
    unique_batches <- unique(batch_info)

    if(length(unique_batches) < 2) {
      warning("Only one batch detected. Consider using Seurat::FindVariableFeatures directly.")
    }

    batch_data_list <- lapply(
      unique_batches,
      function(batch){
        count_matrix[, batch_info == batch, drop = FALSE]})

    message(paste0('Running VST and selecting top ', ngenes_per_batch, ' genes per batch \U0001F92F'))
    batch_hvgs <- lapply(
      batch_data_list,
      function(x){
        if(ncol(x) < 3) {
          warning(paste0("Batch has fewer than 3 cells. Skipping VST for this batch."))
          return(character(0))
        }
        Seurat::FindVariableFeatures(
          as.matrix(x),
          selection.method = 'vst',
          nfeatures = min(ngenes_per_batch, nrow(x)))[['variable']]
      }) %>% as.data.frame()

    message(paste0('Finding intersection between the ', ngenes_per_batch, ' genes per batch \U0001F92F'))
    common_hvgs <- rowSums(batch_hvgs) >= length(unique_batches) - n_not_detected_batch_to_permit
    message(paste0('Found ', sum(common_hvgs), ' HVGs across all batches \U0001F92F'))

    if(sum(common_hvgs) == 0) {
      warning("No common highly variable genes found across batches. Consider increasing ngenes_per_batch or n_not_detected_batch_to_permit")
      return(obj)
    }

    # Subset the Seurat object
    genes_to_keep <- rownames(count_matrix)[common_hvgs]
    filtered_obj <- Seurat::subset(obj, features = genes_to_keep)
    return(filtered_obj)
  })

#' @rdname FilterLowlyExpressedGenes
#' @param assay_name Character string specifying the assay to use (for SingleCellExperiment objects only, default: "counts")
#' @export
setMethod(
  'FilterLowlyExpressedGenes',
  signature = c(obj = 'SingleCellExperiment'),
  function(obj,
           batch_variable,
           assay = 'counts',
           ngenes_per_batch,
           n_not_detected_batch_to_permit){

    # Input validation
    if(missing(batch_variable)) {
      stop("batch_variable must be provided")
    }

    if(!batch_variable %in% colnames(SummarizedExperiment::colData(obj))){
      stop(paste0("batch_variable '", batch_variable, "' not found in colData"))
    }

    if(!assay_name %in% SummarizedExperiment::assayNames(obj)){
      stop(paste0("assay '", assay_name, "' not found in SingleCellExperiment object"))
    }

    if(!is.numeric(ngenes_per_batch) || ngenes_per_batch <= 0) {
      stop("ngenes_per_batch must be a positive integer")
    }

    if(!is.numeric(n_not_detected_batch_to_permit) || n_not_detected_batch_to_permit < 0) {
      stop("n_not_detected_batch_to_permit must be a non-negative integer")
    }

    # Get the expression matrix
    count_matrix <- SummarizedExperiment::assay(obj, assay_name)

    # Check if there's enough data
    if(nrow(count_matrix) == 0 || ncol(count_matrix) == 0) {
      stop("The count matrix in the SingleCellExperiment object is empty")
    }

    # Get batch information
    batch_info <- SummarizedExperiment::colData(obj)[[batch_variable]]

    message('Splitting data by batch \U0001F92F')
    unique_batches <- unique(batch_info)

    if(length(unique_batches) < 2) {
      warning("Only one batch detected. Consider using Seurat::FindVariableFeatures directly.")
    }

    batch_data_list <- lapply(
      unique_batches,
      function(batch){
        count_matrix[, batch_info == batch, drop = FALSE]})

    message(paste0('Running VST and selecting top ', ngenes_per_batch, ' genes per batch \U0001F92F'))
    batch_hvgs <- lapply(
      batch_data_list,
      function(x){
        if(ncol(x) < 3) {
          warning(paste0("Batch has fewer than 3 cells. Skipping VST for this batch."))
          return(character(0))
        }
        Seurat::FindVariableFeatures(
          as.matrix(x),
          selection.method = 'vst',
          nfeatures = min(ngenes_per_batch, nrow(x)))[['variable']]
      }) %>% as.data.frame()

    message(paste0('Finding intersection between the ', ngenes_per_batch, ' genes per batch \U0001F92F'))
    common_hvgs <- rowSums(batch_hvgs) >= length(unique_batches) - n_not_detected_batch_to_permit
    message(paste0('Found ', sum(common_hvgs), ' HVGs across all batches \U0001F92F'))

    if(sum(common_hvgs) == 0) {
      warning("No common highly variable genes found across batches. Consider increasing ngenes_per_batch or n_not_detected_batch_to_permit")
      return(obj)
    }

    # Subset the SingleCellExperiment object
    genes_to_keep <- rownames(count_matrix)[common_hvgs]
    filtered_obj <- obj[genes_to_keep, ]
    return(filtered_obj)
  })
