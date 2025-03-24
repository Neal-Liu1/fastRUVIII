#' Create all Pseudo-replicate sets for a given dataset
#' @description Given a dgCMatrix, Seurat or SingleCellExperiment, compute
#' all prpc sets and outputs a single matrix.
#'
#' @param obj A dgCMatrix, Seurat or SingleCellExperiment
#' @param bio_vars A character indicating the names of the biological groups in your metadata.
#' RUVIII will assume the variation in these groups is biology and will not be removed.
#' If you do not already know some previous biology in your data,
#' labels from unsupervised clustering can also be used. If using multi-omic data,
#' the FindCorrectedMultimodalNeighbours function could be used for joint multi-omic
#' clustering (more accurate)
#' @param uv_vars A character indicating the names of the unwanted variables
#' in your metadata that you want to remove.
#' @param group_by_vars A character indicating if you want to partition each
#' unwanted variable by another variable. For example, make libsize replicates
#' for each batch separately. Must be the same length as uv_vars, and set to NA
#' if you dont want to group the corresponding uv_variable.
#' @param separate_bins_by_biology A logical indicating which continuous uv variable
#' should be binned per biological group instead of globally.
#' Must be the same length as uv_vars, and set to NA if you dont want
#' to separate the corresponding uv_variable.
#' @param assay The assay containing your RAW COUNTS if using Seurat or SingleCellExperiment
#' @param sampling_amount How much to sample for each biological group. Eg if set to 3,
#' then each celltype will have 3 replicates per batch from random sampling
#' @param metadata A DATAFRAME of metadata containing the bio and uv variables
#' (rows as cells and columns as variables)
#' @param continuous_bins Number of bins to bin a continuous uv_variable. Default is 3
#' @return A matrix with pseudo-replicates as columns and features/genes as rows.
#' @export

CreatePRPC <- function(
    obj,
    bio_vars,
    uv_vars,
    group_by_vars,
    separate_bins_by_biology,
    assay = NULL,
    sampling_amount = 3,
    metadata = NULL,
    continuous_bins = 3){

  if(class(obj) == 'dgCMatrix'){matrix <- obj}
  else if(class(obj) == 'Seurat'){
    matrix <- Seurat::GetAssayData(obj, assay = assay, layer = 'counts')
    metadata <- obj@meta.data
  }
  else if(class(obj) == 'SingleCellExperiment'){
    matrix <- SummarizedExperiment::assay(obj, assay)
    metadata <- SummarizedExperiment::colData(obj)
  }
  else{stop('Your main object is not a dgCMatrix, Seurat or SingleCelleExperiment \U1F92F')}

  if(length(bio_vars) > 1)
  {bio_labels <- do.call(paste0, metadata[bio_vars])}
  if(length(bio_vars) == 1)
  {bio_labels <- metadata[[bio_vars]]}
  if(length(unique(sapply(list(uv_vars, group_by_vars, separate_bins_by_biology), length))) != 1)
  {stop("The length of uv_vars, group_by_vars, separate_bins_by_biology must all be the same \U1F92F")}
  if(length(sampling_amount) == 1)
  {sampling_amount <- rep(sampling_amount, length(uv_vars))}
  if(length(sampling_amount) != 1 && length(sampling_amount) != length(uv_vars))
  {stop('sampling_amount must be one integer or a vector the same length as uv_vars \U1F92F')}
  if(class(obj) == 'dgCMatrix' && is.null(metadata))
  {stop("You need to provide metadata if using just a matrix as the main object \U1F92F")}
  if(class(obj) == 'Seurat' && is.null(assay))
  {stop("You inputted a Seurat object but didn't specify an assay \U1F92F")}
  if(class(obj) == 'SingleCellExperiment' && is.null(assay))
  {stop("You inputted a SingleCellExperiment object but didn't specify an assay \U1F92F")}

  # Log2 the matrix
  matrix <- log2_sparse(matrix, pseudocount = 1)

  prpc <- list()
  for(i in 1:length(uv_vars)){
    # Get the name and vector of the uv variable
    variable <- uv_vars[i]
    uv_vector <- metadata[[variable]]

    # Check if we're grouping by another variable
    if(!is.na(group_by_vars[i])){
      message(paste0(
        "Grouping for ", variable, " is set to '", group_by_vars[i],
        "'. Pseudo replicates will be created for each ", group_by_vars[i], ' separately.'))

      # Make indices for subsetting out each group
      group_vector <- metadata[[group_by_vars[i]]]
      group_indices <- lapply(unique(group_vector), function(x){which(group_vector == x)})
      names(group_indices) <- unique(group_vector)

      # Do prpc separately for each group
      for(j in names(group_indices)){
        message(paste0('Creating pseudo-replicates for ', variable, ' and ',
                       group_by_vars[i],' ', j, ' \U1F483'))
        # Give different names for each group
        grouped_variable <- paste0(variable,'_', j)
        indices <- group_indices[[j]]

        prpc[[grouped_variable]] <- createPrPc_default(
          matrix[,indices],
          bio_vector = bio_labels[indices],
          uv_vector = uv_vector[indices],
          sampling = sampling_amount[i],
          continuous_bins = continuous_bins,
          colname_suffix = grouped_variable,
          separate_bins_by_biology = separate_bins_by_biology[i],
          log_transform = FALSE)
      }
    }

    # If not grouping, just use the bio_vector and uv_vector
    if(is.na(group_by_vars[i])){
      message(paste0('Creating pseudo-replicates for ', variable,' \U1F483'))
      prpc[[variable]] <- createPrPc_default(
        matrix,
        bio_vector = bio_labels,
        uv_vector = uv_vector,
        sampling = sampling_amount,
        continuous_bins = continuous_bins,
        colname_suffix = variable,
        separate_bins_by_biology = separate_bins_by_biology[i],
        log_transform = FALSE)
    }
  }

  prpc <- do.call(cbind, prpc)
  return(prpc)
}


# Creates a prpc set for a single unwanted variable
createPrPc_default <- function(
    sparse_matrix, bio_vector, uv_vector,
    sampling = 3,
    continuous_bins = 3,
    colname_suffix = NULL,
    separate_bins_by_biology = TRUE,
    log_transform = TRUE){
  # Check bio and uv vectors
  if(!class(uv_vector) %in% c('integer', 'character', 'factor', 'numeric')){
    stop('Your uv_vector is not a continuous or categorical vector.')}
  if(!class(bio_vector) %in% c('integer', 'character', 'factor')){
    stop('Your bio_vector is not a continuous or categorical vector.')}

  # Do sampling
  sample_ <- sample(1:sampling, length(bio_vector), replace = TRUE)
  sampled_bio_labels <- paste0(bio_vector, '-', sample_)

  # Bin continuous uv
  if(class(uv_vector) == 'numeric' & separate_bins_by_biology == TRUE){
    # Create df to bin per celltype
    data <- data.frame(bio_vector = bio_vector, uv_vector = uv_vector)
    binned_data <- data %>%
      dplyr::group_by(bio_vector) %>%
      dplyr::mutate(uv_vector = dplyr::ntile(uv_vector, continuous_bins)) %>%
      dplyr::ungroup()
    uv_vector <- binned_data$uv_vector
  }
  if(class(uv_vector) == 'numeric' & separate_bins_by_biology == FALSE){
    uv_vector <- dplyr::ntile(uv_vector, continuous_bins)
  }

  if(log_transform){sparse_matrix <- log2_sparse(sparse_matrix, pseudocount = 1)}

  prpc <- average_celltypes_sparse(
    sparse_matrix,
    batch_info = uv_vector,
    celltype_info = sampled_bio_labels, min_cells = 5)
  colnames(prpc) <- sub("\\-.*$", "", colnames(prpc))
  if(!is.null(colname_suffix)){colnames(prpc) <- paste0(colnames(prpc), '_', colname_suffix)}

  return(prpc)
}
