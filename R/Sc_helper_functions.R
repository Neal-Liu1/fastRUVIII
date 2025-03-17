


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




#' Create all Pseudo-replicate sets for a given dataset
#' @description Given a dgCMatrix, Seurat or SingleCellExperiment, compute 
#' all prpc sets and outputs a single matrix.
#' 
#' @param obj A dgCMatrix, Seurat or SingleCellExperiment
#' @param bio_vars A character indicating the names of the biological variables 
#' @param uv_vars A character indicating the names of the unwanted variables you want to remove
#' @param group_by_vars A character indicating if you want to partition each 
#' unwanted variable by another variable. For example, make libsize replicates for each batch separately.
#' Must be the same length as uv_vars, and set to NA if you dont want to group the corresponding uv_variable.
#' @param separate_bins_by_biology A logical indicating which continuous uv variable 
#' should be binned per biological group instead of globally.
#' Must be the same length as uv_vars, and set to NA if you dont want to separate the corresponding uv_variable.
#' @param assay Which assay to use if using Seurat or SingleCellExperiment
#' @param sampling_amount How much to sample for each biological group. Eg if set to 3, 
#' then each celltype will have 3 replicates per batch from random sampling
#' @param metadata Metadata containing the bio and uv variables. This is a dataframe 
#' with rows as cells and columns as variables
#' @param continuous_bins Number of bins to bin a continuous uv_variable. Default is 3
#' @return A BenchmarkMetrics object with the adjusted data, their PCs and runtimes 
#' added to the corresponding slots.
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

  if(class(obj) == 'dgCMatrix'){matrix <- obj}
  else if(class(obj) == 'Seurat'){
    matrix <- Seurat::GetAssayData(obj, assay = assay, layer = 'counts')
    metadata <- obj@meta.data
  }
  else if(class(obj) == 'SingleCellExperiment'){
    matrix <- assay(obj, assay)
    metadata <- colData(obj)
  }
  else{stop('Your main object is not a dgCMatrix, Seurat or SingleCelleExperiment \U1F92F')}
  
  
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
        message(paste0('Creating pseudo-replicates for ', variable, ' and ', group_by_vars[i],' ', j, ' \U1F483'))
        # Give different names for each group
        grouped_variable <- paste0(variable,'_', j)
        indices <- group_indices[[j]]
        
        prpc[[grouped_variable]] <- createPrPc_default(
          matrix[,indices], bio_vector = bio_labels[indices], uv_vector = uv_vector[indices], 
          sampling = sampling_amount[i], continuous_bins = continuous_bins, 
          colname_suffix = grouped_variable, 
          separate_bins_by_biology = separate_bins_by_biology[i],
          log_transform = FALSE)
      }
    }
    
    # If not grouping, just use the bio_vector and uv_vector
    if(is.na(group_by_vars[i])){
      message(paste0('Creating pseudo-replicates for ', variable,' \U1F483'))
      prpc[[variable]] <- createPrPc_default(
        matrix, bio_vector = bio_labels, uv_vector = uv_vector, 
        sampling = sampling_amount[i], continuous_bins = continuous_bins, 
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




#' Run all current popular RNA normalizations
#' @description Runs Seurat logNormalize, Seurat CCA, Seurat SCT, fastMNN and Harmony, 
#' also computes their PCAs and records runtimes. All stored back into the BenchmarkMetrics object. 
#' 
#' @param obj A BenchmarkMetrics object with populated raw_data and metadata
#' @param batch_variable A string indicating which metadata column corresponds to batch information
#' @param num_pcs An integer indicating how many PCs to compute for dim reduction
#' @param ... Additional arguments for other methods (currently not used)
#' @return A BenchmarkMetrics object with the adjusted data, their PCs and runtimes added to the corresponding slots.
#' @export

setGeneric(
  'NormalizeRNA',
  function(obj, batch_variable = 'batch', num_pcs = 30, ...)
  {standardGeneric('NormalizeRNA')})

#' @describeIn NormalizeRNA S4 method for class BenchmarkMetrics
#' @export
setMethod(
  'NormalizeRNA',
  signature = c(obj = 'BenchmarkMetrics'),
  function(obj, batch_variable, num_pcs, ...)
  {
    if(is.null(obj@Raw_data))
      {stop('Your metrics object does not have any raw data \U0001F92F')}
    if(is.null(obj@Metadata))
      {stop('Your metrics object does not have any metadata \U0001F92F')}
    if(!batch_variable %in% colnames(obj@Metadata))
      {stop("The batch variable name you entered doesn't exist in the metadata \U0001F92F")}
    
    message('Computing PCA for raw data \U0001F92F')
    obj@PCs[['Raw_counts']] <- run_PCA(obj@Raw_data, pcs = num_pcs)$u
    message('Computing HVGs \U0001F92F')
    HVGs <- Seurat::FindVariableFeatures(obj@Raw_data, 
                                         selection.method = 'vst',
                                         nfeatures = 2000)[['variable']]
    
    message('Starting Seurat LogNormalize \U0001F92F')
    start <- Sys.time()
    obj@Adj_data[['Seurat_LogNormalize']] <- Seurat::NormalizeData(obj@Raw_data)
    obj@PCs[['Seurat_LogNormalize']] <- run_PCA(obj@Adj_data[['Seurat_LogNormalize']], pcs = num_pcs)$u
    obj@RunningTime[['Seurat_LogNormalize']] <- difftime(Sys.time(), start, units = 'mins')
    
    message('Starting Seurat rPCA')
    
    seurat_object <- Seurat::CreateSeuratObject(
      counts = obj@Raw_data, 
      meta.data = obj@Metadata,
      assay = 'RNA')
    seurat_object@assays$RNA@layers$counts@Dimnames[[1]] <- rownames(obj@Raw_data)
    
    seurat_list <- Seurat::SplitObject(
      seurat_object, split.by = batch_variable)
    seurat_list <- lapply(X = seurat_list, function(x) {
      x <- Seurat::NormalizeData(x) %>% 
        Seurat::FindVariableFeatures() %>%
        Seurat::ScaleData(do.scale = FALSE, do.center = TRUE) %>% 
        Seurat::RunPCA()
    }) 
    features <- Seurat::SelectIntegrationFeatures(seurat_list)
    start <- Sys.time()
    rpca_anchors <- Seurat::FindIntegrationAnchors(
      object.list = seurat_list,
      anchor.features = features,
      reduction = 'rpca')
    
    Integrated_rpca <- Seurat::IntegrateData(anchorset = rpca_anchors)
    obj@RunningTime[['Seurat_rPCA']] <- difftime(Sys.time(), start, units = 'mins')
    obj@Adj_data[['Seurat_rPCA']] <- Integrated_rpca@assays$integrated@data
    
    message('Starting Seurat SCTransform \U0001F92F')
    start <- Sys.time()
    obj@Adj_data[['SCTransform']] <- Seurat::SCTransform(
      obj@Raw_data,
      cell.attr = obj@Metadata,
      return.only.var.genes = F)$y %>% as(Class = 'dgCMatrix')
    obj@RunningTime[['SCTransform']] <- difftime(Sys.time(), start, units = 'mins')
    
    message('Starting fastMNN \U0001F92F')
    matrix <- as.matrix(obj@Adj_data[['Seurat_LogNormalize']])
    start <- Sys.time()
    obj@PCs[['fastMNN']] <- batchelor::fastMNN(
      matrix,
      batch = obj@Metadata[[batch_variable]],
      subset.row = HVGs,
      d = num_pcs)@int_colData$reducedDims$corrected
    obj@RunningTime[['fastMNN']] <- difftime(Sys.time(), start, units = 'mins')
    
    message('FastMNN finished. Starting Harmony \U0001F92F')
    start <- Sys.time()
    obj@PCs[['Harmony_LogNormalize']] <- harmony::RunHarmony(
      obj@PCs[['Seurat_LogNormalize']],
      meta_data = obj@Metadata[[batch_variable]])
    obj@RunningTime[['Harmony_LogNormalize']] <- difftime(Sys.time(), start, units = 'mins')
    
    
    # totalVI, RUVIII
    
    message('Starting PCA for the adjusted data \U0001F92F')
    obj@PCs[['SCTransform']] <- run_PCA(obj@Adj_data[['SCTransform']], pcs = num_pcs)$u
    obj@PCs[['Seurat_rPCA']] <- run_PCA(obj@Adj_data[['Seurat_rPCA']], pcs = num_pcs)$u
    
    return(obj)
    
  })



#' Run all current popular ADT normalizations
#' @description Runs Seurat CLR, DSB, ADTnorm and Harmony, 
#' also computes their PCAs and records runtimes. All stored back into the BenchmarkMetrics object. 
#' 
#' @param obj A BenchmarkMetrics object with populated raw_data and metadata
#' @param batch_variable A string indicating which metadata column corresponds to batch information
#' @param num_pcs An integer indicating how many PCs to compute for dim reduction. 
#' Default is set to 10, but you might want to change depending on how many ADTs you have.
#' @param ... Additional arguments for other methods (currently not used)
#' @return A BenchmarkMetrics object with the adjusted data, their PCs and runtimes added to the corresponding slots.
#' @export


setGeneric(
  'NormalizeADT',
  function(obj, params, batch_variable = 'batch', num_pcs = 10)
  {standardGeneric('NormalizeADT')
  }
)

setMethod(
  'NormalizeADT',
  signature = c(obj = 'BenchmarkMetrics'),
  function(obj, params, batch_variable, num_pcs)
  {
    if(is.null(obj@Raw_data))
    {stop('Your metrics object does not have any raw data')}
    if(is.null(obj@Metadata))
    {stop('Your metrics object does not have any metadata')}
    
    message('Starting DSB')
    start <- Sys.time()
    obj@Adj_data[['DSB']] <- dsb::ModelNegativeADTnorm(
      t(obj@Raw_data), 
      use.isotype.control = F)
    obj@RunningTime[['DSB']] <- difftime(Sys.time(), start, units = 'mins')
    obj@Adj_data[['DSB']] <- obj@Adj_data[['DSB']] %>% t()
    
    message('Starting CLR')
    start <- Sys.time()
    obj@Adj_data[['CLR']] <- NormalizeData(
      obj@Raw_data, 
      normalization.method = 'CLR', 
      margin = 2)
    obj@RunningTime[['CLR']] <- difftime(Sys.time(), start, units = 'mins')
    
    message('Starting PCA for raw counts & CLR')
    obj@PCs[['Raw_counts']] <- run_PCA(obj@Raw_data, pcs = num_pcs)$u
    obj@PCs[['CLR']] <- run_PCA(obj@Adj_data[['CLR']], pcs = num_pcs)$u
    
    message('Starting Harmony on PCA from CLR')
    start <- Sys.time()
    obj@PCs[['Harmony_CLR']] <- harmony::RunHarmony(
      data_mat = obj@PCs[['CLR']],
      meta_data = obj@Metadata[[batch_variable]])
    obj@RunningTime[['Harmony_CLR']] <- difftime(Sys.time(), start, units = 'mins')
    
    message('Starting ADTnorm')
    start <- Sys.time()
    obj@Adj_data[['ADTnorm']] <- ADTnorm::ADTnorm(
      cell_x_adt = as.matrix(t(obj@Raw_data)), 
      cell_x_feature = data.frame(sample = obj@Metadata[[batch_variable]]),
      save_outpath = '/vast/scratch/users/liu.ne/')
    obj@RunningTime[['ADTnorm']] <- difftime(Sys.time(), start, units = 'mins')
    obj@Adj_data[['ADTnorm']] <- obj@Adj_data[['ADTnorm']] %>% t()
    
    # NEED TO ADD TOTALVI
    
    message('Starting PCA for the adjusted data')
    obj@PCs[['DSB']] <- run_PCA(obj@Adj_data[['DSB']], pcs = num_pcs)$u
    # obj@PCs[['CLR']] <- run_PCA(obj@Adj_data[['CLR']], pcs = num_pcs)$u
    obj@PCs[['ADTnorm']] <- run_PCA(obj@Adj_data[['ADTnorm']], pcs = num_pcs)$u
    
    return(obj)
    
  }
)




#' Compute Assessments
#' 
#' @description This function computes various assessments (LISI, Silhouette, ARI) for multiple 
#' variables in a BenchmarkMetrics object.
#' @name ComputeAssessments Computes assessment metrics
#' 
setGeneric('ComputeAssessments',
           function(obj, variables, ...){standardGeneric('ComputeAssessments')})

#' @describeIn ComputeAssessments
#' 
#' @param obj A BenchmarkMetrics object.
#' @param variables A vector of variable names to compute the assessments for.
#' @param ... Additional parameters for other methods (currently not used)
#' @return A BenchmarkMetrics object with the computed assessments added.
#' 
#' @examples
#' \dontrun{
#' # Assuming `bm` is a BenchmarkMetrics object and `vars` is a vector of variable names
#' bm <- ComputeAssessments(bm, vars)
#' }
#' @export
#'
setMethod(
  'ComputeAssessments',
  signature = c(obj = 'BenchmarkMetrics'),
  function(obj,
           variables,
           ari_sampling = 0.1,
           ari_cv = 5,
           ari_resolution = 1,
           ...){
    require(parallel)
    if(!all(variables %in% colnames(obj@Metadata))){
      stop('Some/all variables you entered is not in the metadata \U0001F92F')
    }
    message('Calculating LISI \U0001F92F')
    start <- Sys.time()
    obj <- ComputeMultipleLISI(obj, variables = variables)
    message(paste0('LISI finished in ', 
                   round(difftime(Sys.time(),start, units = 'secs'), 2),
                   ' seconds \U0001F92F\nStarting Silhouette calculation \U0001F92F'))
    start <- Sys.time()
    obj <- ComputeMultipleSilhouette(obj, variables = variables)
    message(paste0('Silhouette finished in ', 
                   round(difftime(Sys.time(),start, units = 'secs'), 2),
                   ' seconds \U0001F92F'))
    "
    start <- Sys.time()
    obj <- ComputeARIs(obj, labels = variables, method = 'graph', clust_resolution = ari_resolution)
    message(paste0('ARI finished in ',
                   round(difftime(Sys.time(),start, units = 'secs'), 2),
                   ' seconds \U0001F92F'))
    "
    return(obj)
  })

#' Compute ARIs
#' 
#' @description Computes Adjusted Rand Index for a specific categorical variable, using the 
#' fastcluster's fast hierarchical clustering implementation.
#' @section Assessment: 
#' 
#' @family ComputeARIs
#' @param obj a BenchmarkMetrics object.
#' @param labels a character vector of variables you want to compute ARI for.
#' @param hclust_method a string indicating the method used for hierarchical clustering. 
#' @param distance_type a string indicating what distance to use for hierarchical clustering
#' @param sample_fraction a number from 0 to 1 indicating how much the data to subset. (Improves speed for large datasets).
#' @param ... Additional params for other methods
#' @param num_cross_validation an integer specifying how fold cross validation to do. Default is 1 (no cross validation).
#' @return a BenchmarkMetrics object with computed ARIs under the ARI slot. 
#' @export
setGeneric('ComputeARIs',
           function(obj,
                    labels, 
                    method = 'graph',
                    clust_resolution = 1,
                    hclust_method = 'complete', 
                    distance_type  = 'euclidean', 
                    sample_fraction = 1, 
                    num_cross_validation = 1,
                    neighbours = NULL,
                    ...)
             {standardGeneric('ComputeARIs')})

#' Compute ARIs for BenchmarkMetrics objects
#' @description Computes Adjusted Rand Index for a specific categorical variable, using the 
#' fastcluster's fast hierarchical clustering implementation.
#' @family ComputeARIs
#' @param obj a BenchmarkMetrics object.
#' @param labels a character vector of variables you want to compute ARI for.
#' @param hclust_method a string indicating the method used for hierarchical clustering. 
#' @param distance_type a string indicating what distance to use for hierarchical clustering
#' @param sample_fraction a number from 0 to 1 indicating how much the data to subset. (Improves speed for large datasets).
#' @param ... Additional params for other methods
#' @param num_cross_validation an integer specifying how fold cross validation to do. Default is 1 (no cross validation).
#' @return a BenchmarkMetrics object with computed ARIs under the ARI slot. 
#' @rdname ComputeARIs_BenchmarkMetrics
#' @export
setMethod('ComputeARIs', 
          signature = c(obj = 'BenchmarkMetrics'),
          function(obj,
                   labels,
                   method,
                   clust_resolution,
                   hclust_method,
                   distance_type,
                   sample_fraction,
                   num_cross_validation,
                   neighbours,
                   ...){
            if(!method %in% c('graph', 'hclust'))
              {stop("Please only specify 'graph' or 'hclust' for the method")}
            
            if(method == 'graph'){
              message("Using 'graph' method. Using Seurat FindNeighbors to construct snn graphs.")
              if(is.null(neighbours)){
              neighbours <- lapply(obj@PCs, function(x){
                Seurat::FindNeighbors(x)[['snn']]})}
              message(paste0("Performing graph partitioning."))
              
              if(length(clust_resolution) == length(neighbours)){
              clusters <- lapply(1:length(neighbours), function(i){
                Seurat::FindClusters(neighbours[[i]], resolution = clust_resolution[[i]], algorithm = 3)})
              names(clusters) <- names(neighbours)
              }
              else if(length(clust_resolution) == 1){
                clusters <- lapply(neighbours, function(x){
                  Seurat::FindClusters(x, resolution = clust_resolution, algorithm = 3)})
              }
              else(stop('clust_resolution must be a vector of same length as neighbours or a single number.'))
              
              message("Computing ARIs")
              for(label in labels){
              label_name <- labels
              labels <- obj@Metadata[[labels]]
              ARIs <- lapply(clusters, function(x){
                mclust::adjustedRandIndex(labels, x[[1]])})
              obj@ARI[[label_name]] <- ARIs}
              return(obj)
              }
            
            else if(method == 'hclust'){
            label_name <- labels
            labels <- obj@Metadata[[labels]]
            ARIs <- vector("list", length = length(obj@PCs))
            names(ARIs) <- names(obj@PCs)
            for (pc_name in names(obj@PCs)) {
              pc <- obj@PCs[[pc_name]]
              ari_values <- numeric(num_cross_validation)
              for (i in 1:num_cross_validation) {
                sample_indices <- sample(
                  seq_len(nrow(pc)),
                  size = ceiling(sample_fraction * nrow(pc)))
                sampled_pc <- pc[sample_indices, , drop = FALSE]
                sampled_labels <- labels[sample_indices]
                clusters <- cutree(
                  fastcluster::hclust(
                    d = dist(
                      sampled_pc, 
                      method = distance_type),
                    method = hclust_method),
                  k = length(unique(sampled_labels)))
                ari_values[i] <- mclust::adjustedRandIndex(
                  clusters, 
                  sampled_labels)}
              ARIs[[pc_name]] <- mean(ari_values)
            }
            obj@ARI[[label_name]] <- ARIs
            return(obj)}
            
          })


#' plot ARIs generic
#' @export
setGeneric('PlotARIs',
           function(obj, variable, title = "ARIs of different methods", ...)
           {standardGeneric('PlotARIs')})

#' plot ARIs for BenchmarkMetrics objects
#' @export
setMethod('PlotARIs',
          signature= c('BenchmarkMetrics'),
          function(obj, variable, title, ...){
            merged_ARIs_df <- data.frame(names = names(obj@ARI[[variable]]), ARI = unlist(obj@ARI[[variable]]))
            merged_ARIs_df$names <- factor(merged_ARIs_df$names, levels = obj@Algorithm)
            ggplot(merged_ARIs_df, aes(x= names , y = ARI, fill = names))+
              geom_bar(stat = "identity", colour = 'black') +
              theme_minimal() +
              labs(title = title,
                   x = NULL,
                   y = "ARI scores")+
              theme(legend.position = 'none',
                    panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
                    aspect.ratio = 1/1.1,
                    axis.line = element_line(colour = "grey45", linewidth = 0.8),
                    panel.grid.major = element_line(color = "grey96"),
                    axis.text.x = element_text(size = 10,angle = 45,hjust = 1))
          })



#' Plot HVG conservation generic
#' @export
setGeneric('PlotHVG_Conservation', 
           function(obj, batch_variable = 'batch', flavour, title = 'HVG conservation plot',
                    n_top_hvgs = 2000)
           {standardGeneric('PlotHVG_Conservation')})

# Take a metrics obj with populated adjusted data slots, 
# compute the percentage of HVGs between the raw data and all the adjusted datasets,
# outputting a bar graph of percentage for each adjusted dataset.
# if raw data has multiple batches, use common HVGs across all batches instead.
#' Plot HVG conservation for BenchmarkMetrics objects
#' @export
setMethod('PlotHVG_Conservation', 
          signature = c(obj = 'BenchmarkMetrics'),
          function(obj, batch_variable, flavour, title, n_top_hvgs)
          {
            # Split raw data by batch, then find common HVGs across all batches
            batch_vector <- obj@Metadata[[batch_variable]]
            unique_batches <- unique(batch_vector)
            batch_data_list <- lapply(unique_batches, function(batch){
              obj@Raw_data[, batch_vector == batch] %>% as.matrix()})
            batch_hvgs <- lapply(batch_data_list, function(x)
            {Seurat::FindVariableFeatures(x, selection.method = flavour, nfeatures = n_top_hvgs)[['variable']]}) %>% 
              as.data.frame()
            
            common_hvgs <- rowSums(batch_hvgs) == length(unique_batches)
            # Find HVGs for each adjusted dataset
            normalized_hvgs <- lapply(obj@Adj_data, function(x)
            {Seurat::FindVariableFeatures(as.matrix(x), selection.method = flavour, nfeatures = n_top_hvgs)[['variable']]})
            
            # Use bitwise AND to find conserved HVGs between adjusted and raw data
            hvg_conserv_scores <- lapply(normalized_hvgs, function(x)
            {mean(x & common_hvgs)/mean(common_hvgs)}) 
            
            # Plot
            df <- data.frame(method = names(hvg_conserv_scores),
                             values = unlist(hvg_conserv_scores))
            df$method <- factor(df$method, levels = obj@Algorithm)
            return(
              ggplot(df, aes(x= method, y = values, fill = method))+
                geom_bar(stat = "identity", colour = 'black') +
                theme_minimal() +
                labs(title = title, x = NULL, y = "% HVG conservation")+
                theme(legend.position = 'none',
                      panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
                      aspect.ratio = 1/1.1,
                      axis.line = element_line(colour = "grey45", linewidth = 0.8),
                      panel.grid.major = element_line(color = "grey96"),
                      axis.text.x = element_text(size = 10,angle = 45,hjust = 1)))
            
          })



# Plot umap
plot_UMAP <- function(matrix, metadata_vector, title = 'UMAP', aspect_ratio = 1/1, run_umap = T, label_is_continuous = F, 
                      continuous_var_upper_lim = NULL, alpha = 1, show_density = T){
  # taking a matrix and a vector of metadata, plot UMAP of the matrix colored by the groups in the metadata vector
  
  if(run_umap){
    umap_result = umap(t(matrix))
    df = data.frame(umap_result$layout)}
  if(!run_umap){df <- as.data.frame(matrix)}
  colnames(df) = c("UMAP1", "UMAP2")
  
  if(!is.null(continuous_var_upper_lim)){
    if(class(continuous_var_upper_lim) != 'numeric'){stop("You didn't enter a numerical value for the continous variable uppe limit. Please only enter numbers.")}
    else(metadata_vector <- lapply(metadata_vector,
                                   function(x) ifelse(x > continuous_var_upper_lim,
                                                      continuous_var_upper_lim,
                                                      x)) %>% unlist())}
  
  df$metadata = metadata_vector 
  
  if(!label_is_continuous){
    centroids <- aggregate(cbind(UMAP1, UMAP2) ~ metadata, df, mean)
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
      geom_point(size = 0.07, alpha = alpha) +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "grey30", linewidth = 0.4),
            panel.border = element_blank(),  #element_rect(colour = "grey90", fill=NA, size=0.7),
            panel.grid.major = element_blank(),  #element_line(color = "grey96"),
            panel.grid.minor = element_blank(),
            aspect.ratio = aspect_ratio,
            legend.position = "none")+
      geom_text(data = centroids, aes(label = metadata), size = 2.2, color = "grey50", hjust = 0.6, vjust = 0.6)
  }
  
  if(label_is_continuous){return(
    ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
      geom_point(size = 0.07, alpha = alpha) +
      scale_fill_viridis() +
      scale_color_viridis() +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "grey30", linewidth = 0.4),
            panel.border = element_blank(),  #element_rect(colour = "grey90", fill=NA, size=0.7),
            panel.grid.major = element_blank(),  #element_line(color = "grey96"),
            panel.grid.minor = element_blank(),
            aspect.ratio = aspect_ratio,
            legend.position = "none")
  )}
  
  if(!label_is_continuous && show_density){
    xdens <- axis_canvas(p, axis = "x")+
      geom_density(df, mapping = aes(x = UMAP1, fill = metadata_vector), color= 'grey55', alpha = 0.50, size = 0.2) +
      theme(legend.position = "none")
    
    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_density(df, mapping = aes(x = UMAP2, fill = metadata_vector), color= 'grey55', alpha = 0.50, size = 0.2) +
      theme(legend.position = "none")+
      coord_flip()
    
    
    p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
    pList <- ggdraw(p2)
    p <- pList
  }
  
  return(p)
}




# FindNCG function

# NOT ROBUST TO LOWLY EXPRESSED GENES
# try doing intersection instead of ranking
# for biology, calculate the corr within batch not globally (try for libsize as well)

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

setGeneric(name = 'FindNCG',
           function(object,
                    unwanted_variables,
                    bio_variables,
                    no.ncg = 1000,
                    apply_log = T,
                    sample_fraction = 0.1)
           {standardGeneric('FindNCG')})


setMethod('FindNCG', 
          signature = c(object = 'Seurat',
                        unwanted_variables = 'character',
                        bio_variables = 'character'),
          function(object,
                   unwanted_variables,
                   bio_variables,
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


setMethod('FindNCG', 
          signature = c(object = 'BenchmarkMetrics',
                        unwanted_variables = 'character',
                        bio_variables = 'character'),
          function(object,
                   unwanted_variables,
                   bio_variables,
                   no.ncg,
                   apply_log,
                   sample_fraction){
            
            sample_num <- sample_fraction * ncol(object@Raw_data)
            sample_ <- sample(ncol(object@Raw_data), sample_num)
            message(paste0('Sampling ', sample_num,' cells from your data'))
            
            # Hardcoding main data matrix to assays > RNA > layers > counts
            if(apply_log){data <- log2_sparse(object@Raw_data, pseudocount = 1)}
            else{data <- object@Raw_data}
            
            data <- data[,sample_] %>% as.matrix()
            metadata <- object@Metadata[sample_,]
            
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
            rownames(corr_data) <- rownames(object@Raw_data)
            ranks <- apply(corr_data, 2, rank) %>% as.data.frame()
            ranks$avg_expr <- Rfast::rowmeans(data)
            prod_ranks <- apply(ranks, 1, prod)
            final_gene_ranks <- prod_ranks[base::order(-prod_ranks)]
            message('FindNCG Completed!')
            
            return(names(final_gene_ranks)[1:no.ncg])
            
          })



# The benchmark metrics object 

setClass(
  'BenchmarkMetrics',
  slots = list(Algorithm = 'character',
               Raw_data = 'dgCMatrix',
               Metadata = 'data.frame',
               Adj_data = 'list',
               PCs = 'list',
               UMAPs = 'list',
               Latent_dims = 'list',
               RunningTime = 'list',
               Silhouette = 'list',
               ARI = 'list',
               LISI = 'list'))

#' Constructor for the BenchmarkMetrics object
BenchmarkMetrics <- function(
    Algorithm = character(),
    Raw_data,
    Metadata,
    Adj_data = list(),
    PCs= list(),
    UMAPs= list(),
    Latent_dims= list(),
    RunningTime= list(),
    Silhouette= list(),
    ARI= list(),
    LISI= list()){
  new('BenchmarkMetrics',
      Algorithm = Algorithm, 
      Raw_data = as(Raw_data, 'dgCMatrix'),
      Metadata = Metadata,
      Adj_data = Adj_data,
      PCs = PCs,
      UMAPs = UMAPs,
      Latent_dims = Latent_dims,
      RunningTime = RunningTime,
      Silhouette = Silhouette,
      ARI = ARI,
      LISI = LISI)
}





#' Plot silhouette for multiple datasets
#' @export
setGeneric("ComputeMultipleSilhouette", 
           function(obj, variables, result_format = 'per_cluster', ...) 
             standardGeneric("ComputeMultipleSilhouette"))

#' Plot silhouette for multiple datasets
#' @export
setMethod('ComputeMultipleSilhouette', 
          signature = c(obj = 'BenchmarkMetrics'),
          function(obj, variables, result_format, ...){
            
            if(!all(variables %in% colnames(obj@Metadata))){
              stop('Some or all variable names you entered is not in the metadata.')}
            
            PCs <- obj@PCs
            num_cores <- detectCores() - 2
            silhouettes <- list()
            for(variable in variables){
              silhouettes <- parallel::mclapply(
                PCs,
                function(x){
                  compute_silhouette(
                    matrix = x,
                    label_vector = obj@Metadata[[variable]],
                    result_format = result_format)
                },
                mc.cores = num_cores)
              obj@Silhouette[[variable]] <- silhouettes
            }
            return(obj)
          })

setMethod('ComputeMultipleSilhouette', 
          signature = c(obj = 'Seurat'), 
          function(obj, 
                   metrics_obj, 
                   reductions, 
                   silhouette_format, 
                   labels){
            
            reductions_list <- list()
            for (name in reductions){
              reductions_list[[name]] <- Embeddings(obj, reduction = name)
            }
            
            label_vector <- obj[[labels]][,1]
            silhouettes <- parallel::mclapply(
              reductions_list,
              function(x){
                compute_silhouette(matrix = x,
                                   label_vector = label_vector,
                                   result_format = silhouette_format)})
            
            for (x in names(silhouettes)) {
              metrics_obj@Silhouette[[x]] <- silhouettes[[x]]
            }
            
            return(metrics_obj)
          })


setMethod('ComputeMultipleSilhouette', 
          signature = c(obj = 'list'), 
          function(obj, 
                   metrics_obj, 
                   reductions, 
                   silhouette_format, 
                   labels){
            # input a list of matrices of dimensionally reduced data (rows are 
            # samples and cols are components), calculate silhouette and store in metrics obj.
            
            reductions_list <- obj
            silhouettes <- parallel::mclapply(
              reductions_list,
              function(x){
                compute_silhouette(matrix = x,
                                   label_vector = labels,
                                   result_format = silhouette_format)})
            
            for (x in names(silhouettes)) {
              metrics_obj@Silhouette[[x]] <- silhouettes[[x]]
            }
            
            return(metrics_obj)
          })
#' Plot silhouettes for multiple datasets
#' @export
setGeneric(
  'PlotMultipleSilhouette',
  function(obj, 
           variable,
           plot_type = 'violin',
           title = 'Silhouette Plot',
           aspect_ratio = 1.3,
           weighted = T,
           ...){
    standardGeneric('PlotMultipleSilhouette')})

#' Plot silhouettes for multiple datasets
#' @export
setMethod(
  'PlotMultipleSilhouette',
  signature = c(obj = 'BenchmarkMetrics'),
  function(
    obj,
    variable,
    plot_type,
    title,
    aspect_ratio,
    weighted = T,
    ...){
    if(variable %in% colnames(obj@Metadata)){
      silhouettes <- obj@Silhouette[[variable]]
      
      if(weighted)
      {silhouettes <- lapply(silhouettes, function(x){
          x$silhouette_score <- (x$silhouette_score * x$n) / sum(x$n); x})}
      
      merged_data <- do.call(rbind, silhouettes)
      merged_data$method <- rep(names(silhouettes),
                                each = length(unique(merged_data$labels)))}
    
    else{stop("Silhouette scores for this variable hasn't been computed yet.")}
    
    merged_data$method <- factor(merged_data$method, levels = obj@Algorithm)
    
    if(plot_type == 'boxplot'){
      plot <- ggplot(merged_data, aes(x=method, y=silhouette_score, fill=method))+
        geom_boxplot(outlier.shape =NA)+
        labs(y = 'Silhouette', x = NULL)+
        ggtitle(title)+
        theme_minimal() +
        theme(panel.border=element_rect(colour = "grey80", fill=NA, size=0.8),
              axis.line = element_line(colour = "grey75", linewidth = 1.1),
              panel.grid.major = element_line(color = "grey96"),
              aspect.ratio = 1/aspect_ratio,
              axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
              legend.position = 'None')}
    
    if(plot_type == 'violin'){
      plot <- ggplot(merged_data, aes(x=method, y=silhouette_score, fill=method))+
        geom_violin(trim = T)+
        labs(y = 'Silhouette', x = NULL)+
        ggtitle(title)+
        theme_minimal() +
        theme(panel.border=element_rect(colour = "grey80", fill=NA, size=0.8),
              axis.line = element_line(colour = "grey75", linewidth = 1.1),
              panel.grid.major = element_line(color = "grey96"),
              aspect.ratio = 1/aspect_ratio,
              axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
              legend.position = 'None')
      plot <- plot + geom_boxplot(width=0.1, alpha=0.5, outlier.shape =NA)
    }
    return(plot)
  })






#' Compute LISI for multiple datasets
#' @export
setGeneric('ComputeMultipleLISI',
           function(obj, reductions , variables, metadata, metrics_obj, ...)
           {standardGeneric('ComputeMultipleLISI')})


# Takes a list of reduction components (eg PCs) as dataframes or matrices (cells as rows and features as cols), 
# and a dataframe of metadata and a vector of metadata variables,
# Compute the LISI scores for every reduction for every variable
# Returns a list of dataframes
setMethod('ComputeMultipleLISI', 
          signature = c(obj = 'list',
                        reductions = NULL,
                        variables = 'character',
                        metadata = 'data.frame',
                        metrics_obj = NULL),
          function(obj, reductions , variables, metadata, metrics_obj){
            LISI_scores <- lapply(obj, function(x){compute_lisi(x,
                                                                meta_data = metadata,
                                                                label_colnames = variables)})
            return(LISI_scores)
          })


# Takes a metrics obj with reductions in the Reductions slot (reductions are list of dataframes of reduction components (eg PCs) (cells as rows and features as cols), 
# and matadata as dataframe in the metadata slot,
# Compute the LISI scores for every reduction for every variable
# Returns the metrics obj with LISI scores in the LISI slots.

#' Compute LISI for mutliple datasets
#' @export
setMethod('ComputeMultipleLISI', 
          signature = c(obj = 'BenchmarkMetrics',
                        reductions = NULL,
                        variables = 'character',
                        metadata = NULL,
                        metrics_obj = NULL),
          function(obj, reductions , variables, metadata, metrics_obj){
            if(length(obj@Metadata) == 0)
            {stop("Your Metrics obj doesn't have any metadata. Please add metadata first.")}
            if(length(obj@PCs) == 0)
            {stop("Your Metrics obj doesn't have any PCs. Please compute PCA first.")}
            
            reduction_matrix_list <- obj@PCs
            metadata <- obj@Metadata
            
            LISI_scores <- lapply(reduction_matrix_list,
                                  function(x){
                                    lisi::compute_lisi(
                                      x[,1:12],
                                      meta_data = metadata,
                                      label_colnames = variables)})
            obj@LISI <- LISI_scores
            return(obj)
          })

#' Plot LISI for multiple datasets
#' @export
setGeneric('PlotMultipleLISI',
           function(obj, variable, reductions, aspect_ratio = 1.3, 
                    title = NULL, levels = NULL)
           {standardGeneric('PlotMultipleLISI')})

#' Plot LISI for multiple datasets
#' @export
setMethod('PlotMultipleLISI',
          signature = c(obj = 'BenchmarkMetrics',
                        reductions = NULL,
                        variable = 'character'),
          function(obj, variable, reductions, aspect_ratio, title, levels){
            
            if(!(variable %in% names(obj@LISI[[1]])))
            {stop('LISI scores have not yet been computed for this variable you entered.')}
            
            scores <- unlist(lapply(obj@LISI, function(x){x[[variable]]}))
            names <- rep(names(obj@LISI), each = ncol(obj@Raw_data))
            data <- data.frame(LISI = scores, method = names)
            
            if(!is.null(levels))
            {data$method <- factor(data$method, levels = obj@Algorithm)}
            
            ggplot(data, aes(x=method, y = LISI, fill=method))+
              geom_boxplot(outlier.shape = NA)+
              labs(x = NULL) +
              ggtitle(title)+
              theme_minimal()+
              theme(legend.position = 'none',
                    panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
                    aspect.ratio = 1/aspect_ratio,
                    axis.line = element_line(colour = "grey50", linewidth = 0.7),
                    panel.grid.major = element_line(color = "grey96"),
                    axis.text.x = element_text(size = 10,angle = 45,hjust = 1))
          })




# Linear and vector correlations as S4 methods

setGeneric('PlotCorrelations', 
           function(obj,
                    variable,
                    reductions = NULL,
                    num_pcs = 10,
                    title = 'Correlation plot')
           {standardGeneric('PlotCorrelations')})

setMethod('PlotCorrelations',
          signature = c(obj = 'BenchmarkMetrics',
                        variable = 'character',
                        reductions = NULL),
          function(obj,
                   variable,
                   reductions,
                   num_pcs,
                   title){
            
            var = obj@Metadata[[variable]]
            reductions = obj@PCs
            
            if(class(var) %in% c('factor', 'character')){
              dummies <- fastDummies::dummy_cols(var)[,-1]
              cancor_scores <- sapply(reductions, function(m) {lapply(1:num_pcs, function(y) {cca <- stats::cancor(x= m[,1:y, drop=F], 
                                                                                                                   y= dummies) 
              1 - prod(1 - cca$cor^2)})}) %>% unlist()
              PCs <- rep(paste0('PC1:', 1:num_pcs), length(reductions))
              Datasets <- rep(names(reductions), each = num_pcs)
              Datasets <- factor(Datasets, levels = obj@Algorithm)
              pc_correlations <- data.frame(PCs, cancor_scores, Datasets)
              pc_correlations$PCs <- factor(pc_correlations$PCs, levels = paste0('PC1:', 1:num_pcs))
              
              return(
                ggplot(pc_correlations, aes(x = PCs, y = cancor_scores, color = Datasets, group = Datasets)) +
                  geom_line(size=0.5,alpha=0.8) +
                  geom_point(alpha=0.8) +
                  labs(x = 'Principal Components', y = 'Correlation', color = 'Method') +
                  ylim(0,1) +
                  theme_minimal() +
                  theme(axis.line = element_line(colour = "grey83", linewidth = 1.1),
                        panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
                        panel.grid.major = element_line(color = "grey96"),
                        aspect.ratio = 1/1.2) +
                  ggtitle(title)
              )}
            
            if(class(var) %in% 'numeric'){
              R_squared <- sapply(reductions, function(matrix, var, num_pcs) {
                sapply(1:num_pcs, function(y) {
                  lm_model <- summary(lm(var ~ matrix[,1:y]))$r.squared})
              }, num_pcs = num_pcs, var = var) %>% as.vector()
              PCs <- rep(paste0('PC1:', 1:num_pcs), length(reductions))
              Datasets <- rep(names(reductions), each = num_pcs) %>% factor(levels = obj@Algorithm)
              pc_correlations <- data.frame(PCs, R_squared, Datasets)
              pc_correlations$PCs <- factor(pc_correlations$PCs, levels = paste0('PC1:', 1:num_pcs))
              
              return(
                ggplot(pc_correlations, aes(x = PCs, y = R_squared, color = Datasets, group = Datasets)) +
                  geom_line(size=0.5,alpha=0.8) +
                  geom_point(alpha=0.8) +
                  labs(x = 'Principal Components', y = 'R-squared', color = 'Dataset') +
                  ylim(0,1) +
                  ggtitle(title)+
                  theme_minimal() +
                  theme(axis.line = element_line(colour = "grey83", linewidth=1.1),
                        panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
                        panel.grid.major = element_line(color = "grey96"),
                        aspect.ratio = 1/1.2)
              )}
          })


#' Plot multiple correlation plots in a panel
#' @export
setGeneric('PlotMultipleCorrelations', function(obj, variables, reductions = 'all', num_pcs = 10, titles = NULL)
{standardGeneric('PlotMultipleCorrelations')})
#' Plot multiple correlation plots in a panel
#' @export
setMethod('PlotMultipleCorrelations',
          signature = c(obj='BenchmarkMetrics',
                        variables = 'character'),
          function(obj, variables, reductions, num_pcs, titles){
            
            theme <- theme(legend.position = 'none',
                           axis.line = element_line(colour = "grey33", linewidth=0.5),
                           panel.border = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text.x = element_text(angle = 45,hjust = 1),
                           axis.title.x  = element_blank(),
                           axis.title.y = element_blank(),
                           plot.title = element_text(hjust = 0.5),
                           axis.line.y = element_blank(),
                           axis.text.y = element_blank())
            
            plots <- lapply(variables, function(x)
            {PlotCorrelations(obj, variable = x, num_pcs = num_pcs, title = x) +
                theme})
            
            if(!(is.null(titles))){
              plots <- lapply(1:length(plots), 
                              function(x){plots[[x]]$labels$title <- titles[[x]]; plots[[x]]})}
            
            plots[[1]] <- plots[[1]] +
              theme(legend.position = 'none',
                    axis.line = element_line(colour = "grey33", linewidth=0.5),
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 45,hjust = 1),
                    plot.title = element_text(hjust = 0.5),
                    axis.line.y = element_line(),
                    axis.text.y = element_text())
            
            return(ggpubr::ggarrange(plotlist = c(plots),
                                     ncol = length(plots),
                                     widths = c(rep(1,length(plots))), 
                                     common.legend = T, 
                                     legend = 'right',
                                     align = 'v')
            )
          }
)


#' Plot runtimes stored in the metrics object
#' @export
PlotRuntime <- function(obj, title = 'Runtime', log = F){
  if(!log){df <- data.frame(time = unlist(obj@RunningTime),
                   method = names(obj@RunningTime))
  ylab = 'Runtime (minutes)'}
  if(log){df <- data.frame(time = log2(unlist(obj@RunningTime)),
                     method = names(obj@RunningTime))
  ylab = 'Runtime (log2 seconds)'}
  df$method <- factor(df$method, levels = obj@Algorithm)
  return(
    ggplot(df, aes(x=method, y= time, fill = time)) +
      geom_bar(stat = "identity", colour = 'black') +
      theme_minimal() +
      labs(title = title, x = NULL, y = ylab)+
      theme(legend.position = 'none',
            panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
            aspect.ratio = 1/1.1,
            axis.line = element_line(colour = "grey45", linewidth = 0.8),
            panel.grid.major = element_line(color = "grey96"),
            axis.text.x = element_text(size = 10,angle = 45,hjust = 1)))
}




#' Filter out lowly expressed genes using Seurat's variance stabilizing transformation

setGeneric('FilterLowlyExpressedGenes', 
           function(obj, 
                    batch_variable, 
                    ngenes_per_batch = 10000, 
                    ...){
             standardGeneric('FilterLowlyExpressedGenes')})

setMethod(
  'FilterLowlyExpressedGenes',
  signature = c(obj = 'BenchmarkMetrics'),
  function(obj, 
           batch_variable,
           ngenes_per_batch, 
           n_not_detected_batch_to_permit = 1,
           ...){
    require('tidyr')
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
    obj@Raw_data <- obj@Raw_data[common_hvgs,]
    return(obj)
  })


#' Compute UMAP
#' @export
setGeneric('ComputeUMAP',
           function(obj, 
                    neighbors = 30, 
                    min_dist = 0.01, 
                    n_components = 2, 
                    parallel = T,
                    nn_method = c('annoy', 'nndescent'),
                    metric = c('cosine', 'euclidean', 'manhattan', 'correlation'),
                    method = c('umap-learn', 'uwot'),
                    n_cores = NULL,
                    ...){
             standardGeneric('ComputeUMAP')
           })

#' Compute UMAP
#' @export
setMethod(
  'ComputeUMAP',
  signature = 'BenchmarkMetrics',
  function(obj, 
           neighbors, 
           min_dist, 
           n_components, 
           parallel, 
           nn_method,
           metric,
           method,
           n_cores,
           ...){
    if(length(nn_method) > 1){nn_method <- nn_method[1]}
    if(length(method) > 1){method <- method[1]}
    if(length(metric) > 1){metric <- metric[1]}
    if(length(obj@PCs) == 0){stop("You haven't computed PCA on any of your data yet.")}
    if(parallel & is.null(n_cores)){n_cores <- parallel::detectCores()-2}
    names_string <- paste(names(obj@PCs), collapse = ", ")
    message('Starting UMAPs for: ', names_string,'.')
    if(parallel){message('Parallel enabled. Using ',n_cores,' cores.')}
    start <- Sys.time()
    
    if(method == 'uwot'){
      if(nn_method == 'nndescent'){require('rnndescent')}
    umaps <- parallel::mclapply(obj@PCs, function(x){
      uwot::umap(
        as.matrix(x), 
        n_neighbors = neighbors, 
        min_dist = min_dist, 
        n_components = n_components, 
        nn_method = nn_method,
        metric = metric)},
      mc.cores = n_cores)}
    
    if(method == 'umap-learn'){
      umaps <- lapply(obj@PCs, function(x){
        py_umap(
          matrix = x, 
          neighbors = neighbors, 
          min_dist = min_dist, 
          n_components = n_components,
          metric = metric)})}
    
    obj@UMAPs <- umaps
    message('UMAP completed in ', round(difftime(Sys.time(), start, 'mins'), digits = 2), ' minutes.')
    return(obj)
  })


# THIS REQUIRES A PYTHON ENVIRONMENT WITH UMAP-LEARN INSTALLED

"
reticulate::virtualenv_create(
  envname = '/home/users/allstaff/liu.ne/scMultiOmics-normalization/fastruv_env',
  python = '/stornext/System/data/apps/python/python-3.11.4/bin/python3.11',
  packages = c('numpy', 'umap-learn'))
reticulate::use_virtualenv(virtualenv='/home/users/allstaff/liu.ne/scMultiOmics-normalization/fastruv_env',
                           required = T)
"

py_umap <- function(matrix, neighbors = 30, min_dist = 0.01, n_components = 2, metric = 'cosine'){
  umap <- reticulate::import("umap")
  data <- reticulate::r_to_py(matrix)
  model <- umap$UMAP(n_neighbors= neighbors, min_dist= min_dist, n_components= n_components, metric = metric)
  embeddings <- model$fit_transform(data)
  embeddings <- reticulate::py_to_r(embeddings)
  return(embeddings)
}



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
    sum_matrix[, i] <- rowSums(expr_sparse[, group_cells, drop = FALSE])
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


# adds totalVI data to the metrics obj
add_totalVI <- function(metrics_obj, singlecellexperiment, assay_name, pca = T, pcs = 30){
  metrics_obj@Adj_data[['totalVI']] <- SummarizedExperiment::assay(
    singlecellexperiment, assay_name) %>% 
    as('dgCMatrix')
  
  if(pca){
    message('Starting PCA')
    metrics_obj@PCs[['totalVI']] <- run_PCA(
      metrics_obj@Adj_data[['totalVI']], pcs = pcs)$u
  }
  
  if("totalVI_runtime" %in% names(singlecellexperiment@metadata)){
    metrics_obj@RunningTime[['totalVI']] <- as.difftime(
      singlecellexperiment@metadata[["totalVI_runtime"]]/60, 
      units = 'mins')}
  return(metrics_obj)
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


