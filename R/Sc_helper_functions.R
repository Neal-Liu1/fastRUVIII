


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


