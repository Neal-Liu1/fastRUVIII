


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


