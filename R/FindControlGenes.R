
# NEED TO FIX ASSAY SPECIFICATION

setGeneric(name = 'FindNCG',
           function(object,
                    unwanted_variables,
                    bio_variables,
                    no.ncg = 1000,
                    apply_log = T,
                    sample_fraction = 0.1)
           {standardGeneric('FindNCG')})

#' @importClassesFrom Seurat Seurat
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
