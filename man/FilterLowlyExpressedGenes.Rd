% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FilterLowlyExpressedGenes.r
\name{FilterLowlyExpressedGenes}
\alias{FilterLowlyExpressedGenes}
\alias{FilterLowlyExpressedGenes,Matrix-method}
\alias{FilterLowlyExpressedGenes,Seurat-method}
\alias{FilterLowlyExpressedGenes,SingleCellExperiment-method}
\title{Filter Lowly Expressed Genes Across Multiple Batches}
\usage{
FilterLowlyExpressedGenes(
  obj,
  batch_variable,
  assay = NULL,
  ngenes_per_batch = 10000,
  n_not_detected_batch_to_permit = 1
)

\S4method{FilterLowlyExpressedGenes}{Matrix}(
  obj,
  batch_variable,
  assay = NULL,
  ngenes_per_batch = 10000,
  n_not_detected_batch_to_permit = 1
)

\S4method{FilterLowlyExpressedGenes}{Seurat}(
  obj,
  batch_variable,
  assay = "RNA",
  ngenes_per_batch = 10000,
  n_not_detected_batch_to_permit = 1
)

\S4method{FilterLowlyExpressedGenes}{SingleCellExperiment}(
  obj,
  batch_variable,
  assay = "counts",
  ngenes_per_batch = 10000,
  n_not_detected_batch_to_permit = 1
)
}
\arguments{
\item{obj}{An object of class 'Matrix', 'Seurat', or 'SingleCellExperiment'
containing gene expression data}

\item{batch_variable}{For 'BenchmarkMetrics' and 'Seurat'/'SingleCellExperiment': A character
string specifying the column name in metadata that identifies batches.
For 'Matrix': A vector of batch identifiers with length equal to the number of columns
in the matrix}

\item{assay}{Character string specifying the assay to use (for Seurat objects only, default: "RNA")}

\item{ngenes_per_batch}{Integer specifying the number of top variable genes to select from
each batch (default: 10000)}

\item{n_not_detected_batch_to_permit}{Integer specifying the number of batches a gene can be
missing from and still be included (default: 1)}

\item{assay_name}{Character string specifying the assay to use (for SingleCellExperiment objects only, default: "counts")}
}
\value{
A filtered version of the input object, containing only the highly variable genes
}
\description{
This function identifies highly variable genes across multiple batches of cells
and filters the data to retain only these genes. It uses Seurat's variance stabilizing
transformation (VST) to identify variable genes separately in each batch, then finds their
intersection.
}
\details{
This function implements a batch-aware approach to gene filtering that helps retain
genes that are consistently variable across experimental batches. It works by:
1. Splitting the data by batch
2. Identifying the top variable genes in each batch using Seurat's VST method
3. Finding genes that are variable in at least (n_batches - n_not_detected_batch_to_permit) batches
4. Filtering the input object to retain only these genes
}
\examples{
\dontrun{

# For Matrix
batch_ids <- rep(c("batch1", "batch2"), each = 500)
filtered_matrix <- FilterLowlyExpressedGenes(
  expression_matrix,
  batch_variable = batch_ids,
  ngenes_per_batch = 3000
)

# For Seurat
filtered_seurat <- FilterLowlyExpressedGenes(
  seurat_obj,
  batch_variable = "orig.ident",
  ngenes_per_batch = 2000
)

# For SingleCellExperiment
filtered_sce <- FilterLowlyExpressedGenes(
  sce_obj,
  batch_variable = "batch",
  ngenes_per_batch = 2000,
  assay_name = "counts"
)
}

}
