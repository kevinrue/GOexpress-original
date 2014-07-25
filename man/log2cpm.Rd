\name{log2cpm}
\alias{log2cpm}
\docType{data}
\title{
Sample expression data from a RNAseq experiment.
}
\description{
A gene-by-sample matrix to test the package with. Subset of data from an
\emph{in vitro} stimulation of bovine macrophages with different mycobacterial
strains. Column names are sample names, and row names are ensembl gene 
identifiers of the \emph{Bos taurus} species. Each cell contains the
log2-transformed normalised expression level of each gene in each sample.
}
\usage{data(log2cpm)}
\format{
A matrix of expression levels for 100 genes measured in 117 samples.
  \describe{
    \item{\code{colnames(log2cpm)}}{Sample names}
    \item{\code{rownames(log2cpm)}}{Ensembl gene identifiers (Bos taurus)}
  }
}
\details{
Gene expression was measured in poly-A purified strand-specific RNA libraries using the RNA-Sequencing Illumina(R) HiSeq(R) 2000 platform as paired-end 2 × 90 nucleotide reads. Raw reads from pooled RNA libraries were first deconvoluted according to sample-specific nucleotide barcodes. Read pairs containing adapter sequence in either read mate were discarded, and similarly read pairs of low overall quality in either mate were also discarded. Paired-end reads from each filtered individual library were aligned to the \emph{Bos taurus} reference genome (\emph{B. taurus} UMD3.1.71 genome release) using the STAR aligner software. For each library, raw counts for each gene based on sense strand data were obtained using the featureCounts software from the Subread package. The featureCounts parameters were set to unambiguously assign uniquely aligned paired-end reads in a stranded manner to the exons of genes within the \emph{Bos taurus} reference genome annotation (\emph{B. taurus} UMD3.1.71 genome annotation). The gene count outputs were further processed using the edgeR Bioconductor package.
 
The gene quantitation pipeline within the edgeR package was customised to: (1) filter out all bovine rRNA genes; (2) filter out genes displaying expression levels below the minimally-set threshold of one count per million [CPM] in at least ten individual libraries (number of biological replicates); (3) calculate normalisation factors for each library using the trimmed mean of M-values method; (4) log2-transform CPM values based on the normalised library size.

To generate this test data subset, we extracted 100 genes from the original dataset of 12,121 genes. All 7 genes associated with the GO term "GO:0034142" (i.e. "toll-like receptor 4 signaling pathway") were kept, plus another random 93 random genes, making a total of 100 genes measured in 117 samples. Samples include all 10 biological replicates collected at four different time-points, see \code{data(targets)}. The TLR4 pathway was found in the full dataset as the top-ranking biological pathway discriminating the different mycobacterial infections (unpublished observations).
}
\source{
  Publication in review process.
}
\examples{
# Load the data
data(log2cpm)

# Structure of the data
str(log2cpm)

# Dimensions (rows, columns) of the data
dim(log2cpm)

# Subset of first 5 rows and columns
log2cpm[1:5, 1:5]
}

\keyword{datasets}
\keyword{GOexpress}
