\name{log2cpm}
\alias{log2cpm}
\docType{data}
\title{
Sample expression data from a RNAseq experiment.
}
\description{
A sample-by-gene matrix to test the package with. Data subset from an \emph{in vitro}
infection of bovine macrophages with different mycobacterial strains. 
Column names are sample names, and row names are ensembl gene identifiers
of the Bos taurus species. Each cell contains the log2-transformed normalised
expression level of each gene in each sample.
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
From the original dataset, all 7 genes associated with the GO term "GO:0034142"
corresponding to "toll-like receptor4 signaling pathway" were kept plus a random
93 random genes, making a total of 100 genes measured in 117 samples.
Samples originate from 10 different animals and gene expression was measured at
4 different time-points, see \code{data(targets)}.
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
