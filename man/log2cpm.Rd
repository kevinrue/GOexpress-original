\name{log2cpm}
\alias{log2cpm}
\docType{data}
\title{
Sample expression data from a RNAseq experiment.
}
\description{
A sample-by-gene matrix to test the package with (although no biological interpretation
may be obtained from this random subset). Column names are sample names, and row names 
are ensembl gene identifiers of the Bos taurus species. Each cell contains the
log2-transformed normalised expression level of each gene  in each sample.
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
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
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
