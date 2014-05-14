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
  The format is:
 num [1:100, 1:117] 6.62 1.31 7.29 4.06 1.8 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:100] "ENSBTAG00000008124" "ENSBTAG00000009649" "ENSBTAG00000000948" "ENSBTAG00000001419" ...
  ..$ : chr [1:117] "N1178_CN_24H" "N1178_CN_2H" "N1178_CN_48H" "N1178_CN_6H" ...
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
data(log2cpm)
str(log2cpm)
log2cpm[1:5, 1:5]
}
\keyword{datasets}
