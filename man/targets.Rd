\name{targets}
\alias{targets}
\docType{data}
\title{
Sample AnnotatedDataFrame from a RNAseq experiment.
}
\description{
An AnnotatedDataFrame to test the package with. Row names are sample names and
columns contain descriptive information about each sample, including
experimental factors (e.g. Treatment, Timepoint, Replicate).
}
\usage{data(targets)}
\format{
A data frame with 117 samples and 7 descriptive fields (e.g. experimental factors)
  \describe{
    \item{\code{rownames(pData(targets))}}{Unique identifier. Here, sample name.}
    \item{\code{File}}{Local filename where the RNAseq counts were obtained from.}
    \item{\code{Sample}}{Sample name.}
    \item{\code{Animal}}{Identifier of animal corresponding to the biological replicate.}
    \item{\code{Treatment}}{Infection status of the sample (CN: Control,
      MB: M. bovis, TB: M. tuberculosis)}
    \item{\code{Time}}{Time of measurement (as a factor).}
    \item{\code{Group}}{Combination of Treatment and Time factors above.}
    \item{\code{Timepoint}}{Time of measurement (as a numeric value). Useful
      for time-series plotting. See function \code{expression_plot()}.}
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
data(targets)
Biobase::pData(targets)
Biobase::pData(targets)$Group
as.character(Biobase::pData(targets)$Group)
}
\keyword{datasets}
