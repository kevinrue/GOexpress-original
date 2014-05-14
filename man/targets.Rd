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
  The format is:
Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  ..@ varMetadata      :'data.frame':	7 obs. of  1 variable:
  .. ..$ labelDescription: chr [1:7] NA NA NA NA ...
  ..@ data             :'data.frame':	117 obs. of  7 variables:
  .. ..$ File     : chr [1:117] "N1178_CN_24H" "N1178_CN_2H" "N1178_CN_48H" "N1178_CN_6H" ...
  .. ..$ Sample   : chr [1:117] "N1178_CN_24H" "N1178_CN_2H" "N1178_CN_48H" "N1178_CN_6H" ...
  .. ..$ Animal   : Factor w/ 10 levels "N1178","N121",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ Treatment: Factor w/ 3 levels "CN","MB","TB": 1 1 1 1 2 2 2 2 3 3 ...
  .. ..$ Time     : Ord.factor w/ 5 levels "0H"<"2H"<"6H"<..: 4 2 5 3 4 2 5 3 4 2 ...
  .. ..$ Group    : Factor w/ 13 levels "CN_0H","CN_2H",..: 4 2 5 3 8 6 9 7 12 10 ...
  .. ..$ Timepoint: num [1:117] 24 2 48 6 24 2 48 6 24 2 ...
  ..@ dimLabels        : chr [1:2] "rowNames" "columnNames"
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 1
  .. .. .. ..$ : int [1:3] 1 1 0
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
