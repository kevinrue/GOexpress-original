\name{plot_design}
\alias{plot_design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Univariate Effects for genes associated with a Gene Ontology
}
\description{
Successively plots univariate effects of one or more
\code{\link[base:factor]{factors}}, typically for a designed experiment
as analyzed by \code{\link[stats:aov]{aov}}().
}
\usage{
plot_design(go_id, expr_data, phenodata, result,
                       factors=colnames(Biobase::pData(phenodata)),
                       main="", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go_id}{
A Gene Ontology (GO) identifier represented by at least one gene
in the dataset.
}
  \item{expr_data}{
A sample-by-gene matrix providing expression levels of each gene
in each sample.
}
  \item{phenodata}{
An AnnotatedDataFrame from package Biobase providing phenotypic information
about the samples and grouping factors with 2 or more levels.
}
  \item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from
\code{subset_scores()}.
}
  \item{factors}{
A set of values from colnames(pData(phenodata)). Each of these values will be
represented on the X-axis to investigate its effect on the average expression
of a given genes for each level of that factor.
}
  \item{main}{
Changes the main title of the plots.
}
  \item{\dots}{
Additional arguments which will be passed on to the \code{plot.design}()
function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The output of the \code{plot.design}() function.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kevin Rue-Albrecht
}
\note{
%%  ~~further notes~~
}
\seealso{
  Method \code{\link[graphics:plot.design]{graphics::plot.design}}
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}

\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
\keyword{expression}

