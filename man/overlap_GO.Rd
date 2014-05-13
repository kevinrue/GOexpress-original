\name{overlap_GO}
\alias{overlap_GO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shared genes between a list of GO terms.
}
\description{
Given a list of two to five GO terms, overlap_GO() will produce a
Venn diagram showing the counts of overlapping genes associated with
those GO terms.
}
\usage{
overlap_GO(go_ids, result, filename, mar = rep(0.1, 4), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{go_ids}{
A character vector of GO term identifiers to compare. For instance,
\code{head(result$scores$go_id, n=5)}.
}
  \item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from
\code{subset_scores()}.
}
  \item{filename}{
The output filename where the Venn diagram will be saved.
}
  \item{mar}{
The margins around the Venn diagram.
}
  \item{\dots}{
Further parameters forwarded to \code{the venn.diagram()} function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the output of the \code{venn.diagram()} function.
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

\section{Warning }{
  An error is returned if the list of GO term identifiers contains less than
  2 elements or more than 5, as the underlying \code{venn.diagram()} method
  does not support values outside that range.
}


\seealso{
  Method \code{\link[VennDiagram:venn.diagram]{VennDiagram::venn.diagram}}.
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
