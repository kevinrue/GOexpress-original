\name{overlap_GO}
\alias{overlap_GO}
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
\arguments{
  \item{go_ids}{
A character vector of GO term identifiers to compare. For instance,
\code{head(result$scores$go_id, n=5)} to compare the first 5 top-
ranked GO terms in the\code{result} variable.
}
  \item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from
\code{subset_scores()}.
}
  \item{filename}{
The output filename where the Venn diagram will be saved.
}
  \item{mar}{
The margins around the Venn diagram. Some Venn diagrams place the GO term
labels outside the visible frame.
}
  \item{\dots}{
Further parameters forwarded to \code{the venn.diagram()} function.
}
}
\value{
  Returns the output of the \code{venn.diagram()} function.
}
\references{
\itemize{
  \item \href{http://cran.r-project.org/web/packages/VennDiagram/index.html}{VennDiagram} package.
  }
}
\author{
Kevin Rue-Albrecht
}
\section{Warning }{
  An error is returned if the list of GO term identifiers contains less than
  2 elements or more than 5, as the underlying \code{venn.diagram()} method
  does not support values outside that range.
}
\seealso{
  Method \code{\link[VennDiagram:venn.diagram]{VennDiagram::venn.diagram}} and
  \href{https://github.com/kevinrue/GOexpress/wiki/Typical-analysis-pipeline}{tutorial} on GitHub.
}
\examples{
# load the sample output data
data(raw_result)

# filter for Biological Processes associated with 10+ genes
filtered_result = subset_scores(result=raw_result, total_count=5, namespace="BP")

# Venn diagram of overlapping genes between top 5 GO terms
overlap_GO(go_ids=head(filtered_result$GO$go_id, n=5), result=filtered_result, filename="VennDiagram.tiff")
}

\keyword{ GOexpress }
\keyword{ anova }
\keyword{ gene }
\keyword{ expression }
\keyword{ Venn }
\keyword{ ontology }
