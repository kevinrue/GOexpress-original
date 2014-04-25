\name{hist_scores}
\alias{hist_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the distribution of scores following an GOexpress analysis.
}
\description{
Plots the an histogram representing the frequencies of scores in the output
variable of the {GO_analyse()} function.

This function can also be used on the output of \code{subset_scores()}
function as it returns a value formatted identically to the output of the
\code{GO_analyse()} function.
}
\usage{
hist_scores(result,
            main=paste("Distribution of average F scores in",
                       deparse(substitute(result))),
            xlab="Average F score", ...)
}
\arguments{
  \item{result}{
The output of the \code{GO_analyse()} function.
}
  \item{main, xlab}{
These arguments to title have useful defaults here.
}
  \item{\dots}{
Additional arguments passed on to \code{hist()}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the output of the \code{graphics::hist()} function.
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
  Method \code{\link[graphics:hist]{graphics::hist}}.
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--  or do  help(data=index)  for the standard data sets.
}

\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
\keyword{expression}
