\name{quantiles_scores}
\alias{quantiles_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the quantiles of scores following an GOexpress analysis.
}
\description{
Returns a set of quantiles in indicating the scores reached by given proportions
of the GO terms.

This function can also be used on the output of \code{subset_scores()}
function as it returns a value formatted identically to the output of the
\code{GO_analyse()} function.
}
\usage{
quantiles_scores(result, probs=c(0.9, 0.95, 0.99, 0.999, 0.9999), quartiles=c(FALSE, TRUE))
}
\arguments{
  \item{result}{
The output of \code{GO_analyse()} or a subset of it obtained from \code{subset_scores()}.
}
  \item{probs}{
numeric vector of probabilities with values in [0,1]. (Values up
to 2e-14 outside that range are accepted and moved to the nearby endpoint.)
See \code{quantile}
}
  \item{quartiles}{
A numeric vector of the percentiles for which the scores are desired. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A named vector of percentiles and corresponding scores.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Method \code{\link[stats:quantile]{stats::quantile}}.
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
