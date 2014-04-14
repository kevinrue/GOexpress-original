\name{quantiles_scores}
\alias{quantiles_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the quantiles of scores in a GO_anova analysis.
}
\description{
Percentiles indicating the scores reached by a proportion of the GO terms tested.
}
\usage{
quantiles_scores(result, probs=c(0.9, 0.95, 0.99, 0.999, 0.9999), quartiles=c(FALSE, TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
The output of \code{GO_anova()} or a subset of it obtained from \code{subset_scores()}.
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
A named vector showing the percentiles and corresponding scores.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
