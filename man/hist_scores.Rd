\name{hist_scores}
\alias{hist_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the distribution of scores in a GO_anova analysis.
}
\description{
An histogram of score frequencies following a GO_anova analysis.
}
\usage{
hist_scores(result,
            main=paste("Distribution of average F scores in",
                       deparse(substitute(result))),
            xlab="Average F score", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
The output of the \code{GO_anova()} function.
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
Same as a call to hist(...)
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
See also \code{\link{hist}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GO }
\keyword{ gene }
\keyword{ ontology }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
