\name{Dbp.merge}
\alias{Dbp.merge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dbp-merge
}
\description{
Merge all sites on an interval, if the sites at the ends of the interval are ``close" to each other 
}
\usage{
Dbp.merge(ordr.vec, thresh.dist, bp.thresh.dist, location.vec, dist.type = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ordr.vec}{
An (n by m) matrix of n measurements of m random variables (methylations). The columns are ordered.
}
  \item{thresh.dist}{
A similarity distance threshold. Sites in the interval are merged to a single cluster if the similarity distance between them is above thresh.dist.}
  \item{bp.thresh.dist}{
A threshold distance in (e.g.) chromosomal location. Any interval to be potentially merged is smaller or equal to bp.thresh.dist.
}
  \item{location.vec}{
An m vector specifying the locations (e.g. chromosomal locations) of the variables measured in the matrix ordr.vec. 
}
  \item{dist.type}{
Type of similarity distance function "spearman", "person" or "euclid". 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An m vector of cluster assignments of the m ordered variables. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
