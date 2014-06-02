\name{calc.dist.clust.point}
\alias{calc.dist.clust.point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL
}
\description{
Calculates the distance between a cluster of vectors and a single vector
}
\usage{
calc.dist.clust.point(clust, pnt, type = "single", dist.type = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust}{
An (n by m) matrix representing $m$ vectors ($n$ realizations for $m$ random variables) in a cluster.
}
  \item{pnt}{
An $n$ vector representing a single random variable, or a single point
}
  \item{type}{
The type composite distance, ``single", ``complete" or ``average"
}
  \item{dist.type}{
distance type ``spearman", "pearson" or ``euclid"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric value, the distance between pnt to clust. 
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
