\name{calc.dist.d.neighbor}
\alias{calc.dist.d.neighbor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL
}
\description{
Calculates distance between variables d neighbors apart
}
\usage{
calc.dist.d.neigbor(ordr.vec, d, dist.type = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ordr.vec}{
An (n by m) ordered matrix of vectors representing realization of random variables. Each column is a variable, and the columns
are ordered, here according to chromosomal location. 
}
  \item{d}{
Calculate distances between vectors  d neighbors apart. 
}
  \item{dist.type}{
The type of distance (``spearman", ``pearson", ``euclid")
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of distances between each vector and the vector d neighbors apart from it. The last d-1 vectors
are given infinity distance, so that returned vector is of length m.
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
