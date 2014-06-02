\name{calc.dist.clusters}
\alias{calc.dist.clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL	
}
\description{
Calculates distance between two clusters
}
\usage{
calc.dist.clusters(clust.1, clust.2, type = "single", dist.type = "spearman")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust.1}{
An (n by $m_1$) matrix of $m_1$ vectors composing a cluster
}
  \item{clust.2}{
An (n by  $m_2$) matrix of $m_2$ vectors composing a cluster
}
  \item{type}{
Clustering by \cor{type}: ``single", ``complete" or ``averge"
}
  \item{dist.type}{
%%     ~~Describe \code{dist.type} here~~
}
}
\details{
The distance between the clusters is
if type = ``single": the minimum of the distances between the vectors in clust.1 and the vectors in clust.2
if type = ``complete": the maximum of the distances between the vectors in clust.1 and the vectors in clust.2 
if type = ``average": the mean of the distances between the vectors in clust.1 and the vectors in clust.2 
}
\value{
A numeric value, the distance between clust.1 and clust.2
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
