\name{update.clust.indicator}
\alias{update.clust.indicator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL	
}
\description{
The function updates the cluster indicator of a set of variables to match to that of another variable (i.e. they are merged to be in the same cluster) 
}
\usage{
update.clust.indicator(which.clust, ind.1, ind.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{which.clust}{
A vector of cluster indicators corresponding to a set of ordered variables 
}
  \item{ind.1}{
A position in the set of ordered variables of one variable
}
  \item{ind.2}{
A position in the set of ordered variables of another variable
}
}
\details{
All variables between ind.2 and ind.1 (including) will be assigned to the cluster of the variable corresponding to ind.2. 
}
\value{
The updated vector of cluster assignments. 
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
