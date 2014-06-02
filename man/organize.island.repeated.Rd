\name{organize.island.repeated}
\alias{organize.island.repeated}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
INTERNAL	
}
\description{
Prepared a data set with multiple outcomes (methylation in multiple CpG sites in multiple rows per subject) to be such that it can be used by the GEE function
}
\usage{
organize.island.repeated(X, covar.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
an (m by n) matrix of methylation values in $m$ sites for $n$ participants
}
  \item{covar.mat}{
an (n by p) matrix of covariates. The rows ordering is the same as the column ordering of the matrix X. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
