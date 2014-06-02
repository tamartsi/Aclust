\name{load.annotation}
\alias{load.annotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load annotation table 
}
\description{
The package is adjusted to Illumina annotation .csv file. 
}
\usage{
load.annotation(annotation.file.name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{annotation.file.name}{
The name of the .csv Illumina annotation file. 
}
}

\value{
A data.table object with the annotation.
}
\examples{

### not run! need to put the appropriate path to the annotation file
# annotation.file.name <- "illumina_450_menifest_v.1.2.csv"
# annot <- load.annotation(annotation.file.name)
# annot[1:5]
}
