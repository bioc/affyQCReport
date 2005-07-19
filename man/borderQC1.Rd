\name{borderQC1}
\alias{borderQC1}
\title{Distribution of intensities of the border elements}
\description{
Creates the forth page of the QC report created by \code{\link[affyQCReport:QCReport]{QCReport}} for an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}   The page contains a boxplots showing the intensities from the border element of arrays in an \code{\link[affy:AffyBatch-class]{AffyBatch}}  object.
}
\usage{
borderQC1(object)
}

\arguments{
  \item{object} {An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}} .}
  
}
\details{
  This creates the forth page of a \code{\link[affyQCReport:QCReport]{QCReport}} report. 
  Boxplots of the positive and negative border elements of arrays in a \code{\link[affy:AffyBatch-class]{AffyBatch}}  object are shown.
  
More details are available in the vignette.

}  
\value{
  TRUE
}

\examples{
    
    library(affydata)
    data(Dilution)

   borderQC1(Dilution)

}
\keyword{hplot}








