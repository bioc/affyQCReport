\name{correlationPlot}
\alias{correlationPlot}
\alias{QCReport-correlationPlot}
\title{correlationPlot}
\description{
Creates the sixth page of the QC report created by   \code{\link[affyQCReport:QCReport]{QCReport}} 
             for an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}   The page displays array-array intensity correlation coefficients.
}

\usage{
correlationPlot(object)
}

\arguments{
 \item{object} {An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  }
\details{
  This creates the sixth page of a \code{\link[affyQCReport:QCReport]{QCReport}} report. 
  The page displays array-array intensity correlation coefficients for all pairs of arrays in an \code{\link[affy:AffyBatch-class]{AffyBatch}} object. Arrays are ordered according to phenotypic data if available. This plot is useful for detecting outlier arrays.
  
More details are available in the vignette.

  
}  
\value{
  TRUE
}

\examples{

   library(affydata)
    
   data(Dilution)

  correlationPlot(Dilution)

}
\keyword{hplot}








