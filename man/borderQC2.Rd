\name{borderQC2}
\alias{borderQC2}

\title{Center of intensity QC plots}
\description{
Creates the fifth page of the  QC report created by \code{\link[affyQCReport:QCReport]{QCReport}} for an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}   The page contains plots showing the center of intensity of the border element of arrays in an \code{\link[affy:AffyBatch-class]{AffyBatch}}  object.
}
\usage{
borderQC2(object)

}

\arguments{
  \item{object} {An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  
}
\details{
  This creates the fifth page of a \code{\link[affyQCReport:QCReport]{QCReport}} report. 
  Plots of the center of intensity of the positive and negative border elements of arrays in a \code{\link[affy:AffyBatch-class]{AffyBatch}} object are shown.  This plot is useful for detecting spatial biases in intensities on an array.

More details are available in the vignette.
 
  
}  
\value{
  TRUE
}

\examples{

   library(affydata)
    data(Dilution)

   borderQC2(Dilution)

}
\keyword{hplot}








