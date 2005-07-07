\name{signalDist}
\alias{signalDist}
\alias{QCReport-signalDist}
\title{signalDist}
\description{
Creates the secondpage of the  QC report created by \code{\link[affyQCReport:QCReport]{QCReport}} for an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}   The page contains a boxplot and intensity distribution plots.
}
\usage{
signalDist(object)

}

\arguments{
  \item{object} {An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  
}
\details{
  This creates the second page of a \code{\link[affyQCReport:QCReport]{QCReport}} report. 
  The first graph is created with the \code{\link[affy:AffyBatch-class]{AffyBatch}}
   method \code{boxplot}. The second graph in created with the 
   \code{\link[affy:AffyBatch-class]{AffyBatch}} method \code{hist}.
 
  More details are available in the vignette.

}  
\value{
  TRUE
}

\examples{
    
    library(affydata)
    data(Dilution)

   signalDist(Dilution)

}
\keyword{hplot}








