\name{titlePage}
\alias{titlePage}
\alias{QCReport-titlePage}
\title{titlePage}
\description{
Creates the titlepage of the  QC report created by \code{\link[affyQCReport:QCReport]{QCReport}} for an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}.
}
\usage{
titlePage(object)
}

\arguments{
  \item{object}{An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  
}
\details{
  This creates the title page of a \code{\link[affyQCReport:QCReport]{QCReport}} report. 
    
  More details are available in the vignette.
  
}  
\value{
  TRUE
}

\examples{
     
    library(affydata)
    data(Dilution)
    
   titlePage(Dilution)

}
\keyword{hplot}








