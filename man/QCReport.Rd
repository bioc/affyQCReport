\name{QCReport}
\alias{QCReport}
\title{QC Report for AffyBatch objects}
\description{
Creates a QC report for an object of class
\code{\link[affy:AffyBatch-class]{AffyBatch}}. 
}
\usage{
QCReport(object,file="AffyQCReport.pdf",...)
}

\arguments{
  \item{object}{An object of class \code{\link[affy:AffyBatch-class]{AffyBatch}}.}
  \item{file}{A valid filename for the output PDF file. }  
  \item{...}{further arguments for \code{\link{plot}} function.}
}
\details{
  This creates a 6 page QC report in PDF file format. 
  \itemize{
  \item The default output file is AffyQCReport.pdf created in the working directory.
  \item The \code{\link[affy:AffyBatch-class]{AffyBatch}} object is not required to contain phenotypic data but that data is used if available.
 
  \item The CDF environment must be available for the method to work fully.
  }
  More details are available in the vignette.
}  
\value{
  TRUE
}

\examples{
  \dontrun{
     
    library(affydata)
    data(Dilution)
    
  

    QCReport(Dilution,file="ExampleQC.pdf")   ##A QC report for the Example data will be generated
  }
}
\keyword{hplot}
