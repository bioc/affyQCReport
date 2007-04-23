\name{affyQAReport}
\alias{affyQAReport}
\alias{openQAReport}
\alias{rmQAReport}
\title{ Generate a QA Report for Affymetrix arrays}
\description{
These functions provide a method for creating, viewing and removing QA reports
on an experiment run on Affymetrix CEL files.
}
\usage{
affyQAReport(affyB, output = "pdf", outdir = file.path(getwd(), "affyQA"), 
       overwrite = FALSE, repName)
rmQAReport(x)
openQAReport(x)
}
\arguments{
  \item{affyB}{ An instance of the \code{AffyBatch} class. }
  \item{output}{Type of output wanted, currently only pdf is supported. }
  \item{outdir}{ The directory (or folder) where output will be directed.}
  \item{overwrite}{Whether an existing report, with the same name, should
      be overwritten.}
  \item{repName}{ A character string indicating the name of the report. }
}
\details{
 A wide variety of tools from different packages are used, in conjunction 
with Sweave tools to produce an integrated document with various statistics
that should help determine whether there are data quality problems with the
data.

 This package only works for Affymetrix data. It requires a working LaTeX
implementation and so is not likely to work on Windows platforms, and
in some cases on OS X unless the user has installed these tools.
}
\value{
 A list with the following components is returned.
  \item{qcStats}{The QC statistics computed by the \pkg{simpleaffy} package.}
  \item{affyPLM}{The QC statistics computed by the \pkg{affyPLM} package.}
  \item{MADS}{The computed MAD statistic on the \sQuote{M} values from
 an M-A decomposition.}
  \item{loc}{The location (directory/folder) containing the report.}
  \item{name}{The name of the report.}
}
\author{R. Gentleman }

\examples{
\dontrun{
   library("affydata")
   data("Dilution")
   affyQAReport(Dilution)
}
}
\keyword{ manip }
