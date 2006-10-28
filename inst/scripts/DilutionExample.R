library("affyQCReport")
library("affydata")
data("Dilution")

xxx = affyQAReport(Dilution)
openQAReport(xxx)

## to delete it
## rmQAReport(xxx)

