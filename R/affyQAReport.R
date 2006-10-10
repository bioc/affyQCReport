##in some ways it would be better to have a .Rnw template, so that one
## could simply process that - but I think this will give us better
## ways to handle the errors

 affyQAReport <- function(affyB, output = "pdf",
     outdir=file.path(tempdir(), "affyQA"), 
     repName=deparse(substitute(affyB)) ) {

   if( !inherits(affyB, "AffyBatch") )
     stop("QA reports can only be generated for AffyBatch instances")

   if( !dir.create(outdir)) 
      stop("could not create output dir")

   cd = setwd(outdir)
   on.exit(setwd(cd))

   outfiles = list(sA = "simpleaffy", hist="affyhist", bxp = "affybxp")

   outf = lapply(outfiles, function(x) paste(x, output, sep=".") )

   numArrays = nrow(phenoData(affyB))

   qcStats = qc(affyB)
   pdf(file=outf$sA)
   plot(qcStats)
   dev.off()

   pdf(file=outf$hist)
   hist(affyB, lty=1)
   dev.off()

   pdf(file=outf$bxp)
   boxplot(affyB)
   dev.off()

   texTemp = system.file("Templates/affyQAtemplate.tex",
		   package="affyQCReport")

   symVals = c(repName=repName, outfiles)

   outFile = file.path(outdir, paste(repName, ".tex", sep=""))

   copySubstitute(texTemp, outFile, symbolValues = symVals)

   syscall = paste("pdflatex", outFile)
   system(syscall)
   return(list(qcStats=qcStats, loc=outdir))   
 }

