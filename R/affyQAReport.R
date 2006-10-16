##in some ways it would be better to have a .Rnw template, so that one
## could simply process that - but I think this will give us better
## ways to handle the errors

## I started thinking that the user might specify the output type -
##  but then hardwired pdf in places - this might be revisited

 rmQAReport = function(x)
    unlink(x$loc, recursive=TRUE)

 openQAReport = function(x)
     openPDF(file.path(x$loc, paste(x$name, ".pdf", sep="")))

 affyQAReport <- function(affyB, output = "pdf",
     outdir=file.path(tempdir(), "affyQA"), 
     repName=deparse(substitute(affyB)) ) {

   if( !inherits(affyB, "AffyBatch") )
     stop("QA reports can only be generated for AffyBatch instances")

   if( !file.exists(outdir) )
     if( !dir.create(outdir)) 
        stop("could not create output directory")
   
   outdir = file.path(outdir, repName)
   if( file.exists(outdir) )
       stop("report already exists")
   if( !dir.create(outdir))
        stop("could not create report directory")

   cd = setwd(outdir)
   on.exit(setwd(cd))

   outfiles = list(sA = "simpleaffy", RNAdeg = "rnadeg", 
       hist="affyhist", bxp = "affybxp", NUSE="NUSE", RLE="RLE")

   outf = lapply(outfiles, function(x) paste(x, output, sep=".") )

   numArrays = nrow(pData(affyB))

   qcStats = qc(affyB)
  
   dfout = data.frame(AvBg = avbg(qcStats), ScaleF=sfs(qcStats),
      PerCPres=percent.present(qcStats))

   tab1 = xtable(dfout, label="table1")
   tcon = textConnection("TAB1", "w", local=TRUE)
   print(tab1, file=tcon)
   close(tcon)
   TAB1 = paste(TAB1, collapse="\n")

   tab2 = xtable(t(ratios(qcStats)), label="table2")
   tcon = textConnection("TAB2", "w", local=TRUE)
   print(tab2, file=tcon)
   close(tcon)
   TAB2 = paste(TAB2, collapse="\n")

   biobOut = qcStats@spikes
   cn = colnames(biobOut)
   cn = gsub("AFFX-", "", cn, fixed=TRUE)
   cn = gsub("-3_at", "", cn, fixed=TRUE)
   colnames(biobOut) = cn
   bbc = qcStats@bioBCalls
   names(bbc) = gsub(".present", "", names(bbc))
   bbO = cbind(BioBCall = bbc, round(biobOut, digits=3))

   tab3 = xtable(bbO, label="table3", caption="BioB and friends")
   tcon = textConnection("TAB3", "w", local=TRUE)
   print(tab3, file=tcon)
   close(tcon)
   TAB3 = paste(TAB3, collapse="\n")


   pdf(file=outf$sA)
   plot(qcStats)
   dev.off()

   pdf(file=outf$hist)
   hist(affyB, lty=1)
   dev.off()

   pdf(file=outf$bxp)
   boxplot(affyB)
   dev.off()

   ##RNA degredation plot
   rnaDeg = AffyRNAdeg(affyB)
   pdf(file=outf$RNAdeg)
   plotAffyRNAdeg(rnaDeg, cols=rep(1, numArrays))
   dev.off()

   ##MA plots - 8 of these per page seems like the right number
   ##normalize and bg correct
   pp1 = preprocess(affyB)
   epp1 = exprs(pp1)
   medArray = rowMedians(epp1)
   M <- sweep(epp1,1,medArray,FUN='-')
   A <- 1/2*sweep(epp1,1,medArray,FUN='+')

   ##FIXME: need a table of the statistics that are often shown in
   ## an ma.plot

   ##if we have 4 or 6 arrays try to make the plots larger
   ##if lots, then 8 arrays per plot
   if( numArrays <= 4 ) 
      app = 4
   else if( numArrays <= 6 )
      app = 6
   else
      app = 8

   nfig = ceiling(numArrays/8)

   ##FIXME: need a more efficient version of ma.plot - 
   ## one that has a smaller graphics file - use hexbin?
   plotNames = paste("MA", 1:nfig, sep="")
   fNames = paste(plotNames, "pdf", sep=".")
   nprint = 1
   for( i in 1:nfig) {
       pdf(file=fNames[i])
       par(mfrow=c(app/2, 2))
       for(j in 1:app ) {
           if(nprint <= numArrays) {
               title <- paste(sampleNames(affyB)[nprint],
                            "vs pseudo-median reference")
               ma.plot(A[,nprint],M[,nprint], main=title, xlab="A",
                        ylab="M", pch=".", show.statistics=FALSE)
               nprint <- nprint + 1
           }
      }
      dev.off()
   } 

  plotnames = paste("\\includegraphics{", plotNames, "}", sep="")

  MALatex = paste("\\begin{figure}[tp]",
                  "\\centering",
		  plotnames,
                  "\\caption{MA plots.}",
                  "\\end{figure}", sep=" \n", collapse="\n\n")

  ##WH's distance plots here

  outM = matrix(0, nrow=numArrays, ncol=numArrays)
  for(i in 1:(numArrays-1)) 
    for(j in (i+1):numArrays ) 
      outM[i,j] = outM[j,i] = mad(epp1[,i] - epp1[,j])

  pdf(file="MADimage.pdf")
  image(1:numArrays, 1:numArrays, outM, xlab="", ylab="",
       main="Distance between arrays\nMeasured by MAD on M-values")
  dev.off()

  ##affyPLM stuff
  ##FIXME: should be able to use pp1, from above here so we
  ## only need to do summarization, and that will make this run
  ## quite a bit faster

  dataPLM = fitPLM(affyB)
  
  #Normalized Unscaled Standard Error (NUSE)

  pdf(file=outf$RLE)
  Mbox(dataPLM, ylim = c(-1, 1), names = NULL, col="lightblue",
   whisklty=0, staplelty=0, main="RLE")
  dev.off()

   pdf(file=outf$NUSE)
   boxplot(dataPLM, ylim = c(0.95, 1.5), names = NULL,
        outline = FALSE, col="lightblue", main="NUSE")
   dev.off()

   ##write the LaTeX
   texTemp = system.file("Templates/affyQAtemplate.tex",
		   package="affyQCReport")

   #get version numbers and sessionInfo
   pkVers =  packageDescription("affyQCReport")$Version
   sessInfo = paste(toLatex(sessionInfo()), collapse="\n")

   symVals = c(repName=repName, outfiles, TABLE1=TAB1, TABLE2=TAB2,
         TABLE3=TAB3, MAPLOTS = MALatex, MADimage="MADimage",
         affyQCVersNO= pkVers, sessionInfo=sessInfo,
         numArrays=as.character(numArrays), chipName = affyB@cdfName )

   outFile = file.path(outdir, paste(repName, ".tex", sep=""))

   copySubstitute(texTemp, outFile, symbolValues = symVals)

   ##set up call to pdflatex and run it twice for x-refs
   syscall = paste("pdflatex", outFile)
   system(syscall)
   system(syscall)

   return(list(qcStats=qcStats, affyPLM=dataPLM, MADS=outM, loc=outdir,
       name=repName))   
 }

