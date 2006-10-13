##in some ways it would be better to have a .Rnw template, so that one
## could simply process that - but I think this will give us better
## ways to handle the errors

 rmQAReport = function(outdir=file.path(tempdir(), "affyQA"))
    unlink(outdir, recursive=TRUE)

 affyQAReport <- function(affyB, output = "pdf",
     outdir=file.path(tempdir(), "affyQA"), 
     repName=deparse(substitute(affyB)) ) {

   if( !inherits(affyB, "AffyBatch") )
     stop("QA reports can only be generated for AffyBatch instances")

   if( !dir.create(outdir)) 
      stop("could not create output dir")

   cd = setwd(outdir)
   on.exit(setwd(cd))

   outfiles = list(sA = "simpleaffy", RNAdeg = "rnadeg", 
       hist="affyhist", bxp = "affybxp", NUSE="NUSE", RLE="RLE")

   outf = lapply(outfiles, function(x) paste(x, output, sep=".") )

   numArrays = nrow(phenoData(affyB))

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
   bbO = cbind(BioBCall = bbc, biobOut)

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
   plotAffyRNAdeg(rnaDeg, cols=1)
   dev.off()

   ##MA plots - 8 of these per page seems like the right number
   ##normalize and bg correct
   pp1 = preprocess(affyB)
   epp1 = exprs(pp1)
   medArray = rowMedians(epp1)
   M <- sweep(epp1,1,medianchip,FUN='-')
   A <- 1/2*sweep(epp1,1,medianchip,FUN='+')

   nfig = ceiling(numArrays/8)
   plotNames = paste("MA", 1:nfig, sep="")
   fNames = paste(plotNames, "pdf", sep=".")
   op = par(mfrow=c(4,2)
   nprint = 1
   for( i in 1:nfig) {
       pdf(file=fNames[i])
       for(j in 1:8) {
           if(nprint <= numArrays) {
               title <- paste(sampleNames(affyB)[nprint],
                            "vs pseudo-median reference chip")
               ma.plot(A[,nprint],M[,nprint], main=title, xlab="A",
                        ylab="M", pch=".")
               nprint <- nprint + 1
           }
      }
   } 

   ##FIXME: add some latex to another variable so we can put
   ## it into the document


   ##affyPLM stuff
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
         TABLE3=TAB3,
         affyQCVersNO= pkVers, sessionInfo=sessInfo )

   outFile = file.path(outdir, paste(repName, ".tex", sep=""))

   copySubstitute(texTemp, outFile, symbolValues = symVals)

   ##set up call to pdflatex and run it twice for x-refs
   syscall = paste("pdflatex", outFile)
   system(syscall)
   system(syscall)

   return(list(qcStats=qcStats, loc=outdir))   
 }

