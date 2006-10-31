#Copyright R. Gentleman, 2006, all rights reserved


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
     outdir=file.path(getwd(), "affyQA"), overwrite = FALSE,
     repName) {
   
   if(missing(repName) )
       repName = deparse(substitute(affyB))

   if( !inherits(affyB, "AffyBatch") )
     stop("QA reports can only be generated for AffyBatch instances")

   sN = sampleNames(affyB)
   lcS = lcsuff(sN)
   if(nchar(lcS) > 0 )
   sN = gsub(paste(lcS, "$", sep=""), "", sN)

   if( !file.exists(outdir) )
     if( !dir.create(outdir)) 
        stop("could not create output directory")
   
   outdir = file.path(outdir, repName)
   if( file.exists(outdir) )
      if( overwrite )
          unlink(outdir, recursive=TRUE)
      else
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

   tab1 = xtable(dfout, label="table1", caption="Background and scale factors")
   tcon = textConnection("TAB1", "w", local=TRUE)
   print(tab1, file=tcon)
   close(tcon)
   TAB1 = paste(TAB1, collapse="\n")

   ##output ratios

   qcratios = ratios(qcStats)
   rnrat = colnames(qcratios)
   rn1 = gsub("^AFFX-", "", rnrat)
   rn2 = strsplit(rn1, "\\.")
   if( any(sapply (rn2, length) != 2) ) 
      colnames(qcratios) = rn1
   else
      colnames(qcratios) = sapply(rn2, paste, collapse="\n")

   tab2 = xtable(qcratios, label="table2", caption="QC ratios")
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

   acol=sample(brewer.pal(8, "Dark2"), numArrays, replace=(8<numArrays))
   argb = sapply(acol, substring, first=c(2,4,6), last=c(3,5,7))
   argb = apply(argb, 2, function(h) as.integer(paste("0x", h, sep=""))/255)

   definecolor = paste(paste("\\definecolor{farbe", seq_along(acol),
     "}{rgb}{",argb[1,], ",",argb[2,],",",argb[3,],"}", sep=""),
     collapse="\n")
   arrayNamesInColors = paste(paste("\\textcolor{farbe",
     seq_along(acol), "}{", sN, "}", sep=""),
     collapse=", ")
   
   pdf(file=outf$sA)
   plot(qcStats)
   dev.off()

   pdf(file=outf$hist)
   hist(affyB, lty=1, col=acol)
   dev.off()

   pdf(file=outf$bxp)
   boxplot(affyB, col=acol, las = 3)
   dev.off()

   ##RNA degredation plot
   rnaDeg = AffyRNAdeg(affyB)
   pdf(file=outf$RNAdeg)
   plotAffyRNAdeg(rnaDeg, cols=acol, lwd=2)
   dev.off()

   ##MA plots - 8 of these per page seems like the right number
   ##normalize and bg correct
   pp1 = preprocess(affyB)
   epp1 = log2(exprs(pp1))
   medArray = rowMedians(epp1)
   M =  epp1-medArray
   A = (epp1+medArray)/2

   ##FIXME: need a table of the statistics that are often shown in
   ## an ma.plot

   ##if we have 4 or 6 arrays try to make the plots larger
   ##if lots, then 8 arrays per plot
   app = 4 + 2*(sum(numArrays>c(4,6)))
   nfig = ceiling(numArrays/8)

   ##FIXME: need a more efficient version of ma.plot - 
   ## one that has a smaller graphics file - use hexbin?
   plotNames = paste("MA", 1:nfig, sep="")
   fNames = paste(plotNames, "pdf", sep=".")
   ## nprint = 1
   xlim = quantile(A, probs=1e-4*c(1,-1)+c(0,1))
   ylim = quantile(M, probs=1e-4*c(1,-1)+c(0,1))


   dummy.df <-
       data.frame(sN = factor(sN, levels = sN),
                  x = seq_along(sN),
                  y = seq_along(sN))

   trobj <-
       xyplot(y ~ x | sN, dummy.df,
              xlim = xlim,
              ylim = ylim,
              xlab = "A",
              ylab = "M",
### possible substitute for xlim and ylim:
              
##               prepanel = function(x, y, ...) {
##                   x <- A[, x]
##                   y <- M[, y]
##                   list(xlim = range(x, na.rm = TRUE),
##                        ylim = range(y, na.rm = TRUE),
##                        dx = diff(x),
##                        dy = diff(y))
##               },

              panel = function(x, y, ...) {
                  x <- A[, x]
                  y <- M[, y]
                  panel.smoothScatter(x, y, ...)
              },

              layout = c(app/2, 2, 1))

   id.firstpage <- seq_len(app)
   
   for(i in seq_len(nfig))
   {
       pdf(file = fNames[i])
       id.thispage <- (i-1) * app + id.firstpage
       id.thispage <- id.thispage[id.thispage <= numArrays]

       ## print(id.thispage)
       ## print(trobj[id.thispage]) # should work, bug?
       print(update(trobj, index.cond = list(id.thispage)))

##      par(mfrow=c(app/2, 2))
##      for(j in seq_len(app)) {
##          if(nprint <= numArrays) {
##              smoothScatter(A[,nprint], M[,nprint],
##                            main=sN[nprint],
##                            xlab="A", ylab="M", xlim=xlim, ylim=ylim)
##              abline(h=0, col="#fe0020")
##              ## ma.plot(A[,nprint],M[,nprint], main=title, xlab="A",
##              ##         ylab="M", pch=".", show.statistics=FALSE)
##              nprint <- nprint + 1
##          }
##      }

       dev.off()
   }

   MALatex = paste("\\begin{figure}[tp]",
   "\\centering",
   paste("\\includegraphics{", plotNames, "}", sep=""),
   "\\caption{MA plots. A \\textit{reference array} array is calculated from",
   "the median across arrays, and for each array $M$ and $A$ values are",
   "calculated for the comparison to that reference.}",
   "\\end{figure}", sep=" \n", collapse="\n\n")

   
  ##WH's distance plots here
  outM = matrix(0, nrow=numArrays, ncol=numArrays)
  for(i in seq_len(numArrays-1))
    for(j in (i+1):numArrays) 
      outM[i,j] = outM[j,i] = mad(epp1[,i] - epp1[,j])

  row.names(outM) = colnames(outM) = sN

  pdf(file="MADimage.pdf", height=6, width=(6-0.7)*1.25+0.7)
  par(mai=c(0.7, 0.7, 0.01, 0.01))
  layout(cbind(1,2), widths=c(4, 1))
  imcol=colorRampPalette(brewer.pal(9, "RdPu"))(256)
  rg=range(outM)
  image(1:numArrays, 1:numArrays, outM, xlab="", ylab="", zlim=rg,
       main="", col=imcol, axes=FALSE)
  axis(1, at=1:numArrays, labels=sN, las=3)
  axis(2, at=1:numArrays, labels=sN, las=2)
  image(1, seq(rg[1], rg[2], length=length(imcol)), rbind(seq_along(imcol)),
        xaxt="n", ylab="", col=imcol)
  text(1, 0, "MAD", xpd=NA)
   
  dev.off()
   
  ##affyPLM stuff
  ## by using pp1 from above we 
  ## only need to do summarization, and that will make this run
  ## quite a bit faster

  dataPLM = fitPLM(pp1, background=FALSE, normalize=FALSE)
  
  #Normalized Unscaled Standard Error (NUSE)

  pdf(file=outf$RLE)
  Mbox(dataPLM, ylim = c(-1, 1), names = sN, col="lightblue",
   whisklty=0, staplelty=0, main="RLE", las=3)
  dev.off()

   pdf(file=outf$NUSE)
   boxplot(dataPLM, ylim = c(0.95, 1.5), names = sN,
        outline = FALSE, col="lightblue", main="NUSE", las=3)
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
         definecolor=definecolor, arrayNamesInColors=arrayNamesInColors,
         numArrays=as.character(numArrays), chipName = affyB@cdfName )

   outFile = file.path(outdir, paste(repName, ".tex", sep=""))

   copySubstitute(texTemp, outFile, symbolValues = symVals)

   ##fix directory structure for latex/win32

   if(.Platform$OS.type == "windows")
     outFile = shortPathName(outFile)

   ##set up call to pdflatex and run it twice for x-refs
   syscall = paste("pdflatex", outFile)
   system(syscall)
   system(syscall)

   return(list(qcStats=qcStats, affyPLM=dataPLM, MADS=outM, loc=outdir,
       name=repName))   
 }

