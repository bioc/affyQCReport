
          correlationPlot <-function(object)
                  {

                   ArrayIndex = as.character(1:length(sampleNames(object)))
			 pmat<-as.matrix(pData(phenoData(object)))
                   phenodepth<-min(ncol(pmat),3)  #allow up to three levels for sorting
                   
                   order<-switch(phenodepth+1, 
                                ArrayIndex,
                                order(pmat[,1]),
                                order(pmat[,1],pmat[,2]),
                                order(pmat[,1],pmat[,2],pmat[,3])
                                )

                   
                   arraypos <- (1:length(ArrayIndex))     *  ( 1/(length(ArrayIndex)-1)) - (1/(length(ArrayIndex)-1))
                  arraypos2=seq(1:length(ArrayIndex)-1)
			for(i in 2:length(ArrayIndex)) 
                  {arraypos2[i-1] <- (arraypos[i]+arraypos[i-1])/2
                  }

                   #win.graph()
                  layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2,3,3,3,4), 4, 4, byrow = TRUE))
                     
                  c<-cor(exprs(object)[ ,order],method = "spearman")
                  
                   image(c ,xaxt="n",yaxt="n",xlab="Array Index"
                        ,ylab="Array Index" , main="Array-Array Intensity Correlation")

                    abline(h = arraypos2, v = arraypos2)
                    


                 axis(1, labels=as.character(order)   , at=arraypos)
                   axis(2, labels=as.character(order)   , at=arraypos)
                    
                     
                    m=matrix(pretty(c,10),nrow=1,ncol=length(pretty(c,10)))

                    image(m,xaxt="n",yaxt="n",ylab="Correlation Coefficient")

                 axis(2, label= as.list(pretty(c,10)),at=seq(0,1,by=  (1/  ( length (pretty(c,10)) -1 )   )) )
                       
                     abline( h= seq( (1/  ( length (pretty(c,10)) -1 )   )/2, 1-(1/  ( length (pretty(c,10)) -1 )   ),by=(1/  ( length (pretty(c,10)) -1 )   )))
 




                   plot(1,1, type="n",xlim=c(0,length(ArrayIndex)),ylim=c(0,phenodepth),
                       xlab="Array Index", ylab="",yaxt="n",xaxt="n",xaxs="i" )
         
              axis(1, labels=as.character(order),at=seq(1-.5,length(ArrayIndex)-.5))
               axis(4, label= as.list(attr(object@phenoData,which="varLabels")[1:phenodepth]),at=seq(1-.5,phenodepth-.5),las=1)
             #  axis(3,label= as.list(substr(sampleNames(object)[order],1,8)), at=seq(1-.5,length(ArrayIndex)-.5),las=3,cex=.2)
  


              abline(h = seq(1,phenodepth), v = seq(1,length(ArrayIndex)), col = "lightgray")
               for(i in 1:phenodepth)
			{

                   text(seq(.5,length(ArrayIndex)-.5,by=1),rep(i-1+.5,length(ArrayIndex)),pmat[order,i],cex=0.7)



                  }

                
                    return(TRUE)
    }
         
        



#correlationPlot(pdata)


