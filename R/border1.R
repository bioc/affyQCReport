        
  borderQC1 <-  function(object)
     {
           
            #get array dimensions 
                 
                 n <- object@nrow
                
            #create indices for sides
                  
                 left  <-  seq(0,n-1)*n+1
                 right <-  seq(1,n)*n

                 top   <-  seq(1,n)
                 bottom<-  seq(1,n)+(n*(n-1)) 
                 
            #the sides will be split into on and off based on intensity
		  #the first array will be used

                 leftmean   <-mean(intensity(object)[left,1])
                 rightmean  <-mean(intensity(object)[right,1])
                 topmean    <-mean(intensity(object)[top,1]) 
                 bottommean <-mean(intensity(object)[bottom,1])

                 lefton     <- left[intensity(object)[left,1]  >1.2*leftmean]
                 righton    <- right[intensity(object)[right,1] >1.2*rightmean]
                 topon      <- top[intensity(object)[top,1]   >1.2*topmean]
                 bottomon   <- bottom[intensity(object)[bottom,1]>1.2*bottommean]
 
		     leftoff     <- left[intensity(object)[left,1] < .8*leftmean]
                 rightoff    <- right[intensity(object)[right,1] <.8*rightmean]
                 topoff      <- top[intensity(object)[top,1]   <.8*topmean]
                 bottomoff   <- bottom[intensity(object)[bottom,1]<.8*bottommean]

                 on<-c(lefton,righton,topon,bottomon)
                 off<-c(leftoff,rightoff,topoff,bottomoff)

             #calculate center of intensity
                  
                 
                 rmon<-mean(data.frame(intensity(object)[righton,]))
                 lmon<-mean(data.frame(intensity(object)[lefton,]))
                 xcmon<- (rmon - lmon)/(rmon+lmon)
                    
                 tmon<-mean(data.frame(intensity(object)[topon,]))
                 bmon<-mean(data.frame(intensity(object)[bottomon,]))
                 ycmon<- (tmon-bmon)/(tmon+bmon)
 


                 rmoff<-mean(data.frame(intensity(object)[rightoff,]))
                 lmoff<-mean(data.frame(intensity(object)[leftoff,]))
                 xcmoff<- (rmoff - lmoff)/(rmoff+lmoff)
                    
                 tmoff<-mean(data.frame(intensity(object)[topoff,]))
                 bmoff<-mean(data.frame(intensity(object)[bottomoff,]))
                 ycmoff<- (tmoff-bmoff)/(tmoff+bmoff)


       #check for out of bounds xcm or ycm

                 ArrayIndex = as.character(1:length(sampleNames(object)))
                 flagcmoff <- ArrayIndex[xcmoff < -.5 | xcmoff > .5 | ycmoff < -.5 | ycmoff > .5 ] 
                 flagcmon <- ArrayIndex[xcmon < -.5 | xcmon > .5 | ycmon < -.5 | ycmon > .5 ] 


                 

                  #win.graph()
                  layout(matrix(c(1,1,2,2,1,1,2,2,3,3,3,3), 3, 4, byrow = TRUE))
                  
                    
                  boxplot(data.frame(intensity(object)[on,]),
                          xlab="Array Index",ylab="Intensity",main="Positive Border Elements",
                          names=ArrayIndex
                          )
                 

                  boxplot(data.frame(intensity(object)[off,]),
                          xlab="Array Index",ylab="Intensity",main="Negative Border Elements",
                           names=ArrayIndex
                          )
                 
                 
                 
                 
                  plot(ArrayIndex,ArrayIndex, type="n",xlim=c(0,length(ArrayIndex)),ylim=c(0,1),
                       xlab="Array Index",ylab="Sample Name",yaxt="n", xaxt="n",xaxs="i" )
              midpoint<-seq(1,length(ArrayIndex)) -0.5
            
               axis(1,label= ArrayIndex, at=seq(1-.5,length(ArrayIndex)-.5),cex=.2)

              abline(v = ArrayIndex)
              for(i in 1:length(ArrayIndex))
			{
                text(midpoint[i],.5,substr(sampleNames(object)[i],1,10),cex=.7)

                }
                 
                 
                 
                 
                 
                 
                 
                 
                    return(TRUE)
                 
         
        }





#borderQC1(pdata)
