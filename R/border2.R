
      
  borderQC2 <-  function(object)
  
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
                  
                 
                 rmon<-colMeans(data.frame(intensity(object)[righton,]))
                 lmon<-colMeans(data.frame(intensity(object)[lefton,]))
                 xcmon<- (rmon - lmon)/(rmon+lmon)
                    
                 tmon<-colMeans(data.frame(intensity(object)[topon,]))
                 bmon<-colMeans(data.frame(intensity(object)[bottomon,]))
                 ycmon<- (tmon-bmon)/(tmon+bmon)
 


                 rmoff<-colMeans(data.frame(intensity(object)[rightoff,]))
                 lmoff<-colMeans(data.frame(intensity(object)[leftoff,]))
                 xcmoff<- (rmoff - lmoff)/(rmoff+lmoff)
                    
                 tmoff<-colMeans(data.frame(intensity(object)[topoff,]))
                 bmoff<-colMeans(data.frame(intensity(object)[bottomoff,]))
                 ycmoff<- (tmoff-bmoff)/(tmoff+bmoff)


       #check for out of bounds xcm or ycm

                 ArrayIndex = as.character(1:length(sampleNames(object)))
                 flagcmoff <- ArrayIndex[xcmoff < -.5 | xcmoff > .5 | ycmoff < -.5 | ycmoff > .5 ] 
                 flagcmon <- ArrayIndex[xcmon < -.5 | xcmon > .5 | ycmon < -.5 | ycmon > .5 ] 


                 

                  #win.graph()
     
                 
                 
                 
                 
                 
                 layout(matrix(c(1,1,2,2,1,1,2,2,3,3,3,3), 3, 4, byrow = TRUE))
                 
                
                  plot(xcmon,ycmon,xlim=c(-1,1),ylim=c(-1,1),
                       xlab="X Center of Intensity position",ylab="Y Center of Intensity position",
                       main="Positive Elements"  
                       )
                     
                   if(length(flagcmon) > 0) {

                                        text(xcmon[as.numeric(flagcmon)],ycmon[as.numeric(flagcmon)],ArrayIndex[as.numeric(flagcmon)],pos=2)
                   
                                          }
                plot(xcmoff,ycmoff,xlim=c(-1,1),ylim=c(-1,1),
                     xlab="X Center of Intensity position",ylab="Y Center of Intensity position",
                       main="Negative Elements" 
                      )
                 

                if(length(flagcmoff) > 0) {
                                  text(xcmoff[as.numeric(flagcmoff)],ycmoff[as.numeric(flagcmoff)],ArrayIndex[as.numeric(flagcmoff)],pos=2)
                                     }


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
         
        





#borderQC2(pdata)









