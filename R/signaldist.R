
  signalDist<- function(object)


{

par(mfrow=c(2, 1))
ArrayIndex = as.character(1:length(sampleNames(object)))
boxplot(object,names=ArrayIndex,ylab="Log2(Intesity)",xlab="Array Index")

hist(x=object,lt=1:length(ArrayIndex),col=1:length(ArrayIndex),which="both")
temppar<-par()

legend(   ((temppar$xaxp[2]-temppar$xaxp[1])/temppar$xaxp[3])*(temppar$xaxp[3]-1) +temppar$xaxp[1]
         ,temppar$yaxp[2]
         ,as.character(ArrayIndex),lt=1:length(ArrayIndex),col=1:length(ArrayIndex)
        ,cex=.5  
        )



}
