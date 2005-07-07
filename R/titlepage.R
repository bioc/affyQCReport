
         titlePage<- function(object)


{
#plot.new()

polygon(c(.0,.0,.9,.9,.0),c(.05,.95,.95,.05,.05))
polygon(c(.45,.45),c(.05,.95))

polygon(c(0,.9),c(.88,.88))
polygon(c(0,.9),c(.87,.87))

text(.5,1,"AffyBatch QC Report",cex=1.5)
text(.4,0,"Produced by AffyQCReport R Package", cex=.5)
text(.4,.03,date(),cex=.5)
text(.45/2, .9, "Array Index")
text(3*.45/2, .9, "Array Name")



n<-length(sampleNames(object))
cexval <- min(10/n,1)
	for( i in 1:n)
	{	

	text(.45/2, .9-i*(.9-.01)/(n+1)  ,as.character(i),cex=cexval)
      
      polygon(c(0,.9),c(.9-(i+.5)*(.9-.01)/(n+1),.9-(i+.5)*(.9-.01)/(n+1))  )

	text(3*.45/2,.9-i*(.9-.01)/(n+1), sampleNames(object)[i],cex=cexval)
     

	}

return(TRUE)
}



