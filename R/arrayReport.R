##  arrayReport.R
##  09-Jun-2005
##
##  Craig Parman craig.parman@bifix.org
##  Conrad Halling conrad.halling@bifx.org

QCReport <-
function(
    object,
    file    = "AffyQCReport.pdf",
    ... )
{
    pdf(
        file    = file,
        width   = 8,
        height  = 11,
        onefile = TRUE )
    plot.window(
        c( 1, 1 ),
        c( 0, 1 ) )
    plot.new()
    titlePage( object )
    signalDist( object )
    plot( qc( object ) )
    borderQC1( object )
    borderQC2( object )
    correlationPlot( object )
    dev.off()
    return( TRUE )
}
