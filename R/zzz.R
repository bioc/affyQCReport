.First.lib <-
function( libname, pkgname )
{
    require( simpleaffy, quietly = TRUE );
    require( affy,       quietly = TRUE );
    cat( "Welcome to 'affyQCReport' V 1.6-1\n" );
    cat( "Further information available at: www.bifix.org\n");
    cat( "      mailto: craig.parman@bifix.org\n" );
}
