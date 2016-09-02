
.onAttach <- function(libname, pkgname) {
	        
	erahV <- utils::packageVersion("erah")     
   
	msg <- paste("Welcome to eRah. This is an early release of eRah (V",erahV,"). For bugs, problems and issues, please use the eRah forum at: http://erah.lefora.com/. Describe your problem and include the output of sessionInfo(). \n \n NOTE OF CAUTION: eRah is currently not compatible with Waters equipment. This is due that Waters chromatograms seem to be downsampled when converted to NCDF, mzXML or mzML files.", sep="")
    
    #packageStartupMessage(msg)            
	
}    
