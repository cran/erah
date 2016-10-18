get.compound.info <- function(k, Spl.List, type)
{	
		#kmain <<-k
		x <- Spl.List[[k]]

		x.split <- unlist(apply(as.matrix(x),1, function(x) {strsplit(x,":")[[1]]}), recursive=T)
	
		Name=" "
		Synon=" "
		RI.VAR5.FAME=0
		RI.VAR5.ALK=0
		RI.MDN35.FAME=0
		RI.MDN35.ALK=0
		CAS=" "
		Formula=" "
		InChi=" "
		KEGG=""
		MW=0
		Comment=" "
		GMD.LINK=" "
		GMD.VERS=" "
		
		Spectra <- ""
		Caps <- apply(as.matrix(x),1, function(x) {strsplit(x,":")[[1]]})
		Sp.Ini <- which(lapply(Caps,function(y)y[1])=="Num Peaks")
		Spectra <- paste(x[(Sp.Ini+1):length(x)],collapse="")
		
		for(j in 1:length(x.split))
		{
			if(x.split[j]==" MST N") Name=x.split[j+1]
			if(x.split[j]==" METB N") {
				if(Synon!=" ") Synon=paste(Synon, "//", x.split[j+1], sep=" ")
				if(Synon==" ") Synon = x.split[j+1]
			}
			if(x.split[j]==" RI")
			{
				if(type=="MDN35.ALK") RI.MDN35.ALK=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="MDN35.FAME") RI.MDN35.FAME=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="VAR5.ALK") RI.VAR5.ALK=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
				if(type=="VAR5.FAME") RI.VAR5.FAME=as.numeric(apply(as.matrix(x.split[j+1]), 2, gsub, patt=",", replace="."))
			}
			if(x.split[j]==" MST SEL MASS") SelMZ <- as.vector(sapply(strsplit(x.split[j+1], "\\|"), as.numeric))
			if(x.split[j]=="CAS#") CAS=x.split[j+1]
			if(x.split[j]=="Formula") Formula=x.split[j+1]
			if(x.split[j]==" METB InChI") InChi= strsplit(x.split[j+1], "=")[[1]][2]
			if(x.split[j]==" METB KEGG") KEGG=x.split[j+1]			
			if(x.split[j]=="MW")
			{
				no.dot <- apply(as.matrix(x.split[j+1]), 2, gsub, patt="\\.", replace="")
				MW=as.numeric(apply(as.matrix(no.dot), 2, gsub, patt=",", replace="."))
			}
			if(x.split[j]=="Comment") Comment=x.split[j+1]
			if(x.split[j]==" GMD LINK") GMD.LINK=paste(x.split[j+1],":",x.split[j+2],sep="")
			if(x.split[j]==" GMD VERS") GMD.VERS= paste(x.split[j+1],":",x.split[j+2],sep="")		
		}
		#cat(Name, "\n")
		s.spl <- strsplit(Name," ")[[1]]; Name <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Synon," ")[[1]]; Synon <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(CAS," ")[[1]]; CAS <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Formula," ")[[1]]; Formula <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(KEGG," ")[[1]]; KEGG <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(Comment," ")[[1]]; Comment <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(GMD.LINK," ")[[1]]; GMD.LINK <- paste(s.spl[-1], collapse=" ")
		s.spl <- strsplit(GMD.VERS," ")[[1]]; GMD.VERS <- paste(s.spl[-1], collapse=" ")
			
		compound.info <- list(Name=Name, Synon=Synon, RI.VAR5.FAME=RI.VAR5.FAME, RI.VAR5.ALK=RI.VAR5.ALK, RI.MDN35.FAME=RI.MDN35.FAME, RI.MDN35.ALK=RI.MDN35.ALK, SelMZ=SelMZ, CAS=CAS, InChi=InChi, Formula=Formula, MW=MW, KEGG=KEGG, Comment=Comment, GMD.LINK=GMD.LINK,	GMD.VERS=GMD.VERS, Spectra=Spectra)
		
	compound.info
}

list.DB <- function(DB.object)
{
		k <- 1
		Cont <- 1
		Spl.List <- list()
		for(i in 1:length(DB.object))
		{
			if(DB.object[i]=="") 
			{	
				Spl.List[Cont] <- list(Compound=DB.object[k:(i-1)])
				k <- i + 1
				Cont <- Cont + 1
			}
		}
		Spl.List
}
	
	
importMSP <- function(filename, DB.name, DB.version, DB.info, type=c("VAR5.ALK","VAR5.FAME","MDN35.ALK", "MDN35.FAME"))
{
	DB.MSP <- readLines(filename)
		
	Spl.List <- list.DB(DB.MSP)
		
	import.database <- list()
	import.database <- lapply(1:length(Spl.List), function(x) get.compound.info(x, Spl.List, type))
	
	final.database <- new("eRah_DB",  name=DB.name, version=DB.version, info=DB.info, database=import.database)
	final.database
}	