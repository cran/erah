compInfo <- function(comp.id, id.database=mslib)
{
	comp.info <- id.database@database[[comp.id]]
	comp.info$Synon <- gsub("//",", ", comp.info$Synon)
	cat("Name: ", comp.info$Name, "\n", "Synonyms: ", comp.info$Synon, "\n", "CAS: ", comp.info$CAS, "\n" ,"Formula: ", comp.info$Formula, "\n" ,"MW: ", comp.info$MW, "\n", "KEGG: ", comp.info$KEGG, "\n", "RI (FAME & Var5): ", comp.info$RI.VAR5.FAME, "\n","RI (ALK & Var5): ", comp.info$RI.VAR5.ALK, "\n","RI (FAME & MDN35): ", comp.info$RI.MDN35.FAME, "\n","RI (ALK & MDN35): ", comp.info$MDN35.ALK, "\n", "------------------------------- \n" ,"Comment: ", comp.info$Comment, sep="")
}

findComp <- function(name=NULL, id.database=mslib, CAS=NULL, chem.form=NULL)
{
	#Only one argument is allowed. If name is introduced, CAS and form, is depreciated, if CAS is introduced, chemical formula is depreciated.
	if(all(is.null(c(name,CAS,chem.form)))) stop("One argument is needed: Name, CAS or formula")
	if(length(which(c(is.null(name),is.null(CAS),is.null(chem.form)))==FALSE)<2) warning("Only one argument will be used! (Name prevails over CAS, and CAS prevails over Chemical Formula)")
	
	db.form <- unlist(lapply(id.database@database, function(x) x$Formula))
	db.cas <- unlist(lapply(id.database@database, function(x) x$CAS))
	db.names <- unlist(lapply(id.database@database, function(x) x$Name))

	if(!is.null(chem.form)) indexes <- grep(chem.form, db.form, ignore.case=T)	
	if(!is.null(CAS)) indexes <- grep(CAS, db.cas, ignore.case=T)	 
	if(!is.null(name)) indexes <- grep(name, db.names, ignore.case=T)

	met.list <- as.data.frame(matrix(c(indexes, db.names[indexes], db.cas[indexes], db.form[indexes]), ncol=4))
	colnames(met.list) <- c("DB.Id","Compound Name","CAS","Formula")
	met.list
}

export2MSP <- function(Experiment, export.id=NULL, id.database = mslib, store.path=getwd())
{
	#Experiment <- ex
	#export.id <- NULL
	#export.id <- c(4,5,10)
	
	#if(dir.exists(paste(c(store.path, "/ExportMSP/"), collapse=""))) stop("Please, delete the following folder before proceeding: ", paste(c(store.path, "/ExportMSP/")))

	
	if(is.null(nrow(Experiment@Results@Identification)) | nrow(Experiment@Results@Identification)==1)
	{
		if(!is.null(export.id)) Experiment@Results@Alignment <- Experiment@Results@Alignment[which(Experiment@Results@Alignment$AlignID %in% export.id),]
			
		SpectList <- sapply(Experiment@Results@Alignment$Spectra, function(x) {
			splitted.spectra.list <- strsplit(as.character(x), split = " ")[[1]]
			splitted.spectra.list <- gsub(",", " ", splitted.spectra.list)
			splitted.spectra.list <- as.character(as.vector(sapply(splitted.spectra.list, function(x) paste(c(x,";"), collapse=""))))
			Npeaks <- length(splitted.spectra.list)
			if(Npeaks>=6)
			{
				sequ <- seq(1, Npeaks, 5)
				if(sequ[length(sequ)]!=Npeaks) sequ <- c(sequ, Npeaks)
				splitted.spectra.list <- paste(unlist(sapply(1:(length(sequ)-1), function(x) paste(c(splitted.spectra.list[sequ[x]:sequ[(x+1)]], "\n")))), collapse=" ")
			}else{
				splitted.spectra.list <- paste(splitted.spectra.list, collapse=" ")
				}
			PeakChar <- paste(gsub(",", " ", splitted.spectra.list), collapse="; ")
			return(list(Npeaks=Npeaks, PeakChar=PeakChar))
		})
		SpectList <- split(SpectList,seq(NROW(SpectList)))

		SpectString <- as.vector(unlist(SpectList[[2]]))
		Npeaks <- as.vector(unlist(SpectList[[1]]))
		
		SpectNames.2 <- as.character(as.vector(Experiment@Results@Alignment$Factor))
		SpectNames.3 <- sapply(as.numeric(as.vector(Experiment@Results@Alignment$tmean)), function(x) paste("Rt:", x))
		
		SpectNames <- apply(cbind(SpectNames.3,SpectNames.2), 1, function(x) paste(x, collapse=" @ "))
	}else{
		
		if(!is.null(export.id)) Experiment@Results@Identification<- Experiment@Results@Identification[which(Experiment@Results@Identification$AlignID %in% export.id),]
		
		SpectList <- sapply(Experiment@Results@Identification$Spectra, function(x) {
			splitted.spectra.list <- strsplit(as.character(x), split = " ")[[1]]
			splitted.spectra.list <- gsub(",", " ", splitted.spectra.list)
			splitted.spectra.list <- as.character(as.vector(sapply(splitted.spectra.list, function(x) paste(c(x,";"), collapse=""))))
			Npeaks <- length(splitted.spectra.list)
			if(Npeaks>=6)
			{
				sequ <- seq(1, Npeaks, 5)
				if(sequ[length(sequ)]!=Npeaks) sequ <- c(sequ, Npeaks)
				splitted.spectra.list <- paste(unlist(sapply(1:(length(sequ)-1), function(x) paste(c(splitted.spectra.list[sequ[x]:sequ[(x+1)]], "\n")))), collapse=" ")
			}else{
				splitted.spectra.list <- paste(splitted.spectra.list, collapse=" ")
				}
			PeakChar <- paste(gsub(",", " ", splitted.spectra.list), collapse="; ")
			return(list(Npeaks=Npeaks, PeakChar=PeakChar))
		})
		SpectList <- split(SpectList,seq(NROW(SpectList)))

		SpectString <- as.vector(unlist(SpectList[[2]]))
		Npeaks <- as.vector(unlist(SpectList[[1]]))
		
        id.found <- as.numeric(as.vector(Experiment@Results@Identification[, "DB.Id.1"]))
		met.name <- unlist(lapply(id.database@database[c(id.found)], function(x) x$Name))
		
		SpectNames.1 <- sapply(as.character(as.vector(Experiment@Results@Identification$AlignID)), function(x) paste("(Align ID: ", x, ")", sep=""))
		SpectNames.2 <- sapply(as.numeric(as.vector(Experiment@Results@Identification$tmean)), function(x) paste("Rt:", x))
		SpectNames.3 <- apply(cbind(SpectNames.2,SpectNames.1), 1, function(x) paste(x, collapse=" "))
		SpectNames <- as.character(as.vector(apply(cbind(SpectNames.3,met.name), 1, function(x) paste(x, collapse=" @ "))))
	}
	
	dir.create(file.path(store.path, "ExportMSP") , showWarnings = FALSE)
	filename <- paste(c(store.path, "/ExportMSP/", "ExportedMSP", ".msp"), collapse="")
		
	fileTag <- character()	
	for(i in 1:length(SpectString))
	{ 
		fileTag[i] <- paste(c("Name: ", SpectNames[i], "\n", "Comments: MSP spectra exported by eRah \n", "Num Peaks: ", Npeaks[i], "\n",SpectString[i]), collapse="")
	}
	fileTagGen <- paste(fileTag, collapse="\n \n")

	writeLines(fileTagGen, filename)	
	cat("Spectra saved at: ", store.path, "/ExportMSP", sep="")
}




# seekSimilar <- function(comp.id, id.database=mslib, n=10)
# {
	# spect.list <- lapply(id.database@database, function(x) x$Spectra)
	# splitted.spectra.list <- lapply(spect.list, function(x) strsplit(as.character(x),split=" ")[[1]])
	# splitted.spectra.list <- lapply(splitted.spectra.list, function(c.bin){ 
		# out <- unlist(strsplit(c.bin,split=":"))
		# mz.ind <- seq(from=1, to=(length(out)-1),by=2)
		# int.ind <- seq(from=2, to=length(out),by=2)
		# list(mz=out[mz.ind], int=out[int.ind])
		# })
	# maxMz <- max(unlist(lapply(splitted.spectra.list, function(x) max(as.numeric(as.vector(x$mz))))))

	
	# spect.mat <- unlist(lapply(splitted.spectra.list, function(x) {
		# out <- rep(0,maxMz)
		# out[as.numeric(as.vector(x$mz))] <- as.numeric(as.vector(x$int))
		# out
	# }))
	# spect.mat <- matrix(spect.mat, nrow=maxMz)
	# cor.mat <- cor(spect.mat[,comp.id], spect.mat[,-comp.id])
	
	# corr.vector <- as.vector(cor.mat)
	# sim.index <- order(corr.vector, decreasing=T)[1:n]
	# sim.corr <- round(corr.vector[sim.index]*100, digits=2)
	
	# sim.list <- as.data.frame(NULL)

	# sim.names <- unlist(lapply(id.database@database[sim.index], function(x){x$Name}))
	# sim.cas <- unlist(lapply(id.database@database[sim.index], function(x){x$CAS}))
	# sim.form <- unlist(lapply(id.database@database[sim.index], function(x){x$Formula}))
		
	# sim.list <- as.data.frame(matrix(c(sim.index,sim.names,sim.corr,sim.cas,sim.form), nrow=n))
	# colnames(sim.list) <- c("DB.Id","Name","MatchFactor","CAS","Formula")
	# sim.list		
# }