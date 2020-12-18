#!/usr/bin/Rscript
# Modified

usage = paste("",
              "Usage: StrainData.R netcdf4File locusId > output.table",
              "       StrainData.R netcdf4File scaffoldId begin end > output.table",
	      "    Given a netcdf4File that contains per-strain fitness values",
	      "    makes a tab-delimited table of the strains within that gene",
	      "    and their fitness values. Includes strains near the edges",
	      "    of the gene or with low coverage that were not used to",
	      "    estimate the gene's fitness (used=FALSE). If locusId is specified,",
              "    only includes strains with an assigned locusId (i.e. strains within",
              "    10-90% of the gene). If the begin and end positions are specified",
              "    instead, then all strains, including those without a locusId, are",
              "    contained.",
	      "", sep="\n");

library("ncdf4");

StrainData = function(args = commandArgs(trailingOnly=TRUE)) {
	source(file.path(GetPath(), "NetCDFInterface_Functions.R"));

	if (length(args)==2){
		filename = args[1];
		locusId = args[2];
	}else if (length(args)==4){
		filename = args[1];
		scaffoldId = args[2];
		begin = as.numeric(args[3]);
		end = as.numeric(args[4]);
		if (is.na(begin) || is.na(end) || end<begin || end<0){
			stop("Illegal coordinates: ", args[2], " ", args[3], " ", args[4])
		}
	}else{
		stop(usage);
	}
	nc_file = nc_open(filename);
	strainCols = dimensionsForVariable(nc_file,"fit/strains")[[3]]$vals;
	strains = ncvar_get(nc_file,"fit/strains");
	lrnCols = dimensionsForVariable(nc_file, "fit/strain_lrn")[[2]]$vals;
	outCols = c(strainCols, lrnCols);
	if (length(args)==2){
		locusIdCol = which(strainCols == "locusId");
		# values are sometimes padded with spaces
		strainLocus = sub("^ +", "", strains[,locusIdCol], perl=T);
		strainLocus = sub(" +$", "", strainLocus, perl=T);
		iRows = which(strainLocus == locusId);
	}else{
		scaffoldIdCol = which(strainCols == "scaffold")
		posCol = which(strainCols == "pos")
		strainScaffold = sub("^ +", "", strains[,scaffoldIdCol], perl=T);
		strainScaffold = sub(" +$", "", strainScaffold, perl=T);
		strainPos = sub("^ +", "", strains[,posCol], perl=T);
		strainPos = sub(" +$", "", strainPos, perl=T);	
		strainPos = suppressWarnings(as.numeric(strainPos))
		iRows = which(strainScaffold == scaffoldId & strainPos>=begin & strainPos<=end)
	}
	if (length(iRows) == 0) {
		# empty table
		write(outCols, stdout(), sep="\t", ncol=length(outCols));
	} else {
		values = sapply(iRows, function(i) ncvar_get(nc_file,"fit/strain_lrn",start=c(i,1),count=c(1,-1)));
		values = t(values) # now, 1 column per fitness experiment, one row per strain
		d = cbind(strains[iRows,], values);
		colnames(d) = outCols;
		for (i in 1:length(strainCols)) {
		    d[,i] = sub("^ +", "", d[,i], perl=T);
		    d[,i] = sub(" +$", "", d[,i], perl=T);
		}
		write.table(d, stdout(), quote=F, row.names=F, sep="\t");
	}
}

GetPath = function() {
    argv = commandArgs(trailingOnly = FALSE)
    dirname(substring(argv[grep("--file=", argv)], 8));
}

# Actually do the work
if(!interactive()) {
	StrainData();
	quit();
}
