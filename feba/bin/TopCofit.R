#!/usr/bin/Rscript
# This is a command-line R script to compute a cofitness table, with "top hits" only, from
# the fitness data.

usage = paste("",
      	"Usage: TopCofit.R fit_logratios_good.tab cofit",
	"    If successful, writes the top 2% of most-cofit genes for each query gene",
	"    to cofit",
	"    Assumes that the fitness data includes locusId, and optionally",
	"    other metadata fields sysName, desc, comb; all other fields should be data fields.",
	"", sep="\n");

RunTopCofit = function(args = commandArgs(trailingOnly=TRUE)) {
	if (length(args) != 2) stop(usage);
	fitfile = args[1];
	outfile = args[2];

	# mode 4 means read permission
	if (file.access(fitfile,mode=4) != 0) stop("Cannot read fitness data file: ",fitfile);
	fit = read.delim(fitfile,as.is=T,check.names=F); # make sure locusId is not converted to a factor
	if (is.null(fit$locusId)) stop("No locusId in fitness data file");
	locusId = fit$locusId;
	# remove metadata
	fit$locusId = NULL;
	fit$sysName = NULL;
	fit$desc = NULL;
	fit$comb = NULL;
	if(ncol(fit) < 5) stop("Cannot run TopCofit.R with fewer than five data columns");
	# all others should be floating point
	for (n in names(fit)) {
	    if(!is.double(fit[[n]])) {
		stop("TopCofit.R -- input column ", n, " is not floating-point");
	    }
            if(any(is.na(fit[[n]]))) stop("TopCofit.R -- input column ", n, " contains missing values");
	}
	writeDelim(TopCofit(locusId, fit), outfile);
}

# g is genes (i.e., locusIds)
# r is a matrix of fitness values with columns as experiments
TopCofit = function(g, r,
		debug = FALSE,
		fraction=0.02,
		n = pmin(pmax(1,round(length(g)*fraction)),length(g)-1)) {
	if (length(g) != nrow(r)) stop("rows and number of genes does not match");
	cofits = cor(t(r)); # correlation between rows; note does not check for NA (is much faster that way)
	if(debug) cat("nTop",n,"\n");
	nOut = length(g) * n;
	if(debug) cat("Making output with",nOut,"rows\n");
        out_hitId = rep("", nOut);
        out_cofit = rep(NA, nOut);
	for (i in 1:length(g)) {
		values = cofits[i,];
		j = order(-values)[2:(n+1)]; # assume self is in position 1
		outi = (i-1)*n + (1:n); # where to put inside out
		out_hitId[outi] = g[j];
		out_cofit[outi] = values[j];
	}
	out = data.frame(locusId=rep(g,each=n), hitId=out_hitId, cofit=out_cofit, rank=rep(1:n,length(g)));
	return(out);
}

writeDelim <- function(table,file,report=FALSE,...) {
	write.table(table,file,sep="\t",quote=FALSE,row.names=FALSE,...);
	if(report) cat("Wrote",nrow(table),"rows","to",file,"\n")
}

if(!interactive()) {
	RunTopCofit();
	quit();
}
