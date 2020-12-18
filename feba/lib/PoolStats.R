#!/usr/bin/Rscript
# This is a command-line R script to analyze pool statistics
# The main line of analysis is in RunPoolStats()

usage = paste("",
        	"Usage: PoolStats.R pool genes.tab totalReads",
		"   If successful, writes metrics to pool.stats and writes",
		"   the set of unhit protein-coding genes (no good insertion) to pool.unhit",
		"", sep="\n");

RunPoolStats = function(args = commandArgs(trailingOnly=TRUE)) {
	if (length(args) != 3) stop(usage);
	poolfile = args[1];
	genesfile = args[2];
	nreadstot = args[3];
	# mode 4 means read permission
	if (file.access(poolfile,mode=4) != 0) stop("Cannot read pool file: ",poolfile);
	if (file.access(genesfile,mode=4) != 0) stop("Cannot read genes file: ",genesfile);
	if (is.na(as.numeric(nreadstot))) stop("nreads not numeric: ",nreadstot);
	nreadstot = as.numeric(nreadstot);

	# Loading and checking input files
	pool = read.delim(poolfile, as.is=T);
	if (is.null(pool) || nrow(pool) == 0) {
	    err_printf("No strains in pool\n");
	    return();
	}
	for (n in words("nTot n scaffold strand pos"))
		if(!n %in% names(pool)) stop("Missing column ", n, " from pool file ", poolfile);
	genes = read.delim(genesfile, as.is=T, quote="");
	for (n in words("scaffoldId begin end strand desc"))
		if(!n %in% names(genes)) stop("Missing column ", n, " from genes file ", genesfile);
	err_printf("%d insertions in genome are at %d different locations\n",
	     sum(pool$scaffold != "pastEnd"),
	     nrow(unique(pool[pool$scaffold != "pastEnd",words("scaffold strand pos")])));
	nSeen = sum(pool$scaffold != "pastEnd" & pool$n >= 2);
	nSeenTwice = sum(pool$scaffold != "pastEnd" & pool$n == 2);
	if (nSeenTwice > 0) diversity_cmp(nreadstot, nSeen, nSeenTwice, sum(pool$n[pool$scaffold != "pastEnd" & pool$n >= 2])/nreadstot);
	poolg = findWithinGrouped(split(pool, pool$scaffold),
				split(without(genes, genes$scaffoldId), genes$scaffoldId),
				"pos", "begin", "end");
	if (is.null(poolg) || nrow(poolg) == 0) {
	    err_printf("No insertions within genes\n");
	} else {
	    poolg$f = (poolg$pos - poolg$begin) / (poolg$end - poolg$begin);
	    PoolReport(poolfile, poolg, genes, nreadstot);
	}
}

# Separated out for convenient debugging
PoolReport = function(poolfile, poolg, genes, nreadstot) {
	if(is.null(genes$type)) genes$type = 1;

	# Insertions within central part of gene
	poolg2 = poolg[poolg$f >= 0.1 & poolg$f <= 0.9,];
	err_printf("Found %d insertions (%d distinct locations) in central 10-90%% of genes\n",
		nrow(poolg2), nrow(unique(poolg2[,words("scaffold pos strand")])));
	err_printf("Found central insertions for %d of %d protein-coding genes\n",
		sum(genes$type==1 & genes$locusId %in% poolg2$locusId),
		sum(genes$type==1));

	genes$expess = genes$type==1 & grepl("ribosomal protein|tRNA synthetase|cell division", genes$desc, ignore.case=T);
	genes$goodhit = genes$locusId %in% poolg2$locusId;
	err_printf("Hit rate in (crude) likely essentials: %.2f other %.2f\n",
			sum(genes$expess & genes$goodhit)/sum(genes$expess),
			sum(!genes$expess & genes$goodhit)/sum(!genes$expess));

	d = without(genes[genes$type==1 & genes$end-genes$begin+1 >= 300 & !genes$goodhit,],words("begin end strand goodhit"));
	d = d[order(d$desc),];
	unhitfile = paste(poolfile,".unhit",sep="");
	writeDelim(d, unhitfile);
	err_printf("Wrote proteins of 300nt or more with no good insertions to %s\n", unhitfile);
	surprisehitfile = paste(poolfile,".surprise",sep="");
	d = sort(genes$desc[genes$expess & genes$goodhit]);
	write(d, surprisehitfile, ncol=1);
	err_printf("Wrote %d genes with surprising insertions in central 10-90%% to %s\n",
			length(d), surprisehitfile);

	# tabulate strains per gene and reads per gene (only genes with hits)
	strainsPerGene = as.data.frame.table(table(poolg2$locusId, dnn="locusId"), responseName="nStrains");
	readsPerGene = aggregate(list(nReads=poolg2$n), list(locusId=poolg2$locusId), sum);
	d = genes[,words("locusId scaffoldId begin")];
	if(!is.null(genes$sysName)) d$sysName = genes$sysName;
	d$desc = genes$desc;
	d = merge(merge(d, strainsPerGene), readsPerGene);
	d = d[order(d$scaffold,d$begin),];
	d$begin = NULL;
	hitfile = paste(poolfile,".hit",sep="");
	writeDelim(d, hitfile);
	err_printf("Wrote read and strain counts for hit genes to %s\n",hitfile);

	nstrain = table(poolg2$locusId[poolg2$locusId %in% genes$locusId[genes$type==1]]);
	err_printf("Strains per hit protein: median %.0f mean %.1f\n", median(nstrain), mean(nstrain));

	plus = poolg2$strand == poolg2$strand.1; # gene and transposon in "same" orientation
	err_printf("Gene and transposon on same strand: %.1f%%\n", mean(plus)*100);

	nreads = aggregate(poolg2[,"n",drop=F], list(locusId=poolg2$locusId), sum);
	err_printf("Reads per hit protein: median %.0f mean %.1f bias (ratio) %.2f\n",
			median(nreads$n), mean(nreads$n), mean(nreads$n)/median(nreads$n));
	# Note reads per million uses total reads from TnSeq, not what is in the pool file
	err_printf("Reads per million for hit proteins: median %.2f mean %.2f\n",
			  median(nreads$n)*1e6/nreadstot, mean(nreads$n)*1e6/nreadstot);
}

diversity_cmp = function(nReads, nStrainsSeen, nSeenTwice, coverage) {
	nStrainsTot = round(nStrainsSeen / coverage);
	fOnce = dbinom(1, size=nReads, prob=1/nStrainsTot);
	fTwice = dbinom(2, size=nReads, prob=1/nStrainsTot);
	err_printf("Naive diversity (under)estimate: %d insertions in genome\n", nStrainsTot);
	err_printf("  coverage: %.3f vs. naive-expected %.3f\n", coverage, 1 - fOnce*nStrainsTot/nReads);
        err_printf("  fraction of strains seen just twice: %.2f vs. naive-expected %.2f\n",
	           nSeenTwice/nStrainsSeen, fTwice/(1-fOnce));
}

# Helper functions

# split a string separated by spaces into a list of word
words = function(s, by=" ", ...) { strsplit(s[1], by, ...)[[1]]; }

without <- function(list,columns=list$without) {
	list2 <- list;
	for (i in columns) { list2[[i]] <- NULL; }
	return(list2);
}

err_printf = function(...) cat(sprintf(...), file=stderr());

# by.x is a column name in x and begin and end are column names in y
# only returns the first matching row in y (as sorted by begin) for each row in x
# but if you specify all="y" then it returns the first matching row in x for
# each row in y
# If you specify unique=T then tries to find cases where y is the only match for that row in x
#	but it still returns some duplicates (e.g. it is fooled if row yA is within row yB)
# If minmax=T, allows begin > end (takes pairwise minima and maxima) and writes "left" and "right"
findWithin = function(x, y, by.x, begin, end, all="x", unique=FALSE, minmax=F) {
	if (nrow(x) == 0 || nrow(y) == 0) return(NULL); # one input is empty

	if(is.null(y[[begin]])) stop("no field named",begin," in y argument");
	if (minmax) {
		y$left = pmin(y[[begin]], y[[end]]);
		y$right = pmax(y[[begin]], y[[end]]);
		begin = "left";
		end = "right";
	}
	if (all=="x") {
		y2 = y[order(y[[begin]]),];
		if (nrow(y) == 1) {
			floorI = rep(1,nrow(x));
		} else {
			floorI = floor(approx(y2[[begin]],1:nrow(y2), xout=x[[by.x]], rule=2, ties=min)$y);
		}
		outy = y2[floorI,];
		keepx = x[[by.x]] >= outy[[begin]] & x[[by.x]] <= outy[[end]];
		if (unique) {
			floorI2 = ceiling(approx(y2[[end]],1:nrow(y2), xout=x[[by.x]], rule=2, ties=max)$y);
			keepx = keepx & floorI==floorI2;
			# Is this where it gets fooled?
		}
		return(data.frame(cbind(x[keepx,], outy[keepx,])));
	} else if (all == "y") {
		if(unique) stop("unique not supported with all=y");
		x2 = x[order(x[[by.x]]),];
		floorI = floor(approx(x2[[by.x]],1:nrow(x2), xout=y[[end]], rule=2,
			ties='ordered')$y);
		outx = x2[floorI,];
		keepy = outx[[by.x]] >= y[[begin]] & outx[[by.x]] <= y[[end]];
		return(data.frame(cbind(outx[keepy,], y[keepy,])));
	} else {
		stop("Unknown value of option all in findWithin",all);
	}
}

# given that data is subgrouped by the same markers, do findWithin on each and merge the results
findWithinGrouped = function(splitx, splity, by.x, begin, end, debug=FALSE, ...) {
	out = NULL;
	for(i in names(splitx)) {
		if (!is.null(splity[[i]])) {
			if(debug) cat("Running group ",i,"\n");
			rows = findWithin(splitx[[i]], splity[[i]], by.x, begin, end, ...);
			if(debug) cat("Ran group ",i,"rows",nrow(rows),"\n");
			if(!is.null(rows) && nrow(rows) > 0) {
				out = if(is.null(out)) rows else rbind(out,rows);
			}
		}
	}
	return(out);
}

writeDelim <- function(table,file,report=FALSE,...) {
	write.table(table,file,sep="\t",quote=FALSE,row.names=FALSE,...);
	if(report) cat("Wrote",nrow(table),"rows","to",file,"\n")
}

# Actually do the work
if(!interactive()) {
	RunPoolStats();
	quit();
}
