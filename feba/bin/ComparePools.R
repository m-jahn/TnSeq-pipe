#!/usr/bin/Rscript
# This is a command-line R script that reports how similar
# two pool files are. The expectation is that these pool files 
# are from independent sequencing runs on the same library
# of mutants.

usage = paste("",
              "Usage: ComparePools.R pool1 pool2",
              "    Creates a short report of how similar two pool files are.",
              "", sep="\n");

ComparePools = function(p1,p2) {
	comb = merge(p1, p2, by=c("barcode","scaffold","pos"));
	all = merge(p1, p2, by="barcode");
        nAgree = nrow(comb);
        cat(sprintf("Shared strains: %d of %d or %d (%.1f%% or %.1f%%)\n",
                    nAgree, nrow(p1), nrow(p2),
                    nAgree*100.0/nrow(p1),
                    nAgree*100.0/nrow(p2) ));
        cat(sprintf("Coverage of shared strains: %.1f%% or %.1f%% of reads\n",
                    100.0*sum(comb$nTot.x)/sum(p1$nTot),
                    100.0*sum(comb$nTot.y)/sum(p2$nTot) ));
        cat(sprintf("Discrepancies: %d barcodes, %.3f%% or %.3f%% of reads\n",
                    nrow(all) - nrow(comb),
                    100.0*(sum(all$nTot.x)-sum(comb$nTot.x))/sum(p1$nTot),
                    100.0*(sum(all$nTot.y)-sum(comb$nTot.y))/sum(p2$nTot) ));
}

if(!interactive()) {
	args = commandArgs(trailingOnly=TRUE);
        if (length(args) != 2) stop(usage);
        for (i in 1:2) {
            # mode 4 means read
	    if (file.access(args[i],mode=4) != 0) stop("Cannot read pool file: ", args[i]);
	}
        p1 = read.delim(args[1], as.is=T);
        p2 = read.delim(args[2], as.is=T);
        fields = c("barcode","nTot","scaffold","strand","pos");
        for (name in fields) {
		stopifnot(!is.null(p1[[name]]));
		stopifnot(!is.null(p2[[name]]));
	}
        cat("Loaded pools\n");
        ComparePools(p1,p2);
}
