#!/usr/bin/Rscript
# A command-line wrapper for running FEBA.R
# Can be invoked from the command-line to analyze the fitness experiments and save the results
# (i.e., FEBA_Fit() and FEBA_Save_Tables())

usage = paste("",
        	"Usage: RunFEBA.R orgname data_directory FEBA_dir",
		"   Compute fitness values from genes, exps, pool, and all.poolcount",
		"   in the data directory.",
		"   If successful, creates a mini web site in the data_directory,",
		"   including tables of per-gene and per-strain fitness values,",
		"   plots and quality metrics to check if the experiments,",
		"   and an R image with everything from the complex fit data structure.",
		"   (RunFeba.R is normally invoked via BarSeqR.pl)",
		"", sep="\n");

RunFEBA = function(args = commandArgs(trailingOnly=TRUE)) {
	if (length(args) != 3) stop(usage);
	org = args[1];
	dir = args[2]
	FEBAdir = args[3];

	allfile = paste(dir,"/all.poolcount",sep="");
	genesfile = paste(dir,"/genes",sep="");
	expsfile = paste(dir,"/exps",sep="");
	poolfile = paste(dir,"/pool",sep="");
	FEBA_R = paste(FEBAdir,"/lib/FEBA.R",sep="");

	# mode 4 means read permission
	if (file.access(allfile, mode=4) != 0) stop("Cannot read all file: ",allfile);
	if (file.access(genesfile, mode=4) != 0) stop("Cannot read genes file: ",genesfile);
	if (file.access(expsfile, mode=4) != 0) stop("Cannot read exps file: ",expsfile);
	if (file.access(poolfile, mode=4) != 0) stop("Cannot read pool file: ",poolfile);
	if (file.access(FEBAdir, mode=4) != 0) stop("Cannot access ",FEBAdir);
	if (file.access(FEBA_R, mode=4) != 0) stop("Cannot find ",FEBA_R);
	source(FEBA_R);

	rules = read.table(paste(FEBAdir,"/lib/desc_short_rules",sep=""),as.is=T);

	genes = read.delim(genesfile,quote="",as.is=T);
        for (n in c("scaffoldId","locusId","sysName","desc"))
          if(!n %in% names(genes)) stop("genes table must include field ", n);
	all = read.delim(allfile,as.is=T,check.names=F);
	exps = read.delim(expsfile,as.is=T);
	d = unique(all$locusId[all$locusId != "" & !is.na(all$locusId)]);
	if (!all(d %in% genes$locusId)) stop("Unknown genes in ",allfile);
	cat(sprintf("Read %d genes, %d exps, and data for %d barcodes\n",
	            nrow(genes), nrow(exps), nrow(all) ));

	# fix up names of all to be shorter
	SetNames = unique(exps$SetName);
	SetNames2 = ShortSetNames(SetNames);
	expNamesNew = paste(SetNames2[match(exps$SetName, SetNames)], exps$Index, sep="");
	
	exps$num = 1:nrow(exps);
	exps$name = paste(exps$SetName,exps$Index,sep=".");
	exps$short = applyRules(rules, exps$Description);
	metacol = 1:7;
	name_list = names(all)[-metacol];
	if(length(name_list) != nrow(exps)) stop("Number of data columns in  ",allfile," does not match number of rows in ",expsfile);
	if(any(name_list != exps$name)) stop("Column names in  ",allfile," do not match names from ",expsfile);

	# remove trailing spaces from Group, Condition_1, Condition_2
	for(n in c("Group","Condition_1","Condition_2")) {
	    if(!is.null(exps[[n]])) exps[[n]] = sub(" +$", "", exps[[n]]);
	}
	names(all)[-metacol] = expNamesNew;
	exps$name = expNamesNew;

        # load strains to use, if this data exists
        strainsUsed = NULL;
        genesUsed = NULL;
        genesUsed12 = NULL;
        fStrainsUsed = paste(dir,"/strainusage.barcodes",sep="");
        fGenesUsed = paste(dir,"/strainusage.genes",sep="");
        fGenesUsed12 = paste(dir,"/strainusage.genes12",sep="");
        # Uses the strain usage files if they exist
	if (file.access(fStrainsUsed, mode=4) == 0) {
        	stopifnot(file.access(fGenesUsed,mode=4)==0);
        	stopifnot(file.access(fGenesUsed12,mode=4)==0);
		barcodesUsed = scan(fStrainsUsed,"");
                strainsUsed = all$barcode %in% barcodesUsed;
		genesUsed = scan(fGenesUsed,"");
		genesUsed12 = scan(fGenesUsed12,"");
                cat(sprintf("Loaded %d strains and %d genes to include in the analysis\n",
                	    length(barcodesUsed), length(genesUsed)));
	}

	options(width=100);
	fit = FEBA_Fit(exps, all, genes, dir=dir,
		strainsUsed=strainsUsed, genesUsed=genesUsed, genesUsed12=genesUsed12);
	FEBA_Save_Tables(fit, genes, org, expsU=exps, dir=dir, FEBAdir=FEBAdir);
	write(date(), paste(dir,"/.FEBA.success", sep="")); # report success to BarSeqR.pl
}

ShortSetNames = function(sets) {
	simple = grepl("(set|test)[0-9A-Z]+[0-9A-Z0-9]*$", sets);
	sets[simple] = sub("^.*(set|test)", "\\1", sets[simple]);
	nleft = sum(!simple);
	candidates = strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]];
	candidates = paste("set",candidates,sep="");
	candidates = setdiff(candidates, sets[simple]);
	if (nleft > length(candidates)) stop("Too many wierd set names");
	oldComplex = sets[!simple];
	sets[!simple] = candidates[1:nleft];
	if(nleft > 0) for (i in 1:nleft) {
	   write(sprintf("Set %s simplified to %s", oldComplex[i], sets[!simple][i]), stderr());
	}
	if (length(sets) != length(unique(sets))) stop("Non-unique sets");
	return(sets);
}

# Actually do the work
if(!interactive()) {
	RunFEBA();
	quit();
}
