# FEBA.R -- analysis scripts for barcode sequencing data
#
# Uses mclapply() from the R parallel package to analyze experiments in parallel --
# set MC_CORES to control the #CPUs used (default is 2).
library(parallel);

# The key routines are:
#
# FEBA_Fit() -- analyze many fitness experiments with AvgStrainFitness() and NormalizeByScaffold()
#      returns a complex data structure
# FEBA_Save_Tables() -- Save the fitness data structure to a mini web-site, tab-delimited files and an R image
#
# Also see:
# AvgStrainFitness -- compute fitness values from counts for the post-experiment counts
#     and the Time0 counts
# NormalizeByScaffold -- normalize the fitness values by scaffold and position
# GeneFitness -- combines the above two, and also computes a t-like test statistic ("t").
#	To do this, it also computes fitness values for the 1st and 2nd half of most genes
#
# Limitations:
# High memory usage (~10GB to process 200K strains x 500 experiments)


# GeneFitness():
# genes -- must include locusId, scaffoldId, and begin
# strainInfo -- must include locusId and (unless use1 is overridden) f, the fraction of the gene
# 	     that the insertion is at
# countCond and countT0 -- counts for each strain
# strainsUsed & genesUsed -- see AvgStrainFitness()
# genesUsed12 -- ditto, for 1st and 2nd half fitness values
# use1 -- which strains are in 1st half (regardless of whether they are usable or not)
# other arguments are passed on to AvgStrainFitness()
# base_se -- likely amount of error in excess of that given by variation within fitness values
# 	for strains in a gene, due to erorrs in normalization or bias in the estimator
#
# Returns a data frame with a row for each gene in genesUsed. It includes
# locusId,
# fit (unnormalized), fitnorm (normalized),
# fit1 or fit2 for 1st- or 2nd-half (unnormalized, may be NA),
# fitnorm1 or fitnorm2 for normalized versions,
# se (estimated standard error of measurement), and t (the test statistic),
# as well as some other values from AvgStrainFitness(), notably sdNaive,
# which is a different (best-case) estimate of the standard error.
GeneFitness = function(genes, strainInfo, countCond, countT0,
	    	    strainsUsed, genesUsed, genesUsed12,
		    use1 = strainInfo$f < 0.5,
		    base_se = 0.1,
		    minGenesPerScaffold=10,
		    ...) {
    d = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
      		         strainsUsed=strainsUsed, genesUsed=genesUsed, ...);
    d$fitnorm = NormalizeByScaffold(d$fit, d$locusId, genes, minToUse=minGenesPerScaffold);

    d1 = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
      		         strainsUsed=strainsUsed & use1, genesUsed=genesUsed12, ...);
    d2 = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
      		         strainsUsed=strainsUsed & !use1, genesUsed=genesUsed12, ...);
    if(any(as.character(d1$locusId) != as.character(d2$locusId))) stop("Non-matching locusId");

    i = match(d$locusId, d1$locusId);

    d$fit1 = d1$fit[i];
    d$fit2 = d2$fit[i];
    d$fitnorm1 = d$fit1 + (d$fitnorm-d$fit);
    d$fitnorm2 = d$fit2 + (d$fitnorm-d$fit);
    d$tot1 = d1$tot[i];
    d$tot1_0 = d1$tot0[i];
    d$tot2 = d2$tot[i];
    d$tot2_0 = d2$tot0[i];

    # for low n, the estimated variance is driven by the overall variance, which can be estimated
    # from the median difference between 1st and 2nd halves via the assumptions
    # Var(fit) = Var((fit1+fit2)/2) ~= Var(fit1-fit2)/4
    # median abs(normal variable) = qnorm(0.75) * sigma = 0.67 * sigma
    # which leads to Var(fit) = Var(fit1-fit2)/4
    # = sigma12**2/4 = median abs diff**2 / (qnorm(0.75)*2)**2
    # The median difference is used because a few genes may have genuine biological differences
    # between the fitness of the two halves.
    # Furthermore, assuming that genes with more reads are less noisy, this
    # pseudovariance should be rescaled based on sdNaive**2
    #
    pseudovar_std = median(abs(d$fit1-d$fit2),na.rm=T)**2 / (2*qnorm(0.75))**2;
    d$pseudovar = pseudovar_std * (d$sdNaive / median(d$sdNaive[!is.na(d$fit1)]))**2;
    # given the variable weighting in sumsq, it is not intuitive that the degrees of freedom is still n-1
    # however, this is the result given the assumption that the weighting is the inverse of the variance
    est_var = (d$pseudovar + d$sumsq)/d$n;
    d$se = sqrt(est_var);
    d$t = d$fitnorm/sqrt(base_se**2 + pmax(d$sdNaive**2, est_var));
    return(d);
}

# AvgStrainFitness():
# strainCounts -- counts at the end of the experiment condition
# strainT0 -- counts for Time0 for each strain
# strainLocus -- which locus the strain is associated with, or NA
#
# If genesUsed (as a list of locusId) and strainsUsed (as boolean vector) are provided,
# then considers only those strains & genes; minimum requirements.
#
# if returnStrainInfo is set, returns a list of two data frames, "genes" and "strains"
# normally returns a per-gene data frame.
#
# debug is for testing purposes.
AvgStrainFitness = function(strainCounts, strainT0, strainLocus,
		 minStrainT0 = 4, minGeneT0 = 40,
		 genesUsed=NULL, strainsUsed=NULL,
		 # maxWeight of N corresponds to having N reads on each side (if perfectly balanced); use 0 for even weighting
		 # 20 on each side corresponds to a standard error of ~0.5; keep maxWeight low because outlier strains
		 # often have higher weights otherwise.
		 # Use maxWeight = 0 for unweighted averages.
		 maxWeight = 20,
		 minGeneFactorNStrains=3,
		 returnStrainInfo=FALSE,
		 debug=FALSE) {
    if (length(strainCounts) < 1 || length(strainT0) < 1 || length(strainLocus) < 1
        || length(strainCounts) != length(strainT0) || length(strainCounts) != length(strainLocus))
        stop("No or misaligned input data");

    if (is.null(strainsUsed)) strainsUsed = strainT0 >= minStrainT0;
    if (is.null(genesUsed)) {
        geneT0 = aggregate(strainT0[strainsUsed], list(locusId=strainLocus[strainsUsed]), sum);
        genesUsed = geneT0$locusId[geneT0$x >= minGeneT0];
    }
    strainsUsed = strainsUsed & strainLocus %in% genesUsed;
    if (!any(strainsUsed)) stop("No usable strains");

    strainT0 = strainT0[strainsUsed];
    strainCounts = strainCounts[strainsUsed];
    strainLocus = strainLocus[strainsUsed];

    readratio = sum(strainCounts) / sum(strainT0);
    # use sqrt(readratio), or its inverse, instead of 1, so that the expectation
    # is about the same regardless of how well sampled the strain or gene is
    strainFit = mednorm(log2(sqrt(readratio) + strainCounts) - log2(1/sqrt(readratio) + strainT0));
    strainFitAdjust = 0;

    # Per-strain "smart" pseudocount to give a less biased per-strain fitness estimate.
    # This is the expected reads ratio, given data for the gene as a whole
    # Arguably, this should be weighted by T0 reads, but right now it isn't.
    # Also, do not do if we have just 1 or 2 strains, as it would just amplify noise

    # note use of as.vector() to remove names -- necessary for speed
    strainLocusF = as.factor(strainLocus);
    nStrains = table(strainLocusF);
    if(!all(names(nStrains)==levels(strainLocusF))) stop("Strain mismatch");
    nStrains = as.vector(nStrains);
    geneFit1 = mednorm(as.vector(tapply(strainFit, strainLocusF, median))); # used to use mean
    i = as.integer(strainLocusF); # from strain index to gene index
    strainPseudoCount = ifelse(nStrains[i] >= minGeneFactorNStrains, 2**geneFit1[i] * readratio, readratio);

    # And apportion the pseudocount equally (in log space) between condition-count and strain-count
    # to minimize the deviations from pseudocount = 1
    condPseudoCount = sqrt(strainPseudoCount);
    t0PseudoCount = 1/sqrt(strainPseudoCount);
    # (or could do some sort of weighted likelihood-based inference of fitness values, might be better)

    # for each strain: fitness, s.d., and weight
    strainFit = log2(condPseudoCount + strainCounts) - log2(t0PseudoCount + strainT0) - strainFitAdjust;
    strainSd = sqrt(1/(1+strainT0) + 1/(1+strainCounts)) / log(2);
    # use harmonic mean for weighting; add as small number to allow maxWeight = 0.
    strainWeight = 0.5 + pmin(maxWeight, 2/( 1/(1+strainT0) + 1/(1+strainCounts) ) );

    fitness = lapply(split(1:length(strainT0), list(locusId=strainLocus)),
     	           function(j) {
		       n = length(j);
                       totw = sum(strainWeight[j]);
		       fitRaw = sum(strainWeight[j] * strainFit[j]) / totw;
		       tot = sum(strainCounts[j]);
		       tot0 = sum(strainT0[j]);
		       sd = sqrt(sum(strainWeight[j]**2 * strainSd[j]))/totw;
		       sumsq = sum(strainWeight[j] * (strainFit[j]-fitRaw)**2)/totw;
		       # high-N estimate of the noise in the log2 ratio of fitNaive
		       # But sdNaive is actually pretty accurate for small n -- e.g.
		       # simulations with E=10 on each side gave slightly light tails
		       # (r.m.s.(z) = 0.94).
		       sdNaive = sqrt( 1/(1+tot) + 1/(1+tot0) ) / log(2);
		       nEff = totw/max(strainWeight[j]);
		       c(fitRaw=fitRaw, sd=sd, sumsq=sumsq, sdNaive=sdNaive, n=n, nEff=nEff,
		         tot=tot, tot0=tot0);
		});
    fitness = data.frame(do.call(rbind, fitness));
    fitness$fit = mednorm(fitness$fit);
    fitness$fitNaive = mednorm(log2(1+fitness$tot) - log2(1+fitness$tot0));
    fitness$locusId = row.names(fitness);
    if (is.integer(strainLocus)) fitness$locusId = as.integer(as.character(fitness$locusId));

    if(returnStrainInfo) return(list(genes=fitness,
        strains=data.frame(strainLocusF,strainCounts,strainT0,strainPseudoCount,strainFit,strainSd,strainWeight)));
    # else
    return(fitness);
}

# NormalizeByScaffold():
# values -- fitness values (as a vector)
# locusId -- the corresponding locusIds
# genes contains locusId, scaffoldId, and begin
# window -- window size for smoothing by medians. Must be odd, default 251. For scaffolds
#     with fewer genes than this, just uses the median.
# minToUse -- if a scaffold has too few genes, cannot correct for possible DNA extraction bias
# 	   so need to remove data for that gene (i.e., returns NA for them).
# returns
# normalized -- data with scaffold and position effects removed
#
NormalizeByScaffold = function(values, locusId, genes, window=251, minToUse=10, debug=FALSE) {
    i = match(locusId, genes$locusId);
    if(any(is.na(i))) stop("Fitness data for loci not in genes");
    beg = genes$begin[i];

    perScaffoldRows = split(1:length(values), genes$scaffoldId[i]);
    for (scaffoldId in names(perScaffoldRows)) {
        rows = perScaffoldRows[[scaffoldId]];
        if (length(rows) < minToUse) {
            if(debug) cat("Removing ",length(rows)," values for ", scaffoldId, "\n");
	    values[rows] = NA;
        } else {
  	    med = median(values[rows]);
	    if(debug) cat("Subtract median for ", scaffoldId, " ", med, "\n");
	    values[rows] = values[rows] - med;

	    if (length(rows) >= window) {
	        # Used to use lowess
                # d = lowess(beg[rows], values[rows]);
	        # if(debug) cat("Subtract loess for ", scaffoldId, " max effect is ", diff(range(d$y)), "\n");
		# values[rows] = values[rows] - approx(d$x,d$y, xout=beg[rows], rule=2, ties="ordered")$y;

		o = order(beg[rows]);
		m = runmed(values[rows[o]], window, endrule="constant");
		if(debug) cat("Subtract smoothed median for ", scaffoldId, " max effect is ",diff(range(m)), "\n");
		values[rows[o]] = values[rows[o]] - m;

	        d = density(values[rows]);
	        mode = d$x[which.max(d$y)];
	        if (debug) cat("Subtract mode for ", scaffoldId, " which is at ", mode, "\n");
	        values[rows] = values[rows] - mode;
            }
        }
    }
    return(values);
}

# simple log-ratio with pseudocount (of 1) and normalized so each scaffold has a median of 0
# note is *not* normalized except to set the total median to 0
StrainFitness = function(count, countT0, scaffolds) {
    fit = mednorm( log2(1+count) - log2(1+countT0) );
    se = sqrt(1/(1+count) + 1/(1+countT0)) / log(2);
    return(data.frame(fit=fit,se=se));
}

# For each strain, find the closest gene, as a row number -- returns a vector
# If there is no gene on that scaffold, returns NA
StrainClosestGenes = function(strains, genes) {
	genes$index = 1:nrow(genes);
	strainSplit = split(strains, strains$scaffold);
	geneSplit = split(genes, genes$scaffold);
	indexSplit = list();
	for (sc in names(strainSplit)) {
		s = strainSplit[[sc]];
		g = geneSplit[[sc]];
		if (is.null(g)) {
		    indexSplit[[sc]] = rep(NA, nrow(s));
		} else if (nrow(g) == 1) {
		    # cannot approx with 1 value so:
		    indexSplit[[sc]] = rep(g$index[1], nrow(s));
		} else {
		    g$pos = (g$begin + g$end) / 2;
		    g = g[order(g$pos),];
		    # rule 2 means use values from extrema
		    i = round(approx(g$pos, 1:nrow(g), xout = s$pos, rule=2)$y);
		    i = pmax(1, pmin(nrow(g), i));
		    indexSplit[[sc]] = g$index[i];
		}
	}
	unsplit(indexSplit, strains$scaffold);
}

getenv_numeric_or_default = function(envname, default) {
	value = Sys.getenv(envname);
	if (value=="") return(default);
	value = as.numeric(value);
	if (is.na(value)) return(default);
	return(value);
}

FEBA_Fit = function(expsUsed, all, genes,
	   		       genesUsed=NULL, strainsUsed=NULL, genesUsed12=NULL,
	   		       minT0Strain=3, minT0Gene=30,
			       minT0GeneSide=minT0Gene/2,
			       minGenesPerScaffold=10,
			       pred=CrudeOp(genes),
			       okDay=TRUE, # OK to use Time0 from another day on the same lane, if necessary?
			       okLane=TRUE, # OK to compare to Time0 from another lane, if necessary?
	   		       metacol=1:7, # for all
			       # names of experiments to ignore; experiments with Drop=TRUE are also ignored
			       ignore=NULL,
			       # ignore those below this threshold, unless ignore is set
			       minSampleReads = getenv_numeric_or_default("FEBA_MIN_SAMPLE_READS", 200*1000),
			       debug=FALSE, computeCofit=TRUE,
                               dir=".",
			       ...) {


	if (is.null(ignore)) {
	    tot = colSums(all[,-metacol]);
	    ignore = names(all)[-metacol][tot < minSampleReads];
        }
	if (!is.null(expsUsed$Drop) && any(expsUsed$Drop, na.rm=TRUE)) {
	    ignore = unique(c(ignore, expsUsed$name[!is.na(expsUsed$Drop) & expsUsed$Drop]));
	}
	if(length(ignore) > 0) {
	    cat("Ignoring ",ignore,"\n");
	    expsUsed = expsUsed[!expsUsed$name %in% ignore,];
	    all = all[, !names(all) %in% ignore,];
	}
        if (nrow(expsUsed)==0) stop("No experiments left to analyze!");

	if(!all(expsUsed$name %in% names(all)))
		stop("names missing from all");
	if(is.null(genes$scaffoldId)) stop("No scaffold for genes");
	if(is.null(genes$begin)) stop("No begin for genes");
	if (is.null(genes$GC)) stop("Warning: no GC field in genes");

	expsUsed$name = as.character(expsUsed$name);

	# if Group = Time0, it is a Time0, even if "short" has a different description
	if(!is.null(expsUsed$Group)) expsUsed$short[expsUsed$Group == "Time0"] = "Time0";
   
        write("Aggregating all over genes",stderr());
	has_gene2 = !is.na(all$f) & all$f >= 0.1 & all$f <= 0.9; # has gene and is central
	all_gN = aggregate(all[has_gene2,-metacol], list(locusId=all$locusId[has_gene2]), sum);

	expsUsed$t0set = paste(expsUsed$Date_pool_expt_started, expsUsed$Set);
	d = expsUsed[expsUsed$short=="Time0",];
	expsT0 = split(d$name, d$t0set);

	# If there is no matched lane/date t0, use the same date from other lane(s)
	for (t0set in setdiff(expsUsed$t0set, expsUsed$t0set[expsUsed$short=="Time0"])) {
	    u = which(expsUsed$t0set == t0set); # affected rows

	    date = unique(expsUsed$Date_pool_expt_started[u]);
	    set = unique(expsUsed$SetName[u]);

	    if (okLane && any(expsUsed$Date_pool_expt_started == date & expsUsed$short == "Time0")) {
		expsT0[[t0set]] = NULL; # no longer used
	        cat("Using Time0 from other lanes instead for ", t0set,"\n");
	        cat("Experiments affected: ", expsUsed$name[u],"\n");
	        expsUsed$t0set[u] = date;
		expsT0[[date]] = expsUsed$name[expsUsed$Date_pool_expt_started == date & expsUsed$short == "Time0"];
	    } else if (okDay && any(expsUsed$SetName == set & expsUsed$short == "Time0")) {
		expsT0[[t0set]] = NULL; # no loner used
		newt0set = unique(expsUsed$t0set[expsUsed$short == "Time0" & expsUsed$SetName == set]);
		newt0set = newt0set[1]; # pick the first t0set arbitrarily
		cat("Warning! Using Time0 from other days instead for ", t0set, "\n");
		cat("Experiments affected: ", expsUsed$name[u], "\n");
		expsUsed$t0set[u] = newt0set;
	    } else {
		stop("No Time0 for", t0set);
	    }
	}

	# note that t0tot does not contain any metadata columns, but t0_gN does
	t0tot = data.frame(lapply(expsT0, function(x) rowSums(all[,x,drop=F])), check.names=F);
        write("Aggregating Time0 totals",stderr());

        # the next few lines are much faster than doing
	# t0_gN = data.frame(locusId = names(t0_gN_list[[1]]), aggregate(t0tot[has_gene2, drop=F], list(locusId=all$locusId[has_gene2]), sum);
        indexBy = as.factor(all$locusId[has_gene2]);
        t0_gN_list = mclapply(names(t0tot), function(n) tapply(t0tot[has_gene2,n], indexBy, sum));
        names(t0_gN_list) = names(t0tot);
        t0_gN = data.frame(locusId=names(t0_gN_list[[1]]), t0_gN_list, check.names=F);
        writeDelim(t0_gN, paste(dir,"/t0_gN",sep="")); # for debugging

	cat("Central Reads per t0set, in millions:\n");
	print(colSums(t0_gN[,-1,drop=F])/1e6, digits=2);

	if(is.null(strainsUsed)) {
	    strainsUsed = has_gene2 & rowMeans(t0tot) >= minT0Strain;
	} else {
	    strainsUsed = strainsUsed & has_gene2;
	}
	if (length(unique(all$locusId[strainsUsed])) < 10) stop("strainsUsed is nearly empty");
	if(is.null(genesUsed)) {
		t0_gN_used = aggregate(t0tot[strainsUsed,], list(locusId=all$locusId[strainsUsed]), sum);
                n0 = rowMeans(t0_gN_used[,-1,drop=F]);
                cat(sprintf("Time0 reads per gene: mean %.1f median %.1f ratio %.2f\n",
			mean(n0), median(n0), mean(n0)/median(n0)));
		genesUsed = t0_gN_used$locusId[ n0 >= minT0Gene ];
		cat("Genes with enough Time0 reads: ",length(genesUsed),"\n");
	}
	genesPerScaffold = table(genes$scaffoldId[genes$locusId %in% genesUsed]);
	smallScaffold = names(genesPerScaffold)[genesPerScaffold < minGenesPerScaffold];
	genesUsed = genesUsed[!genesUsed %in% genes$locusId[genes$scaffoldId %in% smallScaffold]];
	if (length(smallScaffold) > 0) cat("Ignoring genes on small scaffolds ",smallScaffold,"\ngenes left: ",length(genesUsed),"\n");
	if(length(genesUsed) < 100 || !all(genesUsed %in% genes$locusId)) stop("Less than 100 genes left!");

	if(length(strainsUsed) != nrow(all)) stop("Invalid strainsUsed");
	cat("Using ",sum(strainsUsed)," of ",sum(has_gene2)," genic strains\n");

	cat("Using ",length(genesUsed)," of ",length(unique(all$locusId[has_gene2]))," genes with data\n");

	if (is.null(genesUsed12)) {
  	    d1 = aggregate(t0tot[strainsUsed & all$f < 0.5,],
	     		   list(locusId=all$locusId[strainsUsed & all$f < 0.5]), sum);
	    d2 = aggregate(t0tot[strainsUsed & all$f >= 0.5,],
	     		   list(locusId=all$locusId[strainsUsed & all$f >= 0.5]), sum);
	    genesUsed12 = intersect(d1$locusId[ MyRowMin(d1[,-1,drop=F]) >= minT0GeneSide],
		      		    d2$locusId[ MyRowMin(d2[,-1,drop=F]) >= minT0GeneSide]);
	    # Should the counts for each half of the gene (d1,d2) be saved as a diagnostic?
            # t0_gN should be enough for now
	    if (length(genesUsed12) < 100) stop(sprintf("genesUsed12 has just %d entries -- check %s/t0_gN",
            			                        length(genesUsed12), dir));
	}
        cat("For cor12, using ",length(genesUsed12),"genes\n");

	if(!all(expsUsed$t0set %in% names(t0tot))) stop("Illegal t0set ", setdiff(expsUsed$t0set, names(t0tot)));

	results = mclapply(names(all)[-metacol], function(n) {
		to_subtract = expsUsed$short[expsUsed$name==n] == "Time0";
		if(debug) cat("GeneFitness() on", n,"t0set",expsUsed$t0set[expsUsed$name==n],
			  			  if(to_subtract) "subtracted" else "", "\n");
                x = all[[n]];
		t0set = as.character(expsUsed$t0set[expsUsed$name==n]);
		t0 = t0tot[[ t0set ]];
		if(to_subtract) t0 = t0 - x;
		if(any(t0 < 0)) stop("Illegal counts under 0 for ",n);
		if(all(t0 == 0)) {
		    cat("Skipping log ratios for ",n," which has no control counts\n");
		    return(NULL);
		}
		gene_fit = GeneFitness(genes, all[has_gene2,words("locusId f")], x[has_gene2], t0[has_gene2],
				    strainsUsed[has_gene2], genesUsed, genesUsed12, 
				    minGenesPerScaffold=minGenesPerScaffold);
                cntrl = setdiff(expsT0[[ t0set ]], n);
		if(length(cntrl) < 1) stop("No Time0 experiments for ",n," should not be reachable");
		if(debug) cat("StrainFitness() on", n, "versus", cntrl,"\n");
		straindata = StrainFitness(all[[n]], rowSums(all[,cntrl,drop=F]));
                gTest = gene_fit$locusId[abs(gene_fit$fitnorm) > 1];
                u = which(strainsUsed & all$locusId %in% gTest);
                polar = PolarTest(gTest, all[u, c("locusId","strand")], straindata$fit[u], debug=debug);
                if(debug) cat("Polar",n,"succeeded","nrows",nrow(polar),"\n");
                if(!is.null(polar) && nrow(polar) > 0) polar$expName = n;
		return(list(gene_fit=gene_fit,
                            strain_fit=straindata$fit, strain_se=straindata$se,
                            polar=polar));
	});
	names(results) = names(all)[-metacol];
	results = results[!sapply(results,is.null)];
	if(length(results) == 0) stop("All comparisons failed\n");
	if(debug) cat("GeneFitness() succeeded\n");
	fit = list(g = results[[1]]$gene_fit$locusId);
	for(n in setdiff(names(results[[1]]$gene_fit), "locusId"))
	    fit[[n]] = data.frame(lapply(results, function(x) x$gene_fit[[n]]));
	names(fit) = sub("fitnorm","lrn",names(fit));
	names(fit) = sub("fit","lr",names(fit));
	if (debug) cat("Extracted fitness values\n");
        fit$polar = do.call(rbind, lapply(results, function(x) x$polar));
        if(!is.null(fit$polar) && nrow(fit$polar) > 0) {
          if(debug) cat("Initial rows in combined polar:",nrow(fit$polar),"\n");
          fit$polar = merge(fit$polar, genes[,c("locusId","sysName","desc")]);
          fit$polar = merge(fit$polar, data.frame(expName=expsUsed$name, expDesc=expsUsed$Description));
          d = which(abs(fit$lrn) >= 1, arr.ind=T);
          d = data.frame(fit=fit$lrn[d], expName=names(fit$lrn)[d[,2]], locusId=fit$g[d[,1]]);
          fit$polar = merge(fit$polar, d, all.x=T);
          if(debug) cat("Final rows in combined polar:",nrow(fit$polar),"\n");
        } else {
          fit$polar = data.frame();
        }
        # Version 1.2: added fit$polar
	fit$version = "1.2.1";

	q_col = words("name short t0set");
	if(!is.null(expsUsed$num)) q_col = c(q_col, "num");
        fit$q = expsUsed[expsUsed$name %in% names(fit$lrn), q_col];
        qnames = as.character(fit$q$name);
	if(!all(qnames == names(fit$lrn))) stop("Mismatched names in fit");
        if(debug) cat("Running FitReadMetrics() and FitQuality()\n");
        fit$q = cbind(fit$q,
			FitReadMetrics(all, qnames, has_gene2),
                        FitQuality(fit, genes, pred));
	if(debug) cat("Running FEBA_Exp_status\n")
	status = FEBA_Exp_Status(fit$q, ...);
	fit$q$u = (status == "OK");
	fit$q$u[is.na(fit$q$u)] = FALSE;

	print(table(status));
	for(s in c("low_count","high_mad12","low_cor12","high_adj_gc_cor")) {
	    if(sum(status==s) > 0) cat(s, ":", fit$q$name[status==s],"\n");
	}

	fit$genesUsed = genesUsed;
	fit$strainsUsed = strainsUsed;
	fit$genesUsed12 = genesUsed12;
	# these gene totals are based on all strains, not on used strains, and will not match tot or tot0
	fit$gN = all_gN;
	fit$t0_gN = t0_gN;

	# These include all strains, not just those in genes
	fit$strains = cbind(all[,metacol], used=fit$strainsUsed, enoughT0=rowMeans(t0tot) >= minT0Strain);
	fit$strain_lr = data.frame(lapply(results, with, strain_fit));
	fit$strain_se = data.frame(lapply(results, with, strain_se));

	# Normalized per-strain values, based on the closest gene
	strainToGene = StrainClosestGenes(fit$strains, genes[match(fit$g, genes$locusId),]);
	fit$strain_lrn = mapply(function(sfit, gdiff) {
		# Add the relevant gene normalization; or, if NA, normalize the scaffold to a median of 0
		sdiffGene = gdiff[strainToGene];
		sdiffSc = -ave(sfit, fit$strains$scaffold, FUN=median);
		sdiff = ifelse(is.na(sdiffGene), sdiffSc, sdiffGene);
		return(sfit + sdiff);
	}, fit$strain_lr, fit$lrn-fit$lr);
	fit$strain_lrn = data.frame(fit$strain_lrn);

	# Statistics of cofitness on pairs, top cofitness hits, specific phenotypes
	if (computeCofit && sum(fit$q$u) >= 20) {
		cat("Computing cofitness with ", sum(fit$q$u), " experiments\n", file=stderr());
                adj = AdjacentPairs(genes);
                adjDiff = adj[adj$strand1 != adj$strand2,];
		adjDiff$rfit = cor12(adjDiff, fit$g, fit$lrn[,fit$q$u]);
		pred$rfit = cor12(pred, fit$g, fit$lrn[,fit$q$u]);
		fit$pairs = list(adjDiff=adjDiff, pred=pred);
		random = data.frame(Gene1 = sample(fit$g, length(fit$g)*2, replace=T),
		       	            Gene2 = sample(fit$g, length(fit$g)*2, replace=T));
		random = random[as.character(random$Gene1) != as.character(random$Gene2),];
		random$rfit = cor12(random, fit$g, fit$lrn[,fit$q$u]);
		fit$pairs = list(adjDiff=adjDiff, pred=pred, random=random);
		fit$cofit = TopCofit(fit$g, fit$lrn[,fit$q$u]);
		d = merge(fit$q[fit$q$u,], expsUsed, by=words("name short"));
		fit$specphe = SpecificPhenotypes(fit$g, d, fit$lrn[,fit$q$u], fit$t[,fit$q$u]);
	} else {
		cat("Only", sum(fit$q$u),"experiments of", nrow(fit$q)," passed quality filters!\n", file=stderr());
	}
        fit$high = HighFit(fit, genes, expsUsed);
	return(fit);
}

FEBA_Save_Tables = function(fit, genes, org="?",
		 topdir="data/FEBA/html/",
		 dir = paste(topdir,org,sep="/"),
		 writeImage=TRUE,
		 FEBAdir="src/feba",
		 template_file=paste(FEBAdir,"/lib/FEBA_template.html",sep=""),
		 expsU=expsUsed,
		 ... # for FEBA_Quality_Plot
		 ) {
	if(!file.exists(dir)) dir.create(dir);

	for (n in words("q lr lrn lrn1 lrn2 t")) {
	    if (is.null(fit[[n]]) || !is.data.frame(fit[[n]])) {
	        stop("Invalid or missing ",n," entry");
	    }
	}
	if (is.null(fit$genesUsed)) stop("Missing genesUsed");
	if (is.null(fit$g)) stop("Missing g -- versioning issue?");

	if(!all(names(fit$lr) == fit$q$name)) stop("Name mismatch");
	if(!all(names(fit$lrn) == fit$q$name)) stop("Name mismatch");

	nameToPath = function(filename) paste(dir,filename,sep="/");
	wroteName = function(x) cat("Wrote ",nameToPath(x),"\n",file=stderr());

	writeDelim(fit$q, nameToPath("fit_quality.tab"));
	wroteName("fit_quality.tab");

	writeDelim(cbind(genes, used=genes$locusId %in% fit$genesUsed), nameToPath("fit_genes.tab"));
	wroteName("fit_genes.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lr));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_unnormalized.tab"));
	wroteName("fit_logratios_unnormalized.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lrNaive));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_unnormalized_naive.tab"));
	wroteName("fit_logratios_unnormalized_naive.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lrn));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios.tab"));
	wroteName("fit_logratios.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lrn1));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_half1.tab"));
	wroteName("fit_logratios_half1.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lrn2));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_half2.tab"));
	wroteName("fit_logratios_half2.tab");

	d = genes[genes$locusId %in% fit$g, c("locusId","sysName","desc")];
	d$comb = paste(d$sysName, d$desc); # for MeV
	if (sum(fit$q$u) == 0) {
		cat("Warning: 0 OK experiments\n");
                d = d[order(d$locusId),]; # ensure same order as other tables
	} else {
		d = merge(d, cbind(locusId=fit$g,fit$lrn[,fit$q$u]));
		names(d)[-(1:4)] = paste(fit$q$name,fit$q$short)[fit$q$u];
	}
	writeDelim(d, nameToPath("fit_logratios_good.tab"));
	cat("Wrote fitness for ",sum(fit$q$u), " successful experiments to ", nameToPath("fit_logratios_good.tab"),"\n",
	    file=stderr());

	d = genes[genes$locusId %in% fit$g, c("locusId","sysName","desc")];
	d$comb = paste(d$sysName, d$desc); # for MeV
        d = merge(d, cbind(locusId=fit$g, fit$tot));
        names(d)[-(1:4)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("gene_counts.tab"));
        wroteName("gene_counts.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$t));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_t.tab"));
	wroteName("fit_t.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$se));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error_obs.tab"));
	wroteName("fit_standard_error_obs.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$sdNaive));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error_naive.tab"));
	wroteName("fit_standard_error_naive.tab");

	writeDelim(cbind(fit$strains,fit$strain_lrn)[order(fit$strains$scaffold, fit$strains$pos),],
		nameToPath("strain_fit.tab"));
	wroteName("strain_fit.tab");

	FEBA_Quality_Plot(fit$q, nameToPath("fit_quality.pdf"), org, ...);
	wroteName("fit_quality.pdf");

	if(is.null(fit$pairs)) {
		paste("No data for cofitness plot\n");
		unlink(nameToPath("cofitness.pdf"));
	} else {
		FEBA_Cofitness_Plot(fit$pairs, nameToPath("cofitness.pdf"), org);
		wroteName("cofitness.pdf");
	}

	pdf(nameToPath("fit_quality_cor12.pdf"),
		pointsize=10, width=6, height=6,
		title=paste(org,"Fitness Cor12 Plots"));
	for (i in 1:nrow(fit$q)) {
	    n = as.character(fit$q$name[i]);
	    changers = fit$g[abs(fit$t[[n]]) >= 3];
	    plot(fit$lrn1[[n]], fit$lrn2[[n]],
	    		  main=sprintf("%s %s #%d (gMed=%.0f rho12=%.3f)\n%s",
			  	org, n, fit$q$num[i], fit$q$gMed[i], fit$q$cor12[i], fit$q$short[i]),
	    		  xlab="First Half", ylab="Second Half",
			  col=ifelse(fit$g %in% changers, 2, 1));
	    eqline(); hline(0); vline(0);
	}
	dev.off();
	wroteName("fit_quality_cor12.pdf");

	labelAll = sprintf("%s #%d gMed=%.0f rho12=%.2f %30.30s",
		      sub("^set","",fit$q$name), fit$q$num, fit$q$gMed, fit$q$cor12, fit$q$short);
        labelAll = ifelse(fit$q$short=="Time0", paste(labelAll, fit$q$t0set), labelAll);

	use = fit$q$short != "Time0";
	if(sum(use) > 2) {
	    lrClust = hclust(as.dist(1-cor(fit$lrn[,as.character(fit$q$name)[use]], use="p")));
	    pdf(nameToPath("fit_cluster_logratios.pdf"),
		pointsize=8, width=0.25*pmax(8,sum(use)), height=8,
		title=paste(org,"Cluster Logratios"));
	    plot(lrClust, labels=labelAll[use], main="");
	    dev.off();
	    wroteName("fit_cluster_logratios.pdf");
	}

	if (ncol(fit$gN)-1 >= 3) { # at least 3 things to cluster
	    countClust = hclust(as.dist(1-cor(log2(1+fit$gN[fit$gN$locusId %in% fit$genesUsed,-1]))));
	    pdf(nameToPath("fit_cluster_logcounts.pdf"),
		pointsize=8, width=pmax(5,0.25*nrow(fit$q)), height=8,
		title=paste(org,"Cluster Log Counts"));
	    # Some Time0s may be missing from fit$q
	    d = match(names(fit$gN)[-1], fit$q$name);
	    labelAll2 = ifelse(is.na(d), paste("Time0", sub("^set","",names(fit$gN)[-1])), labelAll[d]);
	    plot(countClust, labels=labelAll2, main="");
	    dev.off();
	    wroteName("fit_cluster_logcounts.pdf");
	}

        d = table(genes$scaffoldId[genes$locusId %in% fit$genesUsed]);
	maxSc = names(d)[which.max(d)];
	if (is.null(maxSc)) stop("Invalid scaffoldId?");
	beg = ifelse(fit$g %in% genes$locusId[genes$scaffold==maxSc],
	    genes$begin[match(fit$g, genes$locusId)], NA);

	pdf(nameToPath("fit_chr_bias.pdf"), pointsize=10, width=6, height=6,
	          title=paste(org,"Chromosome Bias"));
	for (i in 1:nrow(fit$q)) {
	    n = as.character(fit$q$name[i]);
	    plot(beg, pmax(-2,pmin(2,fit$lr[[n]])),
	    		  main=sprintf("%s %s #%d (gMed=%.0f rho12=%.3f)\n%s",
			  	org, sub("^set","",n), fit$q$num[i], fit$q$gMed[i], fit$q$cor12[i], fit$q$short[i]),
	    		  xlab="Position on Main Scaffold",
			  ylab="Fitness (Unnormalized)",
			  ylim=c(-2,2), col="darkgrey");
	    o = order(beg);
	    lines(beg[o], (fit$lr[[n]] - fit$lrn[[n]])[o], col="darkgreen", lwd=2);
	    hline(0,lty=1,col=1);
	}
	dev.off();
	wroteName("fit_chr_bias.pdf");

	if (!is.null(expsU)) {
		writeDelim(expsU, nameToPath("expsUsed"));
		wroteName("expsUsed");
	}

	if (is.null(fit$cofit)) {
	    d = data.frame(locusId="",sysName="",desc="",hitId="",cofit=0,rank=0,hitSysName="",hitDesc="");
	} else {
	    d = merge(genes[,words("locusId sysName desc")], fit$cofit, by="locusId");
	    d = merge(d, data.frame(hitId=genes$locusId, hitSysName=genes$sysName, hitDesc=genes$desc));
	    d = d[order(d$locusId,d$rank),];
	}
	writeDelim(d, nameToPath("cofit"));
	wroteName("cofit");

	if (is.null(fit$specphe)) {
	   d = data.frame(locusId="",sysName="",desc="",name="",short="",Group="",Condition_1="",Concentration_1="",Units_1="",
	                  Condition_2="",Concentration_2="",Units_2="");
	   d = d[0,];
	} else {
	   d = merge(genes[,words("locusId sysName desc")], fit$specphe, by="locusId");
	}
	writeDelim(d, nameToPath("specific_phenotypes"));
	wroteName("specific_phenotypes");

        d = which(abs(fit$lrn) > 2 & abs(fit$t) > 5, arr.ind=T);
        if (nrow(d) >= 1) {
	  out = data.frame(locusId=fit$g[d[,1]], name=names(fit$lrn)[d[,2]], lrn=fit$lrn[d], t=fit$t[d]);
	  out = merge(genes[,words("locusId sysName desc")], merge(expsU[,c("name","short")], out));
	  writeDelim(out, nameToPath("strong.tab"));
	  wroteName("strong.tab");
	}

        writeDelim(fit$high, nameToPath("high_fitness.tab"));
        wroteName("high_fitness.tab");

        if(nrow(fit$polar) > 0) {
          writeDelim(fit$polar[,words("locusId sysName desc expName expDesc fit absDiff p")], nameToPath("polar_effects.tab"));
          wroteName("polar_effects.tab");
        }

	if(writeImage) {
	    img = format(Sys.time(),"fit%Y%b%d.image"); # e.g., fit2013Oct24.image
	    expsUsed = expsU;
	    save(fit, genes, expsUsed, file=nameToPath(img));
	    wroteName(img);
	    unlink(nameToPath("fit.image"));
	    file.symlink(img, nameToPath("fit.image"));
	    cat("Created link for ",nameToPath("fit.image"),"\n", file=stderr());
	}

	if(!is.null(template_file)) {
	    FEBA_Save_HTML(nameToPath("index.html"), template_file,
			   list(ORG=org,
			        NEXPS=sum(fit$q$short != "Time0"),
				NSUCCESS=sum(fit$q$u),
				VERSION=fit$version,
				DATE=date()));
	    wroteName("index.html");
	}
}

FEBA_Quality_Plot = function(q, pdfFile, org,
		             min_gMed=50, max_mad12=0.5, max_adjcor = 0.25, max_gccor = 0.2, min_cor12 = 0.1,
			     multiples=TRUE) {
	qCol = ifelse(q$short=="Time0", "grey",
		ifelse(q$u, ifelse(q$maxFit > 5, "blue", "darkgreen"), "red"));
	qLab = q$num;

	if(!is.null(pdfFile)) pdf(pdfFile,
		pointsize=10, width=8, height=8,
		title=paste(org,"Fitness Quality"));
	oldpar = par(mfrow=c(2,2));

	# mad12 vs. gMed
	plotlab(0.1 + q$gMed, pmin(0.75,q$mad12), qLab, col=qCol, ylim=c(0,0.75), log="x",
		            cex=0.8,
			    xlab="Median Reads per Gene",
			    ylab="Median abs. diff.(1st half, 2nd half)",
			    main=paste(org,"mad12"));
	vline(min_gMed); hline(max_mad12);

	# cor12 vs. mad12
	plotlab(pmin(0.75,q$mad12), pmax(0,q$cor12), qLab, col=qCol, xlim=c(0,0.75), ylim=0:1,
			     cex=0.8,
			     xlab="Median abs. diff.(1st half, 2nd half)", ylab="rho(1st half, 2nd half)",
			     main=paste(org,"rho12"));
	vline(max_mad12); hline(min_cor12);

	# opcor vs. adjcor
	plotlab(abs(q$adjcor), pmax(0,q$opcor), qLab, col=qCol, xlim=0:1, ylim=0:1,
			     cex=0.8,
			     xlab="| rho(adj. genes on diff. strands) |",
			     ylab="rho(operon pairs)",
			     main=paste(org,"rho(operons)"));
	eqline();
	vline(max_adjcor);

	# adjcor vs. gccor
	plotlab(abs(q$gccor), abs(q$adjcor), qLab, col=qCol, xlim=0:1, ylim=0:1,
			     cex=0.8,
			     xlab="| cor(gene GC, fitness) |",
			     ylab="rho(adj. genes on diff. strands)",
			     main=paste(org,"GC effects"));
	eqline();
	vline(max_gccor); hline(max_adjcor);

	if (!is.null(pdfFile) && multiples) {
	    for(t0set in unique(q$t0set)) {
	        FEBA_Quality_Plot(q[q$t0set==t0set,], pdfFile=NULL, org=t0set, multiples=FALSE);
	    }
	}
	par(oldpar);
	if(!is.null(pdfFile)) dev.off();
}

FEBA_Cofitness_Plot = function(pairs, pdfFile, org) {
	if(!is.null(pdfFile)) pdf(pdfFile,
		pointsize=10, width=4, height=4,
		title=paste(org,"Cofitness"));
	CompareDensities(list(Operon=withoutNA(pairs$pred$rfit[pairs$pred$bOp]),
			          Adjacent=withoutNA(pairs$adjDiff$rfit),
				  Random=withoutNA(pairs$random$rfit)),
			 legendX="topleft",
			 xlim=c(-1,1),
			 xlab="Cofitness", ylab="Density", lwd=c(1,2,2), col=c(3,2,1), lty=c(1,2,4),
			 main=paste("Cofitness in",org));
	if(!is.null(pdfFile)) dev.off();
}

FEBA_Save_HTML = function(outfile, template, values) {
	lines = readLines(template);
	for (n in names(values)) lines = gsub(n, values[[n]], lines);
	writeLines(lines, outfile);
}

# Utilities

# split a string separated by spaces into a list of word
words = function(s, by=" ", ...) { strsplit(s[1], by, ...)[[1]]; }

# Tab-delimited output
writeDelim = function(table,file,report=FALSE,...) {
	write.table(table,file,sep="\t",quote=FALSE,row.names=FALSE,...);
	if(report) cat("Wrote",nrow(table),"rows","to",file,"\n")
}

# Just like plot.default() but adds on the labels
# By default plotting symbols are off but you can override that
plotlab = function(x,y,labels, cex=1, col=1, pch="", ...) {
	plot.default(x,y,cex=cex,col=col,pch=pch,...);
	text(x,y,labels,cex=cex,col=col);
}

### Graphics utilities
hline <- function(y,col="grey",lty=2,lwd=1) {
	lines(c(-1e20,1e-40,1e20),c(y,y,y),col=col,lty=lty,lwd=lwd);
}

vline <- function(x,col="grey",lty=2,lwd=1) {
	lines(c(x,x,x),c(-1e20,1e-40,1e20),col=col,lty=lty,lwd=lwd);
}

eqline <- function(col="grey",lty=2,lwd=1) {
	x <- 10**(-25:25);
	lines(c(-rev(x),x),c(-rev(x),x),col=col,lty=lty,lwd=lwd);
}

# Crude operon predictions -- pairs of genes that are on the same strand and
# separated by less than the median amount are predicted to be in the same opron
# Input genes is a data frame with locusId, strand, begin, end, with genes in sorted order
# Returns a data frame with Gene1, Gene2, Sep for separation, and bOp (TRUE if predicted operon pair)
CrudeOp = function(genes) {
	d = merge(merge(data.frame(Gene1=genes$locusId[-nrow(genes)],Gene2=genes$locusId[-1]), genes, by.x="Gene1", by.y="locusId"), genes, by.x="Gene2", by.y="locusId",suffixes=1:2);
	d = d[d$strand1==d$strand2,]
	d$Sep = pmin(abs(d$begin1-d$end2),abs(d$end1-d$begin2));
	d$bOp = d$Sep < median(d$Sep);
	return(d);
}

# are values correlated for pairs? Intended for operon pairs but would work with other types too
paircor = function(pairs, genes, values, use="p", method="pearson", names=c("Gene1","Gene2")) {
	d = merge(merge(pairs[,names], data.frame(Gene1=genes, value1=values), by.x=names[1], by.y="Gene1"),
			data.frame(Gene2=genes, value2=values), by.x=names[2], by.y="Gene2");
	return(cor(d$value1, d$value2, use=use, method=method));
}

# median-based normalization
mednorm = function(x) x - median(x);

# replace missing values (NA) with 0s
na0 = function(x) ifelse(is.na(x),0,x);

MyRowMin = function(x) {
	 if (is.vector(x)) return(x);
	 apply(x, 1, min);
}

# For reformatting per-lane pool count tables into the "all" table
prefixName = function(x, prefix) { names(x) = paste(prefix,names(x),sep=""); return(x); }

# For shortening the experiment descriptions
applyRules = function(rules, desc) {
    for (i in 1:nrow(rules)) desc = sub(rules[i,1], rules[i, 2], desc);
    return(desc);
}

# to count number of A, C, and G nucleotides in each barcode
CountACG = function(rcbarcodes) {
	 nA = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "A"));
	 nC = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "C"));
	 nG = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "G"));
	 return(data.frame(nA=nA,nC=nC,nG=nG));
}

without = function(list,columns=list$without) {
	list2 = list;
	for (i in columns) { list2[[i]] = NULL; }
	return(list2);
}

AdjacentPairs = function(genes) {
	genes = genes[order(genes$scaffold, genes$begin),];
	adj = data.frame(Gene1 = genes$locusId, Gene2=c(genes$locusId[-1],genes$locusId[1]));
	# add metadata and only keep pairs with same scaffold
	adj = merge(merge(adj, genes, by.x="Gene1", by.y="locusId"), genes, by.x=c("Gene2","scaffoldId"), by.y=c("locusId","scaffoldId"), suffixes=1:2);
	return(adj);
}

# Given a list of Gene1 Gene2 pairs, and a matrix of data (as genes and data-only matrix),
# compute correlations for each pair or NA
cor12 = function(pairs, genes, data, use="p", method="pearson", names=c("Gene1","Gene2")) {
	i1 = match(pairs[[names[1]]], genes);
	i2 = match(pairs[[names[2]]], genes);
	return(sapply(1:nrow(pairs), function(x) if(is.na(i1[x]) | is.na(i2[x])) NA else
		cor(c(data[i1[x],], recursive=T), c(data[i2[x],], recursive=T), method=method, use=use)));
}

# Compute read metrics -- nMapped, nPastEnd, nGenic, for the given data columns
# The final argument is used to define genic
FitReadMetrics = function(all, cols, in_gene) {
	data.frame(nMapped  = colSums(all[, cols, drop=F]),
                   nPastEnd = colSums(all[all$scaffold=="pastEnd", cols, drop=F]),
                   nGenic = colSums(all[in_gene, cols, drop=F]));
}

# Compute the quality metrics from fitness values, fitness values of halves of genes, or
# counts per gene (for genes or for halves of genes)
FitQuality = function(fit, genes, pred=CrudeOp(genes)) {
	adj = AdjacentPairs(genes);
	adjDiff = adj[adj$strand1 != adj$strand2,];

	data.frame(
		nUsed = colSums(fit$tot),
		gMed = apply(fit$tot, 2, median),
		gMedt0 = apply(fit$tot0, 2, median),
		gMean = apply(fit$tot, 2, mean),
		cor12 = mapply(function(x,y) cor(x,y,method="s",use="p"), fit$lrn1, fit$lrn2),
		mad12 = apply(abs(fit$lrn1-fit$lrn2), 2, median, na.rm=T),
		# consistency of log2 counts for 1st and 2nd half, for sample and for time0
		mad12c = apply(abs(log2(1+fit$tot1) - log2(1+fit$tot2)), 2, median, na.rm=T),
		mad12c_t0 = apply(abs(log2(1+fit$tot1_0) - log2(1+fit$tot2_0)), 2, median, na.rm=T),
		opcor = apply(fit$lrn, 2, function(x) paircor(pred[pred$bOp,], fit$g, x, method="s")),
		adjcor = sapply(names(fit$lrn), function(x) paircor(adjDiff, fit$g, fit$lrn[[x]], method="s")),
		gccor = c( cor(fit$lrn, genes$GC[ match(fit$g, genes$locusId) ], use="p") ),
		maxFit = apply(fit$lrn,2,max,na.rm=T)
	);
}        
        
# Returns status of each experiment -- "OK" is a non-Time0 experiment that passes all quality metrics
# Note -- arguably min_cor12 should be based on linear correlation not Spearman.
# 0.1 threshold was chosen based on Marinobacter set5, in which defined media experiments with cor12 = 0.1-0.2
# clearly worked, and Kang Polymyxin B (set1), with cor12 ~= 0.13 and they barely worked.
FEBA_Exp_Status = function(q, min_gMed = 50, max_mad12 = 0.5, min_cor12 = 0.1,
				 max_gccor = 0.2, max_adjcor = 0.25) {
    with(q, ifelse(short=="Time0", "Time0",
	           ifelse(gMed < min_gMed, "low_count",
		   ifelse(mad12 > max_mad12, "high_mad12",
		   ifelse(cor12 < min_cor12, "low_cor12",
		   ifelse(abs(gccor) > max_gccor | abs(adjcor) > max_adjcor, "high_adj_gc_cor", "OK"))))));
}

withoutNA = function(x) x[!is.na(x)];

CompareDensities <- function(list,labels=names(list),xlim=range(unlist(list)),ylim=c(0,3),
		col=1:length(labels), lty=1:length(labels), lwd=rep(1,length(labels)),
		legendX=mean(xlim),legendY=ylim[2],main="",xlab="",ylab="",showCounts=FALSE,
		showLegend=TRUE, bty="o") {

	for (i in 1:length(labels)) {
		x = list[[ names(list)[i] ]];
		d <- density(x,from=xlim[1],to=xlim[2]);
		if(i==1) {
			plot(d$x,d$y,type="l",col=col[i],lty=lty[i],lwd=lwd[i],
				main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim);
		} else {
			lines(d$x,d$y,col=col[i],lty=lty[i],lwd=lwd[i]);
		}
		if (showCounts) labels[i] <- paste(labels[i]," (n=", sum(!is.na(x)), ")", sep="");
	}
	if(showLegend) {
	       if(is.numeric(legendX)) {
	           legend(legendX,legendY,labels,col=col,lty=lty,lwd=lwd, bty=bty);
	       } else { # e.g. "topleft"
	           legend(legendX, labels,col=col,lty=lty,lwd=lwd, bty=bty);
 	       }
	}
}

# Given pairs of adjacent genes, identify the upstream and downstream genes,
# as "up" and "dn"
OpPairUpDn = function(oppairs, genes) {
	oppairs$strand = genes$strand[match(oppairs$Gene1, genes$locusId)];
	if(any(is.na(oppairs$strand))) stop("Unknown genes in Gene1 in OpPairUpDn()");

	# mixing factors is problematic so make sure they are character or integer
	if(is.factor(oppairs$Gene1)) oppairs$Gene1 = as.character(oppairs$Gene1);
	if(is.factor(oppairs$Gene2)) oppairs$Gene2 = as.character(oppairs$Gene2);

	oppairs$up = ifelse(oppairs$strand=="+", oppairs$Gene1, oppairs$Gene2);
	oppairs$dn = ifelse(oppairs$strand=="-", oppairs$Gene1, oppairs$Gene2);
	return(oppairs);
}

# How often upstream-only is sick vs. downstream-only
# oppairs must include "up", "dn" (as in returned values from OpPairUpDn)
OperonPairSick = function(oppairs, g, lrn, t,
			sick = -1, min_diff = 0.75, max_t = -4) {
	oppairs = oppairs[oppairs$up %in% g & oppairs$dn %in% g,];
	if(is.vector(lrn)) {
		i1 = match(oppairs$up, g);
		i2 = match(oppairs$dn, g);
		uponly = sum(lrn[i1] < sick & lrn[i2] > sick & t[i1] < max_t & t[i2] > max_t
			& abs(lrn[i1]-lrn[i2]) > min_diff);
		dnonly = sum(lrn[i2] < sick & lrn[i1] > sick & t[i2] < max_t & t[i1] > max_t
			& abs(lrn[i1]-lrn[i2]) > min_diff);
		both = sum(lrn[i1] < sick & lrn[i2] < sick & t[i1] < sick & t[i2] < sick);
		return(c(uponly=uponly, dnonly=dnonly, both=both));
	}
	#else make a table, where the row names will be the experiment names
	as.data.frame(t(mapply(function(x,y)
	    OperonPairSick(oppairs, g, x, y, sick=sick, min_diff=min_diff, max_t=max_t), lrn, t)));
}

# Utilities for mapping items within regions
# e.g. for identifying strains that lie within genes
# These routines were originally used to map strains to genes, but that has been
# replaced by BarSeqR.pl
# Limitations: No items will map within regions that wrap around the origin (i.e., begin > end)

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

# given that data is subgrouped by the same markers (such as a scaffold),
# do findWithin on each and merge the results
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

# Identify "specific phenotypes" -- cases where a gene is sick
# in some experiment(s), with |fit| > minFit and |fit| > percentileFit + minDelta and |t| > minT
# percentileFit is defined as the 95th percentile (by default) of |fit| for that gene
#
# exps ideally includes name (the column names of lrn and tval) along with
# short, Group, Condition_1, Concentration_1, Units_1, Condition_2, Concentration_2, Units_2
#
# Returns a data frame of locusId, fit, t, name, short, etc.
SpecificPhenotypes = function(g, exps, lrn, tval,
	   	    minT = 5, minFit = 1.0,
		    percentile = 0.95, percentileFit = 1.0, minDelta = 0.5,
		    expsFields = intersect(names(exps),
		         words("name short Group Condition_1 Concentration_1 Units_1 Condition_2 Concentration_2 Units_2 Condition_3 Concentration_3 Units_3 Condition_4 Concentration_4 Units_4")))
{
	rowHi = apply(abs(lrn), 1, quantile, percentile);
	bool = abs(lrn) > minFit & abs(lrn) > rowHi+minDelta & rowHi < percentileFit & abs(tval) > minT;
	specsick = data.frame(which(bool, arr.in=T));
	specsick$locusId = g[specsick$row];
	specsick$name = names(lrn)[specsick$col];
	specsick$lrn = as.matrix(lrn)[cbind(specsick$row,specsick$col)];
	specsick$t = as.matrix(tval)[cbind(specsick$row,specsick$col)];
	specsick$row = NULL;
	specsick$col = NULL;
	return(merge(specsick, exps[,expsFields]));
}

# Save the which-strains-to-use components of a fit data structure to the given directory
SaveStrainUsage = function(fit, dir=".") {
	stopifnot(is.logical(fit$strains$used));
	stopifnot(is.character(fit$strains$barcode));
	stopifnot(!is.null(fit$genesUsed));
	stopifnot(!is.null(fit$genesUsed12));
	write(fit$strains$barcode[fit$strains$used], paste(dir,"/strainusage.barcodes",sep=""), ncol=1);
	write(fit$genesUsed, paste(dir,"/strainusage.genes",sep=""), ncol=1);
	write(fit$genesUsed12, paste(dir,"/strainusage.genes12",sep=""), ncol=1);
        cat(sprintf("Wrote strain usage to %s/strainusage.*\n", dir));
}

# Note thresholds are different than in high_fit.pl
HighFit = function(fit, genes, expsUsed, min.fit=4, min.t=5, max.se=2, min.gMean=10, max.below=8) {
  wHigh = which(fit$lrn >= min.fit & fit$t >= min.t, arr.ind=T);
  high = data.frame(locusId=fit$g[wHigh[,1]], expName=names(fit$lrn)[wHigh[,2]], fit=fit$lrn[wHigh], t=fit$t[wHigh]);
  # t ~= fit/standard_error, so estimate s.e. = fit/t
  high$se = high$fit/high$t;
  high$sdNaive = fit$sdNaive[wHigh];
  high = subset(high, se <= max.se);

  # which experiments are ok
  fields = words("name Group Condition_1 Concentration_1 Units_1 Media short");
  fields = fields[fields %in% names(expsUsed)];
  exps = expsUsed[, fields];
  exps = merge(exps, fit$q[,words("name u short maxFit gMean")]);
  high = merge(high, exps, by.x="expName", by.y="name");
  high = subset(high, gMean >= min.gMean & fit >= maxFit - max.below);
  names(high)[names(high)=="u"] = "used";
  high = merge(genes[,c("locusId","sysName","desc")], high);
  high = high[order(high$expName, -high$fit),];
  return(high);
}

# genes is a list of locusIds
# strains includes locusId and strand
# returns the difference between the two strands and a p value,
# if it meets the thresholds
PolarTest = function(genes, strains, strain_fit,
                     mindiff=1, maxp=0.01, debug=F) {
  if(length(genes)==0) return(NULL);
  stopifnot(nrow(strains) == length(strain_fit));
  if(debug) cat("PolarTest has",length(genes),"genes",nrow(strains),"strains","\n");
  groups = split(1:nrow(strains), strains$locusId);
  results = lapply(groups, function(i) {
    values = strain_fit[i];
    isPlus = strains$strand[i] == "+";
    if (sum(isPlus) < 2 || sum(!isPlus) < 2) return(NULL);
    diff = abs(mean(values[isPlus]) - mean(values[!isPlus]));
    if (diff < mindiff) return(NULL);
    data.frame(locusId = strains$locusId[i[1]], absDiff=diff,
               # t.test will occasionally fail because all values are equal;
               # return a p-value of 0 in this case (have already verified a large difference)
               p=tryCatch( { t.test(values[isPlus], values[!isPlus])$p.value },
                           error = function(e) { 0 } ));
  });
  results = do.call(rbind, results);
  if (is.null(results) || nrow(results) == 0) return(NULL);
  return(subset(results, p < maxp));
}
