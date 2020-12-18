# Interactive Utilities for viewing fitness data across organisms:

# get_genes(org, list of gene names)
#         given a list or space-delimited set of names, which may be
#	  locus tags (sysName) or VIMSS ids as well as locusIds, return the 
#	  corresponding locusIds, in the same order (or NA for items that do not match).

# gene_info(org, list of gene names)
#	a short report for each of the specified genes, including whether they have specific phenotypes,
#	whether those are conserved, and their top cofitness hits.
#	The gene names are as for get_genes()

# get_exps(org, pattern)
#	returns a list of sucessful experiment ids ("name", i.e. set1H5) that match the pattern

# exp_info(org, list of experiment names)
#	prints information about those experiments

# show_fit(org, list of gene names)
#	an interactive heatmap of the fitness data for those genes
#	The gene names are as for get_genes()
#	Also try show_fit_dot(org, gene) for one gene.

# compare_genes(org, gene1, gene2)
#	interactive scatterplot of fitness data for the two genes
#	The gene names are as for get_genes()

# volcano_fit(org, gene) -- interactive volcano plot (|t| vs. fitness)

# compare_exps(org, exp1, exp2, minT=4)
#	interactive scatterplot of fitness data for the two experiments
#	each name should be like "setH5", or a list of replicates or near-replicates
#	genes with |t1| < minT are in grey
#
# PlotSpecList -- make a specific phenotypes plot

# Programmer Utilities:
# Use LoadOrgs(list of organism nicknames) to create the data structures
#     Then use SetupEssentials() as well to load essentiality information.
#     Use CensorExperiments() to remove experiments from the data set.
# Use names(orgs) to list which FEBA nicknames are included
# Use get_genes(), get_exps(), get_fit(), get_t(), all_cofitness() to get at data items

# Data structures (created by LoadOrgs()):
# orgs -- a separate data structure for each bug
#	genes -- metadata about genes, including VIMSS ids when possible
#	exps -- metadata about experiments (including Time 0s)
#	g -- gene ids in the order that they appear in lrn, t
#	q -- quality metrics for fitness experiments, including success or not (u)
#	lrn -- fitness values, 1 row per g, 1 column per successful experiment
#	t -- similarly for t values (but includes all experiments)
#	specsick -- specific phenotypes [these can be positive, actually, so "sick" is a misnomer]
#	cofit -- most cofit hits for each gene
# If global is TRUE, then these are also computed:
# specsicks -- the union of the specific phenotypes across all organisms
# css -- specific phenotypes that are conserved for orthologs genes (conserved specific sick)
# ccofit -- cofitness that is conserved for orthologous pairs (conserved cofit)
# 	 only has pairs in the order locus1 < locus2
# para -- information about paralogs (and orgs[[org]]$para and $parahits are built from it)


LoadOrgs = function(orgnames,
	base="data/FEBA/",
        html=paste(base,"/html/",sep=""),
        global=TRUE,
	debug=F) {
	orgs = list();
	for(n in orgnames) {
	       cat("Loading ",n,"\n",sep="");
		exps = read.delim(paste(html,n,"/expsUsed",sep=""),as.is=T, quote="");
		genes = read.delim(paste(html,n,"/genes",sep=""), as.is=T, quote="");
		q = read.delim(paste(html,n,"/fit_quality.tab",sep=""), as.is=T, quote="");
		if(!all(q$name %in% exps$name)) stop(n, ": Unmatched names in q: ", setdiff(q$name,exps$name));
		lrn = read.delim(paste(html,n,"/fit_logratios.tab",sep=""),check.names=F,as.is=T,quote="");
		if(debug) cat("Read ", nrow(lrn)," rows from fit_logratios.tab\n");
		names(lrn) = sub(" .*","",names(lrn));
		if(is.null(q$u)) stop("Missing u column in experiment quality table is no longer supported");
		g = lrn$locusId;
		if(!all(q$name %in% names(lrn))) stop(n, ": in q but not in lrn: ", setdiff(q$name, names(lrn)));

		tval = read.delim(paste(html,n,"/fit_t.tab",sep=""),as.is=T,check.names=F,quote="");
		names(tval) = sub(" .*","",names(tval));
		if(!all(q$name %in% names(tval))) stop(n, ": in q but not in t: ", setdiff(q$name, names(tval)));

		# Save Time0s, then remove unsuccessful experiments and Time0s from the main tables
		lrn_t0 = lrn[,q$name[q$short=="Time0"]];
		t_t0 = tval[,names(lrn_t0)];
		lrn = lrn[,q$name[q$u]];
		tval = tval[,names(lrn)];

		file = paste(html,n,"/cofit",sep="");
		if (file.access(file,mode=4) != 0) {
			cat("Warning: cannot read cofitness file for ",n,"\n",sep="");
			cofit = NULL;
		} else {
			cofit = read.delim(file, as.is=T, quote="");
			if (nrow(cofit) > 0 && all(is.na(cofit$locusId))) cofit=data.frame();
		}

		file = paste(html,n,"/gffmo",sep="");
		if (file.access(file,mode=4) == 0) {# readable
		    gffmo = read.delim(file, as.is=T, quote="");
		    genes = merge(genes, gffmo[,words("locusId VIMSS")], all.x=T);
		}

		file = paste(base,"/g/",n,"/pfam.tab", sep="");
		if (file.access(file,mode=4) == 0) {
		    pfam = read.delim(file, as.is=T, quote="");
		} else {
		    cat("Warning: no pfam.tab file for ", n, "\n");
		    pfam = data.frame(locusId=genes$locusId[1], domainId="PF0000.00", domainName="undef")[0,];
		}

		file = paste(base,"/g/",n,"/tigrfam.tab", sep="");
		if (file.access(file,mode=4) == 0) {
		    tigrfam = read.delim(file, as.is=T, quote="");
		} else {
		    cat("Warning: no tigrfam.tab file for ", n, "\n");
		    tigrfam = data.frame(locusId=genes$locusId[1], domainId="TIGR00000", domainName="undef")[0,];
		}

		specsick = read.delim(paste(html,n,"/specific_phenotypes",sep=""), as.is=T, quote="");

		if (nrow(specsick) > 0 && !all(specsick$locusId %in% g)) stop("Unknown genes in specsick");
		if (nrow(cofit) > 0 && !all(cofit$locusId %in% g)) stop("Unknown genes in cofit");
		if (nrow(cofit) > 0 && !all(g %in% cofit$locusId)) stop("Not all genes with data are in cofit");
		orgs[[n]] = list(org=n, genes=genes, exps=exps, g=g, q=q, lrn=lrn, t=tval, lrn_t0=lrn_t0, t_t0=t_t0,
			         specsick=specsick, cofit=cofit,
				 pfam=pfam, tigrfam=tigrfam);
	}

        if (global) {
	  d = lapply(orgs, function(x) x$specsick);
	  for(i in names(d)) {
	    if(nrow(d[[i]]) > 0) {
	        d[[i]]$org = i;
	    } else {
	        d[[i]] = NULL; # remove it;
	    }
            # add Condition/Units/Concentration/_3,4 if not present already
            for (n in c("Condition_3","Units_3","Concentration_3","Condition_4","Units_4","Concentration_4")) {
              if (is.null(d[[i]][[n]])) {
                d[[i]][[n]] = "";
              }
            }
	  }
	  if (!is.null(d) && length(d) > 0) {
          if(debug) cat("Combining specsicks\n");
	   specsicks <<- do.call(rbind,d);
           if(debug) cat("specsicks:", nrow(specsicks), " rows\n");
	  } else {
	   cat("Warning: all specific_phenotypes files were empty\n");
	   specsicks <<- NULL;
 	  }

	  css <<- read.delim(paste(base,"css",sep=""), as.is=T); # conserved specific sick
	  ccofit <<- read.delim(paste(base,"comb_cofit.pairs",sep=""),as.is=T);
	  para <<- read.delim(paste(base,"aaseqs.para",sep=""), as.is=T, header=F, col.names=words("org locusId para bits ratio"));
	  para <<- split(para, para$org);
	  for(org in names(orgs)) {
	    if (!is.null(para[[org]])) {
	        orgs[[org]]$parahits = para[[org]];
	        orgs[[org]]$para = aggregate(para[[org]][,"ratio",drop=F], para[[org]][,"locusId",drop=F], max);
	    }
	  }
        }

	orgs <<- orgs;

	tigrfams = read.delim(paste(base,"/tigrinfo",sep=""), as.is=T, quote="");
	tigrfams$roleId = NULL;
	tigrfams$ignore = NULL;
	rolelink = read.delim(paste(base,"/TIGRFAMS_ROLE_LINK",sep=""), as.is=T, header=F, col.names=c("tigrId","roleId"));
	rolelink = rolelink[is.unique(rolelink$tigrId),]; # ignore any multiple mappings
	roles = read.delim(paste(base,"/tigrroles",sep=""), as.is=T, quote="");
	roles = merge(roles[roles$level=="main",], roles[roles$level=="sub1",], by="roleId", suffixes=c("",".sub"));
	rolelink = merge(rolelink, data.frame(roleId=roles$roleId, toprole=roles$description, subrole=roles$description.sub));
	# as of TIGRFAMs 15, a few roles are lost due to lack of metadata, odd.
	if(!all(is.unique(rolelink$tigrId))) stop("Non-unique role descriptions");
	tigrfams = merge(tigrfams, rolelink, by="tigrId", all.x=T);
	tigrfams <<- tigrfams;
	tigrroles <<- tigrfams[!is.na(tigrfams$toprole),];
}

SetupEssentials = function(base="data/FEBA", debug=TRUE, ...) {
    for (org in names(orgs)) {
	if(debug) cat("Setting up essentials for ", org, "\n");
	ess_genes_table = read.delim(paste(base,"/g/",org,"/essentiality.genes",sep=""), as.is=T);
	orgs[[org]]$esstable <<- Essentials(orgs[[org]]$genes, ess_genes_table, orgs[[org]]$g, debug=debug, ...);
	orgs[[org]]$ess <<- with(orgs[[org]]$esstable, locusId[ess]);
	length(orgs[[org]]$ess);
    }
}

# Given the genes information, including GC content (in genes), and
# statistics on the uniqueness of the gene from Essentiality.pl, and
# whether or not we estimated fitness values, returns a table with
# ess=TRUE for likely essential protein-coding loci.  Strictly
# speaking, these are loci that are probably very important for
# fitness, but we don't know if they are truly essential.
Essentials = function(genes, ess, g, max.fp=0.02, debug=TRUE) {
	if (is.null(genes$type)) genes$type = 1; # workaround for PS
	ess = merge(ess, genes[,words("type locusId GC nTA")]);
	# ignore non-proteins, genes with parts that might not map uniquely, and very short genes
	ess = ess[ess$type==1 & ess$dupScore == 0 & ess$ntLenNoOverlap >= 100,];

	# Choose a minimum gene length such that there is just
	# ~1% odds of no central insertions by bad luck
	candLen = seq(100,1000,25);
	d = ess[ess$ntLenNoOverlap >= 500,];
	medCentral = median(d$nPosCentral);
	medLen = median(d$ntLenNoOverlap);
	p0 = sapply(candLen, function(minLen) ppois(0, medCentral * minLen/medLen));
	if (all(p0 > max.fp)) stop("All probabilities too high in Essentials()");
	minLen = candLen[ min(which(p0 < max.fp)) ];
	if(debug) cat("Chose length ", minLen, "minimum fp rate", p0[candLen==minLen], "\n");

	# density of insertion locations in central 10-90% of gene, normalized
	ess$dens = ess$nPosCentral/ess$ntLenNoOverlap;
	ess$dens = ess$dens / median(ess$dens);

	# rate of reads, normalized for %GC
	ess$reads = ess$nReads / ess$ntLenNoOverlap;
	ess = ess[order(ess$GC),];
	ess$normreads = ess$reads / runmed(ess$reads, 201, endrule="constant");

	ess$ess = ess$ntLenNoOverlap >= minLen & ess$normreads < 0.2 & ess$dens < 0.2 & !ess$locusId %in% g;
	return(ess);
}

# turn list (or space-delimited) systematic names or VIMSS ids or locusIds into a list of locusIds
get_genes = function(org, specs) {
	genes = orgs[[org]]$genes;
	if(is.null(genes)) stop("Invalid org: ",org);
	if(is.character(specs) && length(specs)==1 && grepl(" ",specs)) {
		specs = words(specs);
	}
	# fetch locusId instead of reusing values to convert to right kind
	ids = genes$locusId[match(specs,genes$locusId)];
	if (!is.null(genes$sysName)) ids = ifelse(is.na(ids), genes$locusId[match(specs, genes$sysName)], ids);
	if (!is.null(genes$name)) ids = ifelse(is.na(ids), genes$locusId[match(specs, genes$name)], ids);
	if (!is.null(genes$VIMSS)) ids = ifelse(is.na(ids), genes$locusId[match(specs, genes$VIMSS)], ids);
	return(ids);
}

get_exps = function(org, pattern, fixed=FALSE, ignore.case=!fixed, perl=!fixed) {
	q = orgs[[org]]$q;
	if(is.null(q)) stop("Invalid organism: ",org);
	return(q$name[q$u & grepl(pattern, q$short, ignore.case=ignore.case, perl=perl, fixed=fixed)]);
}

gene_info = function(org, loci, n=5, minStrong=2, minT=4) {
	arg = loci;
	loci = get_genes(org, loci);
	loci = loci[!is.na(loci)];
	genes = orgs[[org]]$genes;
	info = genes[match(loci, genes$locusId),];
	hasStrong = apply(abs(orgs[[org]]$lrn), 1, max) > minStrong;
	hasSig = apply(abs(orgs[[org]]$t), 1, max) > minT;
	# Uses orgs[[org]]$ess if available as a list of locusIds
	# (LoadOrgs() does not set that up.)
	info$class = ifelse(info$locusId %in% orgs[[org]]$ess, "Essential",
			ifelse(!info$locusId %in% orgs[[org]]$g, "No data",
			ifelse(info$locusId %in% orgs[[org]]$g[hasStrong & hasSig], "Strong",
			ifelse(info$locusId %in% orgs[[org]]$g[hasSig], "Significant",
			"No phenotype"))));
	if(is.null(info) || nrow(info) == 0) stop("No such gene in ",org,": ",arg);
	out = info[,words("locusId sysName name desc class")];
	out$desc = sub("(NCBI ptt file)","", out$desc, fixed=T);
	if(!is.null(info$VIMSS)) out$VIMSS = info$VIMSS;
	row.names(out) = 1:nrow(out);
	print(out);
	cofit = orgs[[org]]$cofit;
	for(locusId in info$locusId) {
	    locusShow = locusId;
	    if (!is.null(info$VIMSS)) locusShow = paste(locusShow, "VIMSS", info$VIMSS[info$locusId %in% locusId]);
	    if(locusId %in% info$locusId[!info$class %in% c("Essential","No data")]) {
                out = orgs[[org]]$specsick;
		out = out[out$locusId %in% locusId,words("name short lrn t")];
		# need to go from specsicks$short to a Condition_1
		out = merge(out, orgs[[org]]$exps[,words("Condition_1 name")]);
                if (exists("css")) out$conserved = out$Condition_1 %in% css$cond[css$locusId %in% locusId & css$tax %in% org];
		cat("\nSpecific phenotypes for", locusShow,
			info$sysName[info$locusId %in% locusId],
			info$desc[info$locusId %in% locusId], "\n");
		if(nrow(out) >= 1) print(out,digits=2) else cat("None\n");

		if(!is.null(cofit)) {
		    out = cofit[cofit$locusId %in% locusId & cofit$rank <= n,];
		    out = merge(out, genes[,words("locusId name")], by.x="hitId", by.y="locusId");
		    names(out)[names(out)=="name"] = "hitName";
		    out = out[order(out$rank),];
		    row.names(out) = 1:nrow(out);
		    out2 = out;
		    out2 = out2[,words("cofit hitId hitSysName hitName hitDesc")];
		    if(exists("ccofit")) out2$conserved = out2$hitId %in% ccofit$locus2[ccofit$tax == org & ccofit$locus1 %in% locusId] |
					                 out2$hitId %in% ccofit$locus1[ccofit$tax == org & ccofit$locus2 %in% locusId];
		    out2$hitDesc = sub("(NCBI ptt file)","", out2$hitDesc, fixed=T);
		    if(!is.null(out$VIMSS)) out2$VIMSS = out$VIMSS;
		    cat("\nTop cofit hits for", locusShow,
			info$sysName[info$locusId %in% locusId],
			info$desc[info$locusId %in% locusId], "\n");
		    print(out2,digits=2);
		}
	    }

	    dom = with(orgs[[org]], rbind(pfam[pfam$locusId %in% locusId,], tigrfam[tigrfam$locusId %in% locusId,]));
	    if (nrow(dom) > 0) {
		    cat("\nPfam and/or TIGRfam hits for", locusShow,
			info$sysName[info$locusId %in% locusId],
			info$desc[info$locusId %in% locusId], "\n");
		    dom = dom[order(dom$begin),];
		    print(dom[,words("domainId domainName begin end score evalue")]);
	    }
	}
}

exp_info = function(org, exps) {
	q = orgs[[org]]$q;
	if(is.null(q)) stop("Invalid org: ",org);
	if(!all(exps %in% q$name)) stop("Invalid exps: ", setdiff(exps,q$name));
	print(q[q$name %in% exps,], digits=2);

}

get_fit = function(org, loci, t=FALSE) {
	loci = get_genes(org, loci);
	g = orgs[[org]]$g;
	q = orgs[[org]]$q;
	indexes = match(loci, g);

	if(t == TRUE) {
		return(t(orgs[[org]]$t[indexes,]));
	} else {
		return(t(orgs[[org]]$lrn[indexes,]));
	}
}

get_t = function(org, loci) get_fit(org, loci, t=TRUE);

identify_exps = function(org, x, y, xlim=NULL, save=FALSE, showY=TRUE, col="darkgreen", cex=1) {
	q = metadata_by_exp(org);
	if(is.null(q)) stop("Invalid org: ",org);
	if(length(x) != nrow(q)) stop("Incorrect lengths: x ",length(x), " metadata ",nrow(q));
	if(length(y) != nrow(q)) stop("Incorrect lengths: y ",length(y), " metadata ",nrow(q));
	cat("Click on points, or right click to exit\n");
	n = 0;
	iSaved = c();
	xUse = if (!is.null(xlim)) pmin(xlim[2], pmax(xlim[1], x)) else x;
	while(TRUE) {
	    i = identify(xUse, y, "", n=1);
	    if (length(i) != 1) break;
	    if(!i %in% iSaved) {
	        iSaved = c(iSaved,i);
		text(xUse[i], y[i], length(iSaved),adj=c(-0.5,0), col=col, cex=cex, xpd=T);
	        if (showY) {
	            cat(sprintf("%d %s: %s (x %.1f y %.1f)\n", length(iSaved), q$name[i], q$short[i], x[i], y[i]));
	        } else {
	            cat(sprintf("%d %s: %s (fit %.1f)\n", length(iSaved), q$name[i], q$short[i], x[i]));
	        }
	    }
	}
	if (save && length(iSaved) > 0) return(click=1:n, exp=q$name[iSaved], short=q$short[iSaved], x=x[iSaved], y=y[iSaved]);
	return(NULL);
}

# For a given organism, for each column in lrn or t, a table of metadata
metadata_by_exp = function(org) with(orgs[[org]], {
	d = merge(q[q$u,], exps, by=words("name short num"));
	d = d[match(names(lrn), d$name),];
	if (!identical(names(lrn), d$name)) stop("used experiments does not match lrn!");
	return(d);
});

# By default, color codes C sources in blue, N sources in dark green, and stresses in red (and otherwise in black)
compare_genes = function(org, locus1, locus2=NULL, xlab=NULL, ylab=NULL, eq=TRUE, locate=TRUE,
	highlight = NULL, # list of experiments to highlight (will use red vs. darkgrey)
	col = NULL,
	pch = 20,
	...) {
	if (is.null(locus2)) {
		if(length(locus1) != 2) stop("Invalid length");
		locus2 = locus1[2];
		locus1 = locus1[1];
	}
	if(is.null(xlab)) xlab=locus1;
	if(is.null(ylab)) ylab=locus2;
	locus1 = get_genes(org, locus1);
	locus2 = get_genes(org, locus2);
	if (is.na(locus1)) stop("No data for: ",locus1, " ", xlab);
	if (is.na(locus2)) stop("No data for: ",locus2, " ", ylab);
	mat = get_fit(org, c(locus1,locus2));
	if (is.null(col)) {
	    if (!is.null(highlight)) {
		col = ifelse(row.names(mat) %in% highlight, 2, 8);
		cat("In red: ", row.names(mat)[row.names(mat) %in% highlight], "\n");
	    } else {
		col = sapply(tolower(metadata_by_exp(org)[["Group"]]), switch,
	    	    "carbon source" = "blue",
	    	    "nitrogen source" = "darkgreen",
		    "stress" = "red",
		    "black");
		cat("C sources in blue, N sources in green, stresses in red\n");
	    }
	}
	plot(mat[,1], mat[,2], xlab=xlab, ylab=ylab, col=col, pch=pch, ...);
	if(eq) eqline();
	if(locate) identify_exps(org, mat[,1], mat[,2]);
}

compare_exps = function(org, exp1, exp2, xlab=exp1, ylab=exp2, eq=TRUE, locate=TRUE, minT=4, col=NULL, ...) {
	if (is.null(orgs[[org]]$lrn[,exp1])) stop("No experiment ",exp1," in ",org);
	if (is.null(orgs[[org]]$lrn[,exp2])) stop("No experiment ",exp2," in ",org);
	x = rowMeans(orgs[[org]]$lrn[,exp1,drop=F]);
	y = rowMeans(orgs[[org]]$lrn[,exp2,drop=F]);
	tx = rowMeans(orgs[[org]]$t[,exp1,drop=F]);
	if(is.null(col)) col=ifelse(abs(tx) >= minT,"black","darkgrey");
	plot(x, y, xlab=xlab, ylab=ylab, col=col, ...);
	if(eq) eqline();
	if(locate) identify_genes(org, x, y);
}

identify_genes = function(org, x, y, col=2, newcol=2, ...) {
	genes = orgs[[org]]$genes;
	g = orgs[[org]]$g;
	if(length(x) != length(g)) stop("Incorrect lengths");
	if(length(y) != length(g)) stop("Incorrect lengths");
	cat("Click on points, or right click to exit\n");
	n = 1;
	while(TRUE) {
	    i = identify(x, y, as.character(n), n=1, col=col, ...);
	    if(length(i) != 1) break;
	    row = genes[genes$locusId == g[i],];
	    if(nrow(row) != 1) stop("Illegal index ",i);
	    id = as.character(row$locusId);
	    if(!is.null(row$VIMSS)) id = paste(id, row$VIMSS);
	    cat(sprintf("%d: %.3f %.3f %s %s %s %s\n", n, x[i], y[i], id, row$sysName, row$name, row$desc));
	    if(!is.null(newcol)) points(x[i],y[i],col=newcol);
	    n = n+1;
	}
}

volcano_fit = function(org, locus, xlab="Fitness", ylab="|t|", locate=TRUE, ...) {
	locus = get_genes(org,locus);
	if (length(locus) != 1) stop("must input one locus");
	x = get_fit(org, locus);
	if(is.na(x[1])) stop("No data for: ",locus);
	tval = get_t(org, locus);
	plot(x, abs(tval), xlab=xlab, ylab=ylab, ...);
	if(locate) identify_exps(org, x, abs(tval));
}

eqline <- function(col="grey",lty=2,lwd=1) {
	x <- 10**(-25:25);
	lines(c(-rev(x),x),c(-rev(x),x),col=col,lty=lty,lwd=lwd);
}

andNoNA = function(x) ifelse(is.na(x),FALSE,x);

# around: show the gene neighborhood around the gene of interest
#	Gene order will be flipped if gene is on the - strand
#	Only an option if only one locus specified
# condspec: choose a subset of conditions to show, as a perl-style regular expression against q$short
# scale: by default, colors are for fitness = -2 to +2. To focus on stronger differences use scale=3.
# ylabmax: maximum number of characters for the condition labels
# cex: adjust the font size (e.g. cex=0.5 for tiny text)
# sort: sort by average fitness value, and show a maximum of maxsort experiments -- very handy for larger data sets
#
show_fit = function(org, loci, labels=NULL, locate=TRUE, around=0, condspec=NULL, scale=2, ylabmax=20, cex=1,
	            sort=FALSE, maxsort=50) {
	if (length(loci==1) && grepl(" ",loci)) loci = words(loci);
	if(is.null(labels)) labels = loci;
	loci = get_genes(org, loci);
	if (length(loci)==1 && around > 0) {
	    genes = orgs[[org]]$genes;
	    if(!is.null(genes$scaffold)) genes = genes[genes$scaffold == genes$scaffold[genes$locusId==loci], ];
	    strand = genes$strand[genes$locusId==loci];
	    if(!is.null(genes$begin)) genes = genes[order(genes$begin),];
	    i = match(loci, genes$locusId);
	    i1 = pmax(1, i - around);
	    i2 = pmin(nrow(genes), i + around);
	    indexes = if(strand=="+") i1:i2 else i2:i1;
	    loci = genes$locusId[i1:i2];
	    labels = genes$sysName[i1:i2];
	}
	mat = get_fit(org, loci); # columns as genes, rows as experiments
	tval = get_t(org, loci);

	q = orgs[[org]]$q;
	q = q[q$u,];
	u = rep(TRUE, nrow(q));
	if(!is.null(condspec)) {
	    u = grepl(condspec,q$short,perl=T);
	    if(sum(u) == 0) stop("No conditions matching: ",condspec);
	}
	if (sort && sum(u) > maxsort) {
	    thresh = sort(abs(rowMeans(mat,na.rm=T)), decreasing=T)[maxsort];
	    cat("Threshold |fit| to show: ", thresh, "\n");
	    u = u & abs(rowMeans(mat,na.rm=T)) > thresh;
	}
	q = q[u,];
	mat = as.matrix(mat[u,]);
	tval = as.matrix(tval[u,]);

	if(length(q$short) != nrow(mat)) stop("Wrong number of rows");
	if (!all(q$name == row.names(mat))) stop("Name mismatch");

	if (sort) {
		o = order(rowMeans(mat,na.rm=T), decreasing=T);
	} else {
		metadata = metadata_by_exp(org);
		if (!identical(metadata$name, q$name)) stop("Metadata mismatch");
		o = order(metadata$Group, metadata$short);
	}
	labRows = q$short[o];
	labRows = sub("[0-9.]+ mM$","",labRows);
	labRows = sub("[0-9.]+ mg/ml$","",labRows);
	labRows = sub("[0-9.]+ vol%$","",labRows);
	# note labRows are shown going up so we want the 1st row to be replaced if a duplicate
	labRows = ifelse(andNoNA(labRows == c(labRows[-1],NA)), ".", labRows);

	# mar is bottom,left,top,right
	ylablen = pmin(ylabmax, max(nchar(labRows)));
	oldpar = par(mgp=c(2,1,0), mar=c(cex*(1+max(nchar(labels))),0.5,0.5,cex*ylablen));
	image(t(mat[o,]), col=myHeatColors(), breaks=breaksUse(scale=scale), useRaster=TRUE,
		xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
	# las = 2 (perpendicular to axis)
	mtext(labRows, las=2, side=4, at=seq(0,1,length.out=length(labRows)), cex=cex); # at right
	mtext(labels, las=2, side=1, at=seq(0,1,length.out=length(labels)), cex=cex); # at bottom

	genes = orgs[[org]]$genes;
	if (locate) {
		cat("Click on cells, or right click to exit\n");
		while(TRUE) {
		    at = locator(1);
		    if(is.null(at)) break;
		    x = at$x[1];
		    y = at$y[1];
		    iGene = pmax(1, pmin(length(labels), 1 + round(x * (length(labels)-1))));
		    g = loci[iGene];
		    geneinfo = genes[genes$locusId == g,];
		    if (nrow(geneinfo) != 1) stop("no metadata for gene: ",g);
		    geneinfo$desc = sub("(NCBI ptt file)","", geneinfo$desc, fixed=T);

		    # iOExp is the sorted order shown, iExp is the order in q
		    iOExp = pmax(1, pmin(length(labRows), 1 + round(y * (length(labRows)-1))));
		    iExp = o[iOExp];
		    locusString = as.character(g);
		    if (!is.null(geneinfo$VIMSS)) locusString = sprintf("%s VIMSS %d",locusString,geneinfo$VIMSS);
		    cat(sprintf("Gene: %s %s %s %s\nExperiment: %s %s\nFitness: %.3f (t %.3f)\n\n",
				locusString, geneinfo$sysName, geneinfo$name, geneinfo$desc,
				q$name[iExp], q$short[iExp],
				mat[iExp,iGene], tval[iExp,iGene]));
		}
	}
	par(oldpar);
}

show_fit_dot = function(org, locus, jitterBy=0.25, xlim=c(-4,4),
	     		xlab=paste("Fitness of",locus), main="",
			pch=20, col=1, locate=T, ...) {
	exps = metadata_by_exp(org);
	nGroups = length(unique(exps$Group));
	exps$fit = get_fit(org, locus);
	exps$t = get_t(org, locus);
	exps$iGroup = as.integer(as.factor(exps$Group));
	# The y layout is to jitter by +/- jitterBy, centered at iGroup
	y = jitterBy * (runif(nrow(exps)) - 0.5) * 2 + exps$iGroup;

	plot(xlim, c(1-0.4, nGroups+0.4), bty="n", xlab=xlab, ylab="", main=main, yaxt="n", pch=NA, yaxs="i");
	points(pmin(xlim[2], pmax(xlim[1], exps$fit)), y, pch=pch, col=col);
	d = unique(exps[,c("Group","iGroup")]);
	text(xlim[1], d$iGroup + jitterBy + strheight("A")/2, as.character(d$Group), adj=c(0,0), xpd=T);
	abline(v=0,col="darkgrey");
	if (locate) {
		cat("Click on points, or right click to exit\n");
		identify_exps(org, exps$fit, y, xlim=xlim, showY=FALSE, ...);
	}
}

breaksUse = function(scale=2) 
{
    d = seq(-3, 3, 0.25) * scale/3;
    d[1] = -100;
    d[length(d)] = 100;
    return(d);
}

# blue to yellow
myHeatColors = function () {
    dn = seq(1, 1/12, -1/12);
    up = seq(1/12, 1, 1/12);
    c(rgb(dn*0.2, dn*0.2, dn), rgb(up, up, 0));
}

# all_cofitness(org, genes) -- compute cofitness of average profile of genes with all
# 		     other genes that have fitness.
# Returns a table sorted by correlation (descending order).
# 		     
all_cofitness = function(org, genes) {
	vec = unlist(rowMeans(get_fit(org, genes)));
	out = data.frame(locusId = orgs[[org]]$g, r = c(cor(vec, t(orgs[[org]]$lrn)), recursive=T));
	out = merge(orgs[[org]]$genes, out);
	out = out[order(-out$r),];
	return(out);
}

# split a string separated by spaces into a list of word
words = function(s, by=" ", ...) { strsplit(s[1], by, ...)[[1]]; }

# returns TRUE or FALSE indicating whether the value is repeated
# input should be a list of values
is.unique = function(x) {
	counts = as.data.frame.table(table(x), dnn="x");
	return(x %in% counts$x[counts$Freq==1]);
}

# Remove some experiments from the data set
CensorExperiments = function(org, expRemove) {
  oldexps = orgs[[org]]$exps;
  oldq = orgs[[org]]$q;
  stopifnot(all(expRemove %in% oldexps$name));
  newexp = setdiff(metadata_by_exp(org)$name, expRemove);

  out = orgs[[org]];
  out$lrn = out$lrn[, newexp];
  out$t = out$t[, newexp];

  out$exps = subset(out$exps, !name %in% expRemove);
  out$q = subset(out$q, !name %in% expRemove);

  out$specsick = subset(out$specsick, !name %in% expRemove);

  # keep the column names but not any actual rows
  out$cofit =   out$cofit[0,];

  orgs[[org]] <<- out;
  cat("#Successful experiments for", org, "reduced from", sum(oldq$u), "to", sum(out$q$u), "\n");
}

# like the old PlotCSS
# but with sysName on left, no organism info, and user-specified label on right
# Each subplot is separated by horizontal lines (unless lines=FALSE)
# Sections can be indicated by loci=""
# Expressions are allowed in labels for loci but not labels for sections 
# Each row of condspec specifies which conditions are colored that way and
# must include Group, Condition_1, and col.
#   Optionally, Concentration_1 or other fields from the metadata_by_exp table can be specified in condspec
#   These can set for some rows and missing (empty or NA) for others.
#   These are all matched in a case insensitive way.
# Rows of condspec can also set cex or pch
#
PlotSpecList = function(org, loci, labels=rep("",length(loci)), condspec,
			default.col=1, default.cex=1, default.pch=20,
                        label.cex=1, jitterBy=0.25,
			xlim=c(-4,4),
			xlab="Gene Fitness", main="", lines=TRUE, stripes=FALSE, stripeColor="lightyellow",
                        shortlabels=FALSE, blanklabels=F) {
	# The y layout is to jitter by +/- jitterBy, centered at nrow(locispec):1
        # And the separator lines are at +/- 0.5
	plot(xlim, c(0.5,length(loci)+0.4), bty="n", xlab=xlab, ylab="", yaxt="n", pch=NA, yaxs="i", main=main);
        if(lines) segments(xlim[1], (2:length(loci))-0.5, xlim[2], (2:length(loci))-0.5);
        stopifnot(c("Group","Condition_1","col") %in% names(condspec));

	# Set up per-experiment colors
        meta = metadata_by_exp(org);
	is_simple = IsSimpleCond(meta$Group, meta$Condition_1, meta$Condition_2);
	col=rep(default.col, nrow(meta));
        assigned = rep(FALSE, nrow(meta));
        cex=rep(default.cex, nrow(meta));
        pch=rep(default.pch, nrow(meta));
	for (j in 1:nrow(condspec)) {
            u = !assigned;
            for (n in intersect(names(condspec), names(meta))) {
              if (condspec[[n]][j] != "" && !is.na(condspec[[n]][j]))
                u = u & tolower(meta[[n]]) %in% tolower(condspec[[n]][j]);
            }
            if (sum(u) == 0) cat("Warning: no matches for row", j, "of condspec\n");
            assigned[u] = TRUE;
            col[u] = if(is.character(default.col)) as.character(condspec$col[j]) else condspec$col[j];
	    if(!is.null(condspec$pch)) pch[u] = condspec$pch[j];
	    if(!is.null(condspec$cex)) cex[u] = condspec$cex[j];
	}
	for (i in 1:length(loci)) {
	    if(as.character(loci[i]) == "") next; # to allow deliberate spacing
            if (stripes && (i %% 2) == 1)
              rect(-100, 1+length(loci)-i-0.5, 100, 1+length(loci)-i+0.5, col=stripeColor, border=NA, xpd=T);
	    locusId = get_genes(org, loci[i]);
	    if(length(locusId) != 1) stop("Multiple genes from ", loci[i]);
	    if(is.na(locusId)) stop("Unrecognized gene in org ", org, " : ", loci[i]);
	    x = get_fit(org,locusId);
	    if(any(is.na(x))) stop("No fitness data for ", locusId," ",loci[i]," in org ",org);

	    y = jitterBy * (runif(length(x)) - 0.5) * 2;
	    points(pmax(xlim[1], pmin(xlim[2], x)), y + length(loci) - i + 1,
		cex=cex, pch=pch, col=col);
            # emphasize assigned points by replotting on top?
	    points(pmax(xlim[1], pmin(xlim[2], x)), y + length(loci) - i + 1,
		cex=cex, pch=pch, col=ifelse(assigned, col, NA));
	}
        if (!blanklabels)
	  mtext(if(shortlabels) sub("^.*_", "_", loci, perl=T) else loci,
              at=length(loci):1, side=2, line=0, las=2, cex=label.cex);
	mtext(ifelse(as.character(loci)=="", "", labels), at=length(loci):1, side=4, line=0, las=2, cex=label.cex, xpd=T);
        # Bold font does not show up unless convert labels to character
	text(xlim[1], length(loci):1, ifelse(as.character(loci)=="", as.character(labels), ""), font=2, cex=label.cex, adj=c(0,0.5), xpd=T);
}
