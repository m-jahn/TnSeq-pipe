#!/usr/bin/Rscript
# By Shiori Sagawa and Morgan Price, UC Berkeley / LBL, 2015

usage = paste("",
      		"Usage: StrainCorPlot.R netcdf4File locusId out_pdf",
      		"       StrainCorPlot.R netcdf4File locusId scaffoldId begin end out_pdf",
		"    Given a netcdf4File that contains gene and strain fitness values,",
		"    the gene of interest, and the region of interest, makes a scatterplot",
		"    into out_pdf with how similar each strain's fitness pattern is",
		"    to the gene's fitness pattern. A second track shows the genes in this",
		"    neighborhood.",
		"", sep="\n");

library("ncdf4");

StrainCorPlot = function(args = commandArgs(trailing=TRUE)) {
	source(file.path(GetPath(), "NetCDFInterface_Functions.R"));
	

	if (length(args) == 3) {
		nc_file = args[1];
		locusId = args[2];
		scaffoldId = NULL;
		out_pdf = args[3];
	} else if (length(args) == 6) {
		nc_file = args[1];
		locusId = args[2];
		scaffoldId = args[3];
		begin = as.integer(args[4]);
		end = as.integer(args[5]);
		if (is.na(begin) || is.na(end) || end < begin || end < 0)
			stop("Illegal coordinates ", args[3], " ", args[4]);
		out_pdf = args[6];
	} else {
		stop(usage);
	}
	nc = nc_open(nc_file);
	load_nc_variables(nc);

	i = which( gene_locusId_axis %in% locusId | genes[,gene_col_name]==locusId | genes[,gene_col_sysname]==locusId );
	if (length(i) < 1) stop("No such gene: ", locusId);
	i = i[1];

	load_strains(nc);

	if (is.null(scaffoldId)) {
		scaffoldId = genes[i, gene_col_scaffold];
		gBeg = as.integer(genes[i, gene_col_begin]);
		gEnd = as.integer(genes[i, gene_col_end]);
		begin = pmax(1, gBeg - 3000);
		end = gEnd + 3000;
	}
	scatter_plot(begin, end, locusId, scaffoldId=scaffoldId, file=out_pdf);
}

load_strains = function(nc){
  strains = ncvar_get(nc, "fit/strains")
  new_df = data.frame(matrix(0, ncol=ncol(strains), nrow=nrow(strains)))
  for (i in 1:ncol(strains)){
    new_df[,i] = type.convert(strains[,i])
  }
  colnames(new_df) = ncvar_get(nc, "strain_fieldname_axis")
  s <<- new_df;
  return(NULL);
}

load_nc_variables = function(nc){
  gene_field_axis <<- ncvar_get(nc, "gene_fieldname_axis")
  nc_lrn <<- ncvar_get(nc, "fit/lrn")
  pergene_locusId_axis <<- ncvar_get(nc, "pergene_locusId_axis") # for lrn
  gene_locusId_axis <<- ncvar_get(nc, "locusId_axis") # for genes
  expsmeta_field_axis <<- ncvar_get(nc, "exps_metadata_fieldname_axis")
  expsmeta <<- ncvar_get(nc, "fit/q")
  genes <<- ncvar_get(nc, "genes")
  s_lrn <<- ncvar_get(nc, "fit/strain_lrn")
  gene_col_name <<- which(gene_field_axis=="name")
  gene_col_sysname <<- which(gene_field_axis=="sysName")
  gene_col_scaffold <<- which(gene_field_axis=="scaffoldId")
  gene_col_begin <<- which(gene_field_axis == "begin")
  gene_col_end <<- which(gene_field_axis == "end")
  gene_col_strand <<- which(gene_field_axis == "strand")
  expsmeta_col_u <<- which(expsmeta_field_axis == "u")
}

scatter_plot= function(begin, end, refGeneName, scaffoldId=NULL, file="scatter.pdf", showUnused=TRUE,
	               width=7, height=5, pointsize=10, geneFrac=0.2) {
  refLocusId = gene_locusId_axis[gene_locusId_axis %in% refGeneName | genes[,gene_col_name]==refGeneName | genes[,gene_col_sysname]==refGeneName ];
  ref_gene_names = genes[which(gene_locusId_axis == refLocusId), c(gene_col_name, gene_col_sysname)]
  exp_u = type.convert(expsmeta[,expsmeta_col_u]);
  pergene_vector = nc_lrn[which(pergene_locusId_axis == refLocusId), exp_u]
  if (is.null(scaffoldId)){
    scaffoldId = maxSc_netcdf(genes, gene_col_scaffold)
  }

  index = which(s$pos >= begin & s$pos<=end & s$scaffold==scaffoldId & (showUnused | s$used));
  strand = as.character(s$strand[index]);
  r = apply(s_lrn[index, exp_u], 1, cor, pergene_vector);
  cat("Computed correlations for ", refLocusId, " with strains using ", sum(exp_u), " experiments\n");

  if(!is.null(file)) pdf(file, width=width, height=height, pointsize=pointsize);
  oldpar = par(no.readonly=T);
  layout(matrix(c(1,2),ncol=1), heights=c(1-geneFrac,geneFrac)); 
  par(mar=c(1,4,2,2),oma=c(0,0,0,0));  
  par(xaxs="i", yaxs="i",bty="n"); 

  gene_header = if (as.character(ref_gene_names[1]) == "NA") ref_gene_names[2] else sprintf("%s (%s)", ref_gene_names[1], ref_gene_names[2]);
  if (gene_header == "NA") gene_header = refLocusId;
  plot(s$pos[index]/1000, r,
		ylim=c(-1,1), xlim=c(floor(begin/1000), ceiling(end/1000)),
		ylab="Correlation",
		main=paste("Strain Similarity with", gene_header),
		col=ifelse(s$used[index], ifelse(strand=="+", "darkgreen", "red"), "darkgrey"),
		pch=strand,
		cex=ifelse(s$used[index], 1, 2/3) * ifelse(strand=="+",1,1.25),
		bty="o"); # box around the plot
  abline(h=0, col="grey");
  par(xpd=T);
  mtext("Position (kilobases)", side=1, line=2)

  plot(c(begin, end)/1000, c(-1,1), pch=NA,
       yaxt="n", xaxt="n", ylab="", xlab="Position (kb)");
  # gets genes of the specified scaffold and position range
  genesI = which(as.integer(genes[,gene_col_begin]) <= end & as.integer(genes[,gene_col_end]) >= begin & genes[,gene_col_scaffold]==scaffoldId);
  geneSet = genes[genesI,];
  if (!is.matrix(geneSet)){
    geneSet = t(as.matrix(geneSet))
  }
  y1 = 0
  y2 = -0.75
  if (nrow(geneSet)!=0){
    gx1 = as.integer(ifelse(geneSet[,gene_col_strand]=="+",geneSet[,gene_col_begin], geneSet[,gene_col_end]));
    gx1 = pmin(pmax(gx1, begin-1), end+1);
    gx2 = as.integer(ifelse(geneSet[,gene_col_strand]=="+",geneSet[,gene_col_end], geneSet[,gene_col_begin]));
    gx2 = pmin(pmax(gx2, begin-1), end+1);
    gy = ifelse(geneSet[,gene_col_strand]=="+", y1, y2)
    # draws arrows for each gene and labels with gene name
    arrows(gx1/1000, gy, gx2/1000, gy, code=2, length=1/12);
    geneStrand = geneSet[,gene_col_strand];
    namesAbove = as.character(geneSet[,gene_col_name]);
    namesAbove = ifelse(namesAbove=="NA","",namesAbove);
    text((gx1+gx2)/2000, ifelse(geneStrand=="+",y1, y2), namesAbove, pos=3); # above
    namesBelow = geneSet[,gene_col_sysname];
    namesBelow = sub("^.*_","_",namesBelow,perl=T);
    namesBelow = ifelse(namesBelow=="NA", gene_locusId_axis[genesI], namesBelow);
    text((gx1+gx2)/2000, ifelse(geneStrand=="+",y1, y2), namesBelow, pos=1); # below
  }
  if(is.null(file)) { par(oldpar); } else { dev.off(); cat("Wrote ",file,"\n");
  }
}

GetPath = function() {
    argv = commandArgs(trailingOnly = FALSE)
    dirname(substring(argv[grep("--file=", argv)], 8));
}

# Actually do the work
if (!interactive()) {
	StrainCorPlot();
	quit();
}
