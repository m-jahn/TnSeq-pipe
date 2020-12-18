#!/usr/bin/Rscript
# By Shiori Sagawa, UC Berkeley/LBL, 2015

library("ncdf4")

usage = paste("",
                "Usage: StrainFitPlot.R netcdf4File genename_or_locusId out_pdf experimentIds(comma-delimited)",
                "       StrainFitPlot.R netcdf4File scaffoldId begin end out_pdf experimentIds(comma-delimited)",
                "    Given a netcdf4File that contains gene and strain fitness values, ",
                "    a gene or a region of interest, and experiments of interest, makes",
		"    a scatterplot into out_pdf. The scatterplot plots each strain's average",
		"    fitness values in the specified experiments against the strain's position.",
		"    If locusId is specified, the plot includes all strains with an assigned",
		"    locusId (i.e. all strains with insertions within 10-90% of the gene, ",
		"    including strains that are not used to estimate gene's fitness). If the",
		"    region is specified instead, then the plot includes all strains,",
		"    including strains without a locusId.",
                "",
               sep='\n')

StrainFitPlot = function(args = commandArgs(trailing=TRUE)){
  getData = function(command){
    cat("Running: ",command,"\n");
    text = system(command, intern=T)
    text = paste(text, collapse="\n")
    text_handle = textConnection(text)
    df = read.csv(text_handle, sep='\t', quote="", header=T, stringsAsFactors = F)
    close(text_handle)
    return (df)
  }
  if (length(args) %in% c(4,6)){
    nc_file = args[1]
    nc = nc_open(nc_file)
    geneCols = ncvar_get(nc, "gene_fieldname_axis")
    locusId = ncvar_get(nc, "locusId_axis")
    locusId = sub("^ +", "", locusId, perl=T)
    locusId = sub(" +$", "", locusId, perl=T);  
    genes = ncvar_get(nc, "genes")
    for (i in 1:length(geneCols)){
      genes[,i] = sub("^ +", "", genes[,i], perl=T)
      genes[,i] =sub(" +$", "", genes[,i], perl=T);  
    }
    colnames(genes) = geneCols
    genes = data.frame(cbind(locusId, genes), stringsAsFactors = F)
    genes$begin = as.integer(genes$begin)
    genes$end = as.integer(genes$end)
  }
  straincmd = file.path(GetPath(), "StrainData.R");
  if (length(args)==4){
    gene = args[2]
    out_pdf = args[3]
    exps = unlist(strsplit(args[4], ',', fixed=T))
    locusId = genes$locusId[genes$locusId==gene | genes$name==gene | genes$sysName==gene]
    data = getData(sprintf("%s %s %s", straincmd,  nc_file, locusId));
    xmin = genes$begin[genes$locusId==gene | genes$name==gene | genes$sysName==gene]
    xmax = genes$end[genes$locusId==gene | genes$name==gene | genes$sysName==gene]     
  }else if (length(args)==6){
    scaffoldId = args[2]
    begin = as.integer(args[3])
    end = as.integer(args[4])
    if (is.na(begin) || is.na(end) || end<begin || end<0){
      stop("Illegal coordinates, ", args[2], " ", args[3], " ", args[4])
    }
    out_pdf = args[5]
    exps = unlist(strsplit(args[6], ',', fixed=T))
    data = getData(sprintf("%s %s %s %d %d", straincmd, nc_file, scaffoldId, begin, end))
    xmin = begin
    xmax = end
  }else{
    stop(usage)
  }
  if (nrow(data)==0){
    stop("No data to be plotted. This may be due to lack of strains in the specified gene/region or invalid arguments")
  }
  if (sum(exps %in% colnames(data))!=length(exps)){
    print(sprintf("The NetCDF file does not contain data for %s", paste(exps[!exps %in% colnames(data)], collapse=', ')))
  }
  if (sum(exps %in% colnames(data))==0){
    stop("No data to be plotted. The NetCDF file does not contain data for any of the specified experiments")
  }
  metacols = c("barcode", "rcbarcode", "scaffold", "strand", "pos", "locusId", "f", "used", "name")
  data = data[,colnames(data)  %in% c(metacols, exps)]
  fitness = data[, colnames(data) %in% exps]
  if (is.vector(fitness)){
    fitness = as.matrix(fitness)
  }
  data["avg"] = rowMeans(fitness)
  geneinfo = genes[match(unique(data$locusId), genes$locusId),]
  geneinfo$name = ifelse(geneinfo$name=="", geneinfo$sysName, geneinfo$name)
  geneinfo$begin = as.integer(geneinfo$begin)
  geneinfo$end = as.integer(geneinfo$end)
  pclass = ifelse(data$used, ifelse(data$strand=='+', 1, 2), ifelse(data$strand=='+', 3, 4))
  pcols = c( "dark green", "red", "grey45", "grey45")
  psizes = c(1.25, 1.5, 0.75, 0.75)
  pshapes = c('+', '-', '+', '-') 
  pdf(out_pdf, width=6, height=4, pointsize=10)
  oldmar=par("mar")
  newmar=oldmar
  newmar[1]=newmar[1]*1.25
  par(mar = newmar)
  plot(data$pos, data$avg, col = pcols[pclass], pch=pshapes[pclass], cex=as.numeric(psizes[pclass]), 
       xlab="", main=sprintf("Position vs. Average Fitness for %s", paste(exps[exps %in% colnames(data)], collapse=", ")),
       ylab="Average Strain Fitness", font.lab=2, xlim=c(xmin, xmax), cex.main=0.85, xaxs='i')
  abline(h=0, col="black", lty=2)
  arrows(ifelse(geneinfo$strand=="+", pmax(xmin,geneinfo$begin), pmin(xmax,geneinfo$end)), 
         min(data$avg)-0.3*(max(data$avg)-min(data$avg)), 
         ifelse(geneinfo$strand=="+", pmin(xmax,geneinfo$end), pmax(geneinfo$begin,xmin)),
         xpd=NA, length=0.1)
  text(x=pmin(pmax((geneinfo$begin + geneinfo$end)*0.5, xmin),xmax),
       y=min(data$avg)-0.4*(max(data$avg)-min(data$avg)),
       labels=geneinfo$name,xpd=NA)
  text(x=0.5*(xmin+xmax), labels="Position (KB)", xpd=NA, 
       y=min(data$avg)-0.5*(max(data$avg)-min(data$avg)), font=2)
  par(mar=oldmar)
  dev.off()
}

GetPath = function() {
    argv = commandArgs(trailingOnly = FALSE)
    dirname(substring(argv[grep("--file=", argv)], 8));
}

if (!interactive()) {
  StrainFitPlot();
  quit();
}
