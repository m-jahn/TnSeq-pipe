#!/usr/bin/Rscript

usage = paste("",
              "Usage: PosHeatmap.R netcdf4File scaffoldId begin end out_pdf",
	      "    Given a netcdf4File that contains per-position average fitness values",
	      "    and the region of the genome of interest, makes a heatmap in out_pdf.",
	      "    The heatmap also has a track that shows the genes in this",
	      "    neighborhood.",
	      "", sep="\n");

library("ncdf4");

PosHeatmap = function(args = commandArgs(trailing=TRUE)) {
	source(file.path(GetPath(), "NetCDFInterface_Functions.R"));

	if (length(args) != 5) stop(usage);
	filename = args[1];
	scaffoldId = args[2];
	begin = as.integer(args[3]);
	end = as.integer(args[4]);
	if (is.na(begin) || is.na(end) || end < begin || end < 0) stop("Illegal coordinates ", args[3], " ", args[4]);
	out_pdf = args[5];

	nc_file = nc_open(filename);
	g = ncvar_get(nc_file, "genes");
	pf = ncvar_get(nc_file, "posfit/pos");
	pf_lrn = ncvar_get(nc_file, "posfit/lrn");
	pf_order = ncvar_get(nc_file, "posfit/clust_order");
	gene_field_axis = ncvar_get(nc_file, "gene_fieldname_axis");
	posfit_field_axis = ncvar_get(nc_file, "posfit_fieldname_axis");
	strainHeatmap_netcdf(posfit = pf, posfit_field_axis = posfit_field_axis,
		posfit_lrn = pf_lrn, posfit_order = pf_order,
		genes = g, gene_field_axis = gene_field_axis,
		scaffoldId=scaffoldId, begin=begin, end=end,
		file=out_pdf, width=7, height=5);
}

GetPath = function() {
    argv = commandArgs(trailingOnly = FALSE)
    dirname(substring(argv[grep("--file=", argv)], 8));
}

# Actually do the work
if(!interactive()) {
	PosHeatmap();
	quit();
}
