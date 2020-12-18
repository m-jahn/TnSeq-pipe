library("ncdf4")

# Gets indices that correspond to the specified gene names in the 
# posfit variables in the specified nc_file
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
posfitIndexByGeneName = function(nc_file, gene_names){
  genes = ncvar_get(nc_file, "genes")
  gene_field_axis = ncvar_get(nc_file, "gene_fieldname_axis")
  indices_genes = geneIndexByGeneName(nc_file, gene_names)
  selected_genes = genes[indices_genes,]
  i_scaffold = which(gene_field_axis == "scaffoldId")
  i_begin = which(gene_field_axis == "begin")
  i_end = which(gene_field_axis == "end") 
  if (is.matrix(selected_genes)){
    g_scaffold = selected_genes[,i_scaffold]
    g_begin = type.convert(selected_genes[,i_begin])
    g_end = type.convert(selected_genes[,i_end])
  }else{
    g_scaffold = selected_genes[i_scaffold]
    g_begin = type.convert(selected_genes[i_begin])
    g_end = type.convert(selected_genes[i_end])  
  }
  pos_list = vector("list", length(g_scaffold))
  for (i in 1:length(g_scaffold)){
    pos_list[[i]] = list(g_scaffold[i], g_begin[i], g_end[i])
  }
  return (posfitIndexByPosition(nc_file, pos_list))
}

# Gets portion of the specified pergene variable in the specified NetCDF
# file that corresponds to the specified gene names
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
posfitVariablesByGeneName = function(nc_file, var_name, gene_names){
  possible_vars = c("posfit/pos", "posfit/lrn")  
  if (!is.element(var_name, possible_vars)){
    return (NULL)
  }
  indices = posfitIndexByGeneName(nc_file, gene_names)
  return (ncvar_get(nc_file, var_name)[indices,])
} 

# Gets indices that correspond to the specified positions in the 
# posfit variables in the specified nc_file
# nc_file: NetCDF file of interest
# pos_list: a list of a 3-element list that includes scaffold ID,
#           starting coordinate, and ending coordinate in order.
posfitIndexByPosition = function(nc_file, pos_list){
  posfit_field_axis = ncvar_get(nc_file, "posfit_fieldname_axis")
  pos = ncvar_get(nc_file, "posfit/pos")
  i_scaffold = which(posfit_field_axis == "scaffold")
  i_pos = which(posfit_field_axis == "pos")
  s_scaffold = pos[,i_scaffold]
  s_pos = type.convert(pos[,i_pos])
  indices = NULL
  for (p in pos_list){
    scaffold = p[[1]]
    start = p[[2]]
    end = p[[3]]
    indices = union(indices, which(s_scaffold==scaffold & s_pos<=end & s_pos>=start))
  }
  return (indices)
}

# Gets portion of the specified posfit variable in the specified NetCDF
# file that corresponds to the specified positions
# nc_file: NetCDF file of interest
# pos_list: a list of a 3-element list that includes scaffold ID,
#           starting coordinate, and ending coordinate in order.
posfitVariablesByPosition = function(nc_file, var_name, pos_list){
  possible_vars = c("posfit/pos", "posfit/lrn")
  if (!is.element(var_name, possible_vars)){
    return (NULL)
  }
  indices = posfitIndexByPosition(nc_file, pos_list)
  return (ncvar_get(nc_file, var_name)[indices,])
}

# Gets indices that correspond to the specified gene names in the 
# perstrain variables in the specified nc_file
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
perstrainIndexByGeneName = function(nc_file, gene_names){
  genes = ncvar_get(nc_file, "genes")
  gene_field_axis = ncvar_get(nc_file, "gene_fieldname_axis")
  indices_genes = geneIndexByGeneName(nc_file, gene_names)
  selected_genes = genes[indices_genes,]
  i_scaffold = which(gene_field_axis == "scaffoldId")
  i_begin = which(gene_field_axis == "begin")
  i_end = which(gene_field_axis == "end") 
  if (is.matrix(selected_genes)){
    g_scaffold = selected_genes[,i_scaffold]
    g_begin = type.convert(selected_genes[,i_begin])
    g_end = type.convert(selected_genes[,i_end])
  }else{
    g_scaffold = selected_genes[i_scaffold]
    g_begin = type.convert(selected_genes[i_begin])
    g_end = type.convert(selected_genes[i_end])  
  }
  pos_list = vector("list", length(g_scaffold))
  for (i in 1:length(g_scaffold)){
    pos_list[[i]] = list(g_scaffold[i], g_begin[i], g_end[i])
  }
  return (perstrainIndexByPosition(nc_file, pos_list))
}

# Gets portion of the specified perstrain variable in the specified NetCDF
# file that corresponds to the specified gene names
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
perstrainVariablesByGeneName = function(nc_file, var_name, gene_names){
  possible_vars = c("fit/strains", "fit/strain_lr", "fit/strain_lrn", "fit/strain_se")  
  if (!is.element(var_name, possible_vars)){
    return (NULL)
  }
  indices = perstrainIndexByGeneName(nc_file, gene_names)
  return (ncvar_get(nc_file, var_name)[indices,])
} 

# Gets indices that correspond to the specified positions in the 
# perstrain variables in the specified nc_file
# nc_file: NetCDF file of interest
# pos_list: a list of a 3-element list that includes scaffold ID,
#           starting coordinate, and ending coordinate in order.
perstrainIndexByPosition = function(nc_file, pos_list){
  strain_field_axis = ncvar_get(nc_file, "strain_fieldname_axis")
  strains = ncvar_get(nc_file, "fit/strains")
  i_scaffold = which(strain_field_axis == "scaffold")
  i_pos = which(strain_field_axis == "pos")
  s_scaffold = strains[,i_scaffold]
  s_pos = type.convert(strains[,i_pos])
  indices = NULL
  for (pos in pos_list){
    scaffold = pos[[1]]
    start = pos[[2]]
    end = pos[[3]]
    indices = union(indices, which(s_scaffold==scaffold & s_pos<=end & s_pos>=start))
  }
  return (indices)
}

# Gets portion of the specified perstrain variable in the specified NetCDF
# file that corresponds to the specified positions
# nc_file: NetCDF file of interest
# pos_list: a list of a 3-element list that includes scaffold ID,
#           starting coordinate, and ending coordinate in order.
perstrainVariablesByPosition = function(nc_file, var_name, pos_list){
  possible_vars = c("fit/strain_lr", "fit/strain_lrn", "fit/strain_se")
  if (!is.element(var_name, possible_vars)){
    return (NULL)
  }
  indices = perstainIndexByPosition(nc_file, pos_list)
  return (ncvar_get(nc_file, var_name)[indices,])
}

# Gets indices that correspond to the specified positions in the 
# pergene variables in the specified nc_file
# nc_file: NetCDF file of interest
# pos_list: a list of a 3-element list that includes scaffold ID,
#           starting coordinate, and ending coordinate in order.
pergeneIndexByPosition = function(nc_file, pos_list){
  pergene_locusId = ncvar_get(nc_file, "pergene_locusId_axis")
  genes_locusId = ncvar_get(nc_file, "locusId_axis")
  genes = ncvar_get(nc_file, "genes")
  gene_field_axis = ncvar_get(nc_file, "gene_fieldname_axis")
  i_scaffold = which(gene_field_axis == "scaffoldId")
  i_begin = which(gene_field_axis == "begin")
  i_end = which(gene_field_axis == "end")  
  g_scaffold = genes[,i_scaffold]
  g_begin = type.convert(genes[,i_begin])
  g_end = type.convert(genes[,i_end])
  indices_genes = NULL
  for (pos in pos_list){
    scaffold = pos[[1]]
    start = pos[[2]]
    end = pos[[3]]
    print(sprintf("%s, %i, %i", scaffold, start, end))
    indices_genes = union(indices_genes, which(g_scaffold==scaffold & g_begin<=end & g_end>=start))
  }
  locusId = genes_locusId[indices_genes]
  indices_pergene = which(is.element(pergene_locusId, locusId))
  return (indices_pergene)
}

# Gets portion of the specified pergene variable in the specified NetCDF
# file that corresponds to the specified positions
# nc_file: NetCDF file of interest
# pos_list: a list of a 3-element list that includes scaffold ID,
#           starting coordinate, and ending coordinate in order.
pergeneVariablesByPosition = function(nc_file, var_name, pos_list){
  possible_vars = c("fit/lr", "fit/lrn", "fit/lr1", "fit/lr2", "fit/lrn1", "fit/lrn2", "fit/lrNaive", "fit/lrRaw", "fit/sumsq", "fit/sd", "fit/sdNaive", "fit/se", "fit/n", "fit/t", "fit/nEff", "fit/tot", "fit/tot0", "fit/tot1", "fit/tot2", "fit/tot1_0", "fit/tot2_0")
  if (!is.element(var_name, possible_vars)){
    return (NULL)
  }
  indices = pergeneIndexByPosition(nc_file, pos_list)
  return (ncvar_get(nc_file, var_name)[indices,])
} 

# Gets indices that correspond to the specified gene names in the 
# pergene variables in the specified nc_file
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
pergeneIndexByGeneName = function(nc_file, gene_names){
  pergene_locusId = ncvar_get(nc_file, "pergene_locusId_axis")
  indices_genes = geneIndexByGeneName(nc_file, gene_names)
  locusId = genes_locusId[indices_genes]
  indices_pergene = which(is.element(pergene_locusId, locusId))
  return (indices_pergene)
}

# Gets portion of the specified pergene variable in the specified NetCDF
# file that corresponds to the specified gene names
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
pergeneVariablesByGeneName = function(nc_file, var_name, gene_names){
  possible_vars = c("fit/lr", "fit/lrn", "fit/lr1", "fit/lr2", "fit/lrn1", "fit/lrn2", "fit/lrNaive", "fit/lrRaw", "fit/sumsq", "fit/sd", "fit/sdNaive", "fit/se", "fit/n", "fit/t", "fit/nEff", "fit/tot", "fit/tot0", "fit/tot1", "fit/tot2", "fit/tot1_0", "fit/tot2_0")
  if (!is.element(var_name, possible_vars)){
    return (NULL)
  }
  indices = pergeneIndexByGeneName(nc_file, gene_names)
  print(indices)
  return (ncvar_get(nc_file, var_name)[indices,])
} 

# Gets indices that correspond to the specified gene names in the 
# variable gene in the specified nc_file
# nc_file: NetCDF file of interest
# gene_names: a vector of gene names of interest. Can be a locus ID,
#             systematic name, or name
geneIndexByGeneName = function(nc_file, gene_names){
  genes_locusId = ncvar_get(nc_file, "locusId_axis")
  genes = ncvar_get(nc_file, "genes")
  gene_field_axis = ncvar_get(nc_file, "gene_fieldname_axis")
  i_name = which(gene_field_axis == "name")
  i_sysname = which(gene_field_axis == "sysName")
  g_name = genes[,i_name]
  g_sysname = genes[,i_sysname]
  indices_genes = which(is.element(genes_locusId, gene_names))
  for (gn in gene_names){
    indices_genes = union(indices_genes, which(g_name==gn))
    indices_genes = union(indices_genes, which(g_sysname==gn))
  }
  return (indices_genes)
}

# Returns a character vector of all variable names
variableNames = function(nc_file){names(nc_file$var)}

# Returns a character vector of all dimension names
dimensionNames = function(nc_file){names(nc_file$dim)}

# Returns a ncvar4 object of a given dimension name
getDimension=function(nc_file, dimensionName){
  nc_file$dim[[dimensionIndex(nc_file, dimensionName)]]
}

# Returns the dimension ID for a dimension
dimensionID = function(nc_file, dimensionName){
  getDimension(nc_file, dimensionName)$id
}

# Returns the dimension indices for a dimension
dimensionIndex = function(nc_file, dimensionName){
  match(c(dimensionName), names(nc_file$dim))
}

# Returns a ncvar4 object of a given variable name
getVariable=function(nc_file, variableName){
  nc_file$var[[variableIndex(nc_file, variableName)]]
}

# Returns the variable ID for a variable
variableID = function(nc_file, variableName){
  getVariable(nc_file, variableName)$id
}

# Returns the variable indices for a variable
variableIndex = function(nc_file, variableName){
  match(c(variableName), names(nc_file$var))
}

# Returns a list of ncdim4 objects corresponding to the dimensions of the specified variable
dimensionsForVariable = function(nc_file, variableName){
  variable = getVariable(nc_file, variableName)
  variable$dim
}

# Returns a character vector of dimension names for a specified variable
dimensionNamesForVariable = function(nc_file, variableName){
  dimensions = dimensionsForVariable(nc_file, variableName)
  names = vector(mode="character", length=length(dimensions))
  for (i in 1:length(dimensions)){
    names[i] = dimensions[[i]]$name
  }
  names
}

# Returns the number of dimensions for a variable
dimensionCountForVariable = function(nc_file, variableName){
  getVariable(nc_file, variableName)$ndim
}

################
### HEATMAPS ###
################

# Returns the color scale (blue to yellow)
myHeatColors = function() {
  c(rgb(0,0,seq(1,1/12,-1/12)),rgb(seq(1/12,1,1/12),seq(1/12,1,1/12),0));
}

# Returns a sequence (seq(-3, 3, 0,25)*2/3) with first and last elements replaced with -100 and 100
# Used to divide the fitness data into intervals, which is then mapped to color
breaksUse = function() { d = seq(-3,3,0.25) * (2/3); d[1] = -100; d[length(d)] = 100; return(d); }

maxSc_rimage = function(genes, fit) {
  d = table(genes$scaffoldId[genes$locusId %in% fit$genesUsed]);
  return( names(d)[which.max(d)] );
}

maxSc_netcdf = function(genes, col_scaffold) {
  d = table(genes[,col_scaffold])
  return (names(d)[which.max(d)]);
}

# Makes a heatmap of the specified scaffoldId and beginning and ending coordinates (in bp). 
# scaffoldId: scaffold ID. possible values are: "Psest_Contig47.1", "Psest_Contig45.2",
# "Psest_Contig44.3", and "Psest_Contig40.4"
# begin: starting coordinate in bp
# end: ending coordinate in bp
# file: title of the exported PDF file
# expand: 0.5 * width of the bars in heatmap
# geneFrac: fraction of the pdf file reserved for labels
# linesAt: vector of coordinates (in bp) to which a vertical line will be drawn
# lineCol: color of the line drawn
strainHeatmap_netcdf = function(posfit, posfit_field_axis, posfit_lrn, posfit_order, genes, gene_field_axis, begin, end, scaffoldId=NULL, clustOrder = TRUE, expand=20, geneFrac=0.1, linesAt = NULL, lineCol=2, file="perstrain.pdf", width=6, height=4) {
  gene_col_scaffold = which(gene_field_axis=="scaffoldId")
  gene_col_begin=which(gene_field_axis == "begin")
  gene_col_end = which(gene_field_axis == "end")
  gene_col_strand = which(gene_field_axis == "strand")  
  if (is.null(scaffoldId)){
    scaffoldId = maxSc_netcdf(genes, gene_col_scaffold)
  }
  # gets indices in the posfit/posfit_lrn table that corresponds to scaffoldId, beginning coordinate,
  # and ending coordinate specified.
  posfit_col_scaffold = which(posfit_field_axis=="scaffold")
  posfit_col_pos = which(posfit_field_axis=="pos")
  indices = which(posfit[,posfit_col_scaffold]==scaffoldId & as.integer(posfit[,posfit_col_pos])>=begin & as.integer(posfit[,posfit_col_pos])<=end)
  # gets the position variable of the selected data points in posfit data
  pos = as.integer(posfit[indices,posfit_col_pos])
  order = order(pos)
  # sorts position
  pos = pos[order]
  # gets fitness data of the specified scaffoldId, beginning coordinate, and ending coordinate
  # and sorts by position
  posfitSet = posfit_lrn[indices,][order,posfit_order]
  # defines the boundaries for each vertical strip in heatmaps
  posmin = pmax(pos-expand, c(begin,0.5+pos[-length(pos)])); 
  posmax = pmin(pos+expand, c(pos[-1]-0.5, end));
  
  if(!is.null(file)) pdf(file, width=width, height=height, pointsize=10); 
  oldpar = par(no.readonly=T);
  layout(matrix(c(1,2),ncol=1), heights=c(1-geneFrac,geneFrac)); 
  par(mar=c(0,4,0,0.5),oma=c(4,0,0,0));  
  par(xaxs="i", yaxs="i", bty="n"); 
  
  breaks = breaksUse();
  cols = myHeatColors();
  
  
  plot(c(begin,end)/1000, c(1, ncol(posfitSet)+1), pch=NA, xaxt="n",
       xlab="", ylab=sprintf("%d conditions", ncol(posfitSet)));
  rect(begin/1000, 1, end/1000, ncol(posfitSet)+1, col="grey", border=F); 
  # for each position in the fitness data, maps the fitness value to a color for each condition
  # and draws a strip in the appropriate location.
  for(i in 1:ncol(posfitSet)) { 
    colsUsed = cols[cut(posfitSet[,i], breaks, ordered=T)]; 
    rect(posmin/1000, rep(i,length(pos)), posmax/1000, rep(i+1,length(pos)),
         col=colsUsed, border=NA); 
  }
  # vertical lines at specified positions
  for (i in linesAt) abline(v=i/1000, col=lineCol);
  
  par(xpd=T);
  plot(c(begin, end)/1000, c(-1,1), pch=NA,
       yaxt="n", ylab="", xlab="");
  # gets genes of the specified scaffold and position range
  geneSet = genes[as.integer(genes[,gene_col_begin]) <= end & as.integer(genes[,gene_col_end]) >= begin & genes[,gene_col_scaffold]==scaffoldId,]
  if (!is.matrix(geneSet)){
    geneSet = t(as.matrix(geneSet))
    print(geneSet)
  }
  gx1 = as.integer(ifelse(geneSet[,gene_col_strand]=="+",geneSet[,gene_col_begin], geneSet[,gene_col_end]));
  gx1 = pmin(pmax(gx1, begin-1), end+1);
  gx2 = as.integer(ifelse(geneSet[,gene_col_strand]=="+",geneSet[,gene_col_end], geneSet[,gene_col_begin]));
  gx2 = pmin(pmax(gx2, begin-1), end+1);
  gy = 0
  # draws arrows for each gene and labels with gene name
  arrows(gx1/1000, gy, gx2/1000, gy, code=2, length=1/12);
  text((gx1+gx2)/2000, 0.5*ifelse(1:nrow(geneSet)%%2 == 0, 1, -1), geneSet[,1])
  mtext("Position (kilobases)", side=1, outer=T, line=2); 
  for (i in linesAt) abline(v=i/1000, col=lineCol);
  if(is.null(file)) { par(oldpar); } else { dev.off(); cat("Wrote ",file,"\n");
  }
}

# Makes a heatmap of the specified scaffoldId and beginning and ending coordinates (in bp). 
# scaffoldId: scaffold ID. possible values are: "Psest_Contig47.1", "Psest_Contig45.2",
# "Psest_Contig44.3", and "Psest_Contig40.4"
# begin: starting coordinate in bp
# end: ending coordinate in bp
# file: title of the exported PDF file
# expand: 0.5 * width of the bars in heatmap
# geneFrac: fraction of the pdf file reserved for labels
# linesAt: vector of coordinates (in bp) to which a vertical line will be drawn
# lineCol: color of the line drawn
strainHeatmap_rimage = function(posfit, posfitClust, genes, fit, begin, end, scaffoldId=NULL, clustOrder = TRUE, metacol=1:3, expand=20, geneFrac=0.1,
                                linesAt = NULL, lineCol=2, file="public_html/tmp/out.pdf") {
  if (is.null(scaffoldId)){
    scaffoldId=maxSc_rimage(genes, fit)
  }
  geneSet = genes[genes$begin <= end & genes$end >= begin & genes$scaffoldId==scaffoldId,]; ##gets genes that are between begin and end
  posfitSet = posfit[posfit$pos >= begin & posfit$pos <= end & posfit$scaffold==scaffoldId,]; ##gets fitness for strains between begin and end
  pos = posfitSet$pos; ##gets the positions for the selected fitness values
  if (clustOrder){
    posfitSet = posfitSet[,-metacol][,posfitClust$order]; 
  }else{
    posfitSet = posfitSet[,-metacol]
  }
  # expand each strain by expand, but leave the neighbor strain at least 0.5 space
  posmin = pmax(pos-expand, c(0,0.5+pos[-length(pos)])); 
  posmax = pmin(pos+expand, c(pos[-1]-0.5, end+1));
  
  if(!is.null(file)) pdf(file, width=6, height=4, pointsize=10); ##makes a pdf graphic
  oldpar = par(no.readonly=T);
  layout(matrix(c(1,2),ncol=1), heights=c(1-geneFrac,geneFrac)); ##Splits the pdf graphic in sections. 1 figure at top with height 1-geneFrac and another figre at the bottom with geneFrac
  par(mgp=c(2,1,0)); ##defines how the axis labels are arranged
  par(mar=c(0,4,0,0.5),oma=c(4,0,0,0));  # bottom,left,top,right ## sets the margins and the outermargins
  par(xaxs="i", yaxs="i", bty="n"); ##makes sure x-axis and y-axis correspond to specified ranges. ##bty: no box around the chart
  
  breaks = breaksUse();
  cols = myHeatColors();
  
  plot(c(begin,end)/1000, c(1, ncol(posfitSet)+1), pch=NA, xaxt="n",
       xlab="", ylab=sprintf("%d conditions", ncol(posfitSet))); ## draws the axes
  rect(begin/1000, 1, end/1000, ncol(posfitSet)+1, col="grey", border=F); ## draw a grey rectangle
  for(i in 1:ncol(posfitSet)) { ##for each position fitness data
    colsUsed = cols[cut(posfitSet[,i], breaks, ordered=T)]; ##maps the fitness value to a color for each condition
    rect(posmin/1000, rep(i,length(pos)), posmax/1000, rep(i+1,length(pos)),
         col=colsUsed, border=NA); ##draws a strip from posmin/1000 to posmax/1000 with the appropriate colors
    
  }
  for (i in linesAt) vline(i/1000, col=lineCol); ##adds lines
  
  par(xpd=T); # no cropping
  plot(c(begin,end)/1000, c(-1,1), pch=NA,
       yaxt="n", ylab="", xlab="");
  gx1 = ifelse(geneSet$strand=="+",geneSet$begin,geneSet$end);
  gx2 = ifelse(geneSet$strand=="+",geneSet$end,geneSet$begin);
  gy = ifelse(geneSet$strand=="+",0.1,-0.1);
  arrows(gx1/1000, gy, gx2/1000, gy, code=2, length=1/12);
  text((gx1+gx2)/2000, 0.5*ifelse(1:nrow(geneSet)%%2 == 0, 1, -1),
       sub("Psest_","", sub("Sama_","",geneSet$sysName)));
  mtext("Position (kilobases)", side=1, outer=T, line=2);  	 # bottom
  for (i in linesAt) vline(i/1000, col=lineCol);
  
  if(is.null(file)) { par(oldpar); } else { dev.off(); cat("Wrote ",file,"\n"); }
}


