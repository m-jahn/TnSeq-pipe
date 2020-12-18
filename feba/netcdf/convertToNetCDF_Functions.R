# The maximum tolerated difference between a value stored in R image
# and a value stored in NetCDF files.
Max_diff = 0.001

# Sets up posfit.
# fit: the fit variable loaded from fit.image. If NULL, image specified
#      by fit_image will be loaded.
# posfit_path: If posfit_path is NULL, then the function does not store
#              posfit into an R image. Otherwise, posfit will be stored
#              as an R image in the path. 
# fit_image: path to a fit.image the posfit will be based off of. The
#            image will be used only if fit is NULL
SetupPosFit = function(fit=NULL, posfit_path = NULL, fit_image = "fit.image"){
  if (is.null(fit)){
    load(fit_image)
  }
  u = fit$q$u
  d = cbind(fit$strains[,c("scaffold", "strand", "pos")],
            lapply(split(as.character(fit$q$name[u]), fit$q$short[u]),
                   function(x) rowMeans(fit$strain_lrn[,x,drop=F])));
  d = d[d$scaffold != "pastEnd" & apply(fit$strain_se,1,median) < 1,];
  #cat("Using ", nrow(d), " strains for posfit\n");
  posfit <<- aggregate(d[,-(1:3)], d[,1:3], mean);
  #cat("Returning ", nrow(posfit), " rows in posfit\n");
  posfitClust <<- hclust(as.dist(cor(posfit[,-(1:3)])));
  if (!is.null(posfit_path)){
    save(posfit, posfitClust, file=posfit_path)
  }
  return(NULL);
}

# Converts an R image specified by fit_image_path into a NetCDF file
# and stores the NetCDF file in the nc_path. 
# c: compression level
# chunk: chunking size
# fit_image_path: path to fit.image of interest
# nc_path: path to the NetCDF file to be created
# posfit_path: If posfit_path is NULL, then the function does not store
#              posfit into an R image. Otherwise, posfit will be stored
#              as an R image in the path. 
make_ncdf4 <- function(c, chunk, fit_image_path, nc_path, posfit_path=NULL){
  start = proc.time()
  print("loading image...")
  load(fit_image_path)
  print("setting up posfit...")
  load(posfit_path)
  #SetupPosFit(fit, posfit_path)
  print("defining axes...")
  num_pos = nrow(posfit)
  num_strain = nrow(fit$strains)
  num_set= ncol(fit$strain_lrn)
  num_posfit_set = ncol(posfit)-3
  num_strain_desc = ncol(fit$strains)
  num_char = 300
  num_gene_pergene = nrow(fit$lrn)
  num_exp_metadata = ncol(fit$q)-1
  num_expsUsed = nrow(expsUsed)
  num_expsUsed_info = ncol(expsUsed)
  num_gene = nrow(genes)
  num_gene_info = ncol(genes)-1
  #DIMENSIONS#
  strain_axis = ncdim_def("strainId_axis", "", 1:num_strain)
  set_axis = ncdim_def("setId_axis", "", 1:num_set, create_dimvar=FALSE)
  strain_desc_axis = ncdim_def("strain_fieldname_axis", "", 1:num_strain_desc, create_dimvar=FALSE)
  char_axis = ncdim_def("string_axis", "", 1:num_char)
  gene_pergene_axis = ncdim_def("pergene_locusId_axis", "", 1:num_gene_pergene, create_dimvar=FALSE)
  exp_metadata_axis = ncdim_def("exps_metadata_fieldname_axis", "", 1:num_exp_metadata, create_dimvar=FALSE)
  expsUsed_axis = ncdim_def("expsId_axis", "", 1:num_expsUsed)
  expsUsed_info_axis = ncdim_def("exps_fieldname_axis", "", 1:num_expsUsed_info, create_dimvar=FALSE)
  gene_axis = ncdim_def("locusId_axis", "", 1:num_gene, create_dimvar=FALSE)
  gene_info_axis = ncdim_def("gene_fieldname_axis", "", 1:num_gene_info, create_dimvar=FALSE)
  posfit_desc_axis = ncdim_def("posfit_fieldname_axis", "", 1:3, create_dimvar = FALSE)
  posfit_axis = ncdim_def("posId_axis", "", 1:num_pos)
  posfit_set_axis = ncdim_def("posfit_set_axis", "", 1:num_posfit_set, create_dimvar = FALSE)
  #VARIABLES#
  print("defining variables...")
  gene_axis_var = ncvar_def("locusId_axis", "", list(char_axis, gene_axis), prec="char", compression=c, chunksizes=c(num_char, min(chunk, num_gene)))
  set_axis_var = ncvar_def("setId_axis", "", list(char_axis, set_axis), prec="char", compression=c, chunksizes=c(num_char,min(chunk, num_set)))
  strain_desc_axis_var = ncvar_def("strain_fieldname_axis", "", list(char_axis, strain_desc_axis), prec="char", compression=c, chunksizes=c(num_char,min(chunk, num_strain_desc)))
  gene_pergene_axis_var = ncvar_def("pergene_locusId_axis", "", list(char_axis, gene_pergene_axis), prec="char", compression=c, chunksizes=c(num_char,min(num_gene_pergene,chunk)))
  exp_metadata_axis_var = ncvar_def("exps_metadata_fieldname_axis", "", list(char_axis, exp_metadata_axis), prec="char", compression=c, chunksizes=c(num_char,min(num_exp_metadata, chunk)))
  expsUsed_info_axis_var = ncvar_def("exps_fieldname_axis", "", list(char_axis, expsUsed_info_axis), prec="char",compression=c,  chunksizes=c(num_char,min(num_expsUsed_info, chunk)))
  gene_info_axis_var = ncvar_def("gene_fieldname_axis", "",list(char_axis, gene_info_axis), prec="char", compression=c, chunksizes=c(num_char,min(chunk, num_gene_info)))
  posfit_desc_axis_var = ncvar_def("posfit_fieldname_axis", "", list(char_axis, posfit_desc_axis), prec="char", compression=c, chunksizes=c(num_char, 3))
  gene_var = ncvar_def("genes", "", list(char_axis, gene_axis, gene_info_axis), prec="char", compression=c, chunksizes=c(num_char, chunk, num_gene_info))
  expsUsed_var = ncvar_def("expsUsed", "", list(char_axis, expsUsed_axis, expsUsed_info_axis), prec="char", compression=c, chunksizes=c(num_char, chunk, num_expsUsed_info))
  posfit_set_axis_var = ncvar_def("posfit_set_axis", "", list(char_axis, posfit_set_axis), prec = "char", compression=c, chunksizes=c(num_char, num_posfit_set))
  #experiment metadata
  q_var = ncvar_def("fit/q", "", list(char_axis, set_axis, exp_metadata_axis), prec="char", compression=c, chunksizes=c(num_char,chunk, num_exp_metadata))
  #per-gene data
  lrRaw_var = ncvar_def("fit/lrRaw", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))
  lrn_var = ncvar_def("fit/lrn", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))
  lrn1_var = ncvar_def("fit/lrn1", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))
  lrn2_var = ncvar_def("fit/lrn2", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  lr_var = ncvar_def("fit/lr", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))
  lr1_var = ncvar_def("fit/lr1", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))
  lr2_var = ncvar_def("fit/lr2", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  lrNaive_var = ncvar_def("fit/lrNaive", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  t_var = ncvar_def("fit/t", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  sd_var = ncvar_def("fit/sd", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))
  se_var = ncvar_def("fit/se", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  sdNaive_var = ncvar_def("fit/sdNaive", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  sumsq_var = ncvar_def("fit/sumsq", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  n_var = ncvar_def("fit/n", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk,num_set))
  nEff_var = ncvar_def("fit/nEff", "", list(gene_pergene_axis, set_axis), chunksizes=c(chunk, num_set))  
  tot_var = ncvar_def("fit/tot", "", list(gene_pergene_axis, set_axis), prec="integer", chunksizes=c(chunk, num_set))  
  tot0_var = ncvar_def("fit/tot0", "", list(gene_pergene_axis, set_axis), prec="integer", chunksizes=c(chunk, num_set))
  tot1_var = ncvar_def("fit/tot1", "", list(gene_pergene_axis, set_axis), prec="integer", chunksizes=c(chunk, num_set))
  tot2_var = ncvar_def("fit/tot2", "", list(gene_pergene_axis, set_axis), prec="integer", chunksizes=c(chunk, num_set))
  tot1_0_var = ncvar_def("fit/tot1_0", "", list(gene_pergene_axis, set_axis), prec="integer", chunksizes=c(chunk, num_set))
  tot2_0_var = ncvar_def("fit/tot2_0", "", list(gene_pergene_axis, set_axis), prec="integer", chunksizes=c(chunk, num_set))
  #per-strain data
  strain_var = ncvar_def("fit/strains", "", list(char_axis, strain_axis, strain_desc_axis), prec="char", chunksizes=c(num_char, chunk, num_strain_desc))
  strain_lr_var = ncvar_def("fit/strain_lr", "", list(strain_axis, set_axis), chunksizes = c(chunk, num_set))
  strain_lrn_var = ncvar_def("fit/strain_lrn", "", list(strain_axis, set_axis), chunksizes=c(chunk, num_set))
  strain_se_var = ncvar_def("fit/strain_se", "", list(strain_axis, set_axis), chunksizes=c(chunk, num_set))
  #posfit data
  posfit_var = ncvar_def("posfit/pos", "", list(char_axis, posfit_axis, posfit_desc_axis), prec="char", chunksizes=c(num_char, chunk, 3))
  posfit_lrn_var = ncvar_def("posfit/lrn", "", list(posfit_axis, posfit_set_axis), chunksizes=c(chunk, num_posfit_set))
  posfit_clust_order_var = ncvar_def("posfit/clust_order", "", list(posfit_set_axis), prec="integer")
  #PUTTING DATA IN#
  print("writing data...")
  nc = nc_create(nc_path, list(lrn1_var, lrn2_var, lrNaive_var, sumsq_var, sd_var, lrRaw_var, strain_lr_var, n_var, gene_axis_var, posfit_clust_order_var, posfit_set_axis_var, posfit_var, posfit_lrn_var, posfit_desc_axis_var, set_axis_var, strain_desc_axis_var, gene_pergene_axis_var, exp_metadata_axis_var, expsUsed_info_axis_var, gene_info_axis_var, gene_var, expsUsed_var, strain_var, strain_lrn_var, strain_se_var, q_var, lrn_var, lr_var, lr1_var, lr2_var, t_var, se_var, sdNaive_var, nEff_var, tot_var, tot0_var, tot1_var, tot2_var, tot1_0_var, tot2_0_var))
  locusId_col = which(colnames(genes)=="locusId")
  ncvar_put(nc, gene_axis_var, genes$locusId)
  ncvar_put(nc, set_axis_var, colnames(fit$strain_lrn))
  ncvar_put(nc, strain_desc_axis_var, colnames(fit$strains))
  ncvar_put(nc, gene_pergene_axis_var, fit$g)
  
  ncvar_put(nc, exp_metadata_axis_var, colnames(fit$q)[-1])
  ncvar_put(nc, expsUsed_info_axis_var, colnames(expsUsed))
  ncvar_put(nc, gene_info_axis_var, colnames(genes)[-locusId_col])
  ncvar_put(nc, posfit_desc_axis_var, c("scaffold", "strand", "pos"))
  ncvar_put(nc, posfit_set_axis_var, colnames(posfit)[-1:-3])

  ncvar_put(nc, sd_var, as.vector(as.matrix(fit$sd)))
  ncvar_put(nc, sumsq_var, as.vector(as.matrix(fit$sumsq)))
  ncvar_put(nc, lrRaw_var, as.vector(as.matrix(fit$lrRaw)))
  ncvar_put(nc, n_var, as.vector(as.matrix(fit$n)))
  ncvar_put(nc, gene_var, as.vector(as.matrix(genes[,-locusId_col])))
  ncvar_put(nc, expsUsed_var, as.vector(as.matrix(expsUsed)))
  used_col = which(colnames(fit$strains)=="used")
  fit_strains = as.vector(as.matrix(cbind(cbind(fit$strains[,1:(used_col-1)], as.character(fit$strains[,used_col])),fit$strains[,-(1:used_col)])))
  ncvar_put(nc, strain_var, fit_strains)
  ncvar_put(nc, strain_lr_var, as.vector(as.matrix(fit$strain_lr)))
  ncvar_put(nc, strain_lrn_var, as.vector(as.matrix(fit$strain_lrn)))
  ncvar_put(nc, strain_se_var, as.vector(as.matrix(fit$strain_se)))
  u_col = which(colnames(fit$q)=="u")
  fit_q = as.vector(as.matrix(cbind(cbind(fit$q[,2:(u_col-1)], as.character(fit$q[,u_col])), fit$q[,-(1:u_col)])))
  ncvar_put(nc, q_var, fit_q)
  ncvar_put(nc, lrn_var, as.vector(as.matrix(fit$lrn)))
  ncvar_put(nc, lrn1_var, as.vector(as.matrix(fit$lrn1)))
  ncvar_put(nc, lrn2_var, as.vector(as.matrix(fit$lrn2)))  
  ncvar_put(nc, lr_var, as.vector(as.matrix(fit$lr)))
  ncvar_put(nc, lrNaive_var, as.vector(as.matrix(fit$lrNaive)))  
  ncvar_put(nc, lr1_var, as.vector(as.matrix(fit$lr1)))
  ncvar_put(nc, lr2_var, as.vector(as.matrix(fit$lr2)))
  ncvar_put(nc, t_var, as.vector(as.matrix(fit$t)))
  ncvar_put(nc, se_var, as.vector(as.matrix(fit$se)))
  ncvar_put(nc, sdNaive_var, as.vector(as.matrix(fit$sdNaive)))
  ncvar_put(nc, nEff_var, as.vector(as.matrix(fit$nEff)))
  ncvar_put(nc, tot_var, as.vector(as.matrix(fit$tot)))
  ncvar_put(nc, tot0_var, as.vector(as.matrix(fit$tot0)))
  ncvar_put(nc, tot1_var, as.vector(as.matrix(fit$tot1)))
  ncvar_put(nc, tot2_var, as.vector(as.matrix(fit$tot2)))
  ncvar_put(nc, tot1_0_var, as.vector(as.matrix(fit$tot1_0)))
  ncvar_put(nc, tot2_0_var, as.vector(as.matrix(fit$tot2_0)))
  ncvar_put(nc, posfit_var, as.vector(as.matrix(posfit[,1:3])))
  ncvar_put(nc, posfit_lrn_var, as.vector(as.matrix(posfit[,-(1:3)])))
  ncvar_put(nc, posfit_clust_order_var, as.vector(posfitClust$order))
  nc_close(nc)
  ncdf4_test(nc_path, genes, fit, expsUsed, posfit, posfitClust)
  remove("fit")
  end = proc.time()
  print(end - start)
}

# Checks that the difference between x and y don't exceed Max_diff
# (global variable defined in the beginning of file)
similar_values <- function (x, y){return(abs(x-y)<=Max_diff)}

# Checks that the dimensions of a specified variable in a specified
# NetCDF file matches the expected dimensions, without taking order
# into consideration. Returns TRUE if the dimensions match and returns
# FALSE otherwise.
# expected: a vector of expected dimension names
# nc_file: NetCDF file of interest
# var_name: a NetCDF variable name of interest
correct_dimensions <- function(expected, nc_file, var_name){
  actual = dimensionNamesForVariable(nc_file, var_name)
  if (length(expected)!= length(actual)){
    return (FALSE)
  }else{
    for (d in actual){
      if (!is.element(d, expected)){
        return (FALSE)
      }
    }
    return (TRUE)
  }
}

# Checks whether a specified vector in an R image and that in a NetCDF 
# file have same/similar values. When the vector contains numbers,
# the function checks that the difference betwen corresponding values
# is at most Max_diff. When the vector contains a logical value or
# a string, the function checks that the corresponding values are 
# identical. Returns TRUE if the vectors have same/similar values
# and returns FALSE otherwise. 
# genes, fit, expsUsed, posfit, posfitClust: values stored in fit.image
#                                            and posfit_setup.image
# nc_file: NetCDF file of interest
# image_var_name: variable name in the R image. string.
# image_col: the column/index of interest in the R image variable. If
#            image_col is NULL, then the entire data table/vector is 
#            considered. string unless NULL.
# nc_var_name: variable name in the NetCDF file. string.
# nc_col: the column/index of interest in the NetCDF variable. If nc_col
#         is NULL, then the entire data table/vector is considered. 
#         string unless NULL.
# charmode: whether the nc_file stores the variable of interest as a
#           character vector
correct_values <- function(genes, fit, expsUsed, posfit, posfitClust, nc_file, image_var_name, image_col, nc_var_name, nc_col, charmode=FALSE){
  eval(parse(text = sprintf("expected = %s", image_var_name)))
  if (!is.null(image_col)){
    if (is.vector(expected)){
      eval(parse(text = sprintf("expected = expected[%s]",image_col)))
    }else{
      eval(parse(text = sprintf("expected = expected[,%s]",image_col)))
    }
  }
  actual = ncvar_get(nc_file,nc_var_name)
  if (!is.null(nc_col)){
    if (is.vector(actual)){
      eval(parse(sprintf(text = "actual = actual[%s]", nc_col)))
    }
    eval(parse(sprintf(text = "actual = actual[,%s]", nc_col)))
  }
  if (is.data.frame(expected) && ((ncol(expected)!=ncol(actual) || nrow(expected)!=nrow(actual)))){
    return (FALSE)
  }
  if (!is.data.frame(expected) && (length(expected) != length(actual))){
    return (FALSE)
  }
  if (!is.data.frame(expected)){
    return (!((charmode && !all(as.vector(na.omit(type.convert(actual)))==as.vector(na.omit(expected)))) || (!charmode && !all(as.vector(na.omit(actual))==as.vector(na.omit(expected))))))
  }
  for (i in 1:ncol(expected)){
    if (charmode){
      actual_converted = na.omit(type.convert(actual[,i]))
      expected_converted = as.vector(na.omit(expected[,i]))
      if (is.character(expected_converted)){
        expected_converted = type.convert(expected_converted)
      }
      if (is.numeric(actual_converted) != is.numeric(expected_converted)){
        return (FALSE)
      }
      if (is.numeric(actual_converted) && !all(similar_values(actual_converted, expected_converted))){
        return (FALSE)
      }
      if (!is.numeric(actual_converted) && !all(actual_converted == expected_converted)){
        return (FALSE)
      }
    }else if (!charmode && !all(similar_values(as.vector(na.omit(actual[,i])),as.vector(na.omit(expected[,i]))))){
      return (FALSE)
    }
  }
  return (TRUE)
}

# Runs correct_values and correct_dimensions tests on specified 
# NetCDF file and R image (genes, fit, expsUsed, posfit, posfitClust)
ncdf4_test <- function(nc_file_path, genes, fit, expsUsed, posfit, posfitClust){
 f = nc_open(nc_file_path)
 dimensions_lst = list(list("genes", c("string_axis", "locusId_axis", "gene_fieldname_axis")),
                   list("expsUsed", c("string_axis", "expsId_axis", "exps_fieldname_axis")),
                   list("fit/lr", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lrn", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lr1", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lr2", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lrn1", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lrn2", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lrNaive", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lrRaw", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/sumsq", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/sd", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/sdNaive", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/se", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/n", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/t", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/nEff", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/tot", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/tot0", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/tot1", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/tot2", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/tot1_0", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/tot2_0", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/lr", c("pergene_locusId_axis", "setId_axis")),
                   list("fit/q", c("string_axis", "setId_axis", "exps_metadata_fieldname_axis")),
                   list("fit/strain_lr", c("strainId_axis", "setId_axis")),
                   list("fit/strain_lrn", c("strainId_axis", "setId_axis")),
                   list("fit/strain_se", c("strainId_axis", "setId_axis")),
                   list("fit/strains", c("string_axis", "strainId_axis", "strain_fieldname_axis")),
                   list("posfit/pos", c("string_axis", "posId_axis", "posfit_fieldname_axis")),
                   list("posfit/lrn", c("posId_axis", "posfit_set_axis")),
                   list("posfit/clust_order", c("posfit_set_axis")))
 values_lst = list(list("genes", "-(which(colnames(genes)=='locusId'))", "genes", NULL, TRUE),
                       list("expsUsed", NULL, "expsUsed", NULL, TRUE),
                       list("fit$lrn1", NULL, "fit/lrn1", NULL, FALSE),
                       list("fit$lrn2", NULL, "fit/lrn2", NULL, FALSE),
                       list("fit$lrNaive", NULL, "fit/lrNaive", NULL, FALSE),
                       list("fit$sumsq", NULL, "fit/sumsq", NULL, FALSE),
                       list("fit$sd", NULL, "fit/sd", NULL, FALSE),
                       list("fit$lrRaw", NULL, "fit/lrRaw", NULL, FALSE),
                       list("fit$strain_lr", NULL, "fit/strain_lr", NULL, FALSE),
                       list("fit$n", NULL, "fit/n", NULL, FALSE),
                       list("fit$strains", NULL, "fit/strains", NULL, TRUE),
                       list("fit$strain_lrn", NULL, "fit/strain_lrn", NULL, FALSE),
                       list("fit$strain_se", NULL, "fit/strain_se", NULL, FALSE),
                       list("fit$lr", NULL, "fit/lr", NULL, FALSE),
                       list("fit$t", NULL, "fit/t", NULL, FALSE),
                       list("fit$nEff", NULL, "fit/nEff", NULL, FALSE),
                       list("fit$tot1", NULL, "fit/tot1", NULL, FALSE),
                       list("fit$tot2_0", NULL, "fit/tot2_0", NULL, FALSE),
                       list("posfit", "-1:-3", "posfit/lrn", NULL, FALSE),
                       list("fit$q", "-1", "fit/q", NULL, TRUE),
                       list("fit$lr1", NULL, "fit/lr1", NULL, FALSE),
                       list("fit$se", NULL, "fit/se", NULL, FALSE),
                       list("fit$tot", NULL, "fit/tot", NULL, FALSE),
                       list("fit$tot", NULL, "fit/tot", NULL, FALSE),
                       list("fit$tot2", NULL, "fit/tot2", NULL, FALSE),
                       list("posfitClust$order", NULL, "posfit/clust_order", NULL, FALSE),
                       list("fit$lrn", NULL, "fit/lrn", NULL, FALSE),
                       list("fit$lr2", NULL, "fit/lr2", NULL, FALSE),
                       list("fit$sdNaive", NULL, "fit/sdNaive", NULL, FALSE),
                       list("fit$tot0", NULL, "fit/tot0", NULL, FALSE),
                       list("fit$tot1_0", NULL, "fit/tot1_0", NULL, FALSE),
                       list("posfit", "1:3", "posfit/pos", NULL, TRUE),
                       list("fit$g", NULL, "pergene_locusId_axis", NULL, TRUE),
                       list("genes", "which(colnames(genes)=='locusId')", "locusId_axis", NULL, TRUE),
                       list("colnames(fit$strains)", NULL, "strain_fieldname_axis", NULL, TRUE),
                       list("colnames(expsUsed)", NULL, "exps_fieldname_axis", NULL, TRUE),
                       list("colnames(fit$lr)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lrn)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lr1)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lr2)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lrn1)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lrn2)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lrNaive)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$lrRaw)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$sumsq)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$sd)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$sdNaive)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$se)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$n)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$t)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$nEff)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$tot)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$tot0)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$tot1)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$tot2)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$tot1_0)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$tot2_0)", NULL, "setId_axis", NULL, TRUE),
                       list("fit$q$name", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$strain_lr)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$strain_lrn)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(fit$strain_se)", NULL, "setId_axis", NULL, TRUE),
                       list("colnames(posfit)", "-1:-3", "posfit_set_axis", NULL, TRUE),
                       list("colnames(posfit)", "1:3", "posfit_fieldname_axis", NULL, TRUE),
                       list("colnames(fit$q)", "-1", "exps_metadata_fieldname_axis", NULL, TRUE),
                       list("colnames(genes)", "-(which(colnames(genes)=='locusId'))", "gene_fieldname_axis", NULL, TRUE))
 for (var_dim in dimensions_lst){
   if (!correct_dimensions(var_dim[[2]], f, var_dim[[1]])){
     print(sprintf("Incorrect dimensions for %s. Expected %s, but received %s", var_dim[[1]], paste(var_dim[[2]], collapse = " "), paste(actual, collapse=" ")))
   }
 }
 for (val in values_lst){
   if (!correct_values(genes, fit, expsUsed, posfit, posfitClust, f, val[[1]], val[[2]], val[[3]], val[[4]], val[[5]])){
     print(sprintf("Values of %s in the original R image and %s in the netCDF file don't match", val[[1]], val[[3]]))
   }
 }
 nc_close(f)
}

