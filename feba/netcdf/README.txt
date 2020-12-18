Author: Shiori Sagawa <ssagawa@berkeley.edu>

************
*** GOAL ***
************
- To convert the per-strain data from an R image into a NetCDF format
  file to improve the efficiency of data retrieval

********************
*** REQUIREMENTS ***
********************
- ncdf4 (a CRAN package)

******************
*** BACKGROUND ***
******************
- For documentation on the NetCDF data format, refer to http://www.unidata.ucar.edu/software/netcdf/docs/html_guide/index.html#user_guide
- For documentation on the CRAN package, ncdf4, refer to http://cran.r-project.org/web/packages/ncdf4/ncdf4.pdf

*************
*** FILES ***
*************
- convertToNetCDF_Functions.R: includes all functions used to...
	- convert an R image to NetCDF-4 format 
	- test the new NetCDF-4 file against the original R image
- convertTONetCDF_Script.R: a script that converts R images into NetCDF-4 format and tests the newly created NetCDF-4 file.
- NetCDFInterface_Functions.R: includes functions used to interact with netcdf objects. Those tasks include:
	- making heatmaps
	- getting a data table corresponding to certain genes
	- getting a data table corresponding to certain positions
	- getting information regarding a variable, dimension, or netcdf file
- strainHeatmapTest.r: tests heatmap-related functions
	- tests that heatmap functions for NetCDF files, written in NetCDFInterface_Functions.R, work by comparing the output of the function against the output of the corresponding function for R images (written by Morgan)

***************************
*** .NC FILE COMPONENTS ***
***************************
For definition of each variables, refer to ../lib/FEBA_template.html
- genes (variable) # locusIds are in gene_locusId
- expsUsed (variable)
- fit (group): holds fitness data
	- q (variable)
	- lr (variable)
	- lrn (variable) # gene axis is in pergene_locusId
	- lr1 (variable)
	- lr2 (variable)
	- lrn1 (variable)
	- lrn2 (variable)
	- lrNaive (variable)
	- lrRaw (variable)
	- sumsq (variable)
	- sd (variable)
	- sdNaive (variable)
	- se (variable)
	- n (variable)
	- t (variable)
	- nEff (variable)
	- tot (variable)
	- tot0 (variable)
	- tot1 (variableot2 (variable)
	- tot1_0 (variable)
	- tot2_0 (variable)
	- strains (variable)
	- strain_lr (variable)
	- strain_lrn (variable)
	- strain_se (variable)
- posfit
	- pos (variable)
	- lrn (variamble)
	- clust_order (variable)
