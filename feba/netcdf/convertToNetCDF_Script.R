library("ncdf4")

source("convertToNetCDF_Functions.R")
source("NetCDFInterface_Functions.R")

TODO = list(c("psRCH2/fit.image", "psRCH2/psRCH2.nc", "psRCH2/posfit_setup.image"), 
            c("MR1/fit.image", "MR1/MR1.nc", "MR1/posfit_setup.image"), 
            c("ANA3/fit.image", "ANA3/ANA3.nc", "ANA3/posfit_setup.image"), 
            c("SB2B/fit.image", "SB2B/SB2B.nc", "SB2B/posfit_setup.image"), 
            c("PV4/fit.image", "PV4/PV4.nc", "PV4/posfit_setup.image"),
            c("Dino/fit.image", "Dino/Dino.nc", "Dino/posfit_setup.image"),
            c("Kang/fit.image", "Kang/Kang.nc", "Kang/posfit_setup.image"),            
            c("Keio/fit.image", "Keio/Keio.nc", "Keio/posfit_setup.image"),
            c("Korea/fit.image", "Korea/Korea.nc", "Korea/posfit_setup.image"),
            c("Marino/fit.image", "Marino/Marino.nc", "Marino/posfit_setup.image"),
            c("Miya/fit.image", "Miya/Miya.nc", "Miya/posfit_setup.image"),
            c("Phaeo/fit.image", "Phaeo/Phaeo.nc", "Phaeo/posfit_setup.image"),
            c("PS/fit.image", "PS/PS.nc", "PS/posfit_setup.image"),
            c("pseudo1_N1B4/fit.image", "pseudo1_N1B4/pseudo1_N1B4.nc", "pseudo1_N1B4/posfit_setup.image"),
            c("pseudo3_N2E3/fit.image", "pseudo3_N2E3/pseudo3_N2E3.nc", "pseudo3_N2E3/posfit_setup.image"),
            c("pseudo5_N2C3_1/fit.image", "pseudo5_N2C3_1/pseudo5_N2C3_1.nc", "pseudo5_N2C3_1/posfit_setup.image"),
            c("SynE/fit.image", "SynE/SynE.nc", "SynE/posfit_setup.image")
            ) 

for (f in TODO){
  print(f[1])
  make_ncdf4(1,4, f[1], f[2], f[3])
}
