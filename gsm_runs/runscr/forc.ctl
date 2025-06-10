DSET %EXP_DIR%/forc.nc
DTYPE netcdf
TITLE ROMS grid file
UNDEF 1.e+37
XDEF 129 linear -91.963 0.19
YDEF 238 linear -42.842 0.14
ZDEF 30 linear 0 30 1
TDEF 1 linear 15jan2000 1dy
VARS 6
shflux   1  t,y,x  hraw
swflux   1  t,y,x  h
SST      1  t,y,x  f
SSS      1  t,y,x  mask_rho
dQSST    1  t,y,x  dqsst
swrad    1  t,y,x  swrad
ENDVARS
