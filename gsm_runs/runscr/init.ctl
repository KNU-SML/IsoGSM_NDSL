DSET %EXP_DIR%/init.nc
DTYPE netcdf
TITLE ROMS grid file
UNDEF -1.e+34
XDEF 129 linear -91.963 0.19
YDEF 238 linear -42.842 0.14
ZDEF 30 linear 0 30 1
TDEF 1 linear 15jan2000 1dy
VARS 4
zeta     1  t,y,x  hraw
temp     30  t,z,y,x  h
salt     30  t,z,y,x  f
mask_rho 30  t,z,y,x  mask_rho
ENDVARS
