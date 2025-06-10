DSET %EXP_DIR%/bdry.nc
DTYPE netcdf
TITLE ROMS grid file
UNDEF -1.e+34
XDEF 128 linear -91.868 0.19
YDEF 238 linear -42.842 0.14
ZDEF 30 linear 0 30 1
TDEF 15 linear 15jan2000 1dy
VARS 8
ubar_north  1  t,x  h
ubar_south  1  t,x  h
ubar_east   1  t,y  h
ubar_west   1  t,y  h
u_north    30  t,z,x  h
u_south    30  t,z,x  h
u_east     30  t,z,y  h
u_west     30  t,z,y  h
ENDVARS
