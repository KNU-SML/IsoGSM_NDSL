DSET %EXP_DIR%/bdry.nc
DTYPE netcdf
TITLE ROMS grid file
UNDEF -1.e+34
XDEF 129 linear -91.963 0.19
YDEF 238 linear -42.842 0.14
ZDEF 30 linear 0 30 1
TDEF 12 linear 15jan2000 1dy
VARS 12
zeta_south     1  t,x  hraw
zeta_north     1  t,x  hraw
zeta_west      1  t,y  hraw
zeta_east      1  t,y  hraw
temp_south    30  t,z,x  hraw
temp_north    30  t,z,x  hraw
temp_west     30  t,z,y  hraw
temp_east     30  t,z,y  hraw
salt_south    30  t,z,x  hraw
salt_north    30  t,z,x  hraw
salt_west     30  t,z,y  hraw
salt_east     30  t,z,y  hraw
ENDVARS
