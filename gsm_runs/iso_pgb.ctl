dset ^pgb.ft%f2
index ^pgb.ft00.idx
undef 9.999E+20
title pgb.ft00
*  produced by grib2ctl v0.9.13
* command line options: pgb.ft00
dtype grib 2
options yrev template
ydef 73 linear -90.000000 2.5
xdef 144 linear 0.000000 2.500000
tdef 10000 linear 00Z01jan2004 24hr
*  z has 17 levels, for prs
zdef 17 levels
1000 925 850 700 600 500 400 300 250 200 150 100 70 50 30 20 10
vars 29
no4LFTXsfc  0 132,1,0  ** surface Best (4-layer) lifted index [K]
ABSVprs 17 41,100,0 ** (profile) Absolute vorticity [/s]
TR3V       0 77,200,0 HDO column amount (kg/m**2/SMOW)
TR2     17 153,100,0 H218O mixing ratio (kg/kg/SMOW)
TR3      17 152,100,0 HDO mixing ratio (kg/kg/SMOW)
TR2V       0 76,200,0 H218O column amount (kg/m**2/SMOW)
HGTsfc  0 7,1,0  ** surface Geopotential height [gpm]
HGTprs 17 7,100,0 ** (profile) Geopotential height [gpm]
LFTXsfc  0 131,1,0  ** surface Surface lifted index [K]
PRESsfc  0 1,1,0  ** surface Pressure [Pa]
PRESmwl  0 1,6,0 ** max wind level Pressure [Pa]
PREStrp  0 1,7,0 ** tropopause Pressure [Pa]
PRMSLmsl  0 2,102,0 ** mean-sea level Pressure reduced to MSL [Pa]
PTENDsfc  0 3,1,0  ** surface Pressure tendency [Pa/s]
PWAT       0 54,200,0 Precipitable water (kg/m**2)
RHprs 17 52,100,0 ** (profile) Relative humidity [%]
RHclm  0 52,200,0 ** atmos column Relative humidity [%]
SPFH      17 51,100,0 Specific humidity (kg/kg)
TMPprs 17 11,100,0 ** (profile) Temp. [K]
TMPmwl  0 11,6,0 ** max wind level Temp. [K]
TMPtrp  0 11,7,0 ** tropopause Temp. [K]
UGRDprs 17 33,100,0 ** (profile) u wind [m/s]
UGRDmwl  0 33,6,0 ** max wind level u wind [m/s]
UGRDtrp  0 33,7,0 ** tropopause u wind [m/s]
VGRDprs 17 34,100,0 ** (profile) v wind [m/s]
VGRDmwl  0 34,6,0 ** max wind level v wind [m/s]
VGRDtrp  0 34,7,0 ** tropopause v wind [m/s]
VVELprs 17 39,100,0 ** (profile) Pressure vertical velocity [Pa/s]
VWSHtrp  0 136,7,0 ** tropopause Vertical speed shear [1/s]
ENDVARS

