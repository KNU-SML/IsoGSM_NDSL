DSET   ^rsm2d-6h.dat
OPTIONS BIG_ENDIAN 
TITLE  regrid_rcm2d-3h_snp    SNP ( Snapshot data )
UNDEF -9.99E33
XDEF  @RSM_IGRD@  LINEAR  224.6853   0.1098
YDEF  @RSM_JGRD@  LINEAR  19.50622   0.1034
ZDEF    1  LEVELS  1000.000
TDEF  @NHM_TDIM@ LINEAR  @NHM_HH@Z@NHM_DD@@NHM_MMC@@NHM_YYYY@   @NHM_INCHOUR@HR
VARS      23
PS         0  99  SURFACE  PRESSURE                                Pa               1
TSEAS      0  99  SST (SST: given from external file ICE: )        K                1
AICE       0  99  Area fraction of sea ice                         %                1
YICE       0  99  Mass of Sea Ice                                  kg/m**2          1
YSNW       0  99  Mass of snow on sea ice                          kg/m**2          1
UA         0  99  Surface Zonal Velocity at 10m                    m/s              1
VA         0  99  Surface Merid. Velocity at 10m                   m/s              1
TA         0  99  Surface Air Temperature at 2m                    K                1
QA         0  99  Surface Air Specific Humidity at 2m              kg/kg            1
TMPCNP     0  99  Canopy temperature                               K                0
TGEF       0  99  Effective Ground temperature (Radiation)         K                1
TMPSL1     0  99  TMP SOIL L1                                      K                0
TMPSL2     0  99  TMP SOIL L2                                      K                0
TMPSL3     0  99  TMP SOIL L3                                      K                0
TMPSL4     0  99  TMP SOIL L4                                      K                0
TSEAI      0  99  Temperature on tne sea ice surface               K                1
Z0         0  99  Roughness Length                                 M                1
WETSL1     0  99  fraction of H2O in SOIL L1                       0-1              0
WETSL2     0  99  fraction of H2O in SOIL L2                       0-1              0
WETSL3     0  99  fraction of H2O in SOIL L3                       0-1              0
SWE        0  99  Snow water equivalent                            kg/m**2          1
DEPSNW     0  99  Snow Depth * CVRSNW                              M                0
CVRSNW     0  99  Snow Coverage                                    0-1              0
ENDVARS
