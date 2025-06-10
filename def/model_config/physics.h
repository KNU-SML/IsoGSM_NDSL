
/* major physical parameterization default options */

/* radiation */
#define RADMDC       /* M-D Chou long wave radiation */
                     /* undef falls back to GFDL long wave */
#define SWCHOU       /* M-D Chou short wave radiation */

/* cloudiness */
#define CLDSLINGO    /* Slingo cloud scheme */
#define INSA_CLD     /* Insa's cloudiness adjustment to Slingo scheme */
#undef  CLDADJ       /* empirical cloudiness adjustment */
#define INTERACTIVE_STRATUS  /* cloudiness based on cloud water and RH */
                             /* only applicable for predicted cloud water */

/* radiation cloud water interaction */
#undef  ICECLOUD     /* cloud water feed into radiation */

/* convective parameterization */
#undef  KUO          /* Kuo scheme */
#undef  SAS          /* Simplified Arakawa Schubert scheme */
#define  RAS          /* Relaxed Arakawa Schubert scheme */
#undef  RASC1        /* RAS with detarained cloud water */
#define  RASC2        /* RAS with detarained cloud water like RASV2 */
#undef  RASV2        /* RAS with downdraft */
#undef  CCMCNV       /* Zhang McFarlane scheme */
#undef  ZM2          /* modified Zhang McFarlane scheme */
#undef KF2          /* Kain Fritsch convective parameterization from WRF */
#undef  CCMSCV       /* Hack shallow convection scheme to be used with CCMCNV */

/* land scheme */
#undef  OSULSM2      /* OSU version 2 land scheme with USGS data */
#undef  HYDRO        /* Residual free hydrology with osu2 */
#undef  OSULSM1      /* OSU version 1 land scheme */
#define NOALSM1      /* NOAH land scheme */
#define NOAHYDRO     /* Residual free hydrology with noah */
#define NOAH_JUN05   /* updated version of NOAH by Massimo */
#undef  VICLSM1      /* VIC single tile scheme */
#undef  NOBSFLW      /* remove base flow from osu2 */
#define CRITSNOW 269.66  /* snow threshold temperature at 850 hPa */

/* orography */
#undef  GTOPO30      /* use original gtopo 30 data for mountain in RSM */
#define _gsmtnres_ 8
#ifdef RDELX
#if ( RDELX < 40000 )
#define _rsmtnres_ 0.5
#define GTOPO30
#else
#define _rsmtnres_ 4
#endif
#endif
#define _romtnres_ _rsmtnres_

/* cloud water prediction */
#define CLD1         /* diagnostic                    */
#undef  CLD2         /* Zhao and Carr  */
#undef  CLDT93       /* Somerville and Iacobellis */
#undef  CLDWSMMPS    /* Updated version of CLD3 */
#undef  CLD3         /* Song-You Hong prognostic qc/qi, qr/qs       */
#undef  CLD5         /* Song-You Hong prognostic qc, qi, qr, qs     */
#undef  CLD6         /* Song-You Hong prognostic qc, qi, qr, qs, qg */
#define NEW_MPHYS    /* additional option for CLD3/CLD5/CLDWSMMPS */
#undef  VCI          /* additional option for CLD3/CLD5/CLDWSMMPS */
#undef  T93PR2       /* new precipitation physics for CLDT93 */

/* thermodynamics */
#define ICE          /* compute saturation on ice */

/*diffusion */
#define DIFFQ        /* horizontal diffusion of moisture */
#undef DIFUC	/* strong horizontal diffusion of cloud water variables to stabilize cloud water prediction schemes and allow for larger time step */

/* trace gasses */
#define O3C          /* o3 climatology */
#undef  CO2PPM       /* default co2 concentration of 348 ppm */
#undef  CH4PPM       /* default ch4 concentration of 1.75 ppm */
#undef  N2OPPM       /* default n2o concentration of 0.28 ppm */

/* miscellaneous */
#undef  R2           /* simulate Reanalysis-2 physics */

/* below are untested options */
#undef  KIMGWD       /* Arakawa Kim gravity wave drag */
#undef  SW99UPDATE   /* Short wave radiation update to Lacis and Hansen scheme */

/* roms ocean model vertical number of levels */
#define ROMS_LEVS 30

/* ocean model bathimetry file */
#define BATH_FILE etopo5.nc      /* low resolution bathimetry.  High resolution is etopo2.nc */
