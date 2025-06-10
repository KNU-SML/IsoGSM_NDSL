
/* dependencies */

#define CLDTABLE
#ifdef CLDSLINGO
#undef CLDTABLE
#endif

#define RADFS
#ifdef RADMDC
#undef RADFS
#endif

#if defined CLDT93 || defined T93PR2 || defined ZM2
#define CCMCOM
#endif

#ifdef SCM
#undef CLDSLINGO    /* Slingo cloud scheme */
#define CLDTABLE
#define RADFS
#define INSA_CLD    /* Insa's marine stratus adjustment to Slingo scheme */
#define  SAS        /* Simplified Arakawa Schubert scheme */
#undef   RAS        /* Relaxed Arakawa Schubert scheme */
#undef   RASC1      /* RAS with detarained cloud water */
#undef   RASC2      /* RAS with detarained cloud water like RASV2 */
#undef  NOALSM1     /* NOAH land scheme */
#undef  NOAHYDRO    /* Residual free hydrology with noah */
#define R2          /* simulate Reanalysis-2 physics */
#undef  MRG_POST           /* merge postprocessor into forecast code */
#undef  PGB_PARALLEL       /* pgb parallel io and parallel computation */
#define PGB_SEMIPARALLEL   /* pgb parallel computation but serial io */
#define DELTAT 600  /*time step for scm */
#endif

#ifndef RSM
#if ( _jcap_ > 42)
#define REDUCE_GRID
#endif
#endif

#ifdef OSULSM1
#define _lsoil_ 2
#define _nsoil_ 0
#define _msub_ 0
#define SFCFTYP osu1
#define _lalbd_ 1
#endif

#ifdef OSULSM2
#define _lsoil_ 2
#define _nsoil_ 0
#define _msub_ 0
#define SFCFTYP osu2
#define _lalbd_ 4
#endif

#ifdef NOALSM1
#define _lsoil_ 4
#define _nsoil_ 0
#define _msub_ 0
#define SFCFTYP noa1
#define _lalbd_ 4
#endif

#ifdef VICLSM1
#define _lsoil_ 3
#define _nsoil_ 5
#define _msub_ 2
#define SFCFTYP vic1
#define _lalbd_ 4
#endif

#ifdef VICLSM2
#define _lsoil_ 3
#define _nsoil_ 5
#define _msub_ 12
#define SFCFTYP vic1
#define _lalbd_ 4
#endif

#ifdef R2
#define NOALSM2
#undef RAD_SMOOTH_CLOUD
#undef CLDADJ
#undef RADMDC
#undef SW99UPDATE
#undef CLDSLINGO
#define SAS
#undef  RAS
#undef RASV2
#undef CCMCNV
#undef CCMSCV
#define SFC
#define SFCMRG
#endif

#ifdef RSM
#undef REDUCE_GRID
#endif

#ifdef MRG_POST
#ifndef MP
#undef MPIGRIB
#else
#ifdef PGB_PARALLEL
#undef PGB_SEMIPARALLEL
#undef PGB_NOPARALLEL
#undef MPIGRIB
#undef PGBGATHER
#endif
#ifdef PGB_SEMIPARALLEL
#undef PGB_PARALLEL
#undef PGB_NOPARALLEL
#undef MPIGRIB
#define PGBGATHER
#endif
#ifdef PGB_NOPARALLEL
#undef PGB_PARALLEL
#undef PGB_SEMIPARALLEL
#define MPIGRIB
#undef PGBGATHER
#endif
#endif
#else
#undef PGB_NOPARALLEL
#undef PGB_PARALLEL
#undef PGB_SEMIPARALLEL
#undef MPIGRIB
#undef PGBGATHER
#endif

#ifdef NOALSM1
#undef HYDRO
#endif

#ifdef GDAS
#undef SFCMRG
#define SKIPSFCMRG
#endif

#ifdef OSULSM1
#undef OSULSM2
#undef NOALSM1
#undef VICLSM1
#undef MRGLSM
#undef USGS_SFC
#undef NOAHYDRO
#endif

#ifdef OSULSM2
#undef OSULSM1
#undef NOALSM1
#undef VICLSM1
#undef MRGLSM
#undef NOAHYDRO
#endif

#ifdef NOALSM1
#undef OSULSM2
#define MRGLSM
#undef USGS_SFC
#undef OSULSM1
#undef VICLSM1
#endif

#undef OSU
#ifdef OSULSM1
#define OSU
#endif
#ifdef OSULSM2
#define OSU
#endif

#undef NOA
#ifdef NOALSM1
#define NOA
#endif
#ifdef NOALSM2
#define NOA
#endif

#undef VIC
#ifdef VICLSM1
#define VIC
#endif
#ifdef VICLSM2
#define VIC
#endif

#ifdef VICLSM1
#define MRGLSM
#undef OSULSM1
#undef OSULSM2
#undef NOALSM1
#undef USGS_SFC
#undef NOAHYDRO
#endif

#ifdef USGS_SFC
#define _nstype_ 16
#define _nvtype_ 12
#define USGS
#define OSULSM2
#define STATSGO_SOIL
#define NCAR_EDIR
#else
#define _nstype_ 9
#define _nvtype_ 13
#undef USGS
#undef NCAR_EDIR
#undef STATSGO_SOIL
#endif

/* tracer */
/* number of tracers.  1 is moisture. 2 is o3 */
#ifdef CLDT93
#define _ntrac_ 3
#else
#ifdef ISOTOPE
#define _ntrac_ 4
#else
#define _ntrac_ 1
#endif
#endif

#ifdef CLD1
#define _ncldg_ 0
#undef RASC1
#undef RASC2
#undef ICE
#undef ICECLOUD
#undef INTERACTIVE_STRATUS
#endif

#ifdef CLD2
#define _ncldg_ 1
#endif

#ifdef CLDT93
#define _ncldg_ 1
#undef  INTERACTIVE_STRATUS
#endif

#ifdef CLDWSMMPS
#define _ncldg_ 2
#endif

#ifdef CLD3
#define _ncldg_ 2
#endif

#ifdef CLD5
#define _ncldg_ 4
#endif

#ifdef CLD6
#define _ncldg_ 5
#endif

#ifndef CLDSLINGO
#undef INSA_CLD
#endif

#ifdef CO2PROG
#define _ntrac_ 3
#define _ncldg_ 0
#endif

/* from Microphysics routines */
#undef  HD1990
#undef  NEW_PSACI
#undef  DEN_SNOW_OPT
                                                                                
#ifdef  NEW_MPHYS
#undef  FLETCHER
#define NEW_ICE_GEN
#define NEW_ICE_AUTO
#define NEW_ICE_ACR
#undef  KESSLER
#undef  DEN_SNOW_WATER
#undef  N0S_CONST
#else
#define FLETCHER
#undef  NEW_ICE_GEN
#undef  NEW_ICE_AUTO
#undef  NEW_ICE_ACR
#define KESSLER
#undef  DEN_SNOW_WATER
#define N0S_CONST
#endif

#undef ALB4
#ifndef OSULSM1
#define ALB4
#endif

#ifdef KIMGWD
#undef GTOPO30
#define _mtnvar_ 10
#else
#define _mtnvar_ 1
#endif

#ifndef RSM
#undef PSPLIT
#endif

#ifndef RSM
#undef GTOPO30
#endif

#ifndef RSM
#undef GET
#undef VECSUM
#undef RSMVECTORIZE
#endif

#ifndef MP
#undef GET
#undef VECSUM
#undef RSMVECTORIZE
#endif

#if ( _igrd_ > 256 )
#define VERY_LARGE_DOMAIN
#endif

#ifdef VERY_LARGE_DOMAIN
#define RINPG_SMPL
#endif

#ifdef DG3
#define DG
#define CLR
#endif

#ifdef DGZ
#define DG
#endif

#ifdef KEN
#undef CLR
#endif

#ifndef MP
#undef RMPGP2F_GLB
#endif

#ifdef ES
#define RESET_FHOUR
#endif

#ifdef SAS
#ifdef SAS2005
#undef  SAS_DIAG
#define  NEW_TRIGGER
#undef  UV_DIFF
#define  UV_TRAN
#define  MUL_CLDTOP
#define _ncldtop_ 1
#define  CLM_CWF
#else
#undef  SAS_DIAG
#undef  NEW_TRIGGER
#undef  UV_DIFF
#undef  UV_TRAN
#undef  MUL_CLDTOP
#define _ncldtop_ 1
#undef  CLM_CWF
#endif
#endif


#ifdef SCM
#define SKIPSFCMRG
#undef PSPLIT
#undef SFCMRG
#endif

#ifdef CONDENS
#define DG3
#define DG
#undef  R_CNVCLD
#undef  CWC_CNV
#endif

#ifdef DG3
#define DG
#define CLR
#endif

#ifdef RIV05
#define RIVER
#define _io2_ 720
#define _jo2_ 361
#endif
#ifdef RIV1
#define RIVER
#define _io2_ 360
#define _jo2_ 181
#endif

#ifdef SPNG_WID
#undef SPNG_NRW
#endif

#ifdef SPNG_NRW
#undef SPNG_WID
#define STDAMP 
#define STDAMP3 
#endif

#ifdef OSULSM2
#ifdef RSM
#define VEG_RES 20m
#else
#define VEG_RES 40m
#endif
#else
#define VEG_RES
#endif

#ifndef PGB_PARALLEL
#undef CATPGB
#else
#define CATPGB
#endif

#ifdef PGB_SEMIPARALLEL
#undef DOPOST
#else
#define DOPOST
#endif

#ifdef CATPGB
#undef DOPOST
#endif



#if ( _jcap_ == 20 )
#define DELTAT 4800
#endif
#if ( _jcap_ == 40 )
#define DELTAT 2400
#endif
#if ( _jcap_ == 42 )
#define DELTAT 2400
#endif
#if ( _jcap_ == 62 )
#define DELTAT 1800
#endif
#if ( _jcap_ == 126 )
#define DELTAT 900
#endif
#if ( _jcap_ == 170 )
#define DELTAT 600
#endif
#if ( _jcap_ == 248 )
#define DELTAT 450
#endif
#if ( _jcap_ == 372 )
#define DELTAT 300
#endif
#if ( _jcap_ == 496 )
#define DELTAT 200
#endif
#if ( _jcap_ == 646 )
#define DELTAT 120
#endif

#ifndef MP
#undef GSMDAMP
#endif

