/*
** svn $Id: adria02.h,v 1.2 2009/10/02 00:34:46 kana Exp $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Adriatic Sea Application.
**
** Application flag:   ADRIA02
** Input script:       ocean_adria02.in
**                     floats_adria02.in
**                     sediment_adria02.in
**                     stations_adria02.in
*/

#define UV_ADV
#define UV_COR
#define UV_PSOURCE
#define DJ_GRADPS
#undef  TS_U3HADVECTION
#undef  TS_SVADVECTION
#define TS_MPDATA
#define TS_DIF2
#define MIX_GEO_TS
#define TS_PSOURCE
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define SPLINES
#define STATIONS
#define CURVGRID
#define FLOATS
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS

#undef NOSEDBBL
#ifdef NOSEDBBL
# undef SEDIMENT
# undef SUSPLOAD
# define ANA_SEDIMENT
# undef  ANA_WWAVE
# undef RIVER_SEDIMENT
#else
# define SEDIMENT
# define SUSPLOAD
# undef  ANA_SEDIMENT
# undef  ANA_WWAVE
# define RIVER_SEDIMENT
#endif

#undef  UV_LOGDRAG
#undef  MB_BBL
#undef  SG_BBL
#define SSW_BBL

#ifdef SG_BBL
# define SG_CALC_ZNOT
# undef  SG_LOGINT
#endif
#ifdef MB_BBL
# define MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif

#undef MY25_MIXING
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define CRAIG_BANNER
# define CHARNOK
#endif

#undef ANA_SRFLUX
#undef ALBEDO
#define DIURNAL_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BPFLUX
#define ANA_BTFLUX
#define ANA_SPFLUX

#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# undef SOLAR_SOURCE
# define ANA_RAIN
# undef COOL_SKIN
#endif

#define WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define RADIATION_2D
#define EAST_M3RADIATION
#define EAST_TRADIATION

#define RAMP_TIDES
#define SSH_TIDES
#ifdef SSH_TIDES
# define EAST_FSCHAPMAN
# define ANA_FSOBC
#else
# define EAST_FSGRADIENT
#endif
#define UV_TIDES
#ifdef UV_TIDES
# define EAST_M2FLATHER
# define ANA_M2OBC
#else
# define EAST_M2RADIATION
#endif
