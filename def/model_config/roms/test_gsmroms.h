/*
** svn $Id: test_gsmroms.h,v 1.3 2009/11/05 19:48:11 haiqinli Exp $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for PERU CHILE (VOCALS) REGION 05/26/08
**
** Application flag:   PCR
** Input script:       ocean_pcr.in
*/

/*
**-------------------------------------------------------------
** Options for the Nonlinear, Tangent Linear and Adjoint Model.
**-------------------------------------------------------------
*/
#define NL_MODEL
#ifdef NL_MODEL
# undef ADJOINT
# undef TANGENT
#endif

/*
**-------------------------------------------------------------
** Options for Model Physics.
**-------------------------------------------------------------
*/
#define UV_ADV
#define UV_COR
#define UV_VIS2
#define UV_LDRAG
#undef  MIX_GEO_UV
#define MIX_S_UV

#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define MIX_S_TS

#define DJ_GRADPS
#define SALINITY
#undef QCORRECTION
#undef SCORRECTION
#define NONLIN_EOS
#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES

#define OUT_DOUBLE
#define AVERAGES
#define AVERAGES_FLUXES

/*
**-------------------------------------------------------------
** Options for MIXING in NL_MODEL.
**-------------------------------------------------------------
*/
#ifdef NL_MODEL
#define LMD_MIXING
#undef  MY25_MIXING
# ifdef LMD_MIXING
#  undef DIURNAL_SRFLUX
#  define SOLAR_SOURCE
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_NONLOCAL
#  define LMD_SKPP
# endif
#else
#undef  LMD_MIXING
#endif

/*
**-------------------------------------------------------------
** Options for Nudging to Climatology.
**-------------------------------------------------------------
*/
#undef CLIM_NUDGING
#ifdef CLIM_NUDGING
# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  ZCLIMATOLOGY
# undef  M2CLM_NUDGING
# undef  M3CLM_NUDGING
# undef  TCLM_NUDGING
# undef  ZCLM_NUDGING
#endif
/*
**-------------------------------------------------------------
** Options for Model Boundary Conditions.
**-------------------------------------------------------------
*/
#undef CLAMPED_BC
#ifdef CLAMPED_BC
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3CLAMPED
# define NORTH_TCLAMPED
# define SOUTH_FSCHAPMAN
# define SOUTH_M2FLATHER
# define SOUTH_M3CLAMPED
# define SOUTH_TCLAMPED
# define WEST_FSCHAPMAN
# define WEST_M2FLATHER
# define WEST_M3CLAMPED
# define WEST_TCLAMPED
#else
# define SPONGE
# define WEST_VOLCONS
# define SOUTH_VOLCONS
# define NORTH_VOLCONS
# define RADIATION_2D
# define WEST_FSGRADIENT
# define WEST_M2RADIATION
# define WEST_M2NUDGING
# define WEST_M3RADIATION
# define WEST_M3NUDGING
# define WEST_TRADIATION
# define WEST_TNUDGING
# define NORTH_FSGRADIENT
# define NORTH_M2RADIATION
# define NORTH_M2NUDGING
# define NORTH_M3RADIATION
# define NORTH_M3NUDGING
# define NORTH_TRADIATION
# define NORTH_TNUDGING
# define SOUTH_FSGRADIENT
# define SOUTH_M2RADIATION
# define SOUTH_M2NUDGING
# define SOUTH_M3RADIATION
# define SOUTH_M3NUDGING
# define SOUTH_TRADIATION
# define SOUTH_TNUDGING
#endif

/* Options for Bottom Flux.  */

#define ANA_BSFLUX
#define ANA_BTFLUX

/*
**-------------------------------------------------------------
** Options for Model Opening Boundary Conditions.
**-------------------------------------------------------------
*/
#undef CLOSED_OBC
#ifdef CLOSED_OBC
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EASTERN_WALL
#define WESTERN_WALL
#endif

