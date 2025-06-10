/*
** svn $Id: test_rsmroms.h,v 1.7 2010/08/25 21:19:03 haiqinli Exp $
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

/* output options */
#undef  DIAGNOSTICS_UV      /* define if writing out momentum diagnostics */
#undef  DIAGNOSTICS_TS      /* define if writing out tracer diagnostics */
#define OUT_DOUBLE
#define AVERAGES
#define AVERAGES_FLUXES
/* advection and coriolis */
#define UV_ADV
#define UV_COR
#define TS_U3HADVECTION
#define TS_C4VADVECTION
 
/* pressure gradient */
#define DJ_GRADPS       /* Splines density  Jacobian (Shchepetkin, 2000) */
/* TIDAL FORCING */
#undef TIDAL
#ifdef TIDAL
#define SSH_TIDES       /* turn on computation of tidal elevation */
#define UV_TIDES        /* turn on computation of tidal currents */
#define ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC       /* Add tidal currents  to processed OBC data */
#endif

/* diffusion and viscosities schemes */
#define UV_VIS2
#define TS_DIF2
#define DIFF_GRID       /* Diffusivity coefficient scaled by grid size */
#define VISC_GRID       /* viscosity coefficient scaled by grid size */
#define MIX_GEO_TS
#define MIX_S_UV
 
#define UV_QDRAG        /* turn ON or OFF quadratic bottom friction */
 
/* vertical mixing schemes */
# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
 
/* equation of state - relaxation - salinity - coordinates - etc */
#define NONLIN_EOS
#define SOLAR_SOURCE

#define SALINITY

#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES
/* boundary conditions */
#define BRYFILE
#ifndef BRYFILE
#define M2CLIMATOLOGY
#define M3CLIMATOLOGY
#define TCLIMATOLOGY
#undef  ZCLIMATOLOGY
#define M2CLM_NUDGING
#define M3CLM_NUDGING
#define TCLM_NUDGING
#endif
#ifdef TIDAL
#define SOUTH_M2FLATHER
#define SOUTH_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_FSCHAPMAN
#undef  WEST_VOLCONS
#undef  SOUTH_VOLCONS
#endif

#define RADIATION_2D
#define SPONGE
/* Climats processing and relaxations */
#define QCORRECTION
#define SCORRECTION
#undef  DIURNAL_SRFLUX  /* impose shortwave radiation local diurnal cycle */
#define ANA_BSFLUX
#define ANA_BTFLUX

/*
**-------------------------------------------------------------
** Options for Model Opening Boundary Conditions.
**-------------------------------------------------------------
*/
#undef  CLOSED_OBC
#ifdef  CLOSED_OBC
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EASTERN_WALL
#define WESTERN_WALL
#endif

