/*
** svn $Id: kelvin.h,v 1.2 2009/10/02 00:34:47 kana Exp $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Kelvin Wave Test.
**
** Application flag:   KELVIN
** Input script:       ocean_kelvin.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_DIF2
#define MIX_S_TS
#define SOLVE3D
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define RADIATION_2D
#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#undef  WEST_FSCLAMPED
#undef  WEST_M2CLAMPED
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SRFLUX
#define ANA_BTFLUX
