/*
** svn $Id: basin.h,v 1.2 2009/10/02 00:34:46 kana Exp $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Big Bad Basin.
**
** Application flag:   BASIN
** Input script:       ocean_basin.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS4
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLVE3D
#define SPLINES
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define BODYFORCE
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
