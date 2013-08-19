/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2013 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Mixed Layer Deepening on a rectangular periodic grid with
** thermal stratification and wind stress.
**
** Application flag:   MIXED_LAYER
** Input script:       ocean_mixed_layer.in
*/

#define WRITE_GRID
#define OUT_DOUBLE
#define ANA_GRID
#define SALINITY
#define SOLVE3D
#define SPLINES
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define TS_U3HADVECTION
#define NORTHERN_WALL
#define SOUTHERN_WALL

#undef  MY25_MIXING
#define GLS_MIXING

#if defined GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define CRAIG_BANNER
# define CHARNOK
# undef  ZOS_HSIG
# undef  TKE_WAVEDISS
#endif

#define UV_LOGDRAG
