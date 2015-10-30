/*
** svn $Id: icetest.h 1307 2008-01-10 00:22:36Z bgalton $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Simplified Ice Shelf Ocean Cavity Model Test.
**
** Application flag:   WAOM10
** Input script:       ocean_waom10.in
*/


#define UV_COR
#define UV_VIS2
#define UV_QDRAG

#define TS_C4HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define SALINITY
#define NONLIN_EOS

#define DJ_GRADPS

#define SOLVE3D
#define CURVGRID
#define SPHERICAL

#define SPLINES

#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SRFLUX

#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX

#define ANA_SEAICE

/*
#define BULK_FLUXES
#define LONGWAVE
#define EMINUSP
#define ALBEDO
#define ALBEDO_CLOUD
#define SNOWFALL
#define SNOW_FROM_RAIN
*/

#define MIX_S_UV
#define MIX_S_TS

#define LMD_MIXING

#ifdef LMD_MIXING
#define LMD_CONVEC
#endif

#undef SSH_TIDES
#ifdef SSH_TIDES
#define RAMP_TIDES
#define FSOBC_REDUCED
#define NODAL_TIDES
#endif

#define ICESHELF
#define ICESHELF_3EQN_VBC

#define PERFECT_RESTART
#define MASKING
#undef AVERAGES
#define WET_DRY
