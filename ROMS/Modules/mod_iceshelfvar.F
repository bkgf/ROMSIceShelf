#include "cppdefs.h"
      MODULE mod_iceshelfvar

/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2013 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module declares Ice Shelf Model Kernel Structure  containing  **
** several variables describing ice draft and thickness properties    **
**                                                                    **
** Note that the  *.h file is located in ROMS/Nonlinear/IceShelf and  **
** it is included within <...> to allow the user customize it in the  **
** project directory, while keeping the distributed file intact       **
** (check build scripts for details).                                 **
**                                                                    **
************************************************************************
*/

#if defined ICESHELF
# include <iceshelfvar_mod.h>
#endif

      END MODULE mod_iceshelfvar
