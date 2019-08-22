#include "cppdefs.h"
       MODULE ad_iceshelf_mod
#if defined ICESHELF && defined ADJOINT

!
!==================================================== Hernan G. Arango =
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!=========================================== Benjmamin K. Galton-Fenzi =
!                                                                      !
!  This is the main driver routine for the ice shelf portion of the    !
!  model.                 
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_iceshelf

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_iceshelf (ng, tile)

      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_forces
      USE mod_iceshelf
      USE ad_iceshelf_vbc_mod, only : ad_iceshelf_vbc
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile

      CALL ad_iceshelf_vbc (ng, tile)

      RETURN
      END SUBROUTINE ad_iceshelf
#endif
      END MODULE ad_iceshelf_mod
