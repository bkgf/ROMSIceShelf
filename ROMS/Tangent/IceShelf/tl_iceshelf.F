#include "cppdefs.h"
       MODULE tl_iceshelf_mod
#if defined ICESHELF && defined TANGENT

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
      PUBLIC  :: tl_iceshelf

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_iceshelf (ng, tile)

      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_forces
      USE mod_iceshelf
      USE tl_iceshelf_vbc_mod, only : tl_iceshelf_vbc
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile

      CALL tl_iceshelf_vbc (ng, tile)

      RETURN
      END SUBROUTINE tl_iceshelf
#endif
      END MODULE tl_iceshelf_mod
