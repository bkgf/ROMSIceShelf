       MODULE iceshelf_mod
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
      PUBLIC  :: iceshelf
      CONTAINS
!
!***********************************************************************
      SUBROUTINE iceshelf (ng, tile)
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_forces
      USE mod_iceshelf
      USE iceshelf_vbc_mod, only : iceshelf_vbc
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      CALL iceshelf_vbc (ng, tile)
      RETURN
      END SUBROUTINE iceshelf
      END MODULE iceshelf_mod
