!
!svn $Id$
!================================================== Hernan G. Arango ===
! Copyright (c) 2002-2013 The ROMS/TOMS Group  Benjamin K. Galton-Fenzi!
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Ice Shelf Model Kernel Variables:                                   !
!                                                                      !
      USE mod_kinds
!
      implicit none
!
      TYPE T_ICESHELFVAR
!
!  Nonlinear model state.
!
#if defined ICESHELF && defined ICESHELF_MORPH
        real(r8), pointer :: iceshelf_draft0(:,:)
        real(r8), pointer :: iceshelf_draft(:,:,:)
#endif
      END TYPE T_ICESHELFVAR

      TYPE (T_ICESHELFVAR), allocatable :: ICESHELFVAR(:)
 
      CONTAINS

      SUBROUTINE allocate_iceshelfvar (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_iceshelf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( ICESHELFVAR(Ngrids) )
!
!  Nonlinear model state.
!
#if defined ICESHELF && defined ICESHELF_MORPH
      allocate ( ICESHELFVAR(ng) % iceshelf_draft0(LBi:UBi,LBj:UBj) )
      allocate ( ICESHELFVAR(ng) % iceshelf_draft(LBi:UBi,LBj:UBj,1:2))
#endif

      RETURN
      END SUBROUTINE allocate_iceshelfvar

      SUBROUTINE initialize_iceshelfvar (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_iceshelf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize iceshelf structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
#if defined ICESHELF && defined ICESHELF_MORPH
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            ICESHELFVAR(ng) % iceshelf_draft0(i,j) = IniVal
            ICESHELFVAR(ng) % iceshelf_draft(i,j,1) = IniVal
            ICESHELFVAR(ng) % iceshelf_draft(i,j,2) = IniVal
          END DO
        END DO
#endif
      END IF
      RETURN
      END SUBROUTINE initialize_iceshelfvar










