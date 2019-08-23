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
#if defined ICESHELF 
# if defined ICESHELF_3EQN_VBC
        real(r8), pointer :: gammaT(:,:)
        real(r8), pointer :: gammaS(:,:)
        real(r8), pointer :: Tb(:,:)
        real(r8), pointer :: Sb(:,:)
#  ifdef TANGENT
        real(r8), pointer :: tl_gammaT(:,:)
        real(r8), pointer :: tl_gammaS(:,:)
        real(r8), pointer :: tl_Tb(:,:)
        real(r8), pointer :: tl_Sb(:,:)
#  endif
#  ifdef ADJOINT
        real(r8), pointer :: ad_gammaT(:,:)
        real(r8), pointer :: ad_gammaS(:,:)
        real(r8), pointer :: ad_Tb(:,:)
        real(r8), pointer :: ad_Sb(:,:)
#  endif
# endif
        real(r8), pointer :: m(:,:)
# ifdef TANGENT
        real(r8), pointer :: tl_m(:,:)
# endif
# ifdef ADJOINT
        real(r8), pointer :: ad_m(:,:)
# endif
# if defined ICESHELF_MORPH
        real(r8), pointer :: iceshelf_draft0(:,:)
        real(r8), pointer :: iceshelf_draft(:,:,:)
#  ifdef TANGENT
        real(r8), pointer :: tl_iceshelf_draft0(:,:)
        real(r8), pointer :: tl_iceshelf_draft(:,:,:)
#  endif
#  ifdef ADJOINT
        real(r8), pointer :: ad_iceshelf_draft0(:,:)
        real(r8), pointer :: ad_iceshelf_draft(:,:,:)
#  endif
# endif
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
#if defined ICESHELF
# if defined ICESHELF_3EQN_VBC
      allocate ( ICESHELFVAR(ng) % gammaT(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % gammaS(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % Tb(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % Sb(LBi:UBi,LBj:UBj))
#  ifdef TANGENT
      allocate ( ICESHELFVAR(ng) % tl_gammaT(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % tl_gammaS(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % tl_Tb(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % tl_Sb(LBi:UBi,LBj:UBj))
#  endif
#  ifdef ADJOINT
      allocate ( ICESHELFVAR(ng) % ad_gammaT(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % ad_gammaS(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % ad_Tb(LBi:UBi,LBj:UBj))
      allocate ( ICESHELFVAR(ng) % ad_Sb(LBi:UBi,LBj:UBj))
#  endif
# endif
      allocate ( ICESHELFVAR(ng) % m(LBi:UBi,LBj:UBj))
# ifdef TANGENT
      allocate ( ICESHELFVAR(ng) % tl_m(LBi:UBi,LBj:UBj))
# endif
# ifdef ADJOINT
      allocate ( ICESHELFVAR(ng) % ad_m(LBi:UBi,LBj:UBj))
# endif
# if defined ICESHELF_MORPH
      allocate ( ICESHELFVAR(ng) % iceshelf_draft0(LBi:UBi,LBj:UBj) )
      allocate ( ICESHELFVAR(ng) % iceshelf_draft(LBi:UBi,LBj:UBj,1:2))
#  ifdef TANGENT
      allocate ( ICESHELFVAR(ng) % tl_iceshelf_draft0(LBi:UBi,LBj:UBj) )
      allocate ( ICESHELFVAR(ng) % tl_iceshelf_draft(LBi:UBi,LBj:UBj,1:2))
#  endif
#  ifdef ADJOINT
      allocate ( ICESHELFVAR(ng) % ad_iceshelf_draft0(LBi:UBi,LBj:UBj) )
      allocate ( ICESHELFVAR(ng) % ad_iceshelf_draft(LBi:UBi,LBj:UBj,1:2))
#  endif
# endif
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
#if defined ICESHELF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# if defined ICESHELF_3EQN_VBC 
            ICESHELFVAR(ng) % gammaT(i,j) = IniVal
            ICESHELFVAR(ng) % gammaS(i,j) = IniVal
            ICESHELFVAR(ng) % Tb(i,j) = IniVal
            ICESHELFVAR(ng) % Sb(i,j) = IniVal
# endif
            ICESHELFVAR(ng) % m(i,j) = IniVal
# if defined ICESHELF_MORPH
            ICESHELFVAR(ng) % iceshelf_draft0(i,j) = IniVal
            ICESHELFVAR(ng) % iceshelf_draft(i,j,1) = IniVal
            ICESHELFVAR(ng) % iceshelf_draft(i,j,2) = IniVal
# endif
          END DO
        END DO
#endif
      END IF
#ifdef TANGENT
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM)) THEN
#if defined ICESHELF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# if defined ICESHELF_3EQN_VBC 
            ICESHELFVAR(ng) % tl_gammaT(i,j) = IniVal
            ICESHELFVAR(ng) % tl_gammaS(i,j) = IniVal
            ICESHELFVAR(ng) % tl_Tb(i,j) = IniVal
            ICESHELFVAR(ng) % tl_Sb(i,j) = IniVal
# endif
            ICESHELFVAR(ng) % tl_m(i,j) = IniVal
# if defined ICESHELF_MORPH
            ICESHELFVAR(ng) % tl_iceshelf_draft0(i,j) = IniVal
            ICESHELFVAR(ng) % tl_iceshelf_draft(i,j,1) = IniVal
            ICESHELFVAR(ng) % tl_iceshelf_draft(i,j,2) = IniVal
# endif
          END DO
        END DO
#endif
      END IF
#endif
#ifdef ADJOINT
!
!  Adjoint linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM)) THEN
#if defined ICESHELF
        DO j=Jmin,Jmax
          DO i=Imin,Imax
# if defined ICESHELF_3EQN_VBC 
            ICESHELFVAR(ng) % ad_gammaT(i,j) = IniVal
            ICESHELFVAR(ng) % ad_gammaS(i,j) = IniVal
            ICESHELFVAR(ng) % ad_Tb(i,j) = IniVal
            ICESHELFVAR(ng) % ad_Sb(i,j) = IniVal
# endif
            ICESHELFVAR(ng) % ad_m(i,j) = IniVal
# if defined ICESHELF_MORPH
            ICESHELFVAR(ng) % ad_iceshelf_draft0(i,j) = IniVal
            ICESHELFVAR(ng) % ad_iceshelf_draft(i,j,1) = IniVal
            ICESHELFVAR(ng) % ad_iceshelf_draft(i,j,2) = IniVal
# endif
          END DO
        END DO
#endif
      END IF
#endif
      RETURN
      END SUBROUTINE initialize_iceshelfvar










