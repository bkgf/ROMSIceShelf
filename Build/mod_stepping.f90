      MODULE mod_stepping
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains time stepping indices.                         !
!                                                                      !
!  Lnew      New descent algorithm state solution index.               !
!  Lold      Previous descent algorithm state solution index.          !
!                                                                      !
!  knew      Barotropic (fast) time-step index corresponding to the    !
!              newest values for 2D primitive equation variables.      !
!  krhs      Barotropic (fast) time-step index used to compute the     !
!              right-hand-terms of 2D primitive equation variables.    !
!  kstp      Barotropic (fast) time-step index to which the current    !
!              changes are added to compute new 2D primitive equation  !
!              variables.                                              !
!                                                                      !
!  nfm3      Float index for time level "n-3".                         !
!  nfm2      Float index for time level "n-2".                         !
!  nfm1      Float index for time level "n-1".                         !
!  nf        Float index for time level "n".                           !
!  nfp1      Float index for time level "n+1".                         !
!                                                                      !
!  nnew      Baroclinic (slow) time-step index corresponding to the    !
!              newest values for 3D primitive equation variables.      !
!  nrhs      Baroclinic (slow) time-step index used to compute the     !
!              right-hand-terms of 3D primitive equation variables.    !
!  nstp      Baroclinic (slow) time-step index to which the current    !
!              changes are added to compute new 3D primitive equation  !
!              variables.                                              !
!                                                                      !
!  NTC       Number of tidal components to consider.                   !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
        integer, allocatable :: knew(:)
        integer, allocatable :: krhs(:)
        integer, allocatable :: kstp(:)
!
        integer, allocatable :: nnew(:)
        integer, allocatable :: nrhs(:)
        integer, allocatable :: nstp(:)
!
        integer, allocatable :: Lnew(:)
        integer, allocatable :: Lold(:)
        integer, allocatable :: NTC(:)
!
      CONTAINS
!
      SUBROUTINE allocate_stepping
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate and intialize time indices.
!-----------------------------------------------------------------------
!
      allocate ( knew(Ngrids) )
      knew(1:Ngrids)=1
      allocate ( krhs(Ngrids) )
      krhs(1:Ngrids)=1
      allocate ( kstp(Ngrids) )
      kstp(1:Ngrids)=1
      allocate ( nnew(Ngrids) )
      nnew(1:Ngrids)=1
      allocate ( nrhs(Ngrids) )
      nrhs(1:Ngrids)=1
      allocate ( nstp(Ngrids) )
      nstp(1:Ngrids)=1
      allocate ( Lnew(Ngrids) )
      Lnew(1:Ngrids)=1
      allocate ( Lold(Ngrids) )
      Lold(1:Ngrids)=1
      allocate ( NTC(Ngrids) )
      RETURN
      END SUBROUTINE allocate_stepping
      END MODULE mod_stepping
