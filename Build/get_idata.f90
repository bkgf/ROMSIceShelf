      SUBROUTINE get_idata (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads input data that needs to be obtained only once.  !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_parallel
      USE mod_scalars
      USE mod_stepping
      USE mod_tides
!
      USE nf_fread3d_mod, ONLY : nf_fread3d
      USE nf_fread4d_mod, ONLY : nf_fread4d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
      integer :: LBi, UBi, LBj, UBj
      integer :: itrc, is
      real(r8) :: time_save = 0.0_r8
!
      SourceFile='get_idata.F'
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
!
!-----------------------------------------------------------------------
!  Tide period, amplitude, phase, and currents.  In grid refinement,
!  only the coarser grid (RefineScale(ng)=0) tidal forcing data is
!  processed and needed.
!-----------------------------------------------------------------------
!
!  Tidal Period.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (iic(ng).eq.0) THEN
          CALL get_ngfld (ng, iNLM, idTper, ncFRCid(idTper,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    1, MTC, 1, 1, 1, NTC(ng), 1,                  &
     &                    TIDES(ng) % Tperiod)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
!
!  Tidal elevation amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (iic(ng).eq.0) THEN
          time_save=time(ng)
          time(ng)=8640000.0_r8
          tdays(ng)=time(ng)*sec2day
          CALL get_2dfld (ng, iNLM, idTzam, ncFRCid(idTzam,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
     &                    TIDES(ng) % SSH_Tamp)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_2dfld (ng, iNLM, idTzph, ncFRCid(idTzph,ng),         &
     &                    nFfiles(ng), FRC(1,ng), update(1),            &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
     &                    TIDES(ng) % SSH_Tphase)
          IF (exit_flag.ne.NoError) RETURN
          time(ng)=time_save
          tdays(ng)=time(ng)*sec2day
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
      RETURN
      END SUBROUTINE get_idata
