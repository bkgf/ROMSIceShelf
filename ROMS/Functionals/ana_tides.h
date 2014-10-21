      SUBROUTINE ana_tides (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical astronomical tidal forcing.            !
!                                                                      !
!============================================ Benjamin K. Galton-Fenzi =
!
      USE mod_param
      USE mod_tides
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tides_tile (ng, tile, model,                             &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   TIDES(ng) % Tperiod,                           &
     &                   TIDES(ng) % SSH_Tamp,                          &
     &                   TIDES(ng) % SSH_Tphase)
     
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(28)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tides
!
!***********************************************************************
      SUBROUTINE ana_tides_tile (ng, tile, model,                       &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Tperiod, SSH_Tamp, SSH_Tphase)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: Tperiod(10)
      real(r8), intent(out) :: SSH_Tamp(10,LBi:,LBj:)
      real(r8), intent(out) :: SSH_Tphase(10,LBi:,LBj:)
#else
      real(r8), intent(out) :: Tperiod(10)
      real(r8), intent(out) :: SSH_Tamp(10,LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: SSH_Tphase(10,LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, n
      integer :: NTC = 10

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tidal angular period (hours) for components
!    M2, S2, N2, K2, K1, O1, P1, Q1, Mm, Mf
!-----------------------------------------------------------------------
!
          Tperiod(1) = 12.4206011981605_r8
          Tperiod(2) = 12.0000000048_r8
          Tperiod(3) = 12.6583482145719_r8 
          Tperiod(4) = 11.9672348025225_r8
          Tperiod(5) = 23.934469605045_r8
          Tperiod(6) = 25.819341694366_r8
          Tperiod(7) = 24.0658902318985_r8
          Tperiod(8) = 26.8683566006764_r8
          Tperiod(9) = 661.309268024546_r8
          Tperiod(10) = 327.858984441058_r8
!
!-----------------------------------------------------------------------
!  Set 
!    SSH_Tphase, tidal elevation phase angle (degrees, time of maximum elevation 
!      with respect chosen time origin), and 
!    SSH_Tamp, tidal elevation amplitude (meter)
!      for components M2, S2, N2, K2, K1, O1, P1, Q1, Mm, Mf
!-----------------------------------------------------------------------
!

# if defined ICESHELF_TIDES 
      DO n=1,NTC
      DO j=JstrT,JendT
        DO i=IstrT,IendT
            IF(j.eq.JendT.and.n.eq.1) THEN
              SSH_Tphase(n,i,j)=0.0_r8
              SSH_Tamp(n,i,j)=0.0_r8
            ELSE
              SSH_Tphase(n,i,j)=0.0_r8
              SSH_Tamp(n,i,j)=0.0_r8
            ENDIF
          END DO
        END DO
      END DO
# else
      DO n=1,NTC
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            SSH_Tphase(n,i,j)=0.0_r8
            SSH_Tamp(n,i,j)=0.0_r8
          END DO
        END DO
      END DO
# endif

!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Tperiod,SSH_Tamp,SSH_Tphase)
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Tperiod,SSH_Tamp,SSH_Tphase`)
#endif

      RETURN
      END SUBROUTINE ana_tides_tile
