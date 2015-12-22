      SUBROUTINE ana_tclima (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets analytical tracer climatology fields.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_tclima_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      CLIMA(ng) % tclm)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(33)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_tclima
!
!***********************************************************************
      SUBROUTINE ana_tclima_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            tclm)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange4d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tclm(LBi:,LBj:,:,:)
#else
      real(r8), intent(out) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
      real(r8) :: val1, val2, val3, val4
#ifdef ISOMIP_PLUS
      real(r8) :: Tbot, Sbot, T02, S02, Bmax
#endif
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer climatology.
!-----------------------------------------------------------------------
!
#if defined DOUBLE_GYRE
      val1=(44.69_r8/39.382_r8)**2
      val2=val1*(rho0*100.0_r8/g)*(5.0E-5_r8/((42.689_r8/44.69_r8)**2))
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            val3=T0(ng)+val2*EXP(GRID(ng)%z_r(i,j,k)/100.0_r8)*         &
     &           (10.0_r8-0.4_r8*TANH(GRID(ng)%z_r(i,j,k)/100.0_r8))
            val4=GRID(ng)%yr(i,j)/el(ng)
            tclm(i,j,k,itemp)=val3-3.0_r8*val4
# ifdef SALINITY
            tclm(i,j,k,isalt)=34.5_r8-0.001_r8*GRID(ng)%z_r(i,j,k)-val4
# endif
          END DO
        END DO
      END DO
#elif defined ISOMIP_PLUS
        Tbot = 1.0_r8 !warm
        Sbot = 34.7_r8 !warm
!        Tbot = -1.9_r8 !cold
!        Sbot = 34.55_r8 !cold
        T02 = -1.9_r8
        S02 = 33.8_r8
        Bmax = 720.0_r8
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            tclm(i,j,k,itemp)=T02 + (Tbot - T02)*(-GRID(ng)%z_r(i,j,k)/Bmax)
            tclm(i,j,k,isalt)=S02 + (Sbot - S02)*(-GRID(ng)%z_r(i,j,k)/Bmax)
!      IF(i.eq.20.and.j.eq.235) THEN
!          write(6,*) tclm(i,j,k,itemp)
!      ENDIF
!      IF(i.eq.20.and.j.eq.240) THEN
!          write(6,*) tclm(i,j,k,itemp)
!      ENDIF


          END DO
        END DO
      END DO

#else
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            tclm(i,j,k,itemp)=???
            tclm(i,j,k,isalt)=???
          END DO
        END DO
      END DO
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NAT
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            tclm(:,:,:,itrc))
        END DO
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tclm)
#endif

      RETURN
      END SUBROUTINE ana_tclima_tile
