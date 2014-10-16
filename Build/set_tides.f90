      MODULE set_tides_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group        Robert Hetland   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine adds tidal elevation (m) and tidal currents (m/s) to   !
!  sea surface height and 2D momentum climatologies, respectively.     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_tides
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_tides (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_tides
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 11)
      CALL set_tides_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     NTC(ng),                                     &
     &                     GRID(ng) % angler,                           &
     &                     TIDES(ng) % SSH_Tamp,                        &
     &                     TIDES(ng) % SSH_Tphase,                      &
     &                     TIDES(ng) % Tperiod)
      CALL wclock_off (ng, iNLM, 11)
      RETURN
      END SUBROUTINE set_tides
!
!***********************************************************************
      SUBROUTINE set_tides_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           NTC,                                   &
     &                           angler,                                &
     &                           SSH_Tamp, SSH_Tphase,                  &
     &                           Tperiod)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod
!
!  Imported variables declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: NTC
      real(r8) :: tide_pred
!
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: Tperiod(MTC)
      real(r8), intent(in) :: SSH_Tamp(LBi:,LBj:,:)
      real(r8), intent(in) :: SSH_Tphase(LBi:,LBj:,:)
!
!  Local variables declarations.
!
      logical :: update
      integer :: i, itide, j
      real(r8) :: Cangle, Cphase, Sangle, Sphase
      real(r8) :: angle, cff, phase, omega, ramp
      real(r8) :: bry_cor, bry_pgr, bry_str, bry_val
      integer :: iye, mon, idy, ihr, min
      real(r8) :: JD, tz_t,Tmax, notide
      real(r4) :: sec
      real(r8) :: htid(NTC), gtid(NTC)
      character*6 :: intide(NTC)
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Etide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Utide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vtide
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!  Set time-ramping parameter.
!
      ramp=1.0_r8
!
!-----------------------------------------------------------------------
!  Add tidal elevation (m) to sea surface height climatology.
!-----------------------------------------------------------------------
!
      Etide(:,:)=0.0_r8
      intide(1) = 'm2    '
      intide(2) = 's2    '
      intide(3) = 'n2    '
      intide(4) = 'k2    '
      intide(5) = 'k1    '
      intide(6) = 'o1    '
      intide(7) = 'p1    '
      intide(8) = 'q1    '
      intide(9) = 'mf    '
      intide(10)= 'mm    '
      JD=tdays(ng)+2440000.5d0              ! Julian datys (from B.C.)
      tz_t=0.0d0                            ! Time zone  0 = GMT
      notide=-9999.0d0
      CALL caldat(iye,mon,idy,ihr,min,sec,JD)
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (SSH_Tamp(i,j,1).gt.notide) THEN 
            DO itide=1,NTC
              htid(itide) = SSH_Tamp(i,j,itide)
              gtid(itide) = SSH_Tphase(i,j,itide)
            END DO
              Etide(i,j)=ramp*tide_pred(iye,mon,idy,ihr,min,sec,tz_t,   &
     &                    htid,gtid,                                    &
     &                    intide,                                       &
     &                    NTC)
          END IF
        END DO
      END DO
!
!  If appropriate, load tidal forcing into boundary arrays.  The "zeta"
!  boundary arrays are important for the Flather or reduced physics
!  boundary conditions for 2D momentum. To avoid having two boundary
!  points for these arrays, the values of "zeta_west" and "zeta_east"
!  are averaged at u-points.  Similarly, the values of "zeta_south"
!  and "zeta_north" is averaged at v-points. Noticed that these
!  arrays are also used for the clamped conditions for free-surface.
!  This averaging is less important for that type ob boundary
!  conditions.
!
      IF (LBC(iwest,isFsur,ng)%acquire.or.                              &
     &    LBC(iwest,isUbar,ng)%acquire.or.                              &
     &    LBC(iwest,isVbar,ng)%acquire) THEN
        update=.FALSE.
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=JstrR,JendR
            BOUNDARY(ng)%zeta_west(j)=0.5_r8*(Etide(Istr-1,j)+          &
     &                                        Etide(Istr  ,j))
          END DO
          update=.TRUE.
        END IF
      END IF
!
      IF (LBC(ieast,isFsur,ng)%acquire.or.                              &
     &    LBC(ieast,isUbar,ng)%acquire.or.                              &
     &    LBC(ieast,isVbar,ng)%acquire) THEN
        update=.FALSE.
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrR,JendR
            BOUNDARY(ng)%zeta_east(j)=0.5_r8*(Etide(Iend  ,j)+          &
     &                                        Etide(Iend+1,j))
          END DO
          update=.TRUE.
        END IF
      END IF
!
      IF (LBC(isouth,isFsur,ng)%acquire.or.                             &
     &    LBC(isouth,isUbar,ng)%acquire.or.                             &
     &    LBC(isouth,isVbar,ng)%acquire) THEN
        update=.FALSE.
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrR,IendR
            BOUNDARY(ng)%zeta_south(i)=0.5_r8*(Etide(i,Jstr-1)+         &
     &                                         Etide(i,Jstr  ))
          END DO
          update=.TRUE.
        END IF
      END IF
!
      IF (LBC(inorth,isFsur,ng)%acquire.or.                             &
     &    LBC(inorth,isUbar,ng)%acquire.or.                             &
     &    LBC(inorth,isVbar,ng)%acquire) THEN
        update=.FALSE.
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=IstrR,IendR
            BOUNDARY(ng)%zeta_north(i)=0.5_r8*(Etide(i,Jend  )+         &
     &                                         Etide(i,Jend+1))
          END DO
          update=.TRUE.
        END IF
      END IF
      RETURN
      END SUBROUTINE set_tides_tile
      END MODULE set_tides_mod
