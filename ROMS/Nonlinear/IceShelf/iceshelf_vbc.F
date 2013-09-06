#include "cppdefs.h"
#undef ICESHELF_VBC_HFPROPTOM
#define GAMMAS_Mc1987 
      MODULE iceshelf_vbc_mod
#ifdef ICESHELF
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!========================================== Benjamin K. Galton-Fenzi ==!
!                                                                      !
! This module sets the ice-shelf/ocean vertical boundary conditions    !
!                                                                      !
!                                                                      !
!   Ice shelf         Ti, Si                                           !
!                                                                      !
!  ================================                                    !
!                                                                      !
!   Laminar sublayer  Tb, Sb, ustar                                    !                      
!                                                                      !
!  --------------------------------                                    !
!                                                                      !
!   Ocean             Tm, Sm, u                                        !
!                                                                      !
!                                                                      !
! References:                                                          ! 
!                                                                      !
! Hellmer, H. H., and D. Olbers, A two-dimensional model for the       !
!  thermohaline circulation under an ice shelf, Antarctic Science,     !
!  1, 325–33, 1989.                                                    ! 
!                                                                      !
! Scheduikat, M., and D. J. Olbers, A one-dimensional mixed layer      !
!  model beneath the Ross Ice Shelf with tidally induced vertical      !
!  mixing, Antarctic Science, 2(1), 29–42, 1990.                       !
!                                                                      !  
! Holland, D. M., and A. Jenkins, Modeling thermodynamic ice-ocean     !
!  interactions at the base of an ice shelf, Journal of Physical       !
!  Oceanography, 29, 1787–1800, 1999.                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: iceshelf_vbc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE iceshelf_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_iceshelf
      USE mod_coupling
      USE mod_stepping
#ifdef ICESHELF_MORPH
      USE mod_iceshelfvar
#endif
!
      implicit none
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
# endif
      CALL iceshelf_vbc_tile (ng, tile,                                 &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng), nnew(ng),                            &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % zice,                               &
     &                   GRID(ng) % f,                                  &
#ifdef ICESHELF_MORPH
     &                   ICESHELFVAR(ng) % iceshelf_draft,              &
#endif
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % rho,                               &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
#   endif
     &                   FORCES(ng) % stflx                             &
     &                   )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
# endif
      RETURN
      END SUBROUTINE iceshelf_vbc
!
!***********************************************************************
      SUBROUTINE iceshelf_vbc_tile (ng, tile,                           &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,nnew,                               &
     &                         z_r, z_w,                                &
     &                         zice,f,                                  &
# ifdef ICESHELF_MORPH
     &                         iceshelf_draft,                          &
# endif
     &                         u, v, t, rho,                            &
     &                         sustr, svstr,                            &
#   ifdef SHORTWAVE
     &                         srflx,                                   &
#   endif
     &                         stflx)                                   
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE mod_iceshelf
#ifdef ICESHELF_MORPH
      USE mod_iceshelfvar
#endif
      USE bc_2d_mod
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs,nnew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
#ifdef ICESHELF_MORPH
      real(r8), intent(inout) :: iceshelf_draft(LBi:,LBj:,:)
#endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
#ifdef ICESHELF_MORPH
      real(r8), intent(in) :: iceshelf_draft(LBi:UBi,LBj:UBj,2)
#endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
# endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
!      integer :: IstrR, IENDR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, itrc
      real(r8) :: Pradj,Scadj, Cdrt, cp_i, gammaS, gammaT
!      real(r8) :: a,b,c,Pr,Sc,Cd,small
!      real(r8) :: cp_w,ustar,visc,turb,TFb,rho_i
      real(r8) :: Sm,Tm,Tb,Sb,rhoi_on_rho0,ustar,TFb,turb
#ifdef ICESHELF_VBC_HFPROPTOM
      real(r8) :: Ex1,Ex2,Ex3,Ex4,Ex5,Ex6,Ep1,Ep2,Ep3,Ep31,Ep4,Ep5
      real(r8) :: Sr1,Sr2,Sf1,Sf2,Tf1,Tf2,rho_r
#else
      real(r8) :: mflag,TFi,cff3 
#endif
      real(r8) :: cff1, cff2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: m
#ifdef ICESHELF_TEOS10
      real(r8) :: ct
#endif
#  if (!defined BBL_MODEL) && defined UV_LOGDRAG
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk
#  endif
#  include "set_bounds.h"
!-----------------------------------------------------------------------
!  If ice shelf cavities, zero out for now the surface tracer flux
!  over the ice.
!-----------------------------------------------------------------------
!
# ifdef SALINITY
          Sm=MAX(0.0_r8,t(i,j,N(ng),nrhs,isalt))
# else
          Sm=0.0_r8
# endif

#    ifdef ANA_SEAICE
! Analytical sea-ice forcing in the open water areas
      tyear = MOD(tdays(ng)-DSTART,365.25_r8)
      IF (tyear.le.59.0_r8) THEN
        sfcTemp = Tmax
        sfcSalt = saltMin
      ELSE IF (tyear.le.259.0_r8) THEN
        sfcTemp = Tmin
        sfcSalt = saltMin + sRateInc * (tyear - 59.0_r8)
      ELSE IF (tyear.le.319.0_r8) THEN
        sfcTemp = Tmin
        sfcSalt = saltMax + sRateDec * (259.0_r8 - tyear)
      ELSE
        sfcTemp = Tmax
        sfcSalt = saltMin
      endif
#    endif

#    ifdef ICESHELF_3EQN_VBC
      Pradj = 12.5_r8*(Pr**(2/3))
      Scadj = 12.5_r8*(Sc**(2/3))
      Cdrt =  sqrt(Cd)
      cp_i =  152.5_r8+7.122_r8*(273.15_r8+Ti)
      rhoi_on_rho0 = rho_i/rho0
#   endif

      DO j=Jstr,JEND
        DO i=Istr,IEND
          IF (zice(i,j).ne.0.0_r8) THEN
#if defined ICESHELF_2EQN_VBC
          TFb = a*Sm+b+c*zice(i,j)
# ifdef ICESHELF_TEOS10
          ct = ct_from_theta(Sm,t(i,j,N(ng),nrhs,itemp))
          Tm = t_from_ct(Sm,ct,-zice(i,j))
# else
         CALL potit(Sm,t(i,j,N(ng),nrhs,itemp),-zice(i,j),0.0_r8,Tm,i,j)
# endif
!            write(6,*) gamma,TFb,Tm
            stflx(i,j,itemp)=0.0_r8 !gamma*(TFb-Tm)
            stflx(i,j,isalt)=Cp*stflx(i,j,itemp)*refSalt/L

#  ifdef ICESHELF_MORPH
            m(i,j) = stflx(i,j,isalt)/Sm
#  endif

#elif defined ICESHELF_3EQN_VBC
! Get the insitu salinity and temperature
# ifdef ICESHELF_TEOS10
          ct = ct_from_theta(Sm,t(i,j,N(ng),nrhs,itemp))
          Tm = t_from_ct(Sm,ct,-zice(i,j))
# else
         CALL potit(Sm,t(i,j,N(ng),nrhs,itemp),-zice(i,j),0.0_r8,Tm,i,j)
# endif

!    Calculate exchange coefficients. Note that this assumes log profile 
            ustar = SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+    &
     &                        (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
!#    ifdef MASKING
!            ustar = ustar*rmask(i,j)
!#    endif

           if(ustar.lt.small) ustar = small
#    ifdef GAMMAS_KY1972
            turb = 2.12*LOG((ustar*Hz(i,j,N(ng)))/visc)
            gammaT = ustar/(turb + Pradj - 9.0_r8)
            gammaS = ustar/(turb + Scadj - 9.0_r8)
#    elif defined GAMMAS_Mc1987
! Uses simplified version of McPhee 1987
            IF (ustar.gt.small.and.ABS(f(i,j)).gt.1.0E-8) THEN
              turb = 2.5*LOG(5300.0_r8*ustar*ustar/ABS(f(i,j))) + 7.12  
            ELSE
              turb = 0.0_r8
            END IF

            gammaT = ustar/(turb + Pradj - 6.0_r8)
            gammaS = ustar/(turb + Scadj - 6.0_r8)
#    else
            gammaT = 1.00e-4_r8
            gammaS = 5.05e-7_r8
#    endif

# ifdef ICESHELF_VBC_HFPROPTOM

! Calculates density RHOW in the boundary layer, interface
! pressure PG[dbar], and solving a quadratic equation for the 
! interface salinity (SB) to determine the melTemp_insitug/freezingrate
! (seta).
          rho_r= rho_i/(1000.0_r8+rho(i,j,N(ng)))
          Ep1  = cp_w*gammaT
          Ep2  = cp_i*gammaS
          Ep3  = L*gammaS
          Ep31 = -rho_r*cp_i*dt_i/zice(i,j)
          Ep4  = b+c*z_w(i,j,N(ng))
          Ep5  = gammaS/rho_r
!
! Temperature gradient within the ice ONLY changes with melTemp_insitug
          TFb = a*Sm+Ep4
          IF (Tm.lt.TFb) THEN
            Ex1 = a*(Ep1+Ep31)
            Ex2 = Ep1*(Tm-Ep4)+Ep3+Ep31*(Ti-Ep4)
            Ex3 = Ep3*Sm
            Ex6 = 0.5
          ELSE
! negative heat flux term in the ice (due to -kappa/D)
            Ex1 = a*(Ep1-Ep2)
            Ex2 = Ep1*(Ep4-Tm)+Ep2*(Ti+a*Sm-Ep4)-Ep3
            Ex3 = Sm*(Ep2*(Ep4-Ti)+Ep3)
            Ex6 = -0.5
          endif
          Ex4 = Ex2/Ex1
          Ex5 = Ex3/Ex1
          Sr1 = 0.25*Ex4*Ex4-Ex5
          Sr2 = Ex6*Ex4
          Sf1 = Sr2+sqrt(Sr1)
          Tf1 = a*Sf1+Ep4
          Sf2 = Sr2-sqrt(Sr1)
          Tf2 = a*Sf2+Ep4
!
! Salinities < 0 psu are NOT defined, i.e....
         IF (Sf1.gt.0.) THEN
           Tb = Tf1
           Sb = Sf1
         ELSE
           Tb = Tf2
           Sb = Sf2
         endif
!
         m(i,j) = Ep5*(1.0-Sm/Sb)
!         m(i,j)= Seta * dtfast
!
! Calculates fluxes:
         stflx(i,j,itemp)=(gammaT-rho_r*m(i,j))*(Tb-Tm)
         stflx(i,j,isalt)=(gammaS-rho_r*m(i,j))*(Sb-Sm)

# else
            TFb = a*Sm + b + c*zice(i,j)
            mflag= 0.5*(1 + SIGN(1.0_r8,Tm-TFb))
            TFi = (1 - mflag)*a*Si + b + c*zice(i,j)

! Calculate coefficents in quadratic to be solved:
           cff1 = L/cp_w + mflag*(cp_i/cp_w)*(TFi-Ti)
           cff2 = gammaS*(L/cp_w + mflag*(cp_i/cp_w)*(TFb-Ti))          &
     &          + gammaT*(TFi-Tm)
           cff3 = gammaS*gammaT*(TFb - Tm)

! Calculate melt rate:
           m(i,j) = -(cff2 - SQRT(cff2*cff2 - 4*cff1*cff3))/(2*cff1)

! Calculate basal temperature and salinity:
           Tb = (gammaT*Tm+mflag*(cp_i/cp_w)*m(i,j)*Ti-(L/cp_w)*m(i,j)) &
              / (gammaT + mflag*(cp_i/cp_w)*m(i,j))
           Sb = (Tb - b - c*zice(i,j))/a

           stflx(i,j,itemp)=(gammaT+rhoi_on_rho0*m(i,j))*(Tb-Tm)
           stflx(i,j,isalt)=(gammaS+rhoi_on_rho0*m(i,j))*(Sb-Sm)
!           write(6,*) gammaT,gammaS,Tb,Sb,Tm,Sm,rhoi_on_rho0
# endif

#else
      DO itrc=1,NT(ng)
!        DO j=JstrR,JENDR
!          DO i=IstrR,IENDR
!            IF (zice(i,j).ne.0.0_r8) THEN
              stflx(i,j,itrc)=0.0_r8
      END DO

          ELSE
            m(i,j) = 0.0_r8
# if defined ANA_SEAICE
            stflx(i,j,itemp)=Hz(i,j,N(ng))*                             &
                             (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax
            stflx(i,j,isalt)=Hz(i,j,N(ng))*                             &
                             (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax
# else
          DO itrc=1,NT(ng)
             stflx(i,j,itrc)=0.0_r8
          END DO
# endif

!          END IF
!        END DO
!      END DO

#endif

          END IF
        END DO
      END DO

!-----------------------------------------------------------------------
! Store old ice shelf thickness.
!-----------------------------------------------------------------------
!
# if defined ICESHELF_MORPH
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          iceshelf_draft(i,j,nnew)=zice(i,j)-                           &
     &                             MIN(m(i,j)*dt(ng),1.0e-3_r8) 
          IF(i.eq.40.and.j.eq.40) THEN
          write(6,*) m(i,j),zice(i,j),iceshelf_draft(i,j,nnew) 
          ENDIF
        END DO
      END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            iceshelf_draft(:,:,nnew))
        END IF
# endif

#   ifdef SHORTWAVE
      DO j=JstrR,JENDR
        DO i=IstrR,IENDR
          IF (zice(i,j).ne.0.0_r8) THEN
            srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------

#   if defined UV_LOGDRAG
!
!  Set logarithmic ice shelf cavity stress.
!
      DO j=JstrV-1,JEND
        DO i=IstrU-1,IEND
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/Zob(ng))
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
        END DO
      END DO
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            sustr(i,j)=-0.5_r8*(wrk(i-1,j)+wrk(i,j))*                   &
     &                 u(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            svstr(i,j)=-0.5_r8*(wrk(i,j-1)+wrk(i,j))*                   &
     &                 v(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
#   elIF defined UV_QDRAG
!
!  Set quadratic ice shelf cavity stress.
!
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff2*cff2)
            sustr(i,j)=-rdrg2(ng)*u(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            svstr(i,j)=-rdrg2(ng)*v(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
#   elif defined UV_LDRAG
!
!  Set linear ice shelf cavity stress.
!
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            sustr(i,j)=-rdrg(ng)*u(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            svstr(i,j)=-rdrg(ng)*v(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
#   else
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    sustr,svstr)
#endif
      RETURN
      END SUBROUTINE iceshelf_vbc_tile

# ifdef ICESHELF_TEOS10
      REAL FUNCTION ct_from_theta(s,theta)
!   conservative temperature from potential temperature, as in
!   Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted
!   JAOT
!   s                : salinity                  (psu)
!   theta            : potential temperature     (deg C, ITS-90)
!   ct_from_theta    : conservative temperature  (deg C, ITS-90) 
!   calls            : penthalpy_F
!   check value      : ct_from_theta(20,20) = 20.45274961282756
!   DRJ on 10/12/03
      implicit real*8(a-h,o-z)

      rec_Cp0 = 2.50494524832013d-4; !Cp0 = 3992.10322329649d0
      ct_from_theta = rec_Cp0*penthalpy_F(s,theta)

      RETURN
      END FUNCTION ct_from_theta

      REAL FUNCTION t_from_ct(s,ct,p)
!   in-situ temperature from conservative temperature
!   Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted
!   JAOT
!   s                : salinity                           (psu)
!   ct               : conservative temperature           (deg C,
!   ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   t                : in-situ temperature                (deg C,
!   ITS-90)
!   calls            : theta_from_ct and theta_from_t
!   check value      : t_from_ct(20,20,1000) = 19.72671624220263 with 1
!                                                      loop in
!                                                      theta_from_ct
!                      t_from_ct(20,20,1000) = 19.72671627735695 with 2
!                                                      loops in
!                                                      theta_from_ct
!   DRJ on 11/10/03

      implicit real*8(a-h,o-z)

      pr0 = 0.d0
      theta = theta_from_ct(s,ct)
      t_from_ct = theta_from_t(s,theta,pr0,p)

      RETURN
      END FUNCTION t_from_ct

      REAL FUNCTION theta_from_ct(s,ct)
!   potential temperature from conservative temperature, as in
!   Jackett, McDougall, Feistel, Wright & Griffies (2004), submitted
!   JAOT
!   s                : salinity                  (psu)
!   ct               : conservative temperature  (deg C, ITS-90)
!   theta_from_ct    : potential temperature     (deg C, ITS-90)
!   calls            : cp0_F and ct_from_theta
!   check value      : theta_from_ct(20,20) = 19.55627910604363 
!   DRJ on 03/06/05

      implicit real*8(a-h,o-z)

      a0 = -1.446013646344788d-2
      a1 = -3.305308995852924d-3
      a2 =  1.062415929128982d-4
      a3 =  9.477566673794488d-1
      a4 =  2.166591947736613d-3
      a5 =  3.828842955039902d-3

      b0 =  1.000000000000000d+0
      b1 =  6.506097115635800d-4
      b2 =  3.830289486850898d-3
      b3 =  1.247811760368034d-6

      a5ct = a5*ct; b3ct = b3*ct
      ct_factor = (a3+a4*s+a5ct)
      th0_num = a0+s*(a1+a2*s)+ct*ct_factor
      rec_th0_den = 1.0d0/(b0+b1*s+ct*(b2+b3ct))
      th0 = th0_num*rec_th0_den
      ct0 = ct_from_theta(s,th0)
      dth_dct = (ct_factor+a5ct-(b2+b3ct+b3ct)*th0)*rec_th0_den
      theta = th0-(ct0-ct)*dth_dct

      rec_Cp0 = 2.50494524832013d-4; ! Cp0 = 3992.10322329649;
      dct = ct_from_theta(s,theta)-ct;
      dct_dth = (theta+273.15d0)*de_dt0_F(s,theta)*rec_Cp0

      theta_from_ct = theta-dct/dct_dth
!       NOTE: theta has a maximum error of 6.0x10^-14

      RETURN
      END FUNCTION theta_from_ct

      REAL FUNCTION theta_from_t(s,t,p,pr)
!   potential temperature from in-situ temperature, as in 
!   Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted
!   JAOT
!   s                : salinity                           (psu)
!   t                : in-situ temperature                (deg C,
!   ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   pr               : reference pressure                 (dbar)
!   theta_from_t     : potential temperature              (deg C,
!   ITS-90)
!   calls            : de_dt_F and entropy_diff_F
!   check values     : theta_from_t(35,20,4000,0) = 19.21108374301637
!                                                           (with
!                                                           nloops=1)
!                      theta_from_t(35,20,4000,0) = 19.21108374301171
!                                                           (with
!                                                           nloops=2)
!
!   DRJ on 10/12/03

      implicit real*8(a-h,o-z)
      integer :: nloops, n

      th0 = t+(p-pr)*( 8.65483913395442d-6   -                          &
     &             s*  1.41636299744881d-6   -                          &
     &        (p+pr)*  7.38286467135737d-9   +                          &
     &             t*(-8.38241357039698d-6   +                          &
     &             s*  2.83933368585534d-8   +                          &
     &             t*  1.77803965218656d-8   +                          &
     &          (p+pr)*  1.71155619208233d-10))
      de_dt = 13.6d0
      nloops = 2                                      ! default
!    NOTE: nloops = 1 gives theta with a maximum error of 5.48x10^-06
!          nloops = 2 gives theta with a maximum error of 2.84x10^-14
      n = 1
      do while(n.le.nloops)
       dentropy = entropy_diff_F(s,t,p,th0,pr)
       theta = th0-dentropy/de_dt
       theta = 0.5d0*(theta+th0)
       de_dt = de_dt_F(s,theta,pr)
       theta = th0-dentropy/de_dt
       n = n+1; th0 = theta
      end do
      theta_from_t = th0

      RETURN
      END FUNCTION theta_from_t

      REAL FUNCTION entropy_diff_F(s,t,p,th0,pr)
!   entropy difference from differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114, and taking differences
!   s                : salinity                           (psu)
!   t                : in-situ temperature                (deg C,
!   ITS-90)
!   th0              : potential temperature              (deg C,
!   ITS-90) 
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   pr               : reference pressure                 (dbar)
!   entropy_diff     : entropy(s,th0,pr)-entropy(s,t,p)   (kJ/kgK)
!   DRJ on 10/12/03

      implicit real*8(a-h,o-z)

      real*8 kTp

      x2 = 2.5d-2*s; x = sqrt(x2); y = 2.5d-2*th0; z = 1.0d-4*pr; z1 = z
      fTp = -5.90578348518236d0 + y*(24715.571866078d0 +                &
     &         y*(-2210.2236124548363d0 + y*(592.743745734632d0 +       &
     &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
     &                               21.35571525415769d0*y)*y))))
      gTp = y*(-1858.920033948178d0 + y*(781.281858144429d0 +           &
     &                y*(-388.6250910633612d0 + 87.18719211065d0*y)))
      hTp = y*(317.440355256842d0 + y*(-202.5696441786141d0 +           &
                                         67.5605099586024d0*y))
      if(z.ne.0.d0) then
        yz = y*z
        fTp = fTp + z*(270.983805184062d0 + z*(-776.153611613101d0 +    &
     &               z*(196.51255088122d0 + z*(-28.9796526294175d0 +    &
     &       2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 +      &
     &           z*(1513.116771538718d0 + z*(-546.959324647056d0 +      &
     &         (111.1208127634436d0 - 8.68841343834394d0*z)*z))) +      &
     &        y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 +      &
     &          z*(718.6359919632359d0 + z*(-146.4037555781616d0 +      &
     &   4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 +      &
     &           z*(1207.261522487504d0 + z*(-608.785486935364d0 +      &
     &  105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z +      &
     &           z*(-381.06836198507096d0 + (133.7383902842754d0 -      &
     &       49.023632509086724d0*z)*z)) + z*(973.091553087975d0 +      &
     &              z*(-602.603274510125d0 + (276.361526170076d0 -      &
                                 32.40953340386105d0*z)*z))))))
        gTp = gTp + z*(-729.116529735046d0 + z*(343.956902961561d0 +    &
     &               z*(-124.687671116248d0 + (31.656964386073d0 -      &
     &    7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 +      &
     &            z*(-674.819060538734d0 + z*(356.629112415276d0 +      &
     &           z*(-88.4080716616d0 + 15.84003094423364d0*z)))) +      &
     &      y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 +      &
     &        145.9491676006352d0*z)) + z*(-2082.7344423998043d0 +      &
                 z*(614.668925894709d0 + z*(-340.685093521782d0 +       &
                                      33.3848202979239d0*z)))))
        hTp = hTp + z*(175.292041186547d0 + z*(-83.1923927801819d0 +    &
     &     29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z -       &
     &          938.26075044542d0*y*z) + z*(-766.116132004952d0 +       &
     &          (108.3834525034224d0 - 51.2796974779828d0*z)*z))
      end if
       y = 2.5d-2*t; z = 1.0d-4*p
        fTp = fTp +5.90578348518236d0 - y*(24715.571866078d0 +          &
     &          y*(-2210.2236124548363d0 + y*(592.743745734632d0 +      &
     &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
     &                              21.35571525415769d0*y)*y))))
        gTp = gTp - y*(-1858.920033948178d0 + y*(781.281858144429d0 +   &
     &       y*(-388.6250910633612d0 + 87.18719211065d0*y)))
        hTp = hTp - y*(317.440355256842d0 + y*(-202.5696441786141d0 +   &
     &                                     67.5605099586024d0*y))
        kTp = 22.6683558512829d0*(z1-z)
      if(z.ne.0.d0) then
        yz = y*z
         fTp = fTp - ( z*(270.983805184062d0 + z*(-776.153611613101d0 + &
     &              z*(196.51255088122d0 + z*(-28.9796526294175d0 +     &
     &         2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 +    &
     &             z*(1513.116771538718d0 + z*(-546.959324647056d0 +    &
     &           (111.1208127634436d0 - 8.68841343834394d0*z)*z))) +    &
     &          y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 +    &
     &            z*(718.6359919632359d0 + z*(-146.4037555781616d0 +    &
     &     4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 +    &
     &             z*(1207.261522487504d0 + z*(-608.785486935364d0 +    &
     &    105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z +    &
     &             z*(-381.06836198507096d0 + (133.7383902842754d0 -    &
     &         49.023632509086724d0*z)*z)) + z*(973.091553087975d0 +    &
     &                z*(-602.603274510125d0 + (276.361526170076d0 -    &
                                 32.40953340386105d0*z)*z)))))))
         gTp = gTp - ( z*(-729.116529735046d0 + z*(343.956902961561d0 + &
     &                z*(-124.687671116248d0 + (31.656964386073d0 -     &
     &     7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 +     &
     &             z*(-674.819060538734d0 + z*(356.629112415276d0 +     &
     &            z*(-88.4080716616d0 + 15.84003094423364d0*z)))) +     &
     &       y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 +     &
     &         145.9491676006352d0*z)) + z*(-2082.7344423998043d0 +     &
     &             z*(614.668925894709d0 + z*(-340.685093521782d0 +     &
     &                                 33.3848202979239d0*z))))))
         hTp = hTp - ( z*(175.292041186547d0 + z*(-83.1923927801819d0 + &
     &       29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z -     &
     &            938.26075044542d0*y*z) + z*(-766.116132004952d0 +     &
     &          (108.3834525034224d0 - 51.2796974779828d0*z)*z)))
      end if
      entropy_diff_F = 2.5d-2*(fTp+x2*(gTp+x*hTp+x2*kTp))

      RETURN
      END FUNCTION entropy_diff_F

      REAL FUNCTION de_dt_F(s,t,p)
!   d(entropy)/dt from twice differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity                           (psu)
!   t                : in-situ temperature                (deg C,
!   ITS-90)
!   p                : gauge pressure                     (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   de_dt_F          : d(entropy)/dt                      (J/kgK^2)
!   check value      : de_dt_F(35,20,4000) = 13.35022011370635
!   DRJ on 10/12/03

      implicit real*8(a-h,o-z)

      x2 = 2.5d-2*s; x = sqrt(x2); y = 2.5d-2*t; z = 1.0d-4*p
      de_dt = 24715.571866078d0 + x2*(-1858.920033948178d0 +            &
     &                x*(317.440355256842d0 + y*(-405.1392883572282d0 + &
     &          202.6815298758072d0*y)) + y*(1562.563716288858d0 +      &
     &          y*(-1165.8752731900836d0 + 348.7487684426d0*y))) +      &
     &         y*(-4420.4472249096725d0 + y*(1778.231237203896d0 +      &
     &            y*(-1160.5182516851419d0 + (569.531539542516d0 -      &
     &                              128.13429152494615d0*y)*y)))
      if(z.ne.0.d0) then
        de_dt = de_dt +                                                 &
     &           z*(-2910.0729080936d0 + x2*(1721.528607567954d0 +      &
     &        y*(-4165.4688847996085d0 + 3572.7449038462437d0*y) +      &
     &            x*(-766.116132004952d0 + (2761.9195908075417d0 -      &
     &         2814.78225133626d0*y)*y)) + y*(4035.04669887042d0 +      &
     &            y*(-4775.621344883664d0 + y*(3892.3662123519d0 +      &
     &      y*(-1905.341809925355d0 + 404.50541014508605d0*y)))) +      &
     &          z*(1513.116771538718d0 + x2*(-674.819060538734d0 +      &
     &              108.3834525034224d0*x + (1229.337851789418d0 -      &
     &        896.713693665072d0*y)*y) + y*(-2996.162344914912d0 +      &
     &            y*(3621.784567462512d0 + y*(-2410.4130980405d0 +      &
     &          668.691951421377d0*y))) + z*(-546.959324647056d0 +      &
     &             x2*(356.629112415276d0 - 51.2796974779828d0*x +      &
     &         y*(-681.370187043564d0 + 437.84750280190565d0*y)) +      &
     &         y*(1437.2719839264719d0 + y*(-1826.356460806092d0 +      &
     &        (1105.446104680304d0 - 245.11816254543362d0*y)*y)) +      &
     &             z*(111.1208127634436d0 + x2*(-88.4080716616d0 +      &
     &           66.7696405958478d0*y) + y*(-292.8075111563232d0 +      &
     &         (316.49805267936244d0 - 129.6381336154442d0*y)*y) +      &
     &             (-8.68841343834394d0 + 15.84003094423364d0*x2 +      &
     &                              9.978426372534301d0*y)*z))))
      end if
      de_dt_F = 6.25d-4*de_dt

      RETURN
      END FUNCTION de_dt_F
      REAL FUNCTION penthalpy_F(s,th)
!   potential enthalpy from differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity (psu)
!   th               : potential temperature (deg C, ITS-90)
!   penthalpy_F      : potential enthalpy (J/kg)
!   check value      : penthalpy_F(20,20) = 81649.4876546448
!   DRJ on 10/12/03

      implicit real*8(a-h,o-z)

      x2 = 2.5d-2*s; x = sqrt(x2); y = 2.5d-2*th;
      penthalpy_F = 61.013624165232955d0 + y*(168776.46138048015d0  +   &
     &            y*(-2735.2785605119643d0 + y*(2574.2164453821442d0 +  &
     &              y*(-1536.6644434977545d0 + y*(545.734049793163d0 +  &
     &        (-50.910917284743334d0 - 18.30489878927802d0*y)*y))))) +  &
     &             x2*(416.31512917743896d0 + x*(937.9793807560891d0 +  &
     &              x*(-3140.435779506947d0 + x*(2975.170149976973d0 +  &
     &            x*(-1760.137081144729d0 + x*414.5655751783703d0))) +  &
     &              y*(2167.72082596016d0 + y*(-1224.5772800562902d0 +  &
     &             y*(326.3074029273967d0 + 50.6703824689518d0*y)))) +  &
     &               y*(-12694.10018182362d0 + y*(4405.71847182968d0 +  &
     &            y*(-2132.9690185026416d0 + y*(303.91071982808035d0 +  &
                                            69.74975368852d0*y)))))

      RETURN
      END FUNCTION penthalpy_F

      REAL FUNCTION de_dt0_F(s,th)
!   d(entropy)/dt from twice differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C,
!   ITS-90)
!   de_dt0_F         : d(entropy)/dt                      (J/kgK^2)
!   check value      : de_dt0_F(35,20) = 13.63256369213874
!   DRJ on 30/06/05

      implicit real*8(a-h,o-z)

      x2 = 2.5d-2*s; x = sqrt(x2); y = 2.5d-2*th
      de_dt = 24715.571866078d0 + x2*(-1858.920033948178d0 +            &
     &           x*(317.440355256842d0 + y*(-405.1392883572282d0 +      &
     &          202.6815298758072d0*y)) + y*(1562.563716288858d0 +      &
     &          y*(-1165.8752731900836d0 + 348.7487684426d0*y))) +      &
     &         y*(-4420.4472249096725d0 + y*(1778.231237203896d0 +      &
     &            y*(-1160.5182516851419d0 + (569.531539542516d0 -      &
     &                              128.13429152494615d0*y)*y)))
      de_dt0_F = 6.25d-4*de_dt

      RETURN
      END FUNCTION de_dt0_F

      REAL FUNCTION fp_theta(s,p,sat)
!   potential temperature freezing point of seawater, as in
!   Jackett, McDougall, Feistel, Wright and Griffies (2004), submitted
!   JAOT
!   s                : salinity                               (psu)
!   p                : gauge pressure                         (dbar)
!                      (absolute pressure - 10.1325 dbar)
!   sat              : string variable
!                       'air-sat'  - air saturated water       
!                       'air-free' - air free water
!   fp_theta         : potential freezing temperature         (deg C,
!   ITS-90)
!   check value      : fp_theta(35,200,'air-sat')   = -2.076426227617581
!   deg C
!                      fp_theta(35,200,'air-free') = -2.074408175943127
!                      deg C
!   DRJ on 2/6/04

      implicit real*8(a-h,o-z)
      character*(*) sat
      sqrts = sqrt(s)
      tf_num =                    2.5180516744541290d-03    +           &
     &                        s*(-5.8545863698926184d-02    +           &
     &                    sqrts*( 2.2979985780124325d-03    -           &
     &                   sqrts *  3.0086338218235500d-04))  +           &
     &                        p*(-7.0023530029351803d-04    +           &
     &                        p*( 8.4149607219833806d-09    +           &
     &                       s *  1.1845857563107403d-11));           
      tf_den =                          1.0000000000000000d+00    +     &
     &                              p*(-3.8493266309172074d-05    +     &
     &                             p *  9.1686537446749641d-10)   +     &
     &                      s*s*sqrts*  1.3632481944285909d-06
      fp_theta = tf_num/tf_den;
      if(sat.eq.'air-sat') then
      fp_theta = fp_theta          -2.5180516744541290d-03    +         &
     &                             s *  1.428571428571429d-05
      elseif(sat.eq.'air-free') then
      continue
      else
      print *, '***  Error in fp_theta.f90: invalid third argument  ***'
      print *
      stop
      endif
      return

      END FUNCTION fp_theta
# else
! *********************************************************************
      SUBROUTINE potit(Sal,theta,Pres,RPres,Temp,i,j)
! *********************************************************************
! Calculates from the salinity (sal, psu), potential temperature 
! (theta, degC) and reference pressure (pres, dbar) the in-situ 
! temperaure (Temp_insitu, degC) related to the in-situ pressure 
! (rfpres, dbar) with the help of an iterative method.
      USE mod_kinds

      integer, intent(in)   :: i, j
      real(r8), intent(in)  :: Sal, Pres,theta
      real(r8), intent(out) :: Temp

      integer               :: ind
      real(r8)              :: tpmd, theta1, thetad, epsi, RPres

      data tpmd / 0.001 /

      epsi = 0.
      do ind=1,100
      Temp   = theta+epsi
      thetad  = thetaa(Sal,Temp,Pres,RPres)-theta
      IF(abs(thetad).lt.tpmd) return
       epsi = epsi-thetad
      ENDdo
      write(6,*) ' WARNING!',                                           &
     & ' in-situ temperature calculation has not converged!', i,j
      RETURN
      END SUBROUTINE potit
! *********************************************************************
      REAL FUNCTION thetaa(Sal,Temp,Pres,RPres)
! Calculates from the salinity (sal, psu), the in-situ temperature 
! (Temp, degC) and the in-situ pressure press, dbar) the potential 
! temperature (Theta, degC) converted to the reference pressure
! (RPres, dbar). A Runge-Kutta procedure of the fourth order is used.
!
! Check value: theta   =    36.89073  degC
!         given sal    =    40.0      psu
!               Temp   =    40.0      degC
!               pres   = 10000.000    dbar
!               rfpres =     0.000    dbar
      USE mod_kinds

      real(r8), intent(in) ::  Sal,Temp,Pres,RPres
      real(r8)             ::  p,t,dp,dt,q,ct2,ct3,cq2a,cq2b,cq3a,cq3b

      data ct2 ,ct3  /0.29289322 ,  1.707106781/
      data cq2a,cq2b /0.58578644 ,  0.121320344/
      data cq3a,cq3b /3.414213562, -4.121320344/

      p  = Pres
      t  = Temp
      dp = RPres-Pres
      dt = dp*dTemp(Sal,t,p)
      t  = t +0.5*dt
      q = dt
      p  = p +0.5*dp
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct2*(dt-q)
      q  = cq2a*dt + cq2b*q
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct3*(dt-q)
      q  = cq3a*dt + cq3b*q
      p  = RPres
      dt = dp*dTemp(Sal,t,p)

      thetaa = t + (dt-q-q)/6.0

      END FUNCTION thetaa

! *********************************************************************
! *********************************************************************
      REAL FUNCTION dTemp(Sal,Temp,Pres)
! Calculates from the salinity (Sal,psu), the in-situ Temperature
! (Temp, degC) and the in-situ pressure (Pres, dbar) the adiabatic 
! temperature gradient (dTemp, K Dbar^-1).
!
! Check values: dTemp  =     3.255976E-4 K dbar^-1
!          given Sal    =    40.0         psu
!                Temp   =    40.0         degC
!                Pres   = 10000.000       dbar
      USE mod_kinds

      real(r8), intent(in) :: Sal, Temp, Pres
      real(r8)             :: s0,a0,a1,a2,a3,b0,b1,c0,c1,c2,c3
      real(r8)             :: d0,d1,e0,e1,e2,ds

      data s0 /35.0D0/
      data a0,a1,a2,a3 /3.5803D-5, 8.5258D-6, -6.8360D-8, 6.6228D-10/
      data b0,b1       /1.8932D-6, -4.2393D-8/
      data c0,c1,c2,c3 /1.8741D-8, -6.7795D-10, 8.7330D-12, -5.4481D-14/
      data d0,d1       /-1.1351D-10, 2.7759D-12/
      data e0,e1,e2    /-4.6206D-13,  1.8676D-14, -2.1687D-16/

      ds = Sal-s0
      dTemp = ( ( (e2*Temp + e1)*Temp + e0 )*Pres                       &
     &      + ( (d1*Temp + d0)*ds                                       &
     &      + ( (c3*Temp + c2)*Temp + c1 )*Temp + c0 ) )*Pres           &
     &      + (b1*Temp + b0)*ds +  ( (a3*Temp + a2)*Temp + a1 )*Temp    &
     &      + a0
      RETURN
      END FUNCTION dTemp
# endif

#endif
      END MODULE iceshelf_vbc_mod

