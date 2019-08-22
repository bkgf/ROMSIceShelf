#include "cppdefs.h"
#ifdef ICESHELF2D_TOY
#  define ICESHELF_VBC_HFPROPTOM
#  undef GAMMAS_Mc1987 
#  undef GAMMAS_KY1972
#else
#  undef ICESHELF_VBC_HFPROPTOM
#  define GAMMAS_Mc1987 
#  undef GAMMAS_KY1972
#endif
      MODULE ad_iceshelf_vbc_mod
#ifdef ICESHELF
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!========================================== Benjamin K. Galton-Fenzi ==!
!
! AMM: NOTE - This code is not consistent with the lastest version of
!      set_vbc (e.g. rdrag and rdrag2 are now arrays not constants).
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
      PUBLIC  :: ad_iceshelf_vbc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_iceshelf_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_iceshelf
      USE mod_coupling
      USE mod_stepping
      USE mod_iceshelfvar
!
      implicit none
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
# endif
      CALL ad_iceshelf_vbc_tile (ng, tile,                              &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng), nnew(ng),                            &
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % ad_Hz,                              &
#endif
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % zice,                               &
     &                   GRID(ng) % f,                                  &
     &                   GRID(ng) % ad_z_r,                             &
     &                   GRID(ng) % ad_z_w,                             &
     &                   GRID(ng) % ad_zice,                            &
#ifdef ICESHELF_MORPH
     &                   ICESHELFVAR(ng) % iceshelf_draft,              &
     &                   ICESHELFVAR(ng) % ad_iceshelf_draft,           &
#endif
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
     &                   OCEAN(ng) % ad_t,                              &
#ifdef ICESHELF_3EQN_VBC
     &                   ICESHELFVAR(ng) % gammaT,                      &
     &                   ICESHELFVAR(ng) % gammaS,                      &
     &                   ICESHELFVAR(ng) % Tb,                          &
     &                   ICESHELFVAR(ng) % Sb,                          &
     &                   ICESHELFVAR(ng) % ad_gammaT,                   &
     &                   ICESHELFVAR(ng) % ad_gammaS,                   &
     &                   ICESHELFVAR(ng) % ad_Tb,                       &
     &                   ICESHELFVAR(ng) % ad_Sb,                       &
#endif
     &                   ICESHELFVAR(ng) % m,                           &
     &                   OCEAN(ng) % rho,                               &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   ICESHELFVAR(ng) % ad_m,                        &
     &                   OCEAN(ng) % ad_rho,                            &
     &                   FORCES(ng) % ad_sustr,                         &
     &                   FORCES(ng) % ad_svstr,                         &
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % ad_srflx,                         &
#   endif
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % ad_stflx                          &
     &                   )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
# endif
      RETURN
      END SUBROUTINE ad_iceshelf_vbc
!
!***********************************************************************
      SUBROUTINE ad_iceshelf_vbc_tile (ng, tile,                        &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,nnew,                               &
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
     &                         Hz,                                      &
     &                         ad_Hz,                                   &
#endif
#ifdef MASKING                 
     &                         rmask,                                   &
#endif
     &                         z_r, z_w,                                &
     &                         zice,f,                                  &
     &                         ad_z_r, ad_z_w,                          &
     &                         ad_zice,                                 &
# ifdef ICESHELF_MORPH
     &                         iceshelf_draft,                          &
     &                         ad_iceshelf_draft,                       &
# endif
     &                         u, v, t,                                 &
     &                         ad_u, ad_v, ad_t,                        &
#ifdef ICESHELF_3EQN_VBC
     &                         gammaT, gammaS,                          &
     &                         Tb, Sb,                                  &
     &                         ad_gammaT, ad_gammaS,                    &
     &                         ad_Tb, ad_Sb,                            &
#endif 
     &                         m, rho,                                  &
     &                         sustr, svstr,                            &
     &                         ad_m, ad_rho,                            &
     &                         ad_sustr, ad_svstr,                      &
#   ifdef SHORTWAVE
     &                         srflx,                                   &
     &                         ad_srflx,                                &
#   endif
     &                         stflx, ad_stflx)                         
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE mod_iceshelf
#ifdef ICESHELF_MORPH
      USE mod_iceshelfvar
#endif
#ifdef ICESHELF_DYNAM
      USE mod_iceshelfdyn
#endif
      USE ad_bc_2d_mod, ONLY: ad_bc_r2d_tile
      USE ad_bc_2d_mod, ONLY: ad_bc_u2d_tile
      USE ad_bc_2d_mod, ONLY: ad_bc_v2d_tile
      USE ad_exchange_2d_mod, ONLY : ad_exchange_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
# endif
# if defined ICESHELF_TEOS10
      USE iceshelf_vbc_mod, ONLY : ct_from_theta
      USE iceshelf_vbc_mod, ONLY : t_from_ct
# else
      USE iceshelf_vbc_mod, ONLY : dTemp
      USE iceshelf_vbc_mod, ONLY : thetaa
      USE iceshelf_vbc_mod, ONLY : potit
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
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_rho(LBi:,LBj:,:)
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
#endif
#ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(inout) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_zice(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
#ifdef ICESHELF_MORPH
      real(r8), intent(inout) :: iceshelf_draft(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_iceshelf_draft(LBi:,LBj:,:)
#endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#ifdef ICESHELF_3EQN_VBC
      real(r8), intent(out) :: gammaT(LBi:,LBj:)
      real(r8), intent(out) :: gammaS(LBi:,LBj:)
      real(r8), intent(inout) :: Tb(LBi:,LBj:)
      real(r8), intent(inout) :: Sb(LBi:,LBj:)
      real(r8), intent(inout) :: ad_gammaT(LBi:,LBj:)
      real(r8), intent(inout) :: ad_gammaS(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Tb(LBi:,LBj:)
      real(r8), intent(inout) :: ad_Sb(LBi:,LBj:)
#endif
      real(r8), intent(out) :: m(LBi:,LBj:)
      real(r8), intent(inout) :: ad_m(LBi:,LBj:)
# ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: ad_srflx(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_rho(LBi:UBi,LBj:UBj,N(ng))
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,N(ng))
#endif
#ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: ad_zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
#ifdef ICESHELF_MORPH
      real(r8), intent(in) :: iceshelf_draft(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ad_iceshelf_draft(LBi:UBi,LBj:UBj,2)
#endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#ifdef ICESHELF_3EQN_VBC
      real(r8), intent(out) :: gammaT(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: gammaS(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Tb(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: Sb(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_gammaT(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_gammaS(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Tb(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_Sb(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: m(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_m(LBi:UBi,LBj:UBj)
# ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_srflx(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
# endif
!
!  Local variable declarations.
!
!      integer :: IstrR, IENDR, JstrR, JendR, IstrU, JstrV
      integer :: i, j, itrc
      real(r8) :: Pradj,Scadj, Cdrt, cp_i
!      real(r8) :: a,b,c,Pr,Sc,Cd,small
!      real(r8) :: cp_w,ustar,visc,turb,TFb,rho_i
      real(r8) :: Sm,Tm,rhoi_on_rho0,ustar,TFb,turb
      real(r8) :: ad_Sm,ad_Tm,ad_rhoi_on_rho0,ad_ustar,ad_TFb,ad_turb
      real(r8) :: ad_RPres
#ifdef ICESHELF_VBC_HFPROPTOM
      real(r8) :: Ex1,Ex2,Ex3,Ex4,Ex5,Ex6,Ep1,Ep2,Ep3,Ep31,Ep4,Ep5
      real(r8) :: ad_Ex1,ad_Ex2,ad_Ex3,ad_Ex4,ad_Ex5,ad_Ex6,ad_Ep1,     &
     &                         ad_Ep2,ad_Ep3,ad_Ep31,ad_Ep4,ad_Ep5
      real(r8) :: Sr1,Sr2,Sf1,Sf2,Tf1,Tf2,rho_r
      real(r8) :: ad_Sr1,ad_Sr2,ad_Sf1,ad_Sf2,ad_Tf1,ad_Tf2,ad_rho_r
#else
      real(r8) :: mflag,TFi 
      real(r8) :: ad_TFi
#endif
      real(r8) :: cff1, cff2, cff3
      real(r8) :: ad_cff1, ad_cff2, ad_cff3, fac, adfac, adfac1,ad_fac
      real(r8) :: ad_Ti
#ifdef ICESHELF_TEOS10
      real(r8) :: ct, ad_ct
#endif
#  if (!defined BBL_MODEL) && defined UV_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_wrk
#  endif
#  include "set_bounds.h"

!
!  Clear adjoint variables
!
      ad_RPres=0.0_r8
#ifdef ICESHELF_VBC_HFPROPTOM
      ad_Ex1=0.0_r8
      ad_Ex2=0.0_r8
      ad_Ex3=0.0_r8
      ad_Ex4=0.0_r8
      ad_Ex5=0.0_r8
      ad_Ex6=0.0_r8
      ad_Ep1=0.0_r8
      ad_Ep2=0.0_r8
      ad_Ep3=0.0_r8
      ad_Ep31=0.0_r8
      ad_Ep4=0.0_r8
      ad_Ep5=0.0_r8
      ad_Sr1=0.0_r8
      ad_Sr2=0.0_r8
      ad_Sf1=0.0_r8
      ad_Sf2=0.0_r8
      ad_Tf1=0.0_r8
      ad_Tf2=0.0_r8
      ad_rho_r=0.0_r8
#endif
      ad_cff1=0.0_r8
      ad_cff2=0.0_r8
      ad_cff3=0.0_r8
      ad_Sm=0.0_r8
      ad_Tm=0.0_r8
      ad_rhoi_on_rho0=0.0_r8
      ad_ustar=0.0_r8
      ad_TFb=0.0_r8
      ad_turb=0.0_r8
      ad_TFi=0.0_r8
      ad_Ti=0.0_r8
#ifdef ICESHELF_TEOS10
      ad_ct=0.0_r8
#endif

!
!     Compute basic state stresses for use below.
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
#   elif defined UV_QDRAG
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
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
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

#   ifdef SHORTWAVE
      DO j=JstrR,JENDR
        DO i=IstrR,IENDR
          IF (zice(i,j).ne.0.0_r8) THEN
!           tl_srflx(i,j)=0.0_r8
            ad_srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif

#    ifdef DISTRIBUTE
!     CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!    &                    LBi, UBi, LBj, UBj,                           &
!    &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
!    &                    tl_stflx(:,:,itemp))
      CALL ad_mp_exchange2d (ng, tile, iNLM, 1,                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ad_stflx(:,:,itemp))
!     CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
!    &                    LBi, UBi, LBj, UBj,                           &
!    &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
!    &                    tl_stflx(:,:,isalt))
      CALL ad_mp_exchange2d (ng, tile, iNLM, 1,                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ad_stflx(:,:,isalt))
#    endif
!
!  Apply adjoinr gradient or periodic boundary conditions for the two fluxes
!
!     CALL bc_r2d_tile (ng, tile,                                       &
!    &                  LBi, UBi, LBj, UBj,                             &
!    &                  tl_stflx(:,:,itemp))
      CALL ad_bc_r2d_tile (ng, tile,                                    &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ad_stflx(:,:,itemp))
!     CALL bc_r2d_tile (ng, tile,                                       &
!    &                  LBi, UBi, LBj, UBj,                             &
!    &                  tl_stflx(:,:,isalt))
      CALL ad_bc_r2d_tile (ng, tile,                                    &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  ad_stflx(:,:,isalt))

!-----------------------------------------------------------------------
! Store old ice shelf thickness.
!-----------------------------------------------------------------------
!
# if defined ICESHELF_MORPH
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
!         CALL exchange_r2d_tile (ng, tile,                             &
!    &                            LBi, UBi, LBj, UBj,                   &
!    &                            tl_iceshelf_draft(:,:,nnew))
          CALL ad_exchange_r2d_tile (ng, tile,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ad_iceshelf_draft(:,:,nnew))
        END IF

      DO j=JstrR,JendR
        DO i=IstrR,IendR
!         tl_iceshelf_draft(i,j,nnew)=tl_zice(i,j)+tl_m(i,j)*dt(ng)
          ad_zice(i,j)=ad_zice(i,j)+ad_iceshelf_draft(i,j,nnew)
          ad_m(i,j)=ad_m(i,j)+dt(ng)*ad_iceshelf_draft(i,j,nnew)
          ad_iceshelf_draft(i,j,nnew)=0.0_r8
        END DO
      END DO

# endif

!-----------------------------------------------------------------------
!  If ice shelf cavities, zero out for now the surface tracer flux
!  over the ice.
!-----------------------------------------------------------------------
!
#    ifdef ANA_SEAICE
#     ifdef SEAICE_CLIMA
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
#     elif defined SEAICE_WINTER
        sfcTemp = Tmin
        sfcSalt = SaltMax
#     endif
#    endif

#    ifdef ICESHELF_3EQN_VBC
      Pradj = 12.5_r8*(Pr**(2.0_r8/3.0_r8))
      Scadj = 12.5_r8*(Sc**(2.0_r8/3.0_r8))
      Cdrt =  sqrt(Cd)
      cp_i =  152.5_r8+7.122_r8*(273.15_r8+Ti)
      rhoi_on_rho0 = rho_i/rho0
#   endif

      DO j=Jstr,JEND
        DO i=Istr,IEND
# ifdef SALINITY
          Sm=MAX(0.0_r8,t(i,j,N(ng),nrhs,isalt))
# else
          Sm=0.0_r8
# endif

          IF (zice(i,j).ne.0.0_r8) THEN
#if defined ICESHELF_2EQN_VBC
          TFb = a*Sm+b+c*zice(i,j)
# ifdef ICESHELF_TEOS10
          ct = ct_from_theta(Sm,t(i,j,N(ng),nrhs,itemp))
          Tm = t_from_ct(Sm,ct,-zice(i,j))
# else
         CALL potit(Sm,t(i,j,N(ng),nrhs,itemp),-zice(i,j),0.0_r8,Tm,i,j)
# endif
!
            stflx(i,j,itemp)=gamma*(TFb-Tm)
            stflx(i,j,isalt)=Cp*stflx(i,j,itemp)*refSalt/L
            m(i,j) = stflx(i,j,isalt)/Sm

!           tl_m(i,j) = tl_stflx(i,j,isalt)/Sm-tl_Sm*m(i,j)/Sm

            ad_stflx(i,j,isalt)=ad_stflx(i,j,isalt)+ad_m(i,j)/Sm
            ad_Sm=ad_Sm-m(i,j)*ad_m(i,j)/Sm
            ad_m(i,j)=0.0_r8

!           tl_stflx(i,j,isalt)=Cp*tl_stflx(i,j,itemp)*refSalt/L
            ad_stflx(i,j,itemp)=ad_stflx(i,j,itemp)+                    &
     &                          Cp*ad_stflx(i,j,isalt)*refSalt/L
            ad_stflx(i,j,isalt)=0.0_r8

!           tl_stflx(i,j,itemp)=gamma*(tl_TFb-tl_Tm)
            ad_TFb=ad_TFb+gamma*ad_stflx(i,j,itemp)
            ad_Tm=ad_Tm-gamma*ad_stflx(i,j,itemp)
            ad_stflx(i,j,itemp)=0.0_r8

# ifdef ICESHELF_TEOS10
!         tl_Tm = tl_t_from_ct(Sm,tl_Sm,ct,tl_ct,-zice(i,j),            &
!    &                                           -tl_zice(i,j))
          ad_Sm=ad_Sm+ad_t_from_ct(1,ad_Tm,Sm,ct,-zice(i,j))
          ad_ct=ad_ct+ad_t_from_ct(1,ad_Tm,Sm,ct,-zice(i,j))
          ad_zice(i,j)=ad_zice(i,j)+                                    &
     &                         ad_t_from_ct(3,ad_Tm,Sm,ct,-zice(i,j))
          ad_Tm=0.0_r8

!         tl_ct = tl_ct_from_theta(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),    &
!    &                             tl_t(i,j,N(ng),nrhs,itemp))

          ad_Sm=ad_Sm+                                                  &
     &              ad_ct_from_theta(1,ad_ct,Sm,t(i,j,N(ng),nrhs,itemp))
          ad_t(i,j,N(ng),nrhs,itemp)=ad_t(i,j,N(ng),nrhs,itemp)+        &
     &              ad_ct_from_theta(2,ad_ct,Sm,t(i,j,N(ng),nrhs,itemp))
          ad_ct=0.0_r8
# else
!        CALL tl_potit(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),               &
!    &                          tl_t(i,j,N(ng),nrhs,itemp),            &
!    &             -zice(i,j),-tl_zice(i,j),0.0_r8,tl_RPres,Tm,tl_Tm,i,j)

         fac=-ad_zice(i,j)
         CALL ad_potit(Sm,ad_Sm,t(i,j,N(ng),nrhs,itemp),               &
     &                          ad_t(i,j,N(ng),nrhs,itemp),            &
     &           -zice(i,j),fac,0.0_r8,ad_RPres,Tm,ad_Tm,i,j)
!AMM     &           -zice(i,j),-ad_zice(i,j),0.0_r8,ad_RPres,Tm,ad_Tm,i,j)
         ad_zice(i,j)=fac
# endif

!         tl_TFb = a*tl_Sm+c*tl_zice(i,j)
          ad_Sm=ad_Sm+a*ad_TFb
          ad_zice(i,j)=ad_zice(i,j)+c*ad_TFb
          ad_TFb=0.0_r8



#elif defined ICESHELF_3EQN_VBC


# ifdef ICESHELF_VBC_HFPROPTOM
!
! Calculates density RHOW in the boundary layer, interface
! pressure PG[dbar], and solving a quadratic equation for the
! interface salinity (SB) to determine the melTemp_insitug/freezingrate
! (seta).
          rho_r= rho_i/(1000.0_r8+rho(i,j,N(ng)))
          Ep1  = cp_w*gammaT(i,j)
          Ep2  = cp_i*gammaS(i,j)
          Ep3  = L*gammaS(i,j)
#  ifdef ICESHELF2D_TOY
          Ep31 = 0.0_r8
#  else
          Ep31 = -rho_r*cp_i*dt_i/zice(i,j)
#  endif
          Ep4  = b+c*z_w(i,j,N(ng))
          Ep5  = gammaS(i,j)/rho_r
!
! Temperature gradient within the ice ONLY changes with melTemp_insitug
          TFb = a*Sm+Ep4
          IF (Tm.lt.TFb) THEN
            Ex1 = a*(Ep1+Ep31)
            Ex2 = Ep1*(Tm-Ep4)+Ep3+Ep31*(Ti-Ep4)
            Ex3 = Ep3*Sm
            Ex6 = 0.5_r8
          ELSE
! negative heat flux term in the ice (due to -kappa/D)
            Ex1 = a*(Ep1-Ep2)
            Ex2 = Ep1*(Ep4-Tm)+Ep2*(Ti+a*Sm-Ep4)-Ep3
            Ex3 = Sm*(Ep2*(Ep4-Ti)+Ep3)
            Ex6 = -0.5r_r8
          ENDIF
          Ex4 = Ex2/Ex1
          Ex5 = Ex3/Ex1
          Sr1 = 0.25_r8*Ex4*Ex4-Ex5
          Sr2 = Ex6*Ex4
          Sf1 = Sr2+sqrt(Sr1)
          Tf1 = a*Sf1+Ep4
          Sf2 = Sr2-sqrt(Sr1)
          Tf2 = a*Sf2+Ep4
!
! Salinities < 0 psu are NOT defined, i.e....
         IF (Sf1.gt.0.) THEN
           Tb(i,j) = Tf1
           Sb(i,j) = Sf1
         ELSE
           Tb(i,j) = Tf2
           Sb(i,j) = Sf2
         endif
         m(i,j) = -Ep5*(1.0-Sm/Sb(i,j))

!
! Calculates fluxes:
!        stflx(i,j,isalt)=(gammaS(i,j)-rho_r*m(i,j))*(Sb(i,j)-Sm)
!        tl_stflx(i,j,isalt)=(tl_gammaS(i,j)-tl_rho_r*m(i,j)-              &
!    &                           rho_r*tl_m(i,j))*(Sb(i,j)-Sm)+            &
!    &                    (gammaS(i,j)-rho_r*m(i,j))*(tl_Sb(i,j)-tl_Sm)

         ad_gammaS(i,j)=ad_gammaS(i,j)+ad_stflx(i,j,isalt)*(Sb(i,j)-Sm)
         ad_rho_r=ad_rho_r-m(i,j)*ad_stflx(i,j,isalt)*(Sb(i,j)-Sm)
         ad_m(i,j)=ad_m(i,j)-rho_r*ad_stflx(i,j,isalt)*(Sb(i,j)-Sm)
         ad_Sb(i,j)=ad_Sb(i,j)+(gammaS(i,j)-rho_r*m(i,j))*              &
     &                               ad_stflx(i,j,isalt)
         ad_Sm=ad_Sm-(gammaS(i,j)-rho_r*m(i,j))*ad_stflx(i,j,isalt)
         ad_stflx(i,j,isalt)=0.0_r8

!        stflx(i,j,itemp)=(gammaT(i,j)-rho_r*m(i,j))*(Tb(i,j)-Tm)
!        tl_stflx(i,j,itemp)=(tl_gammaT(i,j)-tl_rho_r*m(i,j)-              &
!    &                           rho_r*tl_m(i,j))*(Tb(i,j)-Tm)+            &
!    &                    (gammaT(i,j)-rho_r*m(i,j))*(tl_Tb(i,j)-tl_Tm)

!        tl_m(i,j) = -tl_Ep5*(1.0-Sm/Sb(i,j))-Ep5*(-tl_Sm/Sb(i,j)+         &
!    &                 tl_Sb(i,j)*Sm/(Sb(i,j)*Sb(i,j)))

         ad_Ep5=ad_Ep5-ad_m(i,j)*(1.0-Sm/Sb(i,j))
         ad_Sm=ad_Sm+Ep5*ad_m(i,j)/Sb(i,j)
         ad_Sb(i,j)=ad_Sb(i,j)-Ep5*Sm*ad_m(i,j)/(Sb(i,j)*Sb(i,j))
         ad_m(i,j)=0.0_r8

!
! Salinities < 0 psu are NOT defined, i.e....
         IF (Sf1.gt.0.) THEN
!          tl_Tb(i,j) = tl_Tf1
           ad_Tf1=ad_Tf1+ad_Tb(i,j)
           ad_Tb(i,j)=0.0_r8
!          tl_Sb(i,j) = tl_Sf1
           ad_Sf1=ad_Sf1+ad_Sb(i,j)
           ad_Sb(i,j)=0.0_r8
         ELSE
!          tl_Tb(i,j) = tl_Tf2
           ad_Tf2=ad_Tf2+ad_Tb(i,j)
           ad_Tb(i,j)=0.0_r8
!          tl_Sb(i,j) = tl_Sf2
           ad_Sf2=ad_Sf2+ad_Sb(i,j)
           ad_Sb(i,j)=0.0_r8
         endif

!         tl_Tf2 = a*tl_Sf2+tl_Ep4
          ad_Sf2=ad_Sf2+a*ad_Tf2
          ad_Ep4=ad_Ep4+ad_Tf2
          ad_Tf2=0.0_r8

!         tl_Sf2 = tl_Sr2-0.5_r8*tl_Sr1/sqrt(Sr1)
          ad_Sr2=ad_Sr2+ad_Sf2
          ad_Sr1=ad_Sr1-0.5_r8*ad_Sf2/sqrt(Sr1)
          ad_Sf2=0.0_r8

!         tl_Tf1 = a*tl_Sf1+tl_Ep4
          ad_Sf1=ad_Sf1+a*ad_Tf1
          ad_Ep4=ad_Ep4+ad_Tf1
          ad_Tf1=0.0_r8

!         tl_Sf1 = tl_Sr2+0.5_r8*tl_Sr1/sqrt(Sr1)
          ad_Sr2=ad_Sr2+ad_Sf1
          ad_Sr1=ad_Sr1+0.5_r8*0.5_r8/sqrt(Sr1)
          ad_Sf1=0.0_r8

!         tl_Sr2 = tl_Ex6*Ex4+Ex6*tl_Ex4
          ad_Ex6=ad_Ex6+ad_Sr2*Ex4
          ad_Ex4=ad_Ex4+Ex6*ad_Sr2
          ad_Sr2=0.0_r8

!         tl_Sr1 = 0.5_r8*tl_Ex4*Ex4-tl_Ex5
          ad_Ex4=ad_Ex4+0.5_r8*Ex4*ad_Sr1
          ad_Ex5=ad_Ex5-ad_Sr1
          ad_Sr1=0.0_r8

!         tl_Ex5 = tl_Ex3/Ex1-tl_Ex1*Ex5/Ex1
          ad_Ex3=ad_Ex3+ad_Ex5/Ex1
          ad_Ex1=ad_Ex1-ad_Ex5*Ex5/Ex1
          ad_Ex5=0.0_r8

!         tl_Ex4 = tl_Ex2/Ex1-tl_Ex1*Ex4/Ex1
          ad_Ex2=ad_Ex2+ad_Ex4/Ex1
          ad_Ex1=ad_Ex1-ad_Ex4*Ex4/Ex1
          ad_Ex4=0.0_r8

          IF (Tm.lt.TFb) THEN
!           tl_Ex6 = 0.0_r8
            ad_Ex6 = 0.0_r8

!           tl_Ex3 = tl_Ep3*Sm+Ep3*tl_Sm
            ad_Ep3=ad_Ep3+ad_Ex3*Sm
            ad_Sm=ad_Sm+Ep3*ad_Ex3
            ad_Ex3=0.0_r8

!           tl_Ex2 = tl_Ep1*(Tm-Ep4)+Ep1*(tl_Tm-tl_Ep4)+tl_Ep3+           &
!    &               tl_Ep31*(Ti-Ep4)-Ep31*tl_Ep4
            ad_Ep1=ad_Ep1+ad_Ex2*(Tm-Ep4)
            ad_Tm=ad_Tm+Ep1*ad_Ex2
            ad_Ep4=ad_Ep4-Ep1*ad_Ex2
            ad_Ep3=ad_Ep3+ad_Ex2
            ad_Ep31=ad_Ep31+ad_Ex2*(Ti-Ep4)
            ad_Ep4=ad_Ep4-Ep31*ad_Ex2
            ad_Ex2=0.0_r8

!           tl_Ex1 = a*(tl_Ep1+tl_Ep31)
            ad_Ep1=ad_Ep1+a*ad_Ex1
            ad_Ep31=ad_Ep31+a*ad_Ex1
            ad_Ex1=0.0_r8
          ELSE

!           tl_Ex6 = 0.0_r8
            ad_Ex6 = 0.0_r8

!           tl_Ex3 = tl_Sm*(Ep2*(Ep4-Ti)+Ep3)+Sm*(tl_Ep2*(Ep4-Ti)+         &
!    &                                              Ep2*tl_Ep4+tl_Ep3)

            ad_Sm=ad_Sm+ad_Ex3*(Ep2*(Ep4-Ti)+Ep3)
            ad_Ep2=ad_Ep2+Sm*(Ep4-Ti)*ad_Ex3
            ad_Ep4=ad_Ep4+Sm*Ep2*ad_Ex3
            ad_Ep3=ad_Ep3+Sm*ad_Ex3
            ad_Ex3=0.0_r8

!           tl_Ex2= tl_Ep1*(Ep4-Tm)+Ep1*(tl_Ep4-tl_Tm)+tl_Ep2*(Ti+a*Sm-Ep4)&
!    &             +Ep2*(a*tl_Sm-tl_Ep4)-tl_Ep3

            ad_Ep1=ad_Ep1+ad_Ex2*(Ep4-Tm)
            ad_Ep4=ad_Ep4+Ep1*ad_Ex2
            ad_Tm=ad_Tm-Ep1*ad_Ex2
            ad_Ep2=ad_Ep2+ad_Ex2*(Ti+a*Sm-Ep4)
            ad_Sm=ad_Sm+Ep2*a*ad_Ex2
            ad_Ep4=ad_Ep4-Ep2*ad_Ex2
            ad_Ep3=ad_Ep3-ad_Ex2
            ad_Ex2=0.0_r8

!           tl_Ex1 = a*(tl_Ep1-tl_Ep2)
            ad_Ep1=ad_Ep1+a*ad_Ex1
            ad_Ep2=ad_Ep2-a*ad_Ex1
            ad_Ex1=0.0_r8

          ENDIF

!         tl_TFb = a*tl_Sm+tl_Ep4
          ad_Sm=ad_Sm+a*ad_TFb
          ad_Ep4=ad_Ep4+ad_TFb
          ad_TFb=0.0_r8

!         tl_Ep5  = tl_gammaS(i,j)/rho_r-tl_rho_r*Ep5/rho_r
          ad_gammaS(i,j)=ad_gammaS(i,j)+ad_Ep5/rho_r
          ad_rho_r=ad_rho_r-ad_Ep5*Ep5/rho_r
          ad_Ep5=0.0_r8

!         tl_Ep4  = c*tl_z_w(i,j,N(ng))
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))+c*ad_Ep4

#  ifdef ICESHELF2D_TOY
!         tl_Ep31 = 0.0_r8
          ad_Ep31 = 0.0_r8
#  else
          Ep31 = -rho_r*cp_i*dt_i/zice(i,j)
!         tl_Ep31 = -tl_rho_r*cp_i*dt_i/zice(i,j)+                       &
!                             tl_zice(i,j)*Ep31/zice(i,j)
          ad_rho_r=ad_rho_r-ad_Ep31*cp_i*dt_i/zice(i,j)
          ad_zice(i,j)=ad_zice(i,j)+ad_Ep31*Ep31/zice(i,j)
#  endif


!         tl_Ep3  = L*tl_gammaS(i,j)
          ad_gammaS(i,j)=ad_gammaS(i,j)+L*ad_Ep3
          ad_Ep3=0.0_r8

!         tl_Ep2  = cp_i*tl_gammaS(i,j)
          ad_gammaS(i,j)=ad_gammaS(i,j)+cp_i*ad_Ep2
          ad_Ep2=0.0_r8

!         tl_Ep1  = cp_w*tl_gammaT(i,j)
          ad_gammaT(i,j)=ad_gammaT(i,j)+cp_w*ad_Ep1
          ad_Ep1=0.0_r8

!         tl_rho_r= -tl_rho(i,j,N(ng))*rho_r/(1000.0_r8+rho(i,j,N(ng)))
          ad_rho(i,j,N(ng))=ad_rho(i,j,N(ng))-ad_rho_r*rho_r/           &
                      (1000.0_r8+rho(i,j,N(ng)))
          ad_rho_r=0.0_r8
# else
            TFb = a*Sm + b + c*zice(i,j)
            mflag= 0.5*(1 + SIGN(1.0_r8,Tm-TFb))
! No TL required for mflag. tl_mflag=0 by definition. See Giering TAMC manual, p32, Table 8.
! The derivative of a SIGN function wrt 2nd argument of zero.

            TFi = (1 - mflag)*a*Si + b + c*zice(i,j)

! Calculate coefficents in quadratic to be solved:
           cff1 = L/cp_w + mflag*(cp_i/cp_w)*(TFi-Ti)
           cff2 = gammaS(i,j)*(L/cp_w + mflag*(cp_i/cp_w)*(TFb-Ti))     &
     &          + gammaT(i,j)*(TFi-Tm)
           cff3 = gammaS(i,j)*gammaT(i,j)*(TFb - Tm)

! Calculate melt rate:
!          m(i,j) = -(cff2 - SQRT(cff2*cff2 - 4.0_r8*cff1*cff3))/          &
!    &                                                (2.0_r8*cff1)
           fac=SQRT(cff2*cff2 - 4.0_r8*cff1*cff3)
           m(i,j) = -(cff2 - fac)/(2.0_r8*cff1)
! Calculate basal temperature and salinity:
           Tb(i,j) = (gammaT(i,j)*Tm+mflag*(cp_i/cp_w)*m(i,j)*Ti        &
     &             - (L/cp_w)*m(i,j))                                   &
     &             /(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))


!        tl_stflx(i,j,isalt)=(tl_gammaS(i,j)+rhoi_on_rho0*tl_m(i,j))*      &
!    &                           (Sb(i,j)-Sm)+                             &
!    &                (gammaS(i,j)+rhoi_on_rho0*m(i,j))*(tl_Sb(i,j)-tl_Sm)

         ad_gammaS(i,j)=ad_gammaS(i,j)+ad_stflx(i,j,isalt)*(Sb(i,j)-Sm)
         ad_m(i,j)=ad_m(i,j)+                                           &
     &             rhoi_on_rho0*ad_stflx(i,j,isalt)*(Sb(i,j)-Sm)
         ad_Sb(i,j)=ad_Sb(i,j)+(gammaS(i,j)+rhoi_on_rho0*m(i,j))*       &
     &                                   ad_stflx(i,j,isalt)
         ad_Sm=ad_Sm-(gammaS(i,j)+rhoi_on_rho0*m(i,j))*                 &
     &                            ad_stflx(i,j,isalt)
         ad_stflx(i,j,isalt)=0.0_r8

!        tl_stflx(i,j,itemp)=(tl_gammaT(i,j)+rhoi_on_rho0*tl_m(i,j))*      &
!    &                           (Tb(i,j)-Tm)+                             &
!    &                (gammaT(i,j)+rhoi_on_rho0*m(i,j))*(tl_Tb(i,j)-tl_Tm)

         ad_gammaT(i,j)=ad_gammaT(i,j)+ad_stflx(i,j,itemp)*(Tb(i,j)-Tm)
         ad_m(i,j)=ad_m(i,j)+rhoi_on_rho0*ad_stflx(i,j,itemp)*          &
     &                                    (Tb(i,j)-Tm)
         ad_Tb(i,j)=ad_Tb(i,j)+(gammaT(i,j)+rhoi_on_rho0*m(i,j))*       &
     &                                    ad_stflx(i,j,itemp)
         ad_Tm=ad_Tm-(gammaT(i,j)+rhoi_on_rho0*m(i,j))*                 &
     &                            ad_stflx(i,j,itemp)
         ad_stflx(i,j,itemp)=0.0_r8

!        tl_Sb(i,j) = (tl_Tb(i,j) - c*tl_zice(i,j))/a

         ad_Tb(i,j)=ad_Tb(i,j)+ad_Sb(i,j)/a
         ad_zice(i,j)=ad_zice(i,j)-c*ad_Sb(i,j)/a
         ad_Sb(i,j)=0.0_r8

! Calculate basal temperature and salinity:
!          tl_Tb(i,j) = (tl_gammaT(i,j)*Tm+gammaT(i,j)*tl_Tm+              &
!    &                   mflag*(cp_i/cp_w)*(tl_m(i,j)*Ti+m(i,j)*tl_Ti)-    &
!    &                   (L/cp_w)*tl_m(i,j))                               &
!    &             /(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))               &
!    &           -(tl_gammaT(i,j) + mflag*(cp_i/cp_w)*tl_m(i,j))*Tb(i,j)   &
!    &             /(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))

         adfac1=ad_Tb(i,j)/(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))
         ad_gammaT(i,j)=ad_gammaT(i,j)+Tm*adfac1
         ad_Tm=ad_Tm+gammaT(i,j)*adfac1
         ad_m(i,j)=ad_m(i,j)+mflag*(cp_i/cp_w)*Ti*adfac1
         ad_Ti=ad_Ti+mflag*(cp_i/cp_w)*m(i,j)*adfac1
         ad_m(i,j)=ad_m(i,j)-(L/cp_w)*adfac1
         ad_gammaT(i,j)=ad_gammaT(i,j)-Tb(i,j)*adfac1
         ad_m(i,j)=ad_m(i,j)-mflag*(cp_i/cp_w)*Tb(i,j)*adfac1
         ad_Tb(i,j)=0.0_r8

! Calculate melt rate:
!          tl_m(i,j) = -(tl_cff2 -tl_fac)/(2.0_r8*cff1)+                   &
!    &                              2.0_r8*tl_cff1*m(i,j)/cff1

           ad_cff2=ad_cff2-ad_m(i,j)/(2.0_r8*cff1)
           ad_fac=ad_fac+ad_m(i,j)/(2.0_r8*cff1)
           ad_cff1=ad_cff1+2.0_r8*ad_m(i,j)*m(i,j)/cff1
           ad_m(i,j)=0.0_r8

!          tl_fac=0.5_r8*(2.0_r8*tl_cff2*cff2 -                          &
!    &                   4.0_r8*(tl_cff1*cff3+cff1*tl_cff3))/fac

           ad_cff2=ad_cff2+0.5_r8*2.0_r8*cff2*ad_fac/fac
           ad_cff1=ad_cff1-0.5_r8*4.0_r8*cff3*ad_fac/fac
           ad_cff3=ad_cff3-0.5_r8*4.0_r8*cff1*ad_fac/fac
           ad_fac=0.0_r8

! Calculate coefficents in quadratic to be solved:
!          tl_cff3 = gammaS(i,j)*gammaT(i,j)*(tl_TFb - tl_Tm)+             &
!                    tl_gammaS(i,j)*gammaT(i,j)*(TFb - Tm)+                &
!                    gammaS(i,j)*tl_gammaT(i,j)*(TFb - Tm)

           ad_TFb=ad_TFb+gammaS(i,j)*gammaT(i,j)*ad_cff3
           ad_Tm=ad_Tm-gammaS(i,j)*gammaT(i,j)*ad_cff3
           ad_gammaS(i,j)=ad_gammaS(i,j)+ad_cff3*gammaT(i,j)*(TFb - Tm)
           ad_gammaT(i,j)=ad_gammaT(i,j)+ad_cff3*gammaS(i,j)*(TFb - Tm)
           ad_cff3=0.0_r8

!          tl_cff2 = gammaS(i,j)*mflag*(cp_i/cp_w)*tl_TFb                  &
!    &          + gammaT(i,j)*(tl_TFi-tl_Tm)+tl_gammaS(i,j)*(L/cp_w +      &
!    &                                  mflag*(cp_i/cp_w)*(TFb-Ti))        &
!    &          +tl_gammaT(i,j)*(TFi-Tm)

           ad_TFb=ad_TFb+gammaS(i,j)*mflag*(cp_i/cp_w)*ad_cff2
           ad_TFi=ad_TFi+gammaT(i,j)*ad_cff2
           ad_Tm=ad_Tm-gammaT(i,j)*ad_cff2
           ad_gammaS(i,j)=ad_gammaS(i,j)+ad_cff2*(L/cp_w +              &
     &                                  mflag*(cp_i/cp_w)*(TFb-Ti))
           ad_cff2=0.0_r8

!          tl_cff1 = mflag*(cp_i/cp_w)*tl_TFi

           ad_TFi=ad_TFi+mflag*(cp_i/cp_w)*ad_cff1
           ad_cff1=0.0_r8

!          tl_TFi = c*tl_zice(i,j)

           ad_zice(i,j)=ad_zice(i,j)+c*ad_TFi
           ad_TFi=0.0_r8

!          tl_TFb = a*tl_Sm  + c*tl_zice(i,j)

           ad_zice(i,j)=ad_zice(i,j)+c*ad_TFb
           ad_Sm=ad_Sm+a*ad_TFb
           ad_TFb=0.0_r8

# endif

!
!    Calculate exchange coefficients. Note that this assumes log profile 
!            ustar = SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+
!            &
!     &                        (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
            cff1=(0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2
            cff2=(0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2
            cff3=SQRT(cff1+cff2)
            ustar = SQRT(cff3)
#    ifdef MASKING
            ustar = ustar*rmask(i,j)
#    endif
           if(ustar.lt.small)then
              ustar = small
           end if

! AMM1
#    ifdef GAMMAS_KY1972
!           turb = 2.12_r8*LOG((ustar*Hz(i,j,N(ng)))/visc)
            cff1=ustar*Hz(i,j,N(ng))/visc
            turb = 2.12_r8*LOG(cff1)
            gammaT(i,j) = ustar/(turb + Pradj - 9.0_r8)
            gammaS(i,j) = ustar/(turb + Scadj - 9.0_r8)

!           tl_gammaS(i,j) = tl_ustar/(turb + Scadj - 9.0_r8)-          &
!    &                    tl_turb*gammaS(i,j)/(turb + Scadj - 9.0_r8)

            ad_ustar=ad_ustar+ad_gammaS(i,j)/(turb + Scadj - 9.0_r8)
            ad_turb=ad_turb-ad_gammaS(i,j)*gammaS(i,j)/                 &
     &                         (turb + Scadj - 9.0_r8)
            ad_gammaS(i,j)=0.0_r8

!           tl_gammaT(i,j) = tl_ustar/(turb + Pradj - 9.0_r8)-          &
!    &                    tl_turb*gammaT(i,j)/(turb + Pradj - 9.0_r8)

            ad_ustar=ad_ustar+ad_gammaT(i,j)/(turb + Pradj - 9.0_r8)
            ad_turb=ad_turb-ad_gammaT(i,j)*gammaT(i,j)/                 &
     &                         (turb + Pradj - 9.0_r8)
            ad_gammaT(i,j)=0.0_r8

!           tl_turb = 2.12_r8*tl_cff1/cff1
            ad_cff1=ad_cff1+2.12_r8*ad_turb/cff1
            ad_turb=0.0_r8

!           tl_cff1=(tl_ustar*Hz(i,j,N(ng))+ustar*tl_Hz(i,j,N(ng)))/visc
            ad_ustar=ad_ustar+Hz(i,j,N(ng)*ad_cff1/visc
            ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+ustar*ad_cff1/visc
            ad_cff1=0.0_r8

#    elif defined GAMMAS_Mc1987
! Uses simplified version of McPhee 1987
            IF (ustar.gt.small.and.ABS(f(i,j)).gt.1.0E-8) THEN
!         turb = 2.5_r8*LOG(5300.0_r8*ustar*ustar/ABS(f(i,j))) + 7.12_r8
              cff1=5300.0_r8*ustar*ustar/ABS(f(i,j))
              turb = 2.5_r8*LOG(cff1) + 7.12_r8
            ELSE
              turb = 0.0_r8
            END IF

            gammaT(i,j) = ustar/(turb + Pradj - 6.0_r8)
            gammaS(i,j) = ustar/(turb + Scadj - 6.0_r8)

!           tl_gammaS(i,j) = tl_ustar/(turb + Scadj - 6.0_r8)-           &
!    &                       tl_turb*gammaS(i,j)/(turb + Scadj - 6.0_r8)
            ad_ustar=ad_ustar+ad_gammaS(i,j)/(turb + Scadj - 6.0_r8)
            ad_turb=ad_turb-ad_gammaS(i,j)*gammaS(i,j)/                 &
     &                               (turb + Scadj - 6.0_r8)
            ad_gammaS(i,j)=0.0_r8

!           tl_gammaT(i,j) = tl_ustar/(turb + Pradj - 6.0_r8)-           &
!    &                       tl_turb*gammaT(i,j)/(turb + Pradj - 6.0_r8)

            ad_ustar=ad_ustar+ad_gammaT(i,j)/(turb + Pradj - 6.0_r8)
            ad_turb=ad_turb-ad_gammaT(i,j)*gammaT(i,j)/                 &
     &                               (turb + Pradj - 6.0_r8)
            ad_gammaT(i,j)=0.0_r8

            IF (ustar.gt.small.and.ABS(f(i,j)).gt.1.0E-8) THEN
!             tl_turb = 2.5_r8*tl_cff1/cff1
              ad_cff1=ad_cff1+2.5_r8*ad_turb/cff1
              ad_turb=0.0_r8
             
!             tl_cff1=2.0_r8*5300.0_r8*tl_ustar*ustar/ABS(f(i,j))
              ad_ustar=ad_ustar+2.0_r8*5300.0_r8*ustar*ad_cff1/         &
     &                          ABS(f(i,j))
              ad_cff1=0.0_r8

            ELSE
!             tl_turb = 0.0_r8
              ad_turb = 0.0_r8
            END IF
#    else
!           tl_gammaT(i,j) = 0.0_r8
            ad_gammaT(i,j) = 0.0_r8

!           tl_gammaS(i,j) = 0.0_r8
            ad_gammaS(i,j) = 0.0_r8

#    endif

           if(ustar.lt.small)then
!             tl_ustar = 0.0_r8
              ad_ustar = 0.0_r8
           end if

#    ifdef MASKING
!           tl_ustar = tl_ustar*rmask(i,j)
            ad_ustar = ad_ustar*rmask(i,j)
#    endif

!           tl_ustar = 0.5_r8*tl_cff3/ustar
            ad_cff3=ad_cff3+0.5_r8*ad_ustar/ustar

!           tl_cff3=0.5_r8*(tl_cff1+tl_cff2)/cff3
            ad_cff1=ad_cff1+0.5_r8*ad_cff3/cff3
            ad_cff2=ad_cff2+0.5_r8*ad_cff3/cff3
            ad_cff3=0.0_r8

!           tl_cff2=(svstr(i,j)+svstr(i,j+1))*                          &
!    &                        (tl_svstr(i,j)+tl_svstr(i,j+1))
            ad_svstr(i,j)=ad_svstr(i,j)+(svstr(i,j)+svstr(i,j+1))*      &
     &                                                ad_cff2
            ad_svstr(i,j+1)=ad_svstr(i,j+1)+(svstr(i,j)+svstr(i,j+1))*  &
     &                                                ad_cff2
            ad_cff2=0.0_r8

!           tl_cff1=(sustr(i,j)+sustr(i+1,j))*                          &
!    &                       (tl_sustr(i,j)+tl_sustr(i+1,j))
            ad_sustr(i,j)=ad_sustr(i,j)+(sustr(i,j)+sustr(i+1,j))*      &
     &                                               ad_cff1
            ad_sustr(i+1,j)=ad_sustr(i+1,j)+(sustr(i,j)+sustr(i+1,j))*  &
     &                                               ad_cff1
            ad_cff1=0.0_r8

! Get the insitu salinity and temperature
# ifdef ICESHELF_TEOS10
          ct = ct_from_theta(Sm,t(i,j,N(ng),nrhs,itemp))

!         tl_Tm = tl_t_from_ct(Sm,tl_Sm,ct,tl_ct,-zice(i,j),            &
!    &                                           -tl_zice(i,j))
          ad_Sm=ad_Sm+ad_t_from_ct(1,ad_Tm,Sm,ct,-zice(i,j))
          ad_ct=ad_ct+ad_t_from_ct(2,ad_Tm,Sm,ct,-zice(i,j))
          ad_zice(i,j)=ad_zice(i,j)+                                    &
     &                          ad_t_from_ct(3,ad_Tm,Sm,ct,-zice(i,j))
          ad_Tm=0.0_r8

!         tl_ct = tl_ct_from_theta(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),    &
!    &                                      tl_t(i,j,N(ng),nrhs,itemp))
          ad_Sm=ad_Sm+                                                  &
     &            ad_ct_from_theta(1,ad_ct,Sm,t(i,j,N(ng),nrhs,itemp))
          ad_t(i,j,N(ng),nrhs,itemp)=ad_t(i,j,N(ng),nrhs,itemp)+        &
     &            ad_ct_from_theta(2,ad_ct,Sm,t(i,j,N(ng),nrhs,itemp))
          ad_ct=0.0_r8
# else
!        CALL tl_potit(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),                &
!    &                 tl_t(i,j,N(ng),nrhs,itemp),-zice(i,j),           &
!    &                 -tl_zice(i,j),0.0_r8,Tm,tl_Tm,i,j)
         fac=-ad_zice(i,j)
         CALL ad_potit(Sm,ad_Sm,t(i,j,N(ng),nrhs,itemp),                &
     &                 ad_t(i,j,N(ng),nrhs,itemp),-zice(i,j),           &
     &                 fac,0.0_r8,ad_RPres,Tm,ad_Tm,i,j)
!AMM     &                 -ad_zice(i,j),0.0_r8,ad_RPres,Tm,ad_Tm,i,j)
         ad_zice(i,j)=fac
# endif

! AMM2
#else
      DO itrc=1,NT(ng)
!             tl_stflx(i,j,itrc)=0.0_r8
              ad_stflx(i,j,itrc)=0.0_r8
      END DO
#endif
          ELSE
# if defined ANA_SEAICE
!           tl_stflx(i,j,isalt)=tl_Hz(i,j,N(ng))*                       &
!                            (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax-  &
!    &                  Hz(i,j,N(ng))*tl_t(i,j,N(ng),nrhs,isalt)/trelax
            ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+ad_stflx(i,j,isalt)*      &
     &                       (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax
            ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)-      &
     &                  Hz(i,j,N(ng))*ad_stflx(i,j,isalt)/trelax
            ad_stflx(i,j,isalt)=0.0_r8

!           tl_stflx(i,j,itemp)=tl_Hz(i,j,N(ng))*                       &
!                            (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax-  &
!    &                  Hz(i,j,N(ng))*tl_t(i,j,N(ng),nrhs,itemp)/trelax

            ad_Hz(i,j,N(ng))=ad_Hz(i,j,N(ng))+ad_stflx(i,j,itemp)*      &
     &                       (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax
            ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)-      &
     &                  Hz(i,j,N(ng))*ad_stflx(i,j,itemp)/trelax
            ad_stflx(i,j,itemp)=0.0_r8

# else
          DO itrc=1,NT(ng)
!            tl_stflx(i,j,itrc)=0.0_r8
             ad_stflx(i,j,itrc)=0.0_r8
          END DO
# endif
!           tl_m(i,j) = 0.0_r8
            ad_m(i,j) = 0.0_r8
          END IF

# ifdef SALINITY
!         tl_Sm=(0.5_r8-SIGN(0.5_r8,-t(i,j,N(ng),nrhs,isalt)))*         &
!    &                                tl_t(i,j,N(ng),nrhs,isalt)
          ad_t(i,j,N(ng),nrhs,isalt)=ad_t(i,j,N(ng),nrhs,isalt)+        &
     &            (0.5_r8-SIGN(0.5_r8,-t(i,j,N(ng),nrhs,isalt)))*ad_Sm
          ad_Sm=0.0_r8
# else
!         tl_Sm=0.0_r8
          ad_Sm=0.0_r8
# endif

        END DO
      END DO
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------

#ifdef DISTRIBUTE
!     CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!    &                    LBi, UBi, LBj, UBj,                           &
!    &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
!    &                    tl_sustr,tl_svstr)
      CALL ad_mp_exchange2d (ng, tile, iNLM, 2,                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    ad_sustr,ad_svstr)
#endif
!
!  Apply boundary conditions.
!
!     CALL bc_u2d_tile (ng, tile,                                       &
!    &                          LBi, UBi, LBj, UBj,                     &
!    &                          tl_sustr)
      CALL ad_bc_u2d_tile (ng, tile,                                    &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ad_sustr)
!     CALL bc_v2d_tile (ng, tile,                                       &
!    &                          LBi, UBi, LBj, UBj,                     &
!    &                          tl_svstr)
      CALL ad_bc_v2d_tile (ng, tile,                                    &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ad_svstr)

#   if defined UV_LOGDRAG
      DO j=JstrV-1,JEND
        DO i=IstrU-1,IEND
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/Zob(ng))
          cff2=vonKar*vonKar*cff1*cff1
          cff3=MAX(Cdb_min,cff2)
          wrk(i,j)=MIN(Cdb_max,cff3)
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

!          tl_svstr(i,j)=-0.5_r8*                                       &
!    &                    ((tl_wrk(i,j-1)+tl_wrk(i,j))*                 &
!    &                     v(i,j,N(ng),nrhs)*cff2+                      &
!    &                     (wrk(i,j-1)+wrk(i,j))*                       &
!    &                     (tl_v(i,j,N(ng),nrhs)*cff2+                  &
!    &                      v(i,j,N(ng),nrhs)*tl_cff2))

            ad_wrk(i,j-1)=ad_wrk(i,j-1)-0.5_r8*ad_svstr(i,j)*           &
     &                                  v(i,j,N(ng),nrhs)*cff2
            ad_wrk(i,j)=ad_wrk(i,j)-0.5_r8*ad_svstr(i,j)*               &
     &                                  v(i,j,N(ng),nrhs)*cff2
            ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)-0.5_r8*           &
     &                   (wrk(i,j-1)+wrk(i,j))*cff2*ad_svstr(i,j)
            ad_cff2=ad_cff2-0.5_r8*(wrk(i,j-1)+wrk(i,j))*               &
     &                         v(i,j,N(ng),nrhs)*ad_svstr(i,j)
            ad_svstr(i,j)=0.0_r8

            IF (cff2.ne.0.0_r8) THEN
!             tl_cff2=(cff1*tl_cff1+                                    &
!    &                  v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
              ad_cff1=ad_cff1+cff1*ad_cff2/cff2
              ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)+                &
     &                          ad_cff2*v(i,j,N(ng),nrhs)/cff2
              ad_cff2=0.0_r8
            ELSE
!             tl_cff2=0.0_r8
              ad_cff2=0.0_r8
            END IF

!           tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
!    &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
!    &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
!    &                       tl_u(i+1,j-1,N(ng),nrhs))

            adfac=0.25_r8*ad_cff1
            ad_u(i  ,j  ,N(ng),nrhs)=ad_u(i  ,j  ,N(ng),nrhs)+adfac
            ad_u(i+1,j  ,N(ng),nrhs)=ad_u(i+1,j  ,N(ng),nrhs)+adfac
            ad_u(i  ,j-1,N(ng),nrhs)=ad_u(i  ,j-1,N(ng),nrhs)+adfac
            ad_u(i+1,j-1,N(ng),nrhs)=ad_u(i+1,j-1,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
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

!           tl_sustr(i,j)=-0.5_r8*                                      &
!    &                    ((tl_wrk(i-1,j)+tl_wrk(i,j))*                 &
!    &                     u(i,j,N(ng),nrhs)*cff2+                      &
!    &                     (wrk(i-1,j)+wrk(i,j))*                       &
!    &                     (tl_u(i,j,N(ng),nrhs)*cff2+                  &
!    &                      u(i,j,N(ng),nrhs)*tl_cff2))

            ad_wrk(i-1,j)=ad_wrk(i-1,j)-0.5_r8*u(i,j,N(ng),nrhs)*cff2*  &
     &                                  ad_sustr(i,j)
            ad_wrk(i,j)=ad_wrk(i,j)-0.5_r8*u(i,j,N(ng),nrhs)*cff2*      &
     &                                  ad_sustr(i,j)
            ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)-0.5_r8*           &
     &                   (wrk(i-1,j)+wrk(i,j))*cff2*ad_sustr(i,j)
            ad_cff2=ad_cff2-0.5_r8*(wrk(i-1,j)+wrk(i,j))*               &
     &                   u(i,j,N(ng),nrhs)*ad_sustr(i,j)
            ad_sustr(i,j)=0.0_r8

            IF (cff2.ne.0.0_r8) THEN
!             tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
!    &                 cff1*tl_cff1)/cff2
              ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)+                &
     &                     u(i,j,N(ng),nrhs)*ad_cff2/cff2
              ad_cff1=ad_cff1+cff1*ad_cff2/cff2
              ad_cff2=0.0_r8
            ELSE
!             tl_cff2=0.0_r8
              ad_cff2=0.0_r8
            END IF
!           tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
!    &                    tl_v(i  ,j+1,N(ng),nrhs)+                     &
!    &                    tl_v(i-1,j  ,N(ng),nrhs)+                     &
!    &                    tl_v(i-1,j+1,N(ng),nrhs))
            
            adfac1=0.25_r8*ad_cff1
            ad_v(i  ,j  ,N(ng),nrhs)=ad_v(i  ,j  ,N(ng),nrhs)+adfac1
            ad_v(i  ,j+1,N(ng),nrhs)=ad_v(i  ,j+1,N(ng),nrhs)+adfac1
            ad_v(i-1,j  ,N(ng),nrhs)=ad_v(i-1,j  ,N(ng),nrhs)+adfac1
            ad_v(i-1,j+1,N(ng),nrhs)=ad_v(i-1,j+1,N(ng),nrhs)+adfac1
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO

!
!  Set logarithmic ice shelf cavity stress.
!
      DO j=JstrV-1,JEND
        DO i=IstrU-1,IEND
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/Zob(ng))
          cff2=vonKar*vonKar*cff1*cff1
          cff3=MAX(Cdb_min,cff2)
          wrk(i,j)=MIN(Cdb_max,cff3)

!         tl_wrk(i,j)=(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*tl_cff3
          ad_cff3=ad_cff3+(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*ad_wrk(i,j)
          ad_wrk(i,j)=0.0

!         tl_cff3=(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*tl_cff2
          ad_cff2=ad_cff2+(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*ad_cff3
          ad_cff3=0.0_r8

!         tl_cff2=vonKar*vonKar*2.0_r8*cff1*tl_cff1
          ad_cff1=ad_cff1+vonKar*vonKar*2.0_r8*cff1*ad_cff2
          ad_cff2=0.0_r8

!         tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng)))/     &
!    &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          ad_z_w(i,j,N(ng))=ad_z_w(i,j,N(ng))-cff1*cff1*ad_cff1/        &
     &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          ad_z_r(i,j,N(ng))=ad_z_r(i,j,N(ng))+cff1*cff1*ad_cff1/        &
     &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          ad_cff1=0.0_r8
        END DO
      END DO

#   elif defined UV_QDRAG
!
!  Set quadratic ice shelf cavity stress.
!
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))

!           tl_svstr(i,j)=-rdrg2(ng)*(tl_v(i,j,N(ng),nrhs)*cff2+        &
!    &                                v(i,j,N(ng),nrhs)*tl_cff2)
            ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)-rdrg2(ng)*cff2*   &
     &                           ad_svstr(i,j)
            ad_cff2=ad_cff2-rdrg2(ng)*v(i,j,N(ng),nrhs)*ad_svstr(i,j)
            ad_svstr(i,j)=0.0_r8

            IF (cff2.ne.0.0_r8) THEN
!             tl_cff2=(cff1*tl_cff1+                                    &
!    &              v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
              ad_cff1=ad_cff1+cff1*ad_cff2/cff2
              ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)+                &
     &                          v(i,j,N(ng),nrhs)*ad_cff2/cff2
              ad_cff2=0.0_r8
            ELSE
!             tl_cff2=0.0_r8
              ad_cff2=0.0_r8
            END IF

!           tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
!    &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
!    &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
!    &                       tl_u(i+1,j-1,N(ng),nrhs))

            adfac=0.25_r8*ad_cff1
            ad_u(i  ,j  ,N(ng),nrhs)=ad_u(i  ,j  ,N(ng),nrhs)+adfac
            ad_u(i+1,j  ,N(ng),nrhs)=ad_u(i+1,j  ,N(ng),nrhs)+adfac
            ad_u(i  ,j-1,N(ng),nrhs)=ad_u(i  ,j-1,N(ng),nrhs)+adfac
            ad_u(i+1,j-1,N(ng),nrhs)=ad_u(i+1,j-1,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
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

!           tl_sustr(i,j)=-rdrg2(ng)*(tl_u(i,j,N(ng),nrhs)*cff2+        &
!    &                                u(i,j,N(ng),nrhs)*tl_cff2)
            ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)-rdrg2(ng)*        &
     &                               cff2*ad_sustr(i,j)
            ad_cff2=ad_cff2-rdrg2(ng)*u(i,j,N(ng),nrhs)*ad_sustr(i,j)
            ad_sustr(i,j)=0.0_r8

            IF (cff2.ne.0.0_r8) THEN
!             tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
!    &                 cff1*tl_cff1)/cff2
              ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)+                &
                          u(i,j,N(ng),nrhs)*ad_cff2/cff2
              ad_cff1=ad_cff1+cff1*ad_cff2/cff2
              ad_cff2=0.0_r8
            ELSE
!             tl_cff2=0.0_r8
              ad_cff2=0.0_r8
            END IF

!           tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
!    &                       tl_v(i  ,j+1,N(ng),nrhs)+                  &
!    &                       tl_v(i-1,j  ,N(ng),nrhs)+                  &
!    &                       tl_v(i-1,j+1,N(ng),nrhs))
            adfac=0.25_r8*ad_cff1
            ad_v(i  ,j  ,N(ng),nrhs)=ad_v(i  ,j  ,N(ng),nrhs)+adfac
            ad_v(i  ,j+1,N(ng),nrhs)=ad_v(i  ,j+1,N(ng),nrhs)+adfac
            ad_v(i-1,j  ,N(ng),nrhs)=ad_v(i-1,j  ,N(ng),nrhs)+adfac
            ad_v(i-1,j+1,N(ng),nrhs)=ad_v(i-1,j+1,N(ng),nrhs)+adfac
            ad_cff1=0.0_r8
          END IF
        END DO
      END DO

#   elif defined UV_LDRAG
!
!  Set linear ice shelf cavity stress.
!
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!           tl_svstr(i,j)=-rdrg(ng)*tl_v(i,j,N(ng),nrhs)
            ad_v(i,j,N(ng),nrhs)=ad_v(i,j,N(ng),nrhs)-rdrg(ng)*         &
     &                          ad_svstr(i,j)
            ad_svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!           tl_sustr(i,j)=-rdrg(ng)*tl_u(i,j,N(ng),nrhs)
            ad_u(i,j,N(ng),nrhs)=ad_u(i,j,N(ng),nrhs)-rdrg(ng)*         &
     &                          ad_sustr(i,j)
            ad_sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO

#   else
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!           tl_sustr(i,j)=0.0_r8
            ad_sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!           tl_svstr(i,j)=0.0_r8
            ad_svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif

      RETURN
      END SUBROUTINE ad_iceshelf_vbc_tile

# ifdef ICESHELF_TEOS10

      REAL FUNCTION ad_theta_from_ct(iflag,ain,s,ct)
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

!
!   Clear adjoint variables
!
      ad_theta=0.0d0
      ad_dct=0.0d0
      ad_ct=0.0d0
      ad_dct_dth=0.0d0
      ad_dth_dct=0.0d0
      ad_th0=0.0d0
      ad_ct0=0.0d0
      ad_a5ct=0.0d0
      ad_b3ct=0.0d0
      ad_ct_factor=0.0d0
      ad_s=0.0d0
      ad_th0_num=0.0d0
      ad_ct_factor=0.0d0
      ad_rec_th0_den=0.0d0

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

!     tl_theta_from_ct = tl_theta-tl_dct/dct_dth+                       &
!    &                   tl_dct_dth*dct/(dct_dth*dct_dth)
      ad_theta=ad_theta+ain
      ad_dct=ad_dct-ain/dct_dth
      ad_dct_dth=ad_dct_dth+ain*dct/(dct_dth*dct_dth)

!     tl_dct_dth = tl_theta*de_dt0_F(s,theta)*rec_Cp0+                  &
!    &       (theta+273.15d0)*tl_de_dt0_F(s,tl_s,theta,tl_theta)*rec_Cp0

      ad_theta=ad_theta+ad_dct_dth*de_dt0_F(s,theta)*rec_Cp0
      ad_s=ad_s+(theta+273.15d0)*rec_Cp0*                               &
     &                             ad_de_dt0_F(1,ad_dct_dth,s,theta)
      ad_theta=ad_theta+(theta+273.15d0)*rec_Cp0*                       &
     &                             ad_de_dt0_F(2,ad_dct_dth,s,theta)
      ad_dct_dth=0.0d0

!     tl_dct = tl_ct_from_theta(s,tl_s,theta,tl_theta)-tl_ct
      ad_ct=ad_ct-ad_dct
      ad_s=ad_s+ad_ct_from_theta(1,ad_dct,s,theta)
      ad_theta=ad_theta+ad_ct_from_theta(2,ad_dct,s,theta)
      ad_dct=0.0d0

!     tl_theta = tl_th0-(tl_ct0-tl_ct)*dth_dct-(ct0-ct)*tl_dth_dct
      ad_th0=ad_th0+ad_theta
      ad_ct0=ad_ct0-ad_theta*dth_dct
      ad_ct=ad_ct+ad_theta*dth_dct
      ad_dth_dct=ad_dth_dct-(ct0-ct)*ad_theta
      ad_theta=0.0d0

!     tl_dth_dct = (tl_ct_factor+tl_a5ct-(b2+b3ct+b3ct)*tl_th0-         &
!    &             (tl_b3ct+tl_b3ct)*th0)*rec_th0_den+                  &
!    &             dth_dct*tl_rec_th0_den

      ad_ct_factor=ad_ct_factor+ad_dth_dct*rec_th0_den
      ad_a5ct=ad_a5ct+ad_dth_dct*rec_th0_den
      ad_th0=ad_th0-(b2+b3ct+b3ct)*ad_dth_dct*rec_th0_den
      ad_b3ct=ad_b3ct-2.0d0*th0*ad_dth_dct*rec_th0_den
      ad_rec_th0_den=ad_rec_th0_den+dth_dct*ad_dth_dct
      ad_dth_dct=0.0d0

!     tl_ct0 = tl_ct_from_theta(s,tl_s,th0,tl_th0)
      ad_s=ad_s+ad_ct_from_theta(1,ad_ct0,s,th0)
      ad_th0=ad_th0+ad_ct_from_theta(2,ad_ct0,s,th0)
      ad_ct0=0.0d0

!     tl_th0 = tl_th0_num*rec_th0_den+th0_num*tl_rec_th0_den
      ad_th0_num=ad_th0_num+rec_th0_den*ad_th0
      ad_rec_th0_den=ad_rec_th0_den+th0_num*ad_th0
      ad_th0=0.0d0

!     tl_rec_th0_den =-(b1*tl_s+tl_ct*(b2+b3ct)+ct*tl_b3ct)*            &
!    &                 rec_th0_den*rec_th0_den

      adfac=rec_th0_den*rec_th0_den*ad_rec_th0_den
      ad_s=ad_s-b1*adfac
      ad_ct=ad_ct+(b2+b3ct)*adfac
      ad_b3ct=ad_b3ct+ct*adfac
      ad_rec_th0_den=0.0d0

!     tl_th0_num = tl_s*(a1+a2*s)+s*a2*tl_s+tl_ct*ct_factor+            &
!    &             ct*tl_ct_factor
      ad_s=ad_s+(a1+a2*s)*ad_th0_num
      ad_s=ad_s+s*a2*ad_th0_num
      ad_ct=ad_ct+ct_factor*ad_th0_num
      ad_ct_factor=ad_ct_factor+ct*ad_th0_num
      ad_th0_num=0.0d0

!     tl_ct_factor = a4*tl_s+tl_a5ct
      ad_s=ad_s+a4*ad_ct_factor
      ad_a5ct=ad_a5ct+ad_ct_factor
      ad_ct_factor=0.0d0

!     tl_b3ct = b3*tl_ct
      ad_ct=ad_ct+b3*ad_b3ct
      ad_b3ct=0.0d0

!     tl_a5ct = a5*tl_ct
      ad_ct=ad_ct+a5*ad_a5ct
      ad_a5ct=0.0d0

      if(iflag.eq.1)ad_theta_from_ct=ad_s
      if(iflag.eq.2)ad_theta_from_ct=ad_ct

      RETURN
      END FUNCTION ad_theta_from_ct

      REAL FUNCTION ad_ct_from_theta(iflag,ain,s,theta)
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
!     tl_ct_from_theta = rec_Cp0*tl_penthalpy_F(s,tl_s,theta,tl_theta)

      if(iflag.eq.1)then
        ad_s=rec_Cp0*ad_penthalpy_F(iflag,ain,s,theta)
        ad_ct_from_theta=ad_s
      endif
      if(iflag.eq.2)then
        ad_theta=rec_Cp0*ad_penthalpy_F(iflag,ain,s,theta)
        ad_ct_from_theta=ad_theta
      endif

      RETURN
      END FUNCTION ad_ct_from_theta

      REAL FUNCTION ad_penthalpy_F(iflag,ain,s,th)
!   potential enthalpy from differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity (psu)
!   th               : potential temperature (deg C, ITS-90)
!   penthalpy_F      : potential enthalpy (J/kg)
!   check value      : penthalpy_F(20,20) = 81649.4876546448
!   DRJ on 10/12/03

      implicit real*8(a-h,o-z)

!
!  Clear adjoint variables
!
      ad_s=0.0d0
      ad_th=0.0d0
      ad_x=0.0d0
      ad_x2=0.0d0
      ad_y=0.0d0
      ad_fac1=0.0d0
      ad_fac2=0.0d0
      ad_fac2a=0.0d0
      ad_fac2b=0.0d0
      ad_fac3=0.0d0
      penthalpy_F=0.0d0

      a0=61.013624165232955d0
      a1=168776.46138048015d0
      a2=-2735.2785605119643d0
      a3=2574.2164453821442d0
      a4=-1536.6644434977545d0
      a5=545.734049793163d0
      a6=-50.910917284743334d0
      a7=-18.30489878927802d0
      a8=416.31512917743896d0
      a9=937.9793807560891d0
      a10=-3140.435779506947d0
      a11=2975.170149976973d0
      a12=-1760.137081144729d0
      a13=414.5655751783703d0
      a14=2167.72082596016d0
      a15=-1224.5772800562902d0
      a16=326.3074029273967d0
      a17=50.6703824689518d0
      a18=-12694.10018182362d0
      a19=4405.71847182968d0
      a20=-2132.9690185026416d0
      a21=303.91071982808035d0
      a22=69.74975368852d0

      x2 = 2.5d-2*s
      x = sqrt(x2)
      y = 2.5d-2*th

      fac1=y*(a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))))
      fac2a=x*(a10+x*(a11+x*(a12+x*a13)))
      fac2b=y*(a14+y*(a15+y*(a16+a17*y)))
      fac2=x*(a9+fac2a+fac2b)
      fac3=y*(a18+y*(a19+y*(a20+y*(a21+a22*y))))


!     tl_penthalpy_F =tl_fac1+tl_x2*(a8+fac2+fac3)+x2*(tl_fac2+tl_fac3)

!     ad_fac1=ad_fac1+ad_penthalpy_F
      ad_fac1=ad_fac1+ain
!     ad_x2=ad_x2+(a8+fac2+fac3)*ad_penthalpy_F
      ad_x2=ad_x2+(a8+fac2+fac3)*ain
!     ad_fac2=ad_fac2+x2*ad_penthalpy_F
      ad_fac2=ad_fac2+x2*ain
!     ad_fac3=ad_fac3+x2*ad_penthalpy_F
      ad_fac3=ad_fac3+x2*ain
!     ad_penthalpy_F=0.0d0
!     ain=0.0d0

!     tl_fac3=tl_y*(a18+y*(a19+y*(a20+y*(a21+a22*y))))+                 &
!    &           y*(tl_y*(a19+y*(a20+y*(a21+a22*y)))+                   &
!    &                 y*(tl_y*(a20+y*(a21+a22*y))+                     &
!    &                       y*(tl_y*(a21+a22*y)+y*a22*tl_y)            &
!    &                   )
!    &             )

      ad_y=ad_y+(a18+y*(a19+y*(a20+y*(a21+a22*y))))*ad_fac3
      ad_y=ad_y+y*(a19+y*(a20+y*(a21+a22*y)))*ad_fac3
      ad_y=ad_y+y*y*(a20+y*(a21+a22*y))*ad_fac3
      ad_y=ad_y+y*y*y*(a21+a22*y)*ad_fac3
      ad_y=ad_y+y*y*y*y*a22*ad_fac3
      ad_fac3=0.0d0

!     tl_fac2=tl_x*(a9+fac2a+fac2b)+x*(tl_fac2a+tl_fac2b)
      ad_x=ad_x+(a9+fac2a+fac2b)*ad_fac2
      ad_fac2a=ad_fac2a+x*ad_fac2
      ad_fac2b=ad_fac2b+x*ad_fac2
      ad_fac2=0.0d0

!     tl_fac2b=tl_y*(a14+y*(a15+y*(a16+a17*y)))+                        &
!    &            y*(tl_y*(a15+y*(a16+a17*y))+                          &
!    &                  y*(tl_y*(a16+a17*y)+y*a17*tl_y)     )

      ad_y=ad_y+(a14+y*(a15+y*(a16+a17*y)))*ad_fac2b
      ad_y=ad_y+y*(a15+y*(a16+a17*y))*ad_fac2b
      ad_y=ad_y+y*y*(a16+a17*y)*ad_fac2b
      ad_y=ad_y+y*y*y*a17*ad_fac2b
      ad_fac2b=0.0d0

!     tl_fac2a=tl_x*(a10+x*(a11+x*(a12+x*a13)))+                        &
!    &            x*(tl_x*(a11+x*(a12+x*a13))+                          &
!    &                  x*(tl_x*(a12+x*a13)+tl_x*x*a13)                 &
!    &              )

      ad_x=ad_x+(a10+x*(a11+x*(a12+x*a13)))*ad_fac2a
      ad_x=ad_x+x*(a11+x*(a12+x*a13))*ad_fac2a
      ad_x=ad_x+x*x*(a12+x*a13)*ad_fac2a
      ad_x=ad_x+x*x*x*a13*ad_fac2a
      ad_fac2a=0.0d0

!     tl_fac1=tl_y*(a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))))+        &
!    &        y*(tl_y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))+            &
!    &              y*(tl_y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))+             &
!    &                    y*(tl_y*(a4+y*(a5+(a6+a7*y)*y))+              &
!    &                          y*(tl_y*(a5+(a6+a7*y)*y)+               &
!    &                                y*((a6+a7*y)*tl_y+a7*y*tl_y)      &
!    &                            )
!    &                      )
!    &                )
!    &          )

      ad_y=ad_y+(a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))))*ad_fac1
      ad_y=ad_y+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))*ad_fac1
      ad_y=ad_y+y*y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))*ad_fac1
      ad_y=ad_y+y*y*y*(a4+y*(a5+(a6+a7*y)*y))*ad_fac1
      ad_y=ad_y+y*y*y*y*(a5+(a6+a7*y)*y)*ad_fac1
      ad_y=ad_y+y*y*y*y*y*(a6+a7*y)*ad_fac1
      ad_y=ad_y+y*y*y*y*y*y*a7*ad_fac1
      ad_fac1=0.0d0

!     tl_y = 2.5d-2*tl_th
      ad_th=ad_th+2.5d-2*ad_y

!     tl_x=0.5d0*tl_x2/x
      ad_x2=ad_x2+0.5d0*ad_x/x

!     tl_x2 = 2.5d-2*tl_s
      ad_s=ad_s+2.5d-2*ad_x2
      ad_x2=0.0d0

      if(iflag.eq.1) ad_penthalpy_F=ad_s
      if(iflag.eq.2) ad_penthalpy_F=ad_th


      RETURN
      END FUNCTION ad_penthalpy_F

      REAL FUNCTION ad_de_dt0_F(iflag,ain,s,th)
!   d(entropy)/dt from twice differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C,
!   ITS-90)
!   de_dt0_F         : d(entropy)/dt                      (J/kgK^2)
!   check value      : de_dt0_F(35,20) = 13.63256369213874
!   DRJ on 30/06/05

      implicit real*8(a-h,o-z)

!
!  Clear adjoint variables
!
      ad_s=0.0d0
      ad_y=0.0d0
      ad_de_dt=0.0d0
      ad_x2=0.0d0
      ad_x=0.0d0
      ad_th=0.0d0
      ad_fac1=0.0d0
      ad_fac2=0.0d0
      ad_fac3=0.0d0

      a0=24715.571866078d0
      a1=-1858.920033948178d0
      a2=317.440355256842d0
      a3=-405.1392883572282d0
      a4=202.6815298758072d0
      a5=1562.563716288858d0
      a6=-1165.8752731900836d0
      a7=348.7487684426d0
      a8=-4420.4472249096725d0
      a9=1778.231237203896d0
      a10=-1160.5182516851419d0
      a11=569.531539542516d0
      a12=-128.13429152494615d0

      x2 = 2.5d-2*s
      x = sqrt(x2)
      y = 2.5d-2*th

      fac1=x*(a2+y*(a3+a4*y))
      fac2=y*(a5+y*(a6+a7*y))
      fac3=y*(a8+y*(a9+y*(a10+(a11+a12*y)*y)))


!     tl_de_dt0_F = 6.25d-4*tl_de_dt
!     ad_de_dt=ad_de_dt+6.25d-4*ad_de_dt0_F
      ad_de_dt=ad_de_dt+6.25d-4*ain

!     tl_de_dt = tl_x2*(a1+fac1+fac2)+x2*(tl_fac1+tl_fac2)+tl_fac3

      ad_x2=ad_x2+(a1+fac1+fac2)*ad_de_dt
      ad_fac1=ad_fac1+x2*ad_de_dt
      ad_fac2=ad_fac2+x2*ad_de_dt
      ad_fac3=ad_fac3+ad_de_dt
      ad_de_dt=0.0d0

!     tl_fac3=tl_y*(a8+y*(a9+y*(a10+(a11+a12*y)*y)))+                   &
!    &           y*(tl_y*(a9+y*(a10+(a11+a12*y)*y))+                    &
!    &                 y*(tl_y*(a10+(a11+a12*y)*y)+                     &
!    &                       y*((a11+a12*y)*tl_y+a12*y*tl_y)            &
!    &                   )
!    &             )

       ad_y=ad_y+(a8+y*(a9+y*(a10+(a11+a12*y)*y)))*ad_fac3
       ad_y=ad_y+y*(a9+y*(a10+(a11+a12*y)*y))*ad_fac3
       ad_y=ad_y+y*y*(a10+(a11+a12*y)*y)*ad_fac3
       ad_y=ad_y+y*y*y*(a11+a12*y)*ad_fac3
       ad_y=ad_y+y*y*y*y*a12*ad_fac3
       ad_fac3=0.0d0

!     tl_fac2=tl_y*(a5+y*(a6+a7*y))+y*(tl_y*(a6+a7*y)+y*a7*tl_y)

      ad_y=ad_y+(a5+y*(a6+a7*y))*ad_fac2
      ad_y=ad_y+y*(a6+a7*y)*ad_fac2
      ad_y=ad_y+y*y*a7*ad_fac2
      ad_fac2=0.0d0

!     tl_fac1=tl_x*(a2+y*(a3+a4*y))+x*(tl_y*(a3+a4*y)+y*a4*tl_y)

      ad_x=ad_x+(a2+y*(a3+a4*y))*ad_fac1
      ad_y=ad_y+x*(a3+a4*y)*ad_fac1
      ad_y=ad_y+x*y*a4*ad_fac1
      ad_fac1=0.0d0

!     tl_y = 2.5d-2*tl_th
      ad_th=ad_th+2.5d-2*ad_y
      ad_y=0.0d0

!     tl_x = 0.5d0*tl_x2/x
      ad_x2=ad_x2+0.5d0*ad_x/x
      ad_x=0.0d0

!     tl_x2 = 2.5d-2*tl_s
      ad_s=ad_s+2.5d-2*ad_x2
      ad_x2=0.0d0

      if(iflag.eq.1)ad_de_dt0_F=ad_s
      if(iflag.eq.2)ad_de_dt0_F=ad_th

      RETURN
      END FUNCTION ad_de_dt0_F

      REAL FUNCTION ad_entropy_diff_F(iflag,ain,s,t,p,th0,pr)
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

!
!    Clear adjoint variables.
!

      ad_fTp=0.0d0
      ad_x2=0.0d0
      ad_x=0.0d0
      ad_gTp=0.0d0
      ad_hTp=0.0d0
      ad_kTp=0.0d0
      ad_fac1=0.0d0
      ad_fac2=0.0d0
      ad_fac3=0.0d0
      ad_fac4=0.0d0
      ad_fac5=0.0d0
      ad_fac6=0.0d0
      ad_y=0.0d0
      ad_z=0.0d0
      ad_z1=0.0d0
      ad_yz=0.0d0
      ad_p=0.0d0
      ad_t=0.0d0
      ad_s=0.0d0
      ad_th0=0.0d0
      ad_pr=0.0d0
!
!    Compute basic state I.
!

      a1=-5.90578348518236d0
      a2=24715.571866078d0
      a3=-2210.2236124548363d0
      a4=592.743745734632d0
      a5=-290.12956292128547d0
      a6=113.90630790850321d0
      a7=-21.35571525415769d0

      x2 = 2.5d-2*s 
      x = dsqrt(x2) 
      y = 2.5d-2*th0 
      z = 1.0d-4*pr 
      z1 = z

!     fTp = -5.90578348518236d0 + y*(24715.571866078d0 +                &
!    &         y*(-2210.2236124548363d0 + y*(592.743745734632d0 +       &
!    &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
!    &                               21.35571525415769d0*y)*y))))

      fTp = a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))

      b1=-1858.920033948178d0
      b2=781.281858144429d0
      b3=-388.6250910633612d0
      b4=87.18719211065d0

!     gTp = y*(-1858.920033948178d0 + y*(781.281858144429d0 +           &
!    &                y*(-388.6250910633612d0 + 87.18719211065d0*y)))

      gTp = y*(b1+y*(b2+y*(b3+b4*y)))

      c1=317.440355256842d0
      c2=-202.5696441786141d0
      c3=67.5605099586024d0

!     hTp = y*(317.440355256842d0 + y*(-202.5696441786141d0 +           &
!                                        67.5605099586024d0*y))

      hTp = y*(c1+y*(c2+c3*y))

      if(z.ne.0.d0) then
        yz = y*z

        d1=270.983805184062d0
        d2=-776.153611613101d0
        d3=196.51255088122d0
        d4=-28.9796526294175d0
        d5=2.13290083518327d0
        d6=-2910.0729080936d0
        d7=1513.116771538718d0
        d8=-546.959324647056d0
        d9=111.1208127634436d0
        d10=-8.68841343834394d0
        d11=2017.52334943521d0
        d12=-1498.081172457456d0
        d13=718.6359919632359d0
        d14=-146.4037555781616d0
        d15=4.9892131862671505d0
        d16=-1591.873781627888d0
        d17=1207.261522487504d0
        d18=-608.785486935364d0
        d19=105.4993508931208d0
        d20=67.41756835751434d0
        d21=-381.06836198507096d0 
        d22=133.7383902842754d0
        d23=-49.023632509086724d0
        d24=973.091553087975d0
        d25=-602.603274510125d0
        d26=276.361526170076d0
        d27=-32.40953340386105d0

!       fTp = fTp + z*(270.983805184062d0 + z*(-776.153611613101d0 +    &
!    &               z*(196.51255088122d0 + z*(-28.9796526294175d0 +    &
!    &       2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 +      &
!    &           z*(1513.116771538718d0 + z*(-546.959324647056d0 +      &
!    &         (111.1208127634436d0 - 8.68841343834394d0*z)*z))) +      &
!    &        y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 +      &
!    &          z*(718.6359919632359d0 + z*(-146.4037555781616d0 +      &
!    &   4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 +      &
!    &           z*(1207.261522487504d0 + z*(-608.785486935364d0 +      &
!    &  105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z +      &
!    &           z*(-381.06836198507096d0 + (133.7383902842754d0 -      &
!    &       49.023632509086724d0*z)*z)) + z*(973.091553087975d0 +      &
!    &              z*(-602.603274510125d0 + (276.361526170076d0 -      &
!    &                           32.40953340386105d0*z)*z))))))

        fTp = fTp+z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))+y*(z*(d6+            &
     &           z*(d7+z*(d8+(d9+d10*z)*z)))+y*(z*(d11+z*(d12+          &
     &          z*(d13+z*(d14+d15*z))))+y*(z*(d16+z*(d17+z*(d18+        &
     &          d19*z)))+y*(y*(d20*y*z+z*(d21+(d22+d23*z)*z))+z*(d24+   &
     &              z*(d25+(d26+d27*z)*z))))))

        e1=-729.116529735046d0
        e2=343.956902961561d0
        e3=-124.687671116248d0
        e4=31.656964386073d0
        e5=7.04658803315449d0
        e6=1721.528607567954d0
        e7=-674.819060538734d0
        e8=356.629112415276d0
        e9=-88.4080716616d0
        e10=15.84003094423364d0
        e11=1190.914967948748d0
        e12=-298.904564555024d0
        e13=145.9491676006352d0
        e14=-2082.7344423998043d0
        e15=614.668925894709d0
        e16=-340.685093521782d0
        e17=33.3848202979239d0


!       gTp = gTp + z*(-729.116529735046d0 + z*(343.956902961561d0 +    &
!    &               z*(-124.687671116248d0 + (31.656964386073d0 -      &
!    &    7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 +      &
!    &            z*(-674.819060538734d0 + z*(356.629112415276d0 +      &
!    &           z*(-88.4080716616d0 + 15.84003094423364d0*z)))) +      &
!    &      y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 +      &
!    &        145.9491676006352d0*z)) + z*(-2082.7344423998043d0 +      &
!                z*(614.668925894709d0 + z*(-340.685093521782d0 +       &
!                                     33.3848202979239d0*z)))))

        gTp = gTp+z*(e1+z*(e2+z*(e3+(e4-e5*z)*z)))+y*(z*(e6+            &
     &            z*(e7+z*(e8+z*(e9+e10*z))))+y*(y*z*(e11+z*(e12+       &
     &            e13*z))+z*(e14+z*(e15+z*(e16+e17*z)))))

        f1=175.292041186547d0
        f2=-83.1923927801819d0
        f3=29.483064349429d0
        f4=1380.9597954037708d0
        f5=-938.26075044542d0
        f6=-766.116132004952d0
        f7=108.3834525034224d0
        f8=-51.2796974779828d0

!       hTp = hTp + z*(175.292041186547d0 + z*(-83.1923927801819d0 +    &
!    &     29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z -       &
!    &          938.26075044542d0*y*z) + z*(-766.116132004952d0 +       &
!    &          (108.3834525034224d0 - 51.2796974779828d0*z)*z))

        hTp =hTp+z*(f1+z*(f2+f3*z))+                                    &
     &                       y*(y*(f4*z+f5*y*z)+z*(f6+(f7+f8*z)*z))

      end if
       y = 2.5d-2*t 
       z = 1.0d-4*p

        g1=5.90578348518236d0
        g2=24715.571866078d0
        g3=-2210.2236124548363d0
        g4=592.743745734632d0
        g5=-290.12956292128547d0
        g6=113.90630790850321d0
        g7=-21.35571525415769d0

!       fTp = fTp +5.90578348518236d0 - y*(24715.571866078d0 +          &
!    &          y*(-2210.2236124548363d0 + y*(592.743745734632d0 +      &
!    &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
!    &                              21.35571525415769d0*y)*y))))

        fTp =fTp+g1-y*(g2+y*(g3+y*(g4+y*(g5+(g6+g7*y)*y))))

        h1=-1858.920033948178d0
        h2=781.281858144429d0
        h3=-388.6250910633612d0
        h4=87.18719211065d0

!       gTp = gTp - y*(-1858.920033948178d0 + y*(781.281858144429d0 +   &
!    &       y*(-388.6250910633612d0 + 87.18719211065d0*y)))

        gTp=gTp- y*(h1+y*(h2+y*(h3+h4*y)))

        p1=317.440355256842d0
        p2=-202.5696441786141d0
        p3=67.5605099586024d0

!       hTp = hTp - y*(317.440355256842d0 + y*(-202.5696441786141d0 +   &
!    &                                     67.5605099586024d0*y))

        hTp=hTp-y*(p1+y*(p2+p3*y))

        kTp = 22.6683558512829d0*(z1-z)

      if(z.ne.0.d0) then
        yz = y*z

         q1=270.983805184062d0
         q2=-776.153611613101d0
         q3=196.51255088122d0
         q4=-28.9796526294175d0
         q5=2.13290083518327d0
         q6=-2910.0729080936d0
         q7=1513.116771538718d0
         q8=-546.959324647056d0
         q9=111.1208127634436d0
         q10=-8.68841343834394d0
         q11=2017.52334943521d0
         q12=-1498.081172457456d0
         q13=718.6359919632359d0
         q14=-146.4037555781616d0
         q15=4.9892131862671505d0
         q16=-1591.873781627888d0
         q17=1207.261522487504d0
         q18=-608.785486935364d0
         q19=105.4993508931208d0
         q20=67.41756835751434d0
         q21=-381.06836198507096d0
         q22=133.7383902842754d0
         q23=-49.023632509086724d0
         q24=973.091553087975d0
         q25=-602.603274510125d0
         q26=276.361526170076d0
         q27=-32.40953340386105d0

!        fTp = fTp - ( z*(270.983805184062d0 + z*(-776.153611613101d0 + &
!    &              z*(196.51255088122d0 + z*(-28.9796526294175d0 +     &
!    &         2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 +    &
!    &             z*(1513.116771538718d0 + z*(-546.959324647056d0 +    &
!    &           (111.1208127634436d0 - 8.68841343834394d0*z)*z))) +    &
!    &          y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 +    &
!    &            z*(718.6359919632359d0 + z*(-146.4037555781616d0 +    &
!    &     4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 +    &
!    &             z*(1207.261522487504d0 + z*(-608.785486935364d0 +    &
!    &    105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z +    &
!    &             z*(-381.06836198507096d0 + (133.7383902842754d0 -    &
!    &         49.023632509086724d0*z)*z)) + z*(973.091553087975d0 +    &
!    &                z*(-602.603274510125d0 + (276.361526170076d0 -    &
!    &                           32.40953340386105d0*z)*z)))))))

         fTp=fTp-(z*(q1+z*(q2+z*(q3+z*(q4+q5*z))))+y*(z*(q6+z*(q7+      &
     &            z*(q8+(q9+q10*z)*z)))+y*(z*(q11+z*(q12+z*(q13+z*(q14+ &
     &            q15*z))))+y*(z*(q16+z*(q17+z*(q18+q19*z)))+           &
     &            y*(y*(q20*y*z+z*(q21+(q22+q23*z)*z))+z*(q24+z*(q25+   &
     &            (q26+q27*z)*z)))))))

         r1=-729.116529735046d0
         r2=343.956902961561d0
         r3=-124.687671116248d0
         r4=31.656964386073d0
         r5=-7.04658803315449d0
         r6=1721.528607567954d0
         r7=-674.819060538734d0
         r8=356.629112415276d0
         r9=-88.4080716616d0
         r10=15.84003094423364d0
         r11=1190.914967948748d0
         r12=-298.904564555024d0
         r13=145.9491676006352d0
         r14=-2082.7344423998043d0
         r15=614.668925894709d0
         r16=-340.685093521782d0
         r17=33.3848202979239d0

!        gTp = gTp - ( z*(-729.116529735046d0 + z*(343.956902961561d0 + &
!    &                z*(-124.687671116248d0 + (31.656964386073d0 -     &
!    &     7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 +     &
!    &             z*(-674.819060538734d0 + z*(356.629112415276d0 +     &
!    &            z*(-88.4080716616d0 + 15.84003094423364d0*z)))) +     &
!    &       y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 +     &
!    &         145.9491676006352d0*z)) + z*(-2082.7344423998043d0 +     &
!    &             z*(614.668925894709d0 + z*(-340.685093521782d0 +     &
!    &                                 33.3848202979239d0*z))))))

         gTp=gTp-(z*(r1+z*(r2+z*(r3+(r4+r5*z)*z)))+y*(z*(r6+z*(r7+      &
     &            z*(r8+z*(r9+r10*z))))+y*(y*z*(r11+z*(r12+r13*z))+     &
     &            z*(r14+z*(r15+z*(r16+r17*z))))))

         s1=175.292041186547d0
         s2=-83.1923927801819d0
         s3=29.483064349429d0
         s4=1380.9597954037708d0
         s5=-938.26075044542d0
         s6=-766.116132004952d0
         s7=108.3834525034224d0
         s8=-51.2796974779828d0


!        hTp = hTp - ( z*(175.292041186547d0 + z*(-83.1923927801819d0 + &
!    &       29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z -     &
!    &            938.26075044542d0*y*z) + z*(-766.116132004952d0 +     &
!    &          (108.3834525034224d0 - 51.2796974779828d0*z)*z)))

         hTp =hTp-(z*(s1+z*(s2+s3*z))+y*(y*(s4*z+s5*y*z)+               &
     &                z*(s6+(s7+s8*z)*z)))

      end if
!
!    End of compute basic state I.
!

!     tl_entropy_diff_F = 2.5d-2*(tl_fTp+tl_x2*(gTp+x*hTp+x2*kTp)+       &
!    &             x2*(tl_gTp+tl_x*hTp+x*tl_hTp+tl_x2*tl_kTp+x2*tl_kTp))

      ad_fTp=ad_fTp+2.5d-2*ain
      ad_x2=ad_x2+2.5d-2*(gTp+x*hTp+x2*kTp)*ain
      ad_gTp=ad_gTp+2.5d-2*x2*ain
      ad_x=ad_x+2.5d-2*x2*hTp*ain
      ad_hTp=ad_hTp+2.5d-2*x2*x*ain
      ad_x2=ad_x2+2.5d-2*x2*kTp*ain
      ad_kTp=ad_kTp+2.5d-2*x2*x2*ain

!
!   Compute basic state II.
!

      x2 = 2.5d-2*s 
      x = dsqrt(x2) 
      y = 2.5d-2*th0 
      z = 1.0d-4*pr 
      z1 = z

!     fTp = -5.90578348518236d0 + y*(24715.571866078d0 +                &
!    &         y*(-2210.2236124548363d0 + y*(592.743745734632d0 +       &
!    &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
!    &                               21.35571525415769d0*y)*y))))

      fTp = a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))


!     gTp = y*(-1858.920033948178d0 + y*(781.281858144429d0 +           &
!    &                y*(-388.6250910633612d0 + 87.18719211065d0*y)))

      gTp = y*(b1+y*(b2+y*(b3+b4*y)))


!     hTp = y*(317.440355256842d0 + y*(-202.5696441786141d0 +           &
!                                        67.5605099586024d0*y))

      hTp = y*(c1+y*(c2+c3*y))

      if(z.ne.0.d0) then
        yz = y*z

!       fTp = fTp + z*(270.983805184062d0 + z*(-776.153611613101d0 +    &
!    &               z*(196.51255088122d0 + z*(-28.9796526294175d0 +    &
!    &       2.13290083518327d0*z)))) + y*(z*(-2910.0729080936d0 +      &
!    &           z*(1513.116771538718d0 + z*(-546.959324647056d0 +      &
!    &         (111.1208127634436d0 - 8.68841343834394d0*z)*z))) +      &
!    &        y*(z*(2017.52334943521d0 + z*(-1498.081172457456d0 +      &
!    &          z*(718.6359919632359d0 + z*(-146.4037555781616d0 +      &
!    &   4.9892131862671505d0*z)))) + y*(z*(-1591.873781627888d0 +      &
!    &           z*(1207.261522487504d0 + z*(-608.785486935364d0 +      &
!    &  105.4993508931208d0*z))) + y*(y*(67.41756835751434d0*y*z +      &
!    &           z*(-381.06836198507096d0 + (133.7383902842754d0 -      &
!    &       49.023632509086724d0*z)*z)) + z*(973.091553087975d0 +      &
!    &              z*(-602.603274510125d0 + (276.361526170076d0 -      &
!    &                           32.40953340386105d0*z)*z))))))

        fTp = fTp+z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))+y*(z*(d6+            &
     &           z*(d7+z*(d8+(d9+d10*z)*z)))+y*(z*(d11+z*(d12+          &
     &          z*(d13+z*(d14+d15*z))))+y*(z*(d16+z*(d17+z*(d18+        &
     &          d19*z)))+y*(y*(d20*y*z+z*(d21+(d22+d23*z)*z))+z*(d24+   &
     &              z*(d25+(d26+d27*z)*z))))))


!       gTp = gTp + z*(-729.116529735046d0 + z*(343.956902961561d0 +    &
!    &               z*(-124.687671116248d0 + (31.656964386073d0 -      &
!    &    7.04658803315449d0*z)*z))) + y*(z*(1721.528607567954d0 +      &
!    &            z*(-674.819060538734d0 + z*(356.629112415276d0 +      &
!    &           z*(-88.4080716616d0 + 15.84003094423364d0*z)))) +      &
!    &      y*(y*z*(1190.914967948748d0 + z*(-298.904564555024d0 +      &
!    &        145.9491676006352d0*z)) + z*(-2082.7344423998043d0 +      &
!                z*(614.668925894709d0 + z*(-340.685093521782d0 +       &
!                                     33.3848202979239d0*z)))))

        gTp = gTp+z*(e1+z*(e2+z*(e3+(e4-e5*z)*z)))+y*(z*(e6+            &
     &            z*(e7+z*(e8+z*(e9+e10*z))))+y*(y*z*(e11+z*(e12+       &
     &            e13*z))+z*(e14+z*(e15+z*(e16+e17*z)))))


!       hTp = hTp + z*(175.292041186547d0 + z*(-83.1923927801819d0 +    &
!    &     29.483064349429d0*z)) + y*(y*(1380.9597954037708d0*z -       &
!    &          938.26075044542d0*y*z) + z*(-766.116132004952d0 +       &
!    &          (108.3834525034224d0 - 51.2796974779828d0*z)*z))

        hTp =hTp+z*(f1+z*(f2+f3*z))+                                    &
     &                       y*(y*(f4*z+f5*y*z)+z*(f6+(f7+f8*z)*z))

      end if
       y = 2.5d-2*t 
       z = 1.0d-4*p


!       fTp = fTp +5.90578348518236d0 - y*(24715.571866078d0 +          &
!    &          y*(-2210.2236124548363d0 + y*(592.743745734632d0 +      &
!    &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
!    &                              21.35571525415769d0*y)*y))))

        fTp =fTp+g1-y*(g2+y*(g3+y*(g4+y*(g5+(g6+g7*y)*y))))

!       gTp = gTp - y*(-1858.920033948178d0 + y*(781.281858144429d0 +   &
!    &       y*(-388.6250910633612d0 + 87.18719211065d0*y)))

        gTp=gTp- y*(h1+y*(h2+y*(h3+h4*y)))


!       hTp = hTp - y*(317.440355256842d0 + y*(-202.5696441786141d0 +   &
!    &                                     67.5605099586024d0*y))

        hTp=hTp-y*(p1+y*(p2+p3*y))

        kTp = 22.6683558512829d0*(z1-z)
!
!  End compute basic state II.

!
      y = 2.5d-2*t
      z = 1.0d-4*p

      if(z.ne.0.d0) then
        yz = y*z

         fac1=z*(s1+z*(s2+s3*z))
         fac2=y*(s4*z+s5*y*z)
         fac3=z*(s6+(s7+s8*z)*z)

!        tl_hTp =tl_hTp-(tl_fac1+tl_y*(fac2+fac3)+y*(tl_fac2+tl_fac3))

         ad_fac1=ad_fac1-ad_hTp
         ad_y=ad_y-(fac2+fac3)*ad_hTp
         ad_fac2=ad_fac2-y*ad_hTp
         ad_fac3=ad_fac3-y*ad_hTp

!        tl_fac3=tl_z*(s6+(s7+s8*z)*z)+z*((s7+s8*z)*tl_z+s8*z*tl_z)

         ad_z=ad_z+(s6+(s7+s8*z)*z)*ad_fac3
         ad_z=ad_z+z*(s7+s8*z)*ad_fac3
         ad_z=ad_z+z*s8*z*ad_fac3
         ad_fac3=0.0d0

!        tl_fac2=tl_y*(s4*z+s5*y*z)+y*(s4*tl_z+s5*(tl_y*z+y*tl_z))

         ad_y=ad_y+(s4*z+s5*y*z)*ad_fac2
         ad_z=ad_z+y*s4*ad_fac2
         ad_y=ad_y+y*s5*z*ad_fac2
         ad_z=ad_z+y*s5*y*ad_fac2
         ad_fac2=0.0d0

!        tl_fac1=tl_z*(s1+z*(s2+s3*z))+
!    &              z*(tl_z*(s2+s3*z)+tl_z*s3*z)

         ad_z=ad_z+(s1+z*(s2+s3*z))*ad_fac1
         ad_z=ad_z+z*(s2+s3*z)*ad_fac1
         ad_z=ad_z+z*s3*z*ad_fac1
         ad_fac1=0.0d0

!

         fac1=z*(r1+z*(r2+z*(r3+(r4+r5*z)*z)))
         fac2=z*(r6+z*(r7+z*(r8+z*(r9+r10*z))))
         fac3=y*z*(r11+z*(r12+r13*z))
         fac4=z*(r14+z*(r15+z*(r16+r17*z)))

!        tl_gTp=tl_gTp-(tl_fac1+tl_y*(fac2+y*(fac3+fac4))+              &
!    &                            y*(tl_fac2+tl_y*(fac3+fac4)+          &
!    &                                          y*(tl_fac3+tl_fac4) ) )

         ad_fac1=ad_fac1-ad_gTp
         ad_y=ad_y-(fac2+y*(fac3+fac4))*ad_gTp
         ad_fac2=ad_fac2-y*ad_gTp
         ad_y=ad_y-y*(fac3+fac4)*ad_gTp
         ad_fac3=ad_fac3-y*y*ad_gTp
         ad_fac4=ad_fac4-y*y*ad_gTp

!        tl_fac4=tl_z*(r14+z*(r15+z*(r16+r17*z)))+                      &
!    &              z*(tl_z*(r15+z*(r16+r17*z))+                        &
!    &                    z*(tl_z*(r16+r17*z)+tl_z*r17*z)   )

         ad_z=ad_z+(r14+z*(r15+z*(r16+r17*z)))*ad_fac4
         ad_z=ad_z+z*(r15+z*(r16+r17*z))*ad_fac4
         ad_z=ad_z+z*z*(r16+r17*z)*ad_fac4
         ad_z=ad_z+z*z*r17*z*ad_fac4
         ad_fac4=0.0d0

!        tl_fac3=(tl_y*z+y*tl_z)*(r11+z*(r12+r13*z))+                   &
!    &            y*z*(tl_z*(r12+r13*z)+tl_z*r13*z)

         ad_y=ad_y+z*(r11+z*(r12+r13*z))*ad_fac3
         ad_z=ad_z+y*(r11+z*(r12+r13*z))*ad_fac3
         ad_z=ad_z+y*z*(r12+r13*z)*ad_fac3
         ad_z=ad_z+y*z*r13*z*ad_fac3
         ad_fac3=0.0d0

!        tl_fac2=tl_z*(r6+z*(r7+z*(r8+z*(r9+r10*z))))+                  &
!    &              z*(tl_z*(r7+z*(r8+z*(r9+r10*z)))+                   &
!    &                    z*(tl_z*(r8+z*(r9+r10*z))+                    &
!    &                          z*(tl_z*(r9+r10*z)+z*r10*tl_z)          &
!    &                      )                                           &
!    &                )

         ad_z=ad_z+(r6+z*(r7+z*(r8+z*(r9+r10*z))))*ad_fac2
         ad_z=ad_z+z*(r7+z*(r8+z*(r9+r10*z)))*ad_fac2
         ad_z=ad_z+z*z*(r8+z*(r9+r10*z))*ad_fac2
         ad_z=ad_z+z*z*z*(r9+r10*z)*ad_fac2
         ad_z=ad_z+z*z*z*z*r10*ad_fac2
         ad_fac2=0.0d0

!        tl_fac1=tl_z*(r1+z*(r2+z*(r3+(r4+r5*z)*z)))+                   &
!    &              z*(tl_z*(r2+z*(r3+(r4+r5*z)*z))+                    &
!    &                    z*(tl_z*(r3+(r4+r5*z)*z)+                     &
!    &                          z*((r4+r5*z)*tl_z+r5*z*tl_z)            &
!    &                      )                                           &
!    &                )

         ad_z=ad_z+(r1+z*(r2+z*(r3+(r4+r5*z)*z)))*ad_fac1
         ad_z=ad_z+z*(r2+z*(r3+(r4+r5*z)*z))*ad_fac1
         ad_z=ad_z+z*z*(r3+(r4+r5*z)*z)*ad_fac1
         ad_z=ad_z+z*z*z*(r4+r5*z)*ad_fac1
         ad_z=ad_z+z*z*z*r5*z*ad_fac1
         ad_fac1=0.0d0
!
         fac1=z*(q1+z*(q2+z*(q3+z*(q4+q5*z))))
         fac2=z*(q6+z*(q7+z*(q8+(q9+q10*z)*z)))
         fac3=z*(q11+z*(q12+z*(q13+z*(q14+q15*z))))
         fac4=z*(q16+z*(q17+z*(q18+q19*z)))
         fac5=y*(q20*y*z+z*(q21+(q22+q23*z)*z))
         fac6=z*(q24+z*(q25+(q26+q27*z)*z))

!        tl_fTp=tl_fTp-                                                 &
!    &             (tl_fac1+tl_y*(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))+&
!    &                  y*(tl_fac2+tl_y*(fac3+y*(fac4+y*(fac5+fac6)))+  &
!    &                  y*(tl_fac3+tl_y*(fac4+y*(fac5+fac6))+           &
!    &               y*(tl_fac4+tl_y*(fac5+fac6)+y*(tl_fac5+tl_fac6)))  &
!    &                               )                                  &
!    &                 )

         ad_fac1=ad_fac1-ad_fTp
         ad_y=ad_y-(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))*ad_fTp
         ad_fac2=ad_fac2-y*ad_fTp
         ad_y=ad_y-y*(fac3+y*(fac4+y*(fac5+fac6)))*ad_fTp
         ad_fac3=ad_fac3-y*y*ad_fTp
         ad_y=ad_y-y*y*(fac4+y*(fac5+fac6))*ad_fTp
         ad_fac4=ad_fac4-y*y*y*ad_fTp
         ad_y=ad_y-y*y*y*(fac5+fac6)*ad_fTp
         ad_fac5=ad_fac5-y*y*y*y*ad_fTp
         ad_fac6=ad_fac6-y*y*y*y*ad_fTp

 
!        tl_fac6=tl_z*(q24+z*(q25+(q26+q27*z)*z))+                      &
!    &              z*(tl_z*(q25+(q26+q27*z)*z)+                        &
!    &                    z*((q26+q27*z)*tl_z+q27*z*tl_z)     )
         
         ad_z=ad_z+(q24+z*(q25+(q26+q27*z)*z))*ad_fac6
         ad_z=ad_z+z*(q25+(q26+q27*z)*z)*ad_fac6
         ad_z=ad_z+z*z*(q26+q27*z)*ad_fac6
         ad_z=ad_z+z*z*q27*z*ad_fac6
         ad_fac6=0.0d0

!        tl_fac5=tl_y*(q20*y*z+z*(q21+(q22+q23*z)*z))+                  &
!    &              y*(q20*(tl_y*z+y*tl_z)+                             &
!    &                     tl_z*(q21+(q22+q23*z)*z)+                    &
!    &                        z*((q22+q23*z)*tl_z+q23*tl_z*z) )
   
         ad_y=ad_y+(q20*y*z+z*(q21+(q22+q23*z)*z))*ad_fac5
         ad_y=ad_y+y*q20*z*ad_fac5
         ad_z=ad_z+y*q20*y*ad_fac5
         ad_z=ad_z+y*(q21+(q22+q23*z)*z)*ad_fac5
         ad_z=ad_z+y*z*(q22+q23*z)*ad_fac5
         ad_z=ad_z+y*z*q23*z*ad_fac5
         ad_fac5=0.0d0

!        tl_fac4=tl_z*(q16+z*(q17+z*(q18+q19*z)))+                      &
!    &              z*(tl_z*(q17+z*(q18+q19*z))+                        &
!    &                    z*(tl_z*(q18+q19*z)+z*q19*tl_z)   )

         ad_z=ad_z+(q16+z*(q17+z*(q18+q19*z)))*ad_fac4
         ad_z=ad_z+z*(q17+z*(q18+q19*z))*ad_fac4
         ad_z=ad_z+z*z*(q18+q19*z)*ad_fac4
         ad_z=ad_z+z*z*z*q19*ad_fac4
         ad_fac4=0.0d0
      

!        tl_fac3=tl_z*(q11+z*(q12+z*(q13+z*(q14+q15*z))))+              &
!    &              z*(tl_z*(q12+z*(q13+z*(q14+q15*z)))+                &
!    &                    z*(tl_z*(q13+z*(q14+q15*z))+                  &
!    &                          z*(tl_z*(q14+q15*z)+z*q15*tl_z)         &
!    &                      )                                           &
!    &                )

         ad_z=ad_z+(q11+z*(q12+z*(q13+z*(q14+q15*z))))*ad_fac3
         ad_z=ad_z+z*(q12+z*(q13+z*(q14+q15*z)))*ad_fac3
         ad_z=ad_z+z*z*(q13+z*(q14+q15*z))*ad_fac3
         ad_z=ad_z+z*z*z*(q14+q15*z)*ad_fac3
         ad_z=ad_z+z*z*z*z*q15*ad_fac3
         ad_fac3=0.0d0

!        tl_fac2=tl_z*(q6+z*(q7+z*(q8+(q9+q10*z)*z)))+                  &
!    &              z*(tl_z*(q7+z*(q8+(q9+q10*z)*z))+                   &
!    &                    z*(tl_z*(q8+(q9+q10*z)*z)+                    &
!    &                          z*((q9+q10*z)*tl_z+q10*z*tl_z)          &
!    &                      )                                           &
!    &                )

         ad_z=ad_z+(q6+z*(q7+z*(q8+(q9+q10*z)*z)))*ad_fac2
         ad_z=ad_z+z*(q7+z*(q8+(q9+q10*z)*z))*ad_fac2
         ad_z=ad_z+z*z*(q8+(q9+q10*z)*z)*ad_fac2
         ad_z=ad_z+z*z*z*(q9+q10*z)*ad_fac2
         ad_z=ad_z+z*z*z*q10*z*ad_fac2
         ad_fac2=0.0d0

!        tl_fac1=tl_z*(q1+z*(q2+z*(q3+z*(q4+q5*z))))+                   &
!    &              z*(tl_z*(q2+z*(q3+z*(q4+q5*z)))+                    &
!    &                    z*(tl_z*(q3+z*(q4+q5*z))+                     &
!    &                          z*(tl_z*(q4+q5*z)+tl_z*q5*z)            &
!    &                      )                                           &
!    &                )
        
         ad_z=ad_z+(q1+z*(q2+z*(q3+z*(q4+q5*z))))*ad_fac1
         ad_z=ad_z+z*(q2+z*(q3+z*(q4+q5*z)))*ad_fac1
         ad_z=ad_z+z*z*(q3+z*(q4+q5*z))*ad_fac1
         ad_z=ad_z+z*z*z*(q4+q5*z)*ad_fac1
         ad_z=ad_z+z*z*z*q5*z*ad_fac1
         ad_fac1=0.0d0

!       tl_yz = tl_y*z+y*tl_z

        ad_z=ad_z+y*ad_yz
        ad_y=ad_y+z*ad_yz
        ad_yz=0.0d0

      endif

!       tl_kTp = 22.6683558512829d0*(tl_z1-tl_z)

        ad_z1=ad_z1+22.6683558512829d0*ad_kTp
        ad_z=ad_z-22.6683558512829d0*ad_kTp
        ad_kTp=0.0d0

!       tl_hTp=tl_hTp-tl_y*(p1+y*(p2+p3*y))-y*(tl_y*(p2+p3*y)+y*p3*tl_y)

        ad_y=ad_y-(p1+y*(p2+p3*y))*ad_hTp
        ad_y=ad_y-y*(p2+p3*y)*ad_hTp
        ad_y=ad_y-y*y*p3*ad_hTp

!       tl_gTp=tl_gTp-tl_y*(h1+y*(h2+y*(h3+h4*y)))-                     &
!    &                   y*(tl_y*(h2+y*(h3+h4*y))+                      &
!    &                         y*(tl_y*(h3+h4*y)+y*h4*tl_y))

        ad_y=ad_y-(h1+y*(h2+y*(h3+h4*y)))*ad_gTp
        ad_y=ad_y-y*(h2+y*(h3+h4*y))*ad_gTp
        ad_y=ad_y-y*y*(h3+h4*y)*ad_gTp
        ad_y=ad_y-y*y*y*h4*ad_gTp

!       tl_fTp =tl_fTp-tl_y*(g2+y*(g3+y*(g4+y*(g5+(g6+g7*y)*y))))-      &
!    &                    y*(tl_y*(g3+y*(g4+y*(g5+(g6+g7*y)*y)))+       &
!    &                          y*(tl_y*(g4+y*(g5+(g6+g7*y)*y))+        &
!    &                                y*(tl_y*(g5+(g6+g7*y)*y)+         &
!    &                                y*((g6+g7*y)*tl_y+g7*y*tl_y)      &
!    &                                  )                               &
!    &                            )                                     &
!    &                      )

        ad_y=ad_y-(g2+y*(g3+y*(g4+y*(g5+(g6+g7*y)*y))))*ad_fTp
        ad_y=ad_y-y*(g3+y*(g4+y*(g5+(g6+g7*y)*y)))*ad_fTp
        ad_y=ad_y-y*y*(g4+y*(g5+(g6+g7*y)*y))*ad_fTp
        ad_y=ad_y-y*y*y*(g5+(g6+g7*y)*y)*ad_fTp
        ad_y=ad_y-y*y*y*y*(g6+g7*y)*ad_fTp
        ad_y=ad_y-y*y*y*y*y*g7*ad_fTp

!      tl_z = 1.0d-4*tl_p

       ad_p=ad_p+1.0d-4*ad_z
       ad_z=0.0d0

!      tl_y = 2.5d-2*tl_t

       ad_t=ad_t+2.5d-2*ad_y
       ad_y=0.0d0

!
      x2 = 2.5d-2*s
      x = dsqrt(x2)
      y = 2.5d-2*th0
      z = 1.0d-4*pr
      z1 = z

      if(z.ne.0.d0) then
        yz = y*z

        fac1=z*(f1+z*(f2+f3*z))
        fac2=y*(f4*z+f5*y*z)
        fac3=z*(f6+(f7+f8*z)*z)

!       tl_hTp =tl_hTp+tl_fac1+tl_y*(fac2+fac3)+y*(tl_fac2+tl_fac3)

        ad_fac1=ad_fac1+ad_hTp
        ad_y=ad_y+(fac2+fac3)*ad_hTp
        ad_fac2=ad_fac2+y*ad_hTp
        ad_fac3=ad_fac3+y*ad_hTp

!       tl_fac3=tl_z*(f6+(f7+f8*z)*z)+z*((f7+f8*z)*tl_z+f8*z*tl_z)

        ad_z=ad_z+(f6+(f7+f8*z)*z)*ad_fac3
        ad_z=ad_z+z*(f7+f8*z)*ad_fac3
        ad_z=ad_z+z*z*f8*ad_fac3
        ad_fac3=0.0d0

!       tl_fac2=tl_y*(f4*z+f5*y*z)+y*(f4*tl_z+f5*(y*tl_z+tl_y*z))

        ad_y=ad_y+(f4*z+f5*y*z)*ad_fac2
        ad_z=ad_z+y*f4*ad_fac2
        ad_z=ad_z+y*y*f5*ad_fac2
        ad_y=ad_y+y*z*f5*ad_fac2
        ad_fac2=0.0d0

!       tl_fac1=tl_z*(f1+z*(f2+f3*z))+z*(tl_z*(f2+f3*z)+z*f3*tl_z)

        ad_z=ad_z+(f1+z*(f2+f3*z))*ad_fac1
        ad_z=ad_z+z*(f2+f3*z)*ad_fac1
        ad_z=ad_z+z*z*f3*ad_fac1
        ad_fac1=0.0d0

!

        fac1=z*(e1+z*(e2+z*(e3+(e4+e5*z)*z)))
        fac2=z*(e6+z*(e7+z*(e8+z*(e9+e10*z))))
        fac3=y*z*(e11+z*(e12+e13*z))
        fac4=z*(e14+z*(e15+z*(e16+e17*z)))

!       tl_gTp = tl_gTp+tl_fac1+tl_y*(fac2+y*(fac3+fac4))+              &
!     &              y*(tl_fac2+tl_y*(fac3+fac4)+y*(tl_fac3+tl_fac4))

        ad_fac1=ad_fac1+ad_gTp
        ad_y=ad_y+(fac2+y*(fac3+fac4))*ad_gTp
        ad_fac2=ad_fac2+y*ad_gTp
        ad_y=ad_y+y*(fac3+fac4)*ad_gTp
        ad_fac3=ad_fac3+y*y*ad_gTp
        ad_fac4=ad_fac4+y*y*ad_gTp

!       tl_fac4=tl_z*(e14+z*(e15+z*(e16+e17*z)))+                       &
!    &             z*(tl_z*(e15+z*(e16+e17*z))+                         &
!    &                   z*(tl_z*(e16+e17*z)+z*e17*tl_z))

        ad_z=ad_z+(e14+z*(e15+z*(e16+e17*z)))*ad_fac4
        ad_z=ad_z+z*(e15+z*(e16+e17*z))*ad_fac4
        ad_z=ad_z+z*z*(e16+e17*z)*ad_fac4
        ad_z=ad_z+z*z*z*e17*ad_fac4
        ad_fac4=0.0d0

!       tl_fac3=(tl_y*z+y*tl_z)*(e11+z*(e12+e13*z))+                    &
!    &           y*z*(tl_z*(e12+e13*z)+z*e13*tl_z)

        ad_y=ad_y+z*(e11+z*(e12+e13*z))*ad_fac3
        ad_z=ad_z+y*(e11+z*(e12+e13*z))*ad_fac3
        ad_z=ad_z+y*z*(e12+e13*z)*ad_fac3
        ad_z=ad_z+y*z*z*e13*ad_fac3
        ad_fac3=0.0d0

!       tl_fac2=tl_z*(e6+z*(e7+z*(e8+z*(e9+e10*z))))+                   &
!    &             z*(tl_z*(e7+z*(e8+z*(e9+e10*z)))+                    &
!    &                   z*(tl_z*(e8+z*(e9+e10*z))+                     &
!    &                         z*(tl_z*(e9+e10*z)+z*e10*tl_z))          &
!    &               )

       ad_z=ad_z+(e6+z*(e7+z*(e8+z*(e9+e10*z))))*ad_fac2
       ad_z=ad_z+z*(e7+z*(e8+z*(e9+e10*z)))*ad_fac2
       ad_z=ad_z+z*z*(e8+z*(e9+e10*z))*ad_fac2
       ad_z=ad_z+z*z*z*(e9+e10*z)*ad_fac2
       ad_z=ad_z+z*z*z*z*e10*ad_fac2
       ad_fac2=0.0d0

!       tl_fac1=tl_z*(e1+z*(e2+z*(e3+(e4+e5*z)*z)))+                    &
!    &             z*(tl_z*(e2+z*(e3+(e4+e5*z)*z))+                     &
!    &                   z*(tl_z*(e3+(e4+e5*z)*z)+                      &
!    &                         z*((e4+e5*z)*tl_z+e5*z*tl_z))            &
!    &               )

        ad_z=ad_z+(e1+z*(e2+z*(e3+(e4+e5*z)*z)))*ad_fac1
        ad_z=ad_z+z*(e2+z*(e3+(e4+e5*z)*z))*ad_fac1
        ad_z=ad_z+z*z*(e3+(e4+e5*z)*z)*ad_fac1
        ad_z=ad_z+z*z*z*(e4+e5*z)*ad_fac1
        ad_z=ad_z+z*z*z*z*e5*ad_fac1
        ad_fac1=0.0d0

!

        fac1=z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))
        fac2=z*(d6+z*(d7+z*(d8+(d9+d10*z)*z)))
        fac3=z*(d11+z*(d12+z*(d13+z*(d14+d15*z))))
        fac4=z*(d16+z*(d17+z*(d18+d19*z)))
        fac5=y*(d20*y*z+z*(d21+(d22+d23*z)*z))
        fac6=z*(d24+z*(d25+(d26+d27*z)*z))

!       tl_fTp = tl_fTp+tl_fac1+                                        &
!    &            tl_y*(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))+          &
!    &               y*(tl_fac2+tl_y*(fac3+y*(fac4+y*(fac5+fac6)))+     &
!    &                    y*(tl_fac3+tl_y*(fac4+y*(fac5+fac6))+         &
!    &                y*(tl_fac4+tl_y*(fac5+fac6)+y*(tl_fac5+tl_fac6))) &
!    &                   )

        ad_fac1=ad_fac1+ad_fTp
        ad_y=ad_y+(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))*ad_fTp
        ad_fac2=ad_fac2+y*ad_fTp
        ad_y=ad_y+y*(fac3+y*(fac4+y*(fac5+fac6)))*ad_fTp
        ad_fac3=ad_fac3+y*y*ad_fTp
        ad_y=ad_y+y*y*(fac4+y*(fac5+fac6))*ad_fTp
        ad_fac4=ad_fac4+y*y*y*ad_fTp
        ad_y=ad_y+y*y*y*(fac5+fac6)*ad_fTp
        ad_fac5=ad_fac5+y*y*y*y*ad_fTp
        ad_fac6=ad_fac6+y*y*y*y*ad_fTp

!       tl_fac6=tl_z*(d24+z*(d25+(d26+d27*z)*z))+                       &
!    &             z*(tl_z*(d25+(d26+d27*z)*z)+                         &
!    &                   z*((d26+d27*z)*tl_z+d27*z*tl_z)                &
!    &               )

        ad_z=ad_z+(d24+z*(d25+(d26+d27*z)*z))*ad_fac6
        ad_z=ad_z+z*(d25+(d26+d27*z)*z)*ad_fac6
        ad_z=ad_z+z*z*(d26+d27*z)*ad_fac6
        ad_z=ad_z+z*z*z*d27*ad_fac6
        ad_fac6=0.0d0

!       tl_fac5=tl_y*(d20*y*z+z*(d21+(d22+d23*z)*z))+                   &
!    &             y*(d20*tl_y*z+d20*y*tl_z+tl_z*(d21+(d22+d23*z)*z)+   &
!    &                      z*((d22+d23*z)*tl_z+d23*z*tl_z)             &
!    &               )

        ad_y=ad_y+(d20*y*z+z*(d21+(d22+d23*z)*z))*ad_fac5
        ad_y=ad_y+y*d20*z*ad_fac5
        ad_z=ad_z+y*d20*y*ad_fac5
        ad_z=ad_z+y*(d21+(d22+d23*z)*z)*ad_fac5
        ad_z=ad_z+y*z*(d22+d23*z)*ad_fac5
        ad_z=ad_z+y*z*z*d23*ad_fac5
        ad_fac5=0.0d0

!       tl_fac4=tl_z*(d16+z*(d17+z*(d18+d19*z)))+                       &
!    &             z*(tl_z*(d17+z*(d18+d19*z))+                         &
!    &                   z*(tl_z*(d18+d19*z)+z*d18*tl_z)                &
!    &               )

       ad_z=ad_z+(d16+z*(d17+z*(d18+d19*z)))*ad_fac4
       ad_z=ad_z+z*(d17+z*(d18+d19*z))*ad_fac4
       ad_z=ad_z+z*z*(d18+d19*z)*ad_fac4
       ad_z=ad_z+z*z*z*d18*ad_fac4
       ad_fac4=0.0d0

!       tl_fac3=tl_z*(d11+z*(d12+z*(d13+z*(d14+d15*z))))+               &
!    &             z*(tl_z*(d12+z*(d13+z*(d14+d15*z)))+                 &
!    &                   z*(tl_z*(d13+z*(d14+d15*z))+                   &
!    &                         z*(tl_z*(d14+d15*z)+z*d15*tl_z)          &
!    &                     )                                            &
!    &               )

        ad_z=ad_z+(d11+z*(d12+z*(d13+z*(d14+d15*z))))*ad_fac3
        ad_z=ad_z+z*(d12+z*(d13+z*(d14+d15*z)))*ad_fac3
        ad_z=ad_z+z*z*(d13+z*(d14+d15*z))*ad_fac3
        ad_z=ad_z+z*z*z*(d14+d15*z)*ad_fac3
        ad_z=ad_z+z*z*z*z*d15*ad_fac3
        ad_fac3=0.0d0

!       tl_fac2=tl_z*(d6+z*(d7+z*(d8+(d9+d10*z)*z)))+                   &
!    &             z*(tl_z*(d7+z*(d8+(d9+d10*z)*z))+                    &
!    &             z*(tl_z*(d8+(d9+d10*z)*z)+                           &
!    &                   z*((d9+d10*z)*tl_z+d10*z*tl_z)) )

        ad_z=ad_z+(d6+z*(d7+z*(d8+(d9+d10*z)*z)))*ad_fac2
        ad_z=ad_z+z*(d7+z*(d8+(d9+d10*z)*z))*ad_fac2
        ad_z=ad_z+z*z*(d8+(d9+d10*z)*z)*ad_fac2
        ad_z=ad_z+z*z*z*(d9+d10*z)*ad_fac2
        ad_z=ad_z+z*z*z*z*d10*ad_fac2
        ad_fac2=0.0d0

!       tl_fac1=tl_z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))+                    &
!    &                     z*(tl_z*(d2+z*(d3+z*(d4+d5*z)))+             &
!    &                           z*(tl_z*(d3+z*(d4+d5*z))+              &
!    &                                 z*(tl_z*(d4+d5*z)+z*d5*tl_z)     &
!    &                             )                                    &
!    &                       )

        ad_z=ad_z+(d1+z*(d2+z*(d3+z*(d4+d5*z))))*ad_fac1
        ad_z=ad_z+z*(d2+z*(d3+z*(d4+d5*z)))*ad_fac1
        ad_z=ad_z+z*z*(d3+z*(d4+d5*z))*ad_fac1
        ad_z=ad_z+z*z*z*(d4+d5*z)*ad_fac1
        ad_z=ad_z+z*z*z*z*d5*ad_fac1
        ad_fac1=0.0d0

!       tl_yz = tl_y*z+y*tl_z

        ad_y=ad_y+z*ad_yz
        ad_z=ad_z+y*ad_yz
        ad_yz=0.0d0

      endif

!     tl_hTp = tl_y*(c1+y*(c2+c3*y))+y*(tl_y*(c2+c3*y)+y*c3*tl_y)

      ad_y=ad_y+(c1+y*(c2+c3*y))*ad_hTp
      ad_y=ad_y+y*(c2+c3*y)*ad_hTp
      ad_y=ad_y+y*y*c3*ad_hTp
      ad_hTp=0.0d0

!     tl_gTp = tl_y*(b1+y*(b2+y*(b3+b4*y)))+                            &
!    &            y*(tl_y*(b2+y*(b3+b4*y))+                             &
!    &                  y*(tl_y*(b3+b4*y)+y*b4*tl_y)                    &
!    &              )

      ad_y=ad_y+(b1+y*(b2+y*(b3+b4*y)))*ad_gTp
      ad_y=ad_y+y*(b2+y*(b3+b4*y))*ad_gTp
      ad_y=ad_y+y*y*(b3+b4*y)*ad_gTp
      ad_y=ad_y+y*y*y*b4*ad_gTp
      ad_gTp=0.0d0

!     tl_fTp = tl_y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))+              &
!    &         y*(a2+tl_y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))+               &
!    &               y*(tl_y*(a4+y*(a5+(a6+a7*y)*y))+                   &
!    &                    y*(tl_y*(a5+(a6+a7*y)*y)+                     &
!                              y*((a6+a7*y)*tl_y+a7*y*tl_y)             &
!    &                      )                                           &
!    &                  )                                               &
!    &           )

      ad_y=ad_y+(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))*ad_fTp
      ad_y=ad_y+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))*ad_fTp
      ad_y=ad_y+y*y*(a4+y*(a5+(a6+a7*y)*y))*ad_fTp
      ad_y=ad_y+y*y*y*(a5+(a6+a7*y)*y)*ad_fTp
      ad_y=ad_y+y*y*y*y*(a6+a7*y)*ad_fTp
      ad_y=ad_y+y*y*y*y*y*a7*ad_fTp
      ad_fTp=0.0d0

!     tl_z1 = tl_z

      ad_z=ad_z+ad_z1
      ad_z1=0.0d0

!     tl_z = 1.0d-4*tl_pr

      ad_pr=ad_pr+1.0d-4*ad_z
      ad_z=0.0d0

!     tl_y = 2.5d-2*tl_th0

      ad_th0=ad_th0+2.5d-2*ad_y
      ad_y=0.0d0

!     tl_x = 0.5d0*tl_x2/x

      ad_x2=ad_x2+0.5d0*ad_x/x
      ad_x=0.0d0

!     tl_x2 = 2.5d-2*tl_s

      ad_s=ad_s+2.5d-2*ad_x2
      ad_x2=0.0d0

      if(iflag.eq.1)ad_entropy_diff_F=ad_s
      if(iflag.eq.2)ad_entropy_diff_F=ad_t
      if(iflag.eq.3)ad_entropy_diff_F=ad_p
      if(iflag.eq.4)ad_entropy_diff_F=ad_th0
      if(iflag.eq.5)ad_entropy_diff_F=ad_pr

      RETURN
      END FUNCTION ad_entropy_diff_F

      REAL FUNCTION ad_de_dt_F(iflag,ain,s,t,p)
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

!
!  Clear adjoint variables.
!
      ad_x2=0.0d0
      ad_x=0.0d0
      ad_y=0.0d0
      ad_z=0.0d0
      ad_p=0.0d0
      ad_t=0.0d0
      ad_s=0.0d0
      ad_fac1=0.0d0
      ad_fac2=0.0d0
      ad_fac3=0.0d0
      ad_fac4=0.0d0
      ad_fac5=0.0d0
      ad_fac6=0.0d0
      ad_fac7=0.0d0
      ad_fac8=0.0d0
      ad_fac9=0.0d0
      ad_de_dt=0.0d0

      x2 = 2.5d-2*s
      x = dsqrt(x2)
      y = 2.5d-2*t
      z = 1.0d-4*p

      a1=24715.571866078d0
      a2=-1858.920033948178d0
      a3=317.440355256842d0
      a4=-405.1392883572282d0
      a5=202.6815298758072d0
      a6=1562.563716288858d0
      a7=-1165.8752731900836d0
      a8=348.7487684426d0
      a9=-4420.4472249096725d0
      a10=1778.231237203896d0
      a11=-1160.5182516851419d0
      a12=569.531539542516d0
      a13=-128.13429152494615d0

!     tl_de_dt_F = 6.25d-4*tl_de_dt
      ad_de_dt=ad_de_dt+6.25d-4*ain

      if(z.ne.0.d0) then

        b1=-2910.0729080936d0
        b2=1721.528607567954d0
        b3=-4165.4688847996085d0
        b4=3572.7449038462437d0
        b5=-766.116132004952d0
        b6=2761.9195908075417d0
        b7=-2814.78225133626d0
        b8=4035.04669887042d0
        b9=-4775.621344883664d0
        b10=3892.3662123519d0
        b11=-1905.341809925355d0
        b12=404.50541014508605d0
        b13=1513.116771538718d0
        b14=-674.819060538734d0
        b15=108.3834525034224d0
        b16=1229.337851789418d0
        b17=-896.713693665072d0
        b18=-2996.162344914912d0
        b19=3621.784567462512d0
        b20=-2410.4130980405d0
        b21=668.691951421377d0
        b22=-546.959324647056d0
        b23=356.629112415276d0
        b24=-51.2796974779828d0
        b25=-681.370187043564d0
        b26=437.84750280190565d0
        b27=1437.2719839264719d0
        b28=-1826.356460806092d0
        b29=1105.446104680304d0
        b30=-245.11816254543362d0
        b31=111.1208127634436d0
        b32=-88.4080716616d0
        b33=66.7696405958478d0
        b34=-292.8075111563232d0
        b35=316.49805267936244d0
        b36=-129.6381336154442d0
        b37=-8.68841343834394d0
        b38=15.84003094423364d0
        b39=9.978426372534301d0

        fac1=x2*(b2+y*(b3+b4*y)+x*(b5+(b6+b7*y)*y))
        fac2=y*(b8+y*(b9+y*(b10+y*(b11+b12*y))))
        fac3=x2*(b14+b15*x+(b16+b17*y)*y)
        fac4=y*(b18+y*(b19+y*(b20+b21*y)))
        fac5=x2*(b23+b24*x+y*(b25+b26*y))
        fac6=y*(b27+y*(b28+(b29+b30*y)*y))
        fac7=x2*(b32+b33*y)
        fac8=y*(b34+(b35+b36*y)*y)
        fac9=(b37+b38*x2+b39*y)*z

!       tl_de_dt = tl_de_dt +                                           &
!    &       tl_z*(b1+fac1+fac2+z*(b13+fac3+fac4+z*(b22+fac5+fac6+      &
!    &                     z*(b31+fac7+fac8+fac9))))+
!    &          z*(tl_fac1+tl_fac2+tl_z*(b13+fac3+fac4+z*(b22+fac5+fac6+&
!    &                     z*(b31+fac7+fac8+fac9)))+                    &
!    &                                z*(tl_fac3+tl_fac4+               &
!    &                  tl_z*(b22+fac5+fac6+z*(b31+fac7+fac8+fac9))+    &
!    &                   z*(tl_fac5+tl_fac6+tl_z*(b31+fac7+fac8+fac9)+  &
!    &                          z*(tl_fac7+tl_fac8+tl_fac9)  )  )  )

        ad_z=ad_z+(b1+fac1+fac2+z*(b13+fac3+fac4+z*(b22+fac5+fac6+      &
     &                     z*(b31+fac7+fac8+fac9))))*ad_de_dt
        ad_fac1=ad_fac1+z*ad_de_dt
        ad_fac2=ad_fac2+z*ad_de_dt
        ad_z=ad_z+z*(b13+fac3+fac4+z*(b22+fac5+fac6+                    &
     &                     z*(b31+fac7+fac8+fac9)))*ad_de_dt
        ad_fac3=ad_fac3+z*z*ad_de_dt
        ad_fac4=ad_fac4+z*z*ad_de_dt
        ad_z=ad_z+z*z*(b22+fac5+fac6+z*(b31+fac7+fac8+fac9))*ad_de_dt
        ad_fac5=ad_fac5+z*z*z*ad_de_dt
        ad_fac6=ad_fac6+z*z*z*ad_de_dt
        ad_z=ad_z+z*z*z*(b31+fac7+fac8+fac9)*ad_de_dt
        ad_fac7=ad_fac7+z*z*z*z*ad_de_dt
        ad_fac8=ad_fac8+z*z*z*z*ad_de_dt
        ad_fac9=ad_fac9+z*z*z*z*ad_de_dt

!       tl_fac9=(b37+b38*x2+b39*y)*tl_z+(b38*tl_x2+b39*tl_y)*z
        ad_z=ad_z+(b37+b38*x2+b39*y)*ad_fac9
        ad_x2=ad_x2+z*b38*ad_fac9
        ad_y=ad_y+z*b39*ad_fac9

!       tl_fac8=tl_y*(b34+(b35+b36*y)*y)+y*((b35+b36*y)*tl_y+b36*y*tl_y)
        ad_y=ad_y+(b34+(b35+b36*y)*y)*ad_fac8
        ad_y=ad_y+y*(b35+b36*y)*ad_fac8
        ad_y=ad_y+y*y*b36*ad_fac8
        ad_fac8=0.0d0

!       tl_fac7=tl_x2*(b32+b33*y)+x2*b33*tl_y
        ad_x2=ad_x2+(b32+b33*y)*ad_fac7
        ad_y=ad_y+x2*b33*ad_fac7
        ad_fac7=0.0d0

!       tl_fac6=tl_y*(b27+y*(b28+(b29+b30*y)*y))+                       &
!    &             y*(tl_y*(b28+(b29+b30*y)*y)+                         &
!    &                   y*((b29+b30*y)*tl_y+b30*y*tl_y)  )

        ad_y=ad_y+(b27+y*(b28+(b29+b30*y)*y))*ad_fac6
        ad_y=ad_y+y*(b28+(b29+b30*y)*y)*ad_fac6
        ad_y=ad_y+y*y*(b29+b30*y)*ad_fac6
        ad_y=ad_y+y*y*y*b30*ad_fac6
        ad_fac6=0.0d0

!       tl_fac5=tl_x2*(b23+b24*x+y*(b25+b26*y))+                        &
!    &             x2*(b24*tl_x+tl_y*(b25+b26*y)+y*b26*tl_y)

        ad_x2=ad_x2+(b23+b24*x+y*(b25+b26*y))*ad_fac5
        ad_x=ad_x+x2*b24*ad_fac5
        ad_y=ad_y+x2*(b25+b26*y)*ad_fac5
        ad_y=ad_y+x2*y*b26*ad_fac5
        ad_fac5=0.0d0
    
!       tl_fac4=tl_y*(b18+y*(b19+y*(b20+b21*y)))+                       &
!    &             y*(tl_y*(b19+y*(b20+b21*y))+                         &
!    &                   y*(tl_y*(b20+b21*y)+y*b21*tl_y) )
      
        ad_y=ad_y+(b18+y*(b19+y*(b20+b21*y)))*ad_fac4
        ad_y=ad_y+y*(b19+y*(b20+b21*y))*ad_fac4
        ad_y=ad_y+y*y*(b20+b21*y)*ad_fac4
        ad_y=ad_y+y*y*y*b21*ad_fac4
        ad_fac4=0.0d0

!       tl_fac3=tl_x2*(b14+b15*x+(b16+b17*y)*y)+                        &
!    &             x2*(b15*tl_x+(b16+b17*y)*tl_y+b17*y*tl_y)

        ad_x2=ad_x2+(b14+b15*x+(b16+b17*y)*y)*ad_fac3
        ad_x=ad_x+x2*b15*ad_fac3
        ad_y=ad_y+x2*(b16+b17*y)*ad_fac3
        ad_y=ad_y+x2*y*b17*ad_fac3
        ad_fac3=0.0d0

!       tl_fac2=tl_y*(b8+y*(b9+y*(b10+y*(b11+b12*y))))+                 &
!    &             y*(tl_y*(b9+y*(b10+y*(b11+b12*y)))+                  &
!    &                   y*(tl_y*(b10+y*(b11+b12*y))+                   &
!    &                         y*(tl_y*(b11+b12*y)+y*b12*tl_y) ) )

        ad_y=ad_y+(b8+y*(b9+y*(b10+y*(b11+b12*y))))*ad_fac2
        ad_y=ad_y+y*(b9+y*(b10+y*(b11+b12*y)))*ad_fac2
        ad_y=ad_y+y*y*(b10+y*(b11+b12*y))*ad_fac2
        ad_y=ad_y+y*y*y*(b11+b12*y)*ad_fac2
        ad_y=ad_y+y*y*y*y*b12*ad_fac2
        ad_fac2=0.0d0

!       tl_fac1=tl_x2*(b2+y*(b3+b4*y)+x*(b5+(b6+b7*y)*y))+              &
!    &             x2*(tl_y*(b3+b4*y)+y*b4*tl_y+tl_x*(b5+(b6+b7*y)*y)+  &
!    &                                    x*((b6+b7*y)*tl_y+b7*y*tl_y))

        ad_x2=ad_x2+(b2+y*(b3+b4*y)+x*(b5+(b6+b7*y)*y))*ad_fac1
        ad_y=ad_y+x2*(b3+b4*y)*ad_fac1
        ad_y=ad_y+x2*y*b4*ad_fac1
        ad_x=ad_x+x2*(b5+(b6+b7*y)*y)*ad_fac1
        ad_y=ad_y+x2*x*(b6+b7*y)*ad_fac1
        ad_y=ad_y+x2*x*y*b7*ad_fac1
        ad_fac1=0.0d0

      end if

      fac1=x*(a3+y*(a4+a5*y))
      fac2=y*(a6+y*(a7+a8*y))
      fac3=y*(a10+y*(a11+(a12+a13*y)*y))

!     tl_de_dt =tl_x2*(a2+fac1+fac2)+x2*(tl_fac1+tl_fac2)+              &
!    &           tl_y*(a9+fac3)+y*tl_fac3

      ad_x2=ad_x2+(a2+fac1+fac2)*ad_de_dt
      ad_fac1=ad_fac1+x2*ad_de_dt
      ad_fac2=ad_fac2+x2*ad_de_dt
      ad_y=ad_y+(a9+fac3)*ad_de_dt
      ad_fac3=ad_fac3+y*ad_de_dt
      ad_de_dt=0.0d0

!     tl_fac3=tl_y*(a10+y*(a11+(a12+a13*y)*y))+                         &
!    &           y*(tl_y*(a11+(a12+a13*y)*y)+                           &
!    &                 y*((a12+a13*y)*tl_y+a13*y*tl_y) )

      ad_y=ad_y+(a10+y*(a11+(a12+a13*y)*y))*ad_fac3
      ad_y=ad_y+y*(a11+(a12+a13*y)*y)*ad_fac3
      ad_y=ad_y+y*y*(a12+a13*y)*ad_fac3
      ad_y=ad_y+y*y*y*a13*ad_fac3
      ad_fac3=0.0d0

!     tl_fac2=tl_y*(a6+y*(a7+a8*y))+y*(tl_y*(a7+a8*y)+y*a8*tl_y)
      ad_y=ad_y+(a6+y*(a7+a8*y))*ad_fac2
      ad_y=ad_y+y*(a7+a8*y)*ad_fac2
      ad_y=ad_y+y*y*a8*ad_fac2
      ad_fac2=0.0d0

!     tl_fac1=tl_x*(a3+y*(a4+a5*y))+x*(tl_y*(a4+a5*y)+y*a5*tl_y)
      ad_x=ad_x+(a3+y*(a4+a5*y))*ad_fac1
      ad_y=ad_y+x*(a4+a5*y)*ad_fac1
      ad_y=ad_y+x*y*a5*ad_fac1
      ad_fac1=0.0d0

!     tl_z = 1.0d-4*tl_p
      ad_p=ad_p+1.0d-4*ad_z
      ad_z=0.0d0

!     tl_y = 2.5d-2*tl_t
      ad_t=ad_t+2.5d-2*ad_y
      ad_y=0.0d0

!     tl_x=0.5d0*tl_x2/x
      ad_x2=ad_x2+0.5d0*ad_x/x
      ad_x=0.0d0

!     tl_x2 = 2.5d-2*tl_s
      ad_s=ad_s+2.5d-2*ad_x2
      ad_x2=0.0d0

      if(iflag.eq.1)ad_de_dt_F=ad_s
      if(iflag.eq.2)ad_de_dt_F=ad_t
      if(iflag.eq.3)ad_de_dt_F=ad_p


      RETURN
      END FUNCTION ad_de_dt_F

      REAL FUNCTION ad_t_from_ct(iflag,ain,s,ct,p)
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

!
!   Clear adjoint variables.
!
      ad_s=0.0d0
      ad_theta=0.0d0
      ad_p=0.0d0
      ad_ct=0.0d0

      pr0 = 0.d0
      theta = theta_from_ct(s,ct)

!AMM I THINK THIS IS WRONG      tl_t_from_ct = tl_theta_from_t(s,tl_s,theta,tl_theta,pr0,tl_p,p)
!     tl_t_from_ct = tl_theta_from_t(s,tl_s,theta,tl_theta,p,tl_p,pr0)
      
      ad_s=ad_s+ad_theta_from_t(1,ain,s,theta,p,pr0)
      ad_theta=ad_theta+ad_theta_from_t(2,ain,s,theta,p,pr0)
      ad_p=ad_p+ad_theta_from_t(3,ain,s,theta,p,pr0)

!     tl_theta = tl_theta_from_ct(s,tl_s,ct,tl_ct)

      ad_s=ad_s+ad_theta_from_ct(1,ad_theta,s,ct)
      ad_ct=ad_ct+ad_theta_from_ct(2,ad_theta,s,ct)
      ad_theta=0.0d0

      if(iflag.eq.1)ad_t_from_ct=ad_s
      if(iflag.eq.2)ad_t_from_ct=ad_ct
      if(iflag.eq.3)ad_t_from_ct=ad_p

      RETURN
      END FUNCTION ad_t_from_ct

      REAL FUNCTION ad_theta_from_t(iflag,ain,s,t,p,pr)
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
      real*8 theta_s(2), th0_s(2), dentropy_s(2)
      real*8 de_dt_s1(2), de_dt_s2(2)
      integer :: nloops, n

!
!  Clear adjoint variables
!
      ad_th0=0.0d0
      ad_t=0.0d0
      ad_p=0.0d0
      ad_pr=0.0d0
      ad_s=0.0d0
      ad_theta=0.0d0
      ad_de_dt=0.0d0
      ad_dentropy=0.0d0

      a1=8.65483913395442d-6
      a2=1.41636299744881d-6
      a3=7.38286467135737d-9
      a4=-8.38241357039698d-6
      a5=2.83933368585534d-8
      a6=1.77803965218656d-8
      a7=1.71155619208233d-10

      th0 = t+(p-pr)*(a1-s*a2-(p+pr)*a3+t*(a4+s*a5+t*a6+(p+pr)*a7))

      de_dt = 13.6d0
      nloops = 2                                      ! default
!
!  Save th0 and theta for later
!
      do n=1,nloops
       th0_s(n) = th0
       de_dt_s1(n)=de_dt
       dentropy = entropy_diff_F(s,t,p,th0,pr)
       dentropy_s(n)=dentropy
       theta = th0-dentropy/de_dt
       theta_s(n) = 0.5d0*(theta+th0)
       de_dt = de_dt_F(s,theta,pr)
       de_dt_s2(n)=de_dt
       theta = th0-dentropy/de_dt
       th0 = theta
      end do

!     tl_theta_from_t = tl_th0

      ad_th0=ad_th0+ain

!    NOTE: nloops = 1 gives theta with a maximum error of 5.48x10^-06
!          nloops = 2 gives theta with a maximum error of 2.84x10^-14
      do n=nloops,1,-1

!      tl_th0 = tl_theta
       ad_theta=ad_theta+ad_th0
       ad_th0=0.0d0

!      tl_theta = tl_th0-tl_dentropy/de_dt+                             &
!    &                            tl_de_dt*dentropy/(de_dt*de_dt)
       ad_th0=ad_th0+ad_theta
       ad_dentropy=ad_dentropy-ad_theta/de_dt_s2(n)
       ad_de_dt=ad_de_dt+dentropy_s(n)*ad_theta/                        &
     &                              (de_dt_s2(n)*de_dt_s2(n))
       ad_theta=0.0d0

!      tl_de_dt = tl_de_dt_F(s,tl_s,theta,tl_theta,pr,tl_pr)
       ad_s=ad_s+ad_de_dt_F(1,ad_de_dt,s,theta_s(n),pr)
       ad_theta=ad_theta+ad_de_dt_F(2,ad_de_dt,s,theta_s(n),pr)
       ad_pr=ad_pr+ad_de_dt_F(3,ad_de_dt,s,theta_s(n),pr)
       ad_de_dt=0.0d0

!      tl_theta = 0.5d0*(tl_theta+tl_th0)
       ad_th0=ad_th0+0.5d0*ad_theta
       ad_theta=0.5d0*ad_theta

!      tl_theta = tl_th0-tl_dentropy/de_dt+                             &
!    &                           tl_de_dt*dentropy/(de_dt*de_dt)
       ad_th0=ad_th0+ad_theta
       ad_dentropy=ad_dentropy-ad_theta/de_dt_s1(n)
       ad_de_dt=ad_de_dt+dentropy_s(n)*ad_theta/                        &
     &                              (de_dt_s1(n)*de_dt_s1(n))
       ad_theta=0.0d0

!      tl_dentropy = tl_entropy_diff_F(s,tl_s,t,tl_t,p,tp_p,th0,tl_th0, &
!                                      pr,tl_pr)
       ad_s=ad_s+ad_entropy_diff_F(1,ad_dentropy,s,t,p,th0_s(n),pr)
       ad_t=ad_t+ad_entropy_diff_F(2,ad_dentropy,s,t,p,th0_s(n),pr)
       ad_p=ad_p+ad_entropy_diff_F(3,ad_dentropy,s,t,p,th0_s(n),pr)
       ad_th0=ad_th0+ad_entropy_diff_F(4,ad_dentropy,s,t,p,th0_s(n),pr)
       ad_pr=ad_pr+ad_entropy_diff_F(5,ad_dentropy,s,t,p,th0_s(n),pr)
       ad_dentropy=0.0d0

      end do

!     tl_de_dt = 0.0d0
      ad_de_dt = 0.0d0

!     tl_th0 = tl_t+tl_p*(a1-s*a2-(p+pr)*a3+t*(a4+s*a5+t*a6+(p+pr)*a7))+&
!    &          (p-pr)*(-tl_s*a2-tl_p*a3+tl_t*(a4+s*a5+t*a6+(p+pr)*a7)+ &
!    &                               t*(tl_s*a5+tl_t*a6+tl_p*a7))

      ad_t=ad_t+ad_th0
      ad_p=ad_p+(a1-s*a2-(p+pr)*a3+t*(a4+s*a5+t*a6+(p+pr)*a7))*ad_th0
      ad_s=ad_s-(p-pr)*a2*ad_th0
      ad_p=ad_p-(p-pr)*a3*ad_th0
      ad_t=ad_t+(p-pr)*(a4+s*a5+t*a6+(p+pr)*a7)*ad_th0
      ad_s=ad_s+(p-pr)*t*a5*ad_th0
      ad_t=ad_t+(p-pr)*t*a6*ad_th0
      ad_p=ad_p+(p-pr)*t*a7*ad_th0
      ad_th0=0.0d0

      if(iflag.eq.1)ad_theta_from_t=ad_s
      if(iflag.eq.2)ad_theta_from_t=ad_t
      if(iflag.eq.3)ad_theta_from_t=ad_p

      RETURN
      END FUNCTION ad_theta_from_t

      REAL FUNCTION ad_fp_theta(iflag,ain,s,p,sat)
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

!
!   Clear adjoint variables.
!
      ad_s=0.0d0
      ad_p=0.0d0
      ad_tf_num=0.0d0
      ad_tf_den=0.0d0
      ad_sqrts=0.0d0

      a1=2.5180516744541290d-03
      a2=-5.8545863698926184d-02
      a3=2.2979985780124325d-03
      a4=3.0086338218235500d-04
      a5=-7.0023530029351803d-04
      a6=8.4149607219833806d-09
      a7=1.1845857563107403d-11

      b1=1.0000000000000000d+00
      b2=-3.8493266309172074d-05
      b3=9.1686537446749641d-10
      b4=1.3632481944285909d-06

      if(sat.eq.'air-sat') then
!        tl_fp_theta=tl_fp_theta+tl_s*1.428571428571429d-05
         ad_s=ad_s+1.428571428571429d-05*ain
      elseif(sat.eq.'air-free') then
      continue
      else
      print *, '***  Error in fp_theta.f90: invalid third argument  ***'
      print *
      stop
      endif

      sqrts = dsqrt(s)
      tf_num =a1+s*(a2+sqrts*(a3-sqrts*a4))+p*(a5+p*(a6+s*a7))
      tf_den =b1+p*(b2+p*b3)+s*s*sqrts*b4
      fp_theta = tf_num/tf_den

!     tl_fp_theta = tl_tf_num/tf_den-tl_tf_den*fp_theta/tf_den
      ad_tf_num=ad_tf_num+ain/tf_den
      ad_tf_den=ad_tf_den-ain*fp_theta/tf_den

!     tl_tf_den =tl_p*(b2+p*b3)+p*tl_p*b3+2.0d0*tl_s*s*sqrts*b4+       &
!    &                                    s*s*tl_sqrts*b4
      ad_p=ad_p+(b2+p*b3)*ad_tf_den
      ad_p=ad_p+p*b3*ad_tf_den
      ad_s=ad_s+2.0d0*s*sqrts*b4*ad_tf_den
      ad_sqrts=ad_sqrts+s*s*b4*ad_tf_den
      ad_tf_den=0.0d0

!     tl_tf_num=tl_s*(a2+sqrts*(a3-sqrts*a4))+                          &
!    &             s*(tl_sqrts*(a3-sqrts*a4)-tl_sqrts*sqrts*a4)+        &
!    &          tl_p*(a5+p*(a6+s*a7))+p*(tl_p*(a6+s*a7)+p*tl_s*a7)

      ad_s=ad_s+(a2+sqrts*(a3-sqrts*a4))*ad_tf_num
      ad_sqrts=ad_sqrts+s*(a3-sqrts*a4)*ad_tf_num
      ad_sqrts=ad_sqrts-s*sqrts*a4*ad_tf_num
      ad_p=ad_p+(a5+p*(a6+s*a7))*ad_tf_num
      ad_p=ad_p+p*(a6+s*a7)*ad_tf_num
      ad_s=ad_s+p*p*a7*ad_tf_num
      ad_tf_num=0.0d0

!     tl_sqrts = 0.5d0*tl_s/sqrts
      ad_s=ad_s+0.5d0*ad_sqrts/sqrts
      ad_sqrts=0.0d0

      if(iflag.eq.1)ad_fp_theta=ad_s
      if(iflag.eq.2)ad_fp_theta=ad_p

      return

      END FUNCTION ad_fp_theta

# else

      REAL FUNCTION ad_dTemp(iflag,ain,Sal,Temp,Pres)
! Calculates from the salinity (Sal,psu), the in-situ Temperature
! (Temp, degC) and the in-situ pressure (Pres, dbar) the adiabatic 
! temperature gradient (dTemp, K Dbar^-1).
!
! Check values: dTemp  =     3.255976E-4 K dbar^-1
!          given Sal    =    40.0         psu
!                Temp   =    40.0         degC
!                Pres   = 10000.000       dbar

      USE mod_kinds

      real(r8), intent(in) :: Sal, Temp, Pres, ain
      real(r8)             :: ad_Sal, ad_Temp, ad_Pres
      real(r8)             :: s0,a0,a1,a2,a3,b0,b1,c0,c1,c2,c3
      real(r8)             :: d0,d1,e0,e1,e2,ds,ad_ds
      real(r8)             :: fac1, fac2, ad_fac1, ad_fac2

      integer :: iflag

      data s0 /35.0D0/
      data a0,a1,a2,a3 /3.5803D-5, 8.5258D-6, -6.8360D-8, 6.6228D-10/
      data b0,b1       /1.8932D-6, -4.2393D-8/
      data c0,c1,c2,c3 /1.8741D-8, -6.7795D-10, 8.7330D-12, -5.4481D-14/
      data d0,d1       /-1.1351D-10, 2.7759D-12/
      data e0,e1,e2    /-4.6206D-13,  1.8676D-14, -2.1687D-16/

!
!  Clear adjoint variables
!
      ad_fac1=0.0_r8
      ad_fac2=0.0_r8
      ad_Temp=0.0_r8
      ad_Sal=0.0_r8
      ad_Pres=0.0_r8
      ad_ds=0.0_r8

      ds = Sal-s0

      fac1=((e2*Temp+e1)*Temp+e0)*Pres
      fac2=(d1*Temp + d0)*ds + ( (c3*Temp + c2)*Temp + c1 )*Temp + c0

!     tl_dTemp = (tl_fac1+tl_fac2)*Pres+(fac1+fac2)*tl_Pres+            &
!    &           b1*tl_Temp*ds+(b1*Temp+b0)*tl_ds+                      &
!    &             ((a3*Temp+a2)*Temp+a1)*tl_Temp+                      &
!    &             ((a3*Temp+a2)*tl_Temp+a3*Temp*tl_Temp)*Temp

      ad_fac1=ad_fac1+Pres*ain
      ad_fac2=ad_fac2+Pres*ain
      ad_Pres=ad_Pres+(fac1+fac2)*ain
      ad_Temp=ad_Temp+b1*ds*ain
      ad_ds=ad_ds+(b1*Temp+b0)*ain
      ad_Temp=ad_Temp+((a3*Temp+a2)*Temp+a1)*ain
      ad_Temp=ad_Temp+Temp*(a3*Temp+a2)*ain
      ad_Temp=ad_Temp+Temp*a3*Temp*ain

!     tl_fac2=d1*tl_Temp*ds+(d1*Temp+d0)*tl_ds+                         &
!    &                  ((c3*Temp+c2)*tl_Temp+c3*Temp*tl_Temp)*Temp+    &
!    &                  ((c3*Temp+c2)*Temp+c1)*tl_Temp
      ad_Temp=ad_Temp+d1*ds*ad_fac2
      ad_ds=ad_ds+(d1*Temp+d0)*ad_fac2
      ad_Temp=ad_Temp+Temp*(c3*Temp+c2)*ad_fac2
      ad_Temp=ad_Temp+Temp*Temp*c3*ad_fac2
      ad_Temp=ad_Temp+((c3*Temp+c2)*Temp+c1)*ad_fac2
      ad_fac2=0.0_r8

!     tl_fac1=((e2*Temp+e1)*Temp+e0)*tl_Pres+((e2*Temp+e1)*tl_Temp+     &
!    &                                         e2*Temp*tl_Temp)*Pres
      ad_Pres=ad_Pres+((e2*Temp+e1)*Temp+e0)*ad_fac1
      ad_Temp=ad_Temp+Pres*(e2*Temp+e1)*ad_fac1
      ad_Temp=ad_Temp+Pres*Temp*e2*ad_fac1
      ad_fac1=0.0_r8

!     tl_ds = tl_Sal
      ad_Sal=ad_Sal+ad_ds
      ad_ds=0.0_r8

      if(iflag.eq.1)ad_dTemp=ad_Sal
      if(iflag.eq.2)ad_dTemp=ad_Temp
      if(iflag.eq.3)ad_dTemp=ad_Pres

      RETURN
      END FUNCTION ad_dTemp

      REAL FUNCTION ad_thetaa(iflag,ain,Sal,Temp,Pres,RPres)
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
      USE iceshelf_vbc_mod, ONLY : dTemp

      real(r8), intent(in) ::  Sal,Temp,Pres,RPres,ain
      real(r8)             ::  ad_Sal,ad_Temp,ad_Pres,ad_RPres
      real(r8)             ::  p,t,dp,dt,q,ct2,ct3,cq2a,cq2b,cq3a,cq3b
      real(r8)             ::  ad_p,ad_t,ad_dp,ad_dt,ad_q

      integer :: iflag

      data ct2 ,ct3  /0.29289322 ,  1.707106781/
      data cq2a,cq2b /0.58578644 ,  0.121320344/
      data cq3a,cq3b /3.414213562, -4.121320344/

!
!  Clear adjoint variables
!
      ad_t=0.0_r8
      ad_p=0.0_r8
      ad_dp=0.0_r8
      ad_dt=0.0_r8
      ad_q=0.0_r8
      ad_Sal=0.0_r8
      ad_Temp=0.0_r8
      ad_Pres=0.0_r8
      ad_RPres=0.0_r8
!
! Recompute p, etc

      p  = Pres
      t  = Temp
      dp = RPres-Pres
      dt = dp*dTemp(Sal,t,p)
      t  = t +0.5_r8*dt
      q = dt
      p  = p +0.5_r8*dp
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct2*(dt-q)
      q  = cq2a*dt + cq2b*q
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct3*(dt-q)
      q  = cq3a*dt + cq3b*q
      p  = RPres
      dt = dp*dTemp(Sal,t,p)

!     tl_thetaa = tl_t + (tl_dt-tl_q-tl_q)/6.0
      ad_t=ad_t+ain
      ad_dt=ad_dt+ain/6.0_r8
      ad_q=ad_q-2.0_r8*ain/6.0_r8

!     tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)
      ad_dp=ad_dp+dTemp(Sal,t,p)*ad_dt
      ad_Sal=ad_Sal+dp*ad_dTemp(1,ad_dt,Sal,t,p)
      ad_t=ad_t+dp*ad_dTemp(2,ad_dt,Sal,t,p)
      ad_p=ad_p+dp*ad_dTemp(3,ad_dt,Sal,t,p)
      ad_dt=0.0_r8

!     tl_p  = tl_RPres
      ad_RPres=ad_RPres+ad_p
      ad_p=0.0_r8

!     tl_q  = cq3a*tl_dt + cq3b*tl_q
      ad_dt=ad_dt+cq3a*ad_q
      ad_q=cq3b*ad_q

!     tl_t  = tl_t + ct3*(tl_dt-tl_q)
      ad_dt=ad_dt+ct3*ad_t
      ad_q=ad_q-ct3*ad_t
!
!  Recompute p, etc
!
      p  = Pres
      t  = Temp
      dp = RPres-Pres
      dt = dp*dTemp(Sal,t,p)
      t  = t +0.5_r8*dt
      q = dt
      p  = p +0.5_r8*dp
      dt = dp*dTemp(Sal,t,p)
      t  = t + ct2*(dt-q)
      q  = cq2a*dt + cq2b*q

!     tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)
      ad_dp=ad_dp+dTemp(Sal,t,p)*ad_dt
      ad_Sal=ad_Sal+dp*ad_dTemp(1,ad_dt,Sal,t,p)
      ad_t=ad_t+dp*ad_dTemp(2,ad_dt,Sal,t,p)
      ad_p=ad_p+dp*ad_dTemp(3,ad_dt,Sal,t,p)
      ad_dt=0.0_r8

!     tl_q  = cq2a*tl_dt + cq2b*tl_q
      ad_dt=ad_dt+cq2a*ad_q
      ad_q=cq2b*ad_q

!     tl_t  = tl_t + ct2*(tl_dt-tl_q)
      ad_dt=ad_dt+ct2*ad_t
      ad_q=ad_q-ct2*ad_t

!
!  Recompute p, etc
!
      p  = Pres
      t  = Temp
      dp = RPres-Pres
      dt = dp*dTemp(Sal,t,p)
      t  = t +0.5_r8*dt
      q = dt
      p  = p +0.5_r8*dp

!     tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)
      ad_dp=ad_dp+dTemp(Sal,t,p)*ad_dt
      ad_Sal=ad_Sal+dp*ad_dTemp(1,ad_dt,Sal,t,p)
      ad_t=ad_t+dp*ad_dTemp(2,ad_dt,Sal,t,p)
      ad_p=ad_p+dp*ad_dTemp(3,ad_dt,Sal,t,p)
      ad_dt=0.0_r8

!     tl_p  = tl_p +0.5_r8*tl_dp
      ad_dp=ad_dp+0.5_r8*ad_p

!     tl_q = tl_dt
      ad_dt=ad_dt+ad_q
      ad_q=0.0_r8

!     tl_t  = tl_t +0.5_r8*tl_dt
      ad_dt=ad_dt+0.5_r8*ad_t

!
!  Recompute p, etc
!
      p  = Pres
      t  = Temp
      dp = RPres-Pres

!     tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)
      ad_dp=ad_dp+dTemp(Sal,t,p)*ad_dt
      ad_Sal=ad_Sal+dp*ad_dTemp(1,ad_dt,Sal,t,p)
      ad_t=ad_t+dp*ad_dTemp(2,ad_dt,Sal,t,p)
      ad_p=ad_p+dp*ad_dTemp(3,ad_dt,Sal,t,p)
      ad_dt=0.0_r8

!     tl_dp = tl_RPres-tl_Pres
      ad_Pres=ad_Pres-ad_dp
      ad_RPres=ad_RPres+ad_dp
      ad_dp=0.0_r8

!     tl_t  = tl_Temp
      ad_Temp=ad_Temp+ad_t
      ad_t=0.0_r8

!     tl_p  = tl_Pres
      ad_Pres=ad_Pres+ad_p
      ad_p=0.0_r8

      if(iflag.eq.1)ad_thetaa=ad_Sal
      if(iflag.eq.2)ad_thetaa=ad_Temp
      if(iflag.eq.3)ad_thetaa=ad_Pres
      if(iflag.eq.4)ad_thetaa=ad_RPres


      END FUNCTION ad_thetaa

! *********************************************************************
      SUBROUTINE ad_potit(Sal,ad_Sal,theta,ad_theta,Pres,ad_Pres,RPres, &
     &                                ad_RPres,Temp,ad_Temp,i,j)
! *********************************************************************
! Calculates from the salinity (sal, psu), potential temperature 
! (theta, degC) and reference pressure (pres, dbar) the in-situ 
! temperaure (Temp_insitu, degC) related to the in-situ pressure 
! (rfpres, dbar) with the help of an iterative method.

      USE mod_kinds
      USE iceshelf_vbc_mod, ONLY : thetaa

      integer, intent(in)   :: i, j
      real(r8), intent(in)  :: Sal, Pres,theta, RPres
      real(r8), intent(inout)  :: ad_Sal, ad_Pres,ad_theta, ad_RPres
      real(r8), intent(out) :: Temp
      real(r8), intent(inout) :: ad_Temp

      integer               :: ind, ind1
      real(r8)              :: tpmd, theta1, thetad, epsi
      real(r8)              :: ad_thetad, ad_epsi

      data tpmd / 0.001_r8 /

!
      ad_epsi=0.0_r8
      ad_thetad=0.0_r8

      do ind=100,1,-1

        epsi = 0.0_r8
        do ind1=1,ind
           Temp   = theta+epsi
           thetad  = thetaa(Sal,Temp,Pres,RPres)-theta
           epsi = epsi-thetad
        enddo

!   Perform 100 iterations in AD routine always.
!     IF(abs(thetad).lt.tpmd) return
!      tl_epsi = tl_epsi-tl_thetad
       ad_thetad=ad_thetad-ad_epsi

!     tl_thetad  = tl_thetaa(Sal,tl_Sal,Temp,tl_Temp,Pres,tl_Pres,      &
!    &                                       RPres,tl_RPres)-tl_theta
      ad_theta=ad_theta-ad_thetad
      ad_Sal=ad_Sal+ad_thetaa(1,ad_thetad,Sal,Temp,Pres,RPres)
      ad_Temp=ad_Temp+ad_thetaa(2,ad_thetad,Sal,Temp,Pres,RPres)
      ad_Pres=ad_Pres+ad_thetaa(3,ad_thetad,Sal,Temp,Pres,RPres)
      ad_RPres=ad_RPres+ad_thetaa(4,ad_thetad,Sal,Temp,Pres,RPres)
      ad_thetad=0.0_r8

!     tl_Temp   = tl_theta+tl_epsi
      ad_theta=ad_theta+ad_Temp
      ad_epsi=ad_epsi+ad_Temp
      ad_Temp=0.0_r8

      end do

!     tl_epsi = 0.0_r8
      ad_epsi = 0.0_r8

      RETURN
      END SUBROUTINE ad_potit
! *********************************************************************

# endif

#endif
      END MODULE ad_iceshelf_vbc_mod

