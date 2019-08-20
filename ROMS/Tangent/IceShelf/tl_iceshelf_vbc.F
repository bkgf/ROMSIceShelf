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
      MODULE tl_iceshelf_vbc_mod
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
      PUBLIC  :: tl_iceshelf_vbc

      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_iceshelf_vbc (ng, tile)
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
      CALL tl_iceshelf_vbc_tile (ng, tile,                              &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng), nnew(ng),                            &
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % tl_Hz,                              &
#endif
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % zice,                               &
     &                   GRID(ng) % f,                                  &
     &                   GRID(ng) % tl_z_r,                             &
     &                   GRID(ng) % tl_z_w,                             &
     &                   GRID(ng) % tl_zice,                            &
#ifdef ICESHELF_MORPH
     &                   ICESHELFVAR(ng) % iceshelf_draft,              &
     &                   ICESHELFVAR(ng) % tl_iceshelf_draft,           &
#endif
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % tl_u,                              &
     &                   OCEAN(ng) % tl_v,                              &
     &                   OCEAN(ng) % tl_t,                              &
#ifdef ICESHELF_3EQN_VBC
     &                   ICESHELFVAR(ng) % gammaT,                      &
     &                   ICESHELFVAR(ng) % gammaS,                      &
     &                   ICESHELFVAR(ng) % Tb,                          &
     &                   ICESHELFVAR(ng) % Sb,                          &
     &                   ICESHELFVAR(ng) % tl_gammaT,                   &
     &                   ICESHELFVAR(ng) % tl_gammaS,                   &
     &                   ICESHELFVAR(ng) % tl_Tb,                       &
     &                   ICESHELFVAR(ng) % tl_Sb,                       &
#endif
     &                   ICESHELFVAR(ng) % m,                           &
     &                   OCEAN(ng) % rho,                               &
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                   ICESHELFVAR(ng) % tl_m,                        &
     &                   OCEAN(ng) % tl_rho,                            &
     &                   FORCES(ng) % tl_sustr,                         &
     &                   FORCES(ng) % tl_svstr,                         &
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % tl_srflx,                         &
#   endif
     &                   FORCES(ng) % stflx                             &
     &                   FORCES(ng) % tl_stflx                          &
     &                   )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
# endif
      RETURN
      END SUBROUTINE tl_iceshelf_vbc
!
!***********************************************************************
      SUBROUTINE tl_iceshelf_vbc_tile (ng, tile,                        &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,nnew,                               &
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
     &                         Hz,                                      &
     &                         tl_Hz,                                   &
#endif
#ifdef MASKING                 
     &                         rmask,                                   &
#endif
     &                         z_r, z_w,                                &
     &                         zice,f,                                  &
     &                         tl_z_r, tl_z_w,                          &
     &                         tl_zice,                                 &
# ifdef ICESHELF_MORPH
     &                         iceshelf_draft,                          &
     &                         tl_iceshelf_draft,                       &
# endif
     &                         u, v, t,                                 &
     &                         tl_u, tl_v, tl_t,                        &
#ifdef ICESHELF_3EQN_VBC
     &                         gammaT, gammaS,                          &
     &                         Tb, Sb,                                  &
     &                         tl_gammaT, tl_gammaS,                    &
     &                         tl_Tb, tl_Sb,                            &
#endif 
     &                         m, rho,                                  &
     &                         sustr, svstr,                            &
     &                         tl_m, tl_rho,                            &
     &                         tl_sustr, tl_svstr,                      &
#   ifdef SHORTWAVE
     &                         srflx,                                   &
     &                         tl_srflx,                                &
#   endif
     &                         stflx, tl_stflx)                                   
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
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(out) :: tl_svstr(LBi:,LBj:)
      real(r8), intent(in) :: tl_rho(LBi:,LBj:,:)
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
#endif
#ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_zice(LBi:,LBj:)
      real(r8), intent(in) :: f(LBi:,LBj:)
#ifdef ICESHELF_MORPH
      real(r8), intent(inout) :: iceshelf_draft(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_iceshelf_draft(LBi:,LBj:,:)
#endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
#ifdef ICESHELF_3EQN_VBC
      real(r8), intent(out) :: gammaT(LBi:,LBj:)
      real(r8), intent(out) :: gammaS(LBi:,LBj:)
      real(r8), intent(out) :: Tb(LBi:,LBj:)
      real(r8), intent(out) :: Sb(LBi:,LBj:)
      real(r8), intent(out) :: tl_gammaT(LBi:,LBj:)
      real(r8), intent(out) :: tl_gammaS(LBi:,LBj:)
      real(r8), intent(out) :: tl_Tb(LBi:,LBj:)
      real(r8), intent(out) :: tl_Sb(LBi:,LBj:)
#endif
      real(r8), intent(out) :: m(LBi:,LBj:)
      real(r8), intent(out) :: tl_m(LBi:,LBj:)
# ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: tl_srflx(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(out) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
#if defined ICESHELF_3EQN_VBC || defined ANA_SEAICE
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
#endif
#ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: tl_zice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: f(LBi:UBi,LBj:UBj)
#ifdef ICESHELF_MORPH
      real(r8), intent(in) :: iceshelf_draft(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tl_iceshelf_draft(LBi:UBi,LBj:UBj,2)
#endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#ifdef ICESHELF_3EQN_VBC
      real(r8), intent(out) :: gammaT(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: gammaS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Tb(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Sb(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_gammaT(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_gammaS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_Tb(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_Sb(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: m(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: tl_m(LBi:UBi,LBj:UBj)
# ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_srflx(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
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
      real(r8) :: tl_Sm,tl_Tm,tl_rhoi_on_rho0,tl_ustar,tl_TFb,tl_turb
      real(r8) :: tl_RPres
#ifdef ICESHELF_VBC_HFPROPTOM
      real(r8) :: Ex1,Ex2,Ex3,Ex4,Ex5,Ex6,Ep1,Ep2,Ep3,Ep31,Ep4,Ep5
      real(r8) :: tl_Ex1,tl_Ex2,tl_Ex3,tl_Ex4,tl_Ex5,tl_Ex6,tl_Ep1,      &
     &                         tl_Ep2,tl_Ep3,tl_Ep31,tl_Ep4,tl_Ep5
      real(r8) :: Sr1,Sr2,Sf1,Sf2,Tf1,Tf2,rho_r
      real(r8) :: tl_Sr1,tl_Sr2,tl_Sf1,tl_Sf2,tl_Tf1,tl_Tf2,tl_rho_r
#else
      real(r8) :: mflag,TFi,cff3 
      real(r8) :: tl_TFi,tl_cff3 
#endif
      real(r8) :: cff1, cff2
      real(r8) :: tl_cff1, tl_cff2
#ifdef ICESHELF_TEOS10
      real(r8) :: ct, tl_ct
#endif
#  if (!defined BBL_MODEL) && defined UV_LOGDRAG
      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: wrk
#  endif
#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------

      tl_RPres=0.0_r8

#   if defined UV_LOGDRAG
!
!  Set logarithmic ice shelf cavity stress.
!
      DO j=JstrV-1,JEND
        DO i=IstrU-1,IEND
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/Zob(ng))
          tl_cff1=-cff1*cff1*(tl_z_w(i,j,N(ng))-tl_z_r(i,j,N(ng)))/     &
     &                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
          cff2=vonKar*vonKar*cff1*cff1
          tl_cff2=vonKar*vonKar*2.0_r8*cff1*tl_cff1
!>AMM      wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
          cff3=MAX(Cdb_min,cff2)
          tl_cff3=(0.5_r8-SIGN(0.5_r8,Cdb_min-cff2))*tl_cff2
          wrk(i,j)=MIN(Cdb_max,cff3)
          tl_wrk(i,j)=(0.5_r8-SIGN(0.5_r8,cff3-Cdb_max))*tl_cff3
        END DO
      END DO
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
     &                    tl_v(i  ,j+1,N(ng),nrhs)+                     &
     &                    tl_v(i-1,j  ,N(ng),nrhs)+                     &
     &                    tl_v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
     &                 cff1*tl_cff1)/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>         sustr(i,j)=-0.5_r8*(wrk(i-1,j)+wrk(i,j))*                   &
!>  &                 u(i,j,N(ng),nrhs)*cff2
            tl_sustr(i,j)=-0.5_r8*                                      &
     &                    ((tl_wrk(i-1,j)+tl_wrk(i,j))*                 &
     &                     u(i,j,N(ng),nrhs)*cff2+                      &
     &                     (wrk(i-1,j)+wrk(i,j))*                       &
     &                     (tl_u(i,j,N(ng),nrhs)*cff2+                  &
     &                      u(i,j,N(ng),nrhs)*tl_cff2))
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
            tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(cff1*tl_cff1+                                    &
     &              v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          svstr(i,j)=-0.5_r8*(wrk(i,j-1)+wrk(i,j))*                   &
!>   &                 v(i,j,N(ng),nrhs)*cff2
            tl_svstr(i,j)=-0.5_r8*                                      &
     &                    ((tl_wrk(i,j-1)+tl_wrk(i,j))*                 &
     &                     v(i,j,N(ng),nrhs)*cff2+                      &
     &                     (wrk(i,j-1)+wrk(i,j))*                       &
     &                     (tl_v(i,j,N(ng),nrhs)*cff2+                  &
     &                      v(i,j,N(ng),nrhs)*tl_cff2))
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
            tl_cff1=0.25_r8*(tl_v(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_v(i  ,j+1,N(ng),nrhs)+                  &
     &                       tl_v(i-1,j  ,N(ng),nrhs)+                  &
     &                       tl_v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(u(i,j,N(ng),nrhs)*tl_u(i,j,N(ng),nrhs)+          &
     &                 cff1*tl_cff1)/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          sustr(i,j)=-rdrg2(ng)*u(i,j,N(ng),nrhs)*cff2
            tl_sustr(i,j)=-rdrg2(ng)*(tl_u(i,j,N(ng),nrhs)*cff2+        &
     &                                u(i,j,N(ng),nrhs)*tl_cff2)
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
            tl_cff1=0.25_r8*(tl_u(i  ,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j  ,N(ng),nrhs)+                  &
     &                       tl_u(i  ,j-1,N(ng),nrhs)+                  &
     &                       tl_u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            IF (cff2.ne.0.0_r8) THEN
              tl_cff2=(cff1*tl_cff1+                                    &
     &               v(i,j,N(ng),nrhs)*tl_v(i,j,N(ng),nrhs))/cff2
            ELSE
              tl_cff2=0.0_r8
            END IF
!>          svstr(i,j)=-rdrg2(ng)*v(i,j,N(ng),nrhs)*cff2
            tl_svstr(i,j)=-rdrg2(ng)*(tl_v(i,j,N(ng),nrhs)*cff2+        &
     &                                v(i,j,N(ng),nrhs)*tl_cff2)
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
!>          sustr(i,j)=-rdrg(ng)*u(i,j,N(ng),nrhs)
            tl_sustr(i,j)=-rdrg(ng)*tl_u(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!>          svstr(i,j)=-rdrg(ng)*v(i,j,N(ng),nrhs)
            tl_svstr(i,j)=-rdrg(ng)*tl_v(i,j,N(ng),nrhs)
          END IF
        END DO
      END DO
#   else
      DO j=Jstr,JEND
        DO i=IstrU,IEND
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
!>          sustr(i,j)=0.0_r8
            tl_sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV,JEND
        DO i=Istr,IEND
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
!>          svstr(i,j)=0.0_r8
            tl_svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!  Apply boundary conditions.
!
!>    CALL bc_u2d_tile (ng, tile,                                       &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          sustr)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_sustr)
!>    CALL bc_v2d_tile (ng, tile,                                       &
!>   &                          LBi, UBi, LBj, UBj,                     &
!>   &                          svstr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          tl_svstr)
#ifdef DISTRIBUTE
!>    CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
!>   &                    LBi, UBi, LBj, UBj,                           &
!>   &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
!>   &                    sustr,svstr)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    tl_sustr,tl_svstr)
#endif
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
          tl_Sm=(0.5_r8-SIGN(0.5_r8,-t(i,j,N(ng),nrhs,isalt)))*         &
     &                                tl_t(i,j,N(ng),nrhs,isalt)
!          write(6,*)
# else
          Sm=0.0_r8
          tl_Sm=0.0_r8
# endif
          IF (zice(i,j).ne.0.0_r8) THEN
#if defined ICESHELF_2EQN_VBC
          TFb = a*Sm+b+c*zice(i,j)
          tl_TFb = a*tl_Sm+c*tl_zice(i,j)
# ifdef ICESHELF_TEOS10
          ct = ct_from_theta(Sm,t(i,j,N(ng),nrhs,itemp))
          tl_ct = tl_ct_from_theta(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),    &
     &                             tl_t(i,j,N(ng),nrhs,itemp))
          Tm = t_from_ct(Sm,ct,-zice(i,j))
          tl_Tm = tl_t_from_ct(Sm,tl_Sm,ct,tl_ct,-zice(i,j),            &
     &                                           -tl_zice(i,j))
# else
         CALL potit(Sm,t(i,j,N(ng),nrhs,itemp),-zice(i,j),0.0_r8,Tm,i,j)
         CALL tl_potit(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),               &
     &                          tl_t(i,j,N(ng),nrhs,itemp),            &
     &           -zice(i,j),-tl_zice(i,j),0.0_r8,tl_RPres,Tm,tl_Tm,i,j)
# endif
            stflx(i,j,itemp)=gamma*(TFb-Tm)
            tl_stflx(i,j,itemp)=gamma*(tl_TFb-tl_Tm)
            stflx(i,j,isalt)=Cp*stflx(i,j,itemp)*refSalt/L
            tl_stflx(i,j,isalt)=Cp*tl_stflx(i,j,itemp)*refSalt/L
            m(i,j) = stflx(i,j,isalt)/Sm
            tl_m(i,j) = tl_stflx(i,j,isalt)/Sm-tl_Sm*m(i,j)/Sm
#elif defined ICESHELF_3EQN_VBC
! Get the insitu salinity and temperature
# ifdef ICESHELF_TEOS10
          ct = ct_from_theta(Sm,t(i,j,N(ng),nrhs,itemp))
          tl_ct = tl_ct_from_theta(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),    &
     &                                      tl_t(i,j,N(ng),nrhs,itemp))
          Tm = t_from_ct(Sm,ct,-zice(i,j))
          tl_Tm = tl_t_from_ct(Sm,tl_Sm,ct,tl_ct,-zice(i,j),            &
     &                                           -tl_zice(i,j))
# else
         CALL potit(Sm,t(i,j,N(ng),nrhs,itemp),-zice(i,j),0.0_r8,Tm,i,j)
         CALL tl_potit(Sm,tl_Sm,t(i,j,N(ng),nrhs,itemp),                &
     &                 tl_t(i,j,N(ng),nrhs,itemp),-zice(i,j),           &
     &                 -tl_zice(i,j),0.0_r8,tl_RPres,Tm,tl_Tm,i,j)
# endif
!
!    Calculate exchange coefficients. Note that this assumes log profile 
!            ustar = SQRT(SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+   &
!     &                        (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2))
            cff1=(0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2
            tl_cff1=(sustr(i,j)+sustr(i+1,j))*                          &
     &                       (tl_sustr(i,j)+tl_sustr(i+1,j))

            cff2=(0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2
            tl_cff2=(svstr(i,j)+svstr(i,j+1))*                          &
     &                        (tl_svstr(i,j)+tl_svstr(i,j+1))

            cff3=SQRT(cff1+cff2)
            tl_cff3=0.5_r8*(tl_cff1+tl_cff2)/cff3

            ustar = SQRT(cff3)
            tl_ustar = 0.5_r8*tl_cff3/ustar

#    ifdef MASKING
            ustar = ustar*rmask(i,j)
            tl_ustar = tl_ustar*rmask(i,j)
#    endif
           if(ustar.lt.small)then
              ustar = small
              tl_ustar = 0.0_r8
           end
#    ifdef GAMMAS_KY1972
!           turb = 2.12_r8*LOG((ustar*Hz(i,j,N(ng)))/visc)
            cff1=ustar*Hz(i,j,N(ng))/visc
            tl_cff1=(tl_ustar*Hz(i,j,N(ng))+ustar*tl_Hz(i,j,N(ng)))/visc
            turb = 2.12_r8*LOG(cff1)
            tl_turb = 2.12_r8*tl_cff1/cff1

            gammaT(i,j) = ustar/(turb + Pradj - 9.0_r8)
            tl_gammaT(i,j) = tl_ustar/(turb + Pradj - 9.0_r8)-          &
     &                    tl_turb*gammaT(i,j)/(turb + Pradj - 9.0_r8)

            gammaS(i,j) = ustar/(turb + Scadj - 9.0_r8)
            tl_gammaS(i,j) = tl_ustar/(turb + Scadj - 9.0_r8)-          &
     &                    tl_turb*gammaS(i,j)/(turb + Scadj - 9.0_r8)

#    elif defined GAMMAS_Mc1987
! Uses simplified version of McPhee 1987
            IF (ustar.gt.small.and.ABS(f(i,j)).gt.1.0E-8) THEN
!         turb = 2.5_r8*LOG(5300.0_r8*ustar*ustar/ABS(f(i,j))) + 7.12_r8  
              cff1=5300.0_r8*ustar*ustar/ABS(f(i,j))
              tl_cff1=2.0_r8*5300.0_r8*tl_ustar*ustar/ABS(f(i,j))
              turb = 2.5_r8*LOG(cff1) + 7.12_r8
              tl_turb = 2.5_r8*tl_cff1/cff1
            ELSE
              turb = 0.0_r8
              tl_turb = 0.0_r8
            END IF

            gammaT(i,j) = ustar/(turb + Pradj - 6.0_r8)
            tl_gammaT(i,j) = tl_ustar/(turb + Pradj - 6.0_r8)-           &
     &                       tl_turb*gammaT(i,j)/(turb + Pradj - 6.0_r8)

            gammaS(i,j) = ustar/(turb + Scadj - 6.0_r8)
            tl_gammaS(i,j) = tl_ustar/(turb + Scadj - 6.0_r8)-           &
     &                       tl_turb*gammaS(i,j)/(turb + Scadj - 6.0_r8)
#    else
            gammaT(i,j) = 1.00e-4_r8
            tl_gammaT(i,j) = 0.0_r8
            gammaS(i,j) = 5.05e-7_r8
            tl_gammaS(i,j) = 0.0_r8
#    endif

# ifdef ICESHELF_VBC_HFPROPTOM

! Calculates density RHOW in the boundary layer, interface
! pressure PG[dbar], and solving a quadratic equation for the 
! interface salinity (SB) to determine the melTemp_insitug/freezingrate
! (seta).
          rho_r= rho_i/(1000.0_r8+rho(i,j,N(ng)))
          tl_rho_r= -tl_rho(i,j,N(ng))*rho_r/(1000.0_r8+rho(i,j,N(ng)))

          Ep1  = cp_w*gammaT(i,j)
          tl_Ep1  = cp_w*tl_gammaT(i,j)

          Ep2  = cp_i*gammaS(i,j)
          tl_Ep2  = cp_i*tl_gammaS(i,j)

          Ep3  = L*gammaS(i,j)
          tl_Ep3  = L*tl_gammaS(i,j)

#ifdef ICESHELF2D_TOY
          Ep31 = 0.0_r8
          tl_Ep31 = 0.0_r8
#else
          Ep31 = -rho_r*cp_i*dt_i/zice(i,j)
          tl_Ep31 = -tl_rho_r*cp_i*dt_i/zice(i,j)+                       &
                              tl_zice(i,j)*Ep31/zice(i,j)
#endif
          Ep4  = b+c*z_w(i,j,N(ng))
          tl_Ep4  = c*tl_z_w(i,j,N(ng))

          Ep5  = gammaS(i,j)/rho_r
          tl_Ep5  = tl_gammaS(i,j)/rho_r-tl_rho_r*Ep5/rho_r
!
! Temperature gradient within the ice ONLY changes with melTemp_insitug
          TFb = a*Sm+Ep4
          tl_TFb = a*tl_Sm+tl_Ep4
          IF (Tm.lt.TFb) THEN
            Ex1 = a*(Ep1+Ep31)
            tl_Ex1 = a*(tl_Ep1+tl_Ep31)
            Ex2 = Ep1*(Tm-Ep4)+Ep3+Ep31*(Ti-Ep4)
            tl_Ex2 = tl_Ep1*(Tm-Ep4)+Ep1*(tl_Tm-tl_Ep4)+tl_Ep3+           &
     &               tl_Ep31*(Ti-Ep4)-Ep31*tl_Ep4
            Ex3 = Ep3*Sm
            tl_Ex3 = tl_Ep3*Sm+Ep3*tl_Sm
            Ex6 = 0.5_r8
            tl_Ex6 = 0.0_r8
          ELSE
! negative heat flux term in the ice (due to -kappa/D)
            Ex1 = a*(Ep1-Ep2)
            tl_Ex1 = a*(tl_Ep1-tl_Ep2)
            Ex2 = Ep1*(Ep4-Tm)+Ep2*(Ti+a*Sm-Ep4)-Ep3
            tl_Ex2= tl_Ep1*(Ep4-Tm)+Ep1*(tl_Ep4-tl_Tm)+tl_Ep2*(Ti+a*Sm-Ep4)&
     &             +Ep2*(a*tl_Sm-tl_Ep4)-tl_Ep3
            Ex3 = Sm*(Ep2*(Ep4-Ti)+Ep3)
            tl_Ex3 = tl_Sm*(Ep2*(Ep4-Ti)+Ep3)+Sm*(tl_Ep2*(Ep4-Ti)+         &
     &                                              Ep2*tl_Ep4+tl_Ep3)
            Ex6 = -0.5r_r8
            tl_Ex6 = 0.0_r8
          endif
          Ex4 = Ex2/Ex1
          tl_Ex4 = tl_Ex2/Ex1-tl_Ex1*Ex4/Ex1
          Ex5 = Ex3/Ex1
          tl_Ex5 = tl_Ex3/Ex1-tl_Ex1*Ex5/Ex1
          Sr1 = 0.25_r8*Ex4*Ex4-Ex5
          tl_Sr1 = 0.5_r8*tl_Ex4*Ex4-tl_Ex5
          Sr2 = Ex6*Ex4
          tl_Sr2 = tl_Ex6*Ex4+Ex6*tl_Ex4
          Sf1 = Sr2+sqrt(Sr1)
          tl_Sf1 = tl_Sr2+0.5_r8*tl_Sr1/sqrt(Sr1)
          Tf1 = a*Sf1+Ep4
          tl_Tf1 = a*tl_Sf1+tl_Ep4
          Sf2 = Sr2-sqrt(Sr1)
          tl_Sf2 = tl_Sr2-0.5_r8*tl_Sr1/sqrt(Sr1)
          Tf2 = a*Sf2+Ep4
          tl_Tf2 = a*tl_Sf2+tl_Ep4
!
! Salinities < 0 psu are NOT defined, i.e....
         IF (Sf1.gt.0.) THEN
           Tb(i,j) = Tf1
           tl_Tb(i,j) = tl_Tf1
           Sb(i,j) = Sf1
           tl_Sb(i,j) = tl_Sf1
         ELSE
           Tb(i,j) = Tf2
           tl_Tb(i,j) = tl_Tf2
           Sb(i,j) = Sf2
           tl_Sb(i,j) = tl_Sf2
         endif
!
         m(i,j) = -Ep5*(1.0-Sm/Sb(i,j))
         tl_m(i,j) = -tl_Ep5*(1.0-Sm/Sb(i,j))-Ep5*(-tl_Sm/Sb(i,j)+         &
     &                 tl_Sb(i,j)*Sm/(Sb(i,j)*Sb(i,j)))
!         m(i,j)= Seta * dtfast
!
! Calculates fluxes:
!        stflx(i,j,itemp)=(gammaT(i,j)-rho_r*m(i,j))*(Tb(i,j)-Tm)
         tl_stflx(i,j,itemp)=(tl_gammaT(i,j)-tl_rho_r*m(i,j)-              &
     &                           rho_r*tl_m(i,j))*(Tb(i,j)-Tm)+            &
     &                    (gammaT(i,j)-rho_r*m(i,j))*(tl_Tb(i,j)-tl_Tm)
!        stflx(i,j,isalt)=(gammaS(i,j)-rho_r*m(i,j))*(Sb(i,j)-Sm)
         tl_stflx(i,j,isalt)=(tl_gammaS(i,j)-tl_rho_r*m(i,j)-              &
     &                           rho_r*tl_m(i,j))*(Sb(i,j)-Sm)+            &
     &                    (gammaS(i,j)-rho_r*m(i,j))*(tl_Sb(i,j)-tl_Sm)

# else
            TFb = a*Sm + b + c*zice(i,j)
            tl_TFb = a*tl_Sm  + c*tl_zice(i,j)
            mflag= 0.5*(1 + SIGN(1.0_r8,Tm-TFb))
! No TL required for mflag. tl_mflag=0 by definition. See Giering TAMC manual, p32, Table 8.
! The derivative of a SIGN function wrt 2nd argument of zero.

            TFi = (1 - mflag)*a*Si + b + c*zice(i,j)
            tl_TFi = c*tl_zice(i,j)

! Calculate coefficents in quadratic to be solved:
           cff1 = L/cp_w + mflag*(cp_i/cp_w)*(TFi-Ti)
           tl_cff1 = mflag*(cp_i/cp_w)*tl_TFi
           cff2 = gammaS(i,j)*(L/cp_w + mflag*(cp_i/cp_w)*(TFb-Ti))        &
     &          + gammaT(i,j)*(TFi-Tm)
           tl_cff2 = gammaS(i,j)*mflag*(cp_i/cp_w)*tl_TFb                  &
     &          + gammaT(i,j)*(tl_TFi-tl_Tm)+tl_gammaS(i,j)*(L/cp_w +      &
     &                                  mflag*(cp_i/cp_w)*(TFb-Ti))        &
     &          +tl_gammaT(i,j)*(TFi-Tm)
           cff3 = gammaS(i,j)*gammaT(i,j)*(TFb - Tm)
           tl_cff3 = gammaS(i,j)*gammaT(i,j)*(tl_TFb - tl_Tm)+             &
                     tl_gammaS(i,j)*gammaT(i,j)*(TFb - Tm)+                &
                     gammaS(i,j)*tl_gammaT(i,j)*(TFb - Tm)

! Calculate melt rate:
!          m(i,j) = -(cff2 - SQRT(cff2*cff2 - 4.0_r8*cff1*cff3))/          &
!    &                                                (2.0_r8*cff1)
           fac=SQRT(cff2*cff2 - 4.0_r8*cff1*cff3)
           tl_fac=0.5d0_r8*(2.0_r8*tl_cff2*cff2 -                          & 
     &                   4.0_r8*(tl_cff1*cff3+cff1*tl_cff3))/fac
           m(i,j) = -(cff2 - fac)/(2.0_r8*cff1)
           tl_m(i,j) = -(tl_cff2 -tl_fac)/(2.0_r8*cff1)+                   &
     &                                  tl_cff1*m(i,j)/cff1

! Calculate basal temperature and salinity:
           Tb(i,j) = (gammaT(i,j)*Tm+mflag*(cp_i/cp_w)*m(i,j)*Ti           &
     &             - (L/cp_w)*m(i,j))                                      &
     &             /(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))
           tl_Tb(i,j) = (tl_gammaT(i,j)*Tm+gammaT(i,j)*tl_Tm+              &
     &                   mflag*(cp_i/cp_w)*(tl_m(i,j)*Ti+m(i,j)*tl_Ti)-    &
     &                   (L/cp_w)*tl_m(i,j))                               &
     &             /(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))               &
     &           -(tl_gammaT(i,j) + mflag*(cp_i/cp_w)*tl_m(i,j))*Tb(i,j)   &
     &             /(gammaT(i,j) + mflag*(cp_i/cp_w)*m(i,j))
           Sb(i,j) = (Tb(i,j) - b - c*zice(i,j))/a
           tl_Sb(i,j) = (tl_Tb(i,j) - c*tl_zice(i,j))/a

         stflx(i,j,itemp)=(gammaT(i,j)+rhoi_on_rho0*m(i,j))*(Tb(i,j)-Tm)
         tl_stflx(i,j,itemp)=(tl_gammaT(i,j)+rhoi_on_rho0*tl_m(i,j))*      &
     &                           (Tb(i,j)-Tm)+                             &
     &                (gammaT(i,j)+rhoi_on_rho0*m(i,j))*(tl_Tb(i,j)-tl_Tm)

         stflx(i,j,isalt)=(gammaS(i,j)+rhoi_on_rho0*m(i,j))*(Sb(i,j)-Sm)
         tl_stflx(i,j,isalt)=(tl_gammaS(i,j)+rhoi_on_rho0*tl_m(i,j))*      &
     &                           (Sb(i,j)-Sm)+                             &
     &                (gammaS(i,j)+rhoi_on_rho0*m(i,j))*(tl_Sb(i,j)-tl_Sm)

!           write(6,*) gammaT,gammaS,Tb,Sb,Tm,Sm,rhoi_on_rho0
# endif

#else
      DO itrc=1,NT(ng)
              stflx(i,j,itrc)=0.0_r8
              tl_stflx(i,j,itrc)=0.0_r8
      END DO
#endif
          ELSE
            m(i,j) = 0.0_r8
            tl_m(i,j) = 0.0_r8
# if defined ANA_SEAICE
            stflx(i,j,itemp)=Hz(i,j,N(ng))*                             &
                             (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax
            tl_stflx(i,j,itemp)=tl_Hz(i,j,N(ng))*                       &
     &                       (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax-  &
     &                  Hz(i,j,N(ng))*tl_t(i,j,N(ng),nrhs,itemp)/trelax
            stflx(i,j,isalt)=Hz(i,j,N(ng))*                             &
                             (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax
            tl_stflx(i,j,isalt)=tl_Hz(i,j,N(ng))*                       &
     &                       (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax-  &
     &                  Hz(i,j,N(ng))*tl_t(i,j,N(ng),nrhs,isalt)/trelax
# else
          DO itrc=1,NT(ng)
             stflx(i,j,itrc)=0.0_r8
             tl_stflx(i,j,itrc)=0.0_r8
          END DO
# endif
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
          iceshelf_draft(i,j,nnew)=zice(i,j)+(m(i,j)*dt(ng))!           &
!     &                             MIN(m(i,j)*dt(ng),1.0e-8_r8) 
          tl_iceshelf_draft(i,j,nnew)=tl_zice(i,j)+tl_m(i,j)*dt(ng)
          IF(i.eq.40.and.j.eq.40) THEN
          write(6,*) m(i,j),zice(i,j),iceshelf_draft(i,j,nnew) 
          ENDIF
        END DO
      END DO

        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            iceshelf_draft(:,:,nnew))
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            tl_iceshelf_draft(:,:,nnew))
        END IF
# endif

!
!  Apply gradient or periodic boundary conditions for the two fluxes
!
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,itemp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_stflx(:,:,itemp))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  stflx(:,:,isalt))
      CALL bc_r2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  tl_stflx(:,:,isalt))
#    ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    tl_stflx(:,:,itemp))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    stflx(:,:,isalt))
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic(ng), NSperiodic(ng), &
     &                    tl_stflx(:,:,isalt))
#    endif
#   ifdef SHORTWAVE
      DO j=JstrR,JENDR
        DO i=IstrR,IENDR
          IF (zice(i,j).ne.0.0_r8) THEN
            srflx(i,j)=0.0_r8
            tl_srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif

      RETURN
      END SUBROUTINE tl_iceshelf_vbc_tile

# ifdef ICESHELF_TEOS10
      REAL FUNCTION tl_theta_from_ct(s,tl_s,ct,tl_ct)
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

      a5ct = a5*ct
      tl_a5ct = a5*tl_ct

      b3ct = b3*ct
      tl_b3ct = b3*tl_ct

      ct_factor = (a3+a4*s+a5ct)
      tl_ct_factor = a4*tl_s+tl_a5ct

      th0_num = a0+s*(a1+a2*s)+ct*ct_factor
      tl_th0_num = tl_s*(a1+a2*s)+s*a2*tl_s+tl_ct*ct_factor+            &
     &             ct*tl_ct_factor

      rec_th0_den = 1.0d0/(b0+b1*s+ct*(b2+b3ct))
      tl_rec_th0_den =-(b1*tl_s+tl_ct*(b2+b3ct)+ct*tl_b3ct)*            &
     &                 rec_th0_den*rec_th0_den

      th0 = th0_num*rec_th0_den
      tl_th0 = tl_th0_num*rec_th0_den+th0_num*tl_rec_th0_den

      ct0 = ct_from_theta(s,th0)
      tl_ct0 = tl_ct_from_theta(s,tl_s,th0,tl_th0)

      dth_dct = (ct_factor+a5ct-(b2+b3ct+b3ct)*th0)*rec_th0_den
      tl_dth_dct = (tl_ct_factor+tl_a5ct-(b2+b3ct+b3ct)*tl_th0-         &
     &             (tl_b3ct+tl_b3ct)*th0)*rec_th0_den+                  &
     &             dth_dct*tl_rec_th0_den

      theta = th0-(ct0-ct)*dth_dct
      tl_theta = tl_th0-(tl_ct0-tl_ct)*dth_dct-(ct0-ct)*tl_dth_dct

      rec_Cp0 = 2.50494524832013d-4; ! Cp0 = 3992.10322329649;
      dct = ct_from_theta(s,theta)-ct
      tl_dct = tl_ct_from_theta(s,tl_s,theta,tl_theta)-tl_ct

      dct_dth = (theta+273.15d0)*de_dt0_F(s,theta)*rec_Cp0
      tl_dct_dth = tl_theta*de_dt0_F(s,theta)*rec_Cp0+                  &
     &       (theta+273.15d0)*tl_de_dt0_F(s,tl_s,theta,tl_theta)*rec_Cp0

      theta_from_ct = theta-dct/dct_dth
      tl_theta_from_ct = tl_theta-tl_dct/dct_dth+                       &
     &                   tl_dct_dth*dct/(dct_dth*dct_dth)
!       NOTE: theta has a maximum error of 6.0x10^-14

      RETURN
      END FUNCTION tl_theta_from_ct

      REAL FUNCTION tl_ct_from_theta(s,tl_s,theta,tl_theta)
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
!>    ct_from_theta = rec_Cp0*penthalpy_F(s,theta)
      tl_ct_from_theta = rec_Cp0*tl_penthalpy_F(s,tl_s,theta,tl_theta)

      RETURN
      END FUNCTION tl_ct_from_theta

      REAL FUNCTION tl_penthalpy_F(s,tl_s,th,tl_th)
!   potential enthalpy from differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity (psu)
!   th               : potential temperature (deg C, ITS-90)
!   penthalpy_F      : potential enthalpy (J/kg)
!   check value      : penthalpy_F(20,20) = 81649.4876546448
!   DRJ on 10/12/03

      implicit real*8(a-h,o-z)

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
      tl_x2 = 2.5d-2*tl_s

      x = sqrt(x2)
!R8     tl_x=0.5d0*tl_x2/x
      tl_x=0.5d0*tl_x2/x

      y = 2.5d-2*th
      tl_y = 2.5d-2*tl_th

!     penthalpy_F = 61.013624165232955d0 + y*(168776.46138048015d0  +   &
!    &            y*(-2735.2785605119643d0 + y*(2574.2164453821442d0 +  &
!    &              y*(-1536.6644434977545d0 + y*(545.734049793163d0 +  &
!    &        (-50.910917284743334d0 - 18.30489878927802d0*y)*y))))) +  &
!    &             x2*(416.31512917743896d0 + x*(937.9793807560891d0 +  &
!    &              x*(-3140.435779506947d0 + x*(2975.170149976973d0 +  &
!    &            x*(-1760.137081144729d0 + x*414.5655751783703d0))) +  &
!    &              y*(2167.72082596016d0 + y*(-1224.5772800562902d0 +  &
!    &             y*(326.3074029273967d0 + 50.6703824689518d0*y)))) +  &
!    &               y*(-12694.10018182362d0 + y*(4405.71847182968d0 +  &
!    &            y*(-2132.9690185026416d0 + y*(303.91071982808035d0 +  &
!    &                                      69.74975368852d0*y)))))

      penthalpy_F = a0 + y*(a1  + y*(a2 + y*(a3 + y*(a4 + y*(a5 +       &
     &             (a6 + a7*y)*y))))) + x2*(a8 + x*(a9 + x*(a10 +       &
     &              x*(a11 + x*(a12 + x*a13))) + y*(a14 + y*(a15 +      &
     &              y*(a16 + a17*y)))) + y*(a18 + y*(a19 + y*(a20 +     &
     &              y*(a21 + a22*y)))))

      fac1=y*(a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))))

      tl_fac1=tl_y*(a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))))+        &
     &        y*(tl_y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))+            &
     &              y*(tl_y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))+             &
     &                    y*(tl_y*(a4+y*(a5+(a6+a7*y)*y))+              &
     &                          y*(tl_y*(a5+(a6+a7*y)*y)+               &
     &                                y*((a6+a7*y)*tl_y+a7*y*tl_y)      &
     &                            )
     &                      )
     &                )
     &          )


      fac2a=x*(a10+x*(a11+x*(a12+x*a13)))

      tl_fac2a=tl_x*(a10+x*(a11+x*(a12+x*a13)))+                        &
     &            x*(tl_x*(a11+x*(a12+x*a13))+                          &
     &                  x*(tl_x*(a12+x*a13)+tl_x*x*a13)                 &
     &              )

      fac2b=y*(a14+y*(a15+y*(a16+a17*y)))

      tl_fac2b=tl_y*(a14+y*(a15+y*(a16+a17*y)))+                        &
     &            y*(tl_y*(a15+y*(a16+a17*y))+                          &
     &                  y*(tl_y*(a16+a17*y)+y*a17*tl_y)     )

      fac2=x*(a9+fac2a+fac2b)

      tl_fac2=tl_x*(a9+fac2a+fac2b)+x*(tl_fac2a+tl_fac2b)

      fac3=y*(a18+y*(a19+y*(a20+y*(a21+a22*y))))

      tl_fac3=tl_y*(a18+y*(a19+y*(a20+y*(a21+a22*y))))+                 &
     &           y*(tl_y*(a19+y*(a20+y*(a21+a22*y)))+                   &
     &                 y*(tl_y*(a20+y*(a21+a22*y))+                     &
     &                       y*(tl_y*(a21+a22*y)+y*a22*tl_y)            &
     &                   )
     &             )

!     penthalpy_F =a0+fac1+x2*(a8+fac2+fac3)

      tl_penthalpy_F =tl_fac1+tl_x2*(a8+fac2+fac3)+x2*(tl_fac2+tl_fac3)

      RETURN
      END FUNCTION tl_penthalpy_F

      REAL FUNCTION tl_de_dt0_F(s,tl_s,th,tl_th)
!   d(entropy)/dt from twice differentiating the Gibbs potential in
!   Feistel (2003), Prog. Ocean. 58, 43-114
!   s                : salinity                           (psu)
!   th               : potential temperature              (deg C,
!   ITS-90)
!   de_dt0_F         : d(entropy)/dt                      (J/kgK^2)
!   check value      : de_dt0_F(35,20) = 13.63256369213874
!   DRJ on 30/06/05

      implicit real*8(a-h,o-z)

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
      tl_x2 = 2.5d-2*tl_s

      x = sqrt(x2)
      tl_x = 0.5d0*tl_x2/x

      y = 2.5d-2*th
      tl_y = 2.5d-2*tl_th

!     de_dt = 24715.571866078d0 + x2*(-1858.920033948178d0 +            &
!    &           x*(317.440355256842d0 + y*(-405.1392883572282d0 +      &
!    &          202.6815298758072d0*y)) + y*(1562.563716288858d0 +      &
!    &          y*(-1165.8752731900836d0 + 348.7487684426d0*y))) +      &
!    &         y*(-4420.4472249096725d0 + y*(1778.231237203896d0 +      &
!    &            y*(-1160.5182516851419d0 + (569.531539542516d0 -      &
!    &                              128.13429152494615d0*y)*y)))

!     de_dt = a0+x2*(a1+x*(a2+y*(a3+a4*y))+y*(a5+y*(a6+a7*y)))+         &
!    &            y*(a8+y*(a9+y*(a10+(a11+a12*y)*y)))

      fac1=x*(a2+y*(a3+a4*y))

      tl_fac1=tl_x*(a2+y*(a3+a4*y))+x*(tl_y*(a3+a4*y)+y*a4*tl_y)

      fac2=y*(a5+y*(a6+a7*y))

      tl_fac2=tl_y*(a5+y*(a6+a7*y))+y*(tl_y*(a6+a7*y)+y*a7*tl_y)

      fac3=y*(a8+y*(a9+y*(a10+(a11+a12*y)*y)))

      tl_fac3=tl_y*(a8+y*(a9+y*(a10+(a11+a12*y)*y)))+                   &
     &           y*(tl_y*(a9+y*(a10+(a11+a12*y)*y))+                    &
     &                 y*(tl_y*(a10+(a11+a12*y)*y)+                     &
     &                       y*((a11+a12*y)*tl_y+a12*y*tl_y)            &
     &                   )
     &             )

!     de_dt = a0+x2*(a1+fac1+fac2)+fac3

      tl_de_dt = tl_x2*(a1+fac1+fac2)+x2*(tl_fac1+tl_fac2)+tl_fac3

!     de_dt0_F = 6.25d-4*de_dt

      tl_de_dt0_F = 6.25d-4*tl_de_dt

      RETURN
      END FUNCTION tl_de_dt0_F

!
!  COMMENTS:
!     (1) Do we need basic state II?
!     (2) Introduce variables z2=z*z, z3=z*z*z, etc to save computations.
!

      REAL FUNCTION tl_entropy_diff_F(s,tl_s,t,tl_t,p,tl_p,th0,tl_th0,&
     &                                          pr,tl_pr)
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

      a1=-5.90578348518236d0
      a2=24715.571866078d0
      a3=-2210.2236124548363d0
      a4=592.743745734632d0
      a5=-290.12956292128547d0
      a6=113.90630790850321d0
      a7=-21.35571525415769d0

      x2 = 2.5d-2*s 
      tl_x2 = 2.5d-2*tl_s 

      x = dsqrt(x2) 
      tl_x = 0.5d0*tl_x2/x 

      y = 2.5d-2*th0
      tl_y = 2.5d-2*tl_th0

      z = 1.0d-4*pr 
      tl_z = 1.0d-4*tl_pr 

      z1 = z
      tl_z1 = tl_z

!     fTp = -5.90578348518236d0 + y*(24715.571866078d0 +                &
!    &         y*(-2210.2236124548363d0 + y*(592.743745734632d0 +       &
!    &         y*(-290.12956292128547d0 + (113.90630790850321d0 -       &
!    &                               21.35571525415769d0*y)*y))))

      fTp = a1+y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))

      tl_fTp = tl_y*(a2+y*(a3+y*(a4+y*(a5+(a6+a7*y)*y))))+              &
     &            y*(tl_y*(a3+y*(a4+y*(a5+(a6+a7*y)*y)))+               &
     &               y*(tl_y*(a4+y*(a5+(a6+a7*y)*y))+                   &
     &                    y*(tl_y*(a5+(a6+a7*y)*y)+                     &
     &                         y*((a6+a7*y)*tl_y+a7*y*tl_y)             &
     &                      )                                           &
     &                  )                                               &
     &           )

      b1=-1858.920033948178d0
      b2=781.281858144429d0
      b3=-388.6250910633612d0
      b4=87.18719211065d0

!     gTp = y*(-1858.920033948178d0 + y*(781.281858144429d0 +           &
!    &                y*(-388.6250910633612d0 + 87.18719211065d0*y)))

      gTp = y*(b1+y*(b2+y*(b3+b4*y)))

      tl_gTp = tl_y*(b1+y*(b2+y*(b3+b4*y)))+                            &
     &            y*(tl_y*(b2+y*(b3+b4*y))+                             &
     &                  y*(tl_y*(b3+b4*y)+y*b4*tl_y)                    &
     &              )

      c1=317.440355256842d0
      c2=-202.5696441786141d0
      c3=67.5605099586024d0

!     hTp = y*(317.440355256842d0 + y*(-202.5696441786141d0 +           &
!                                        67.5605099586024d0*y))

      hTp = y*(c1+y*(c2+c3*y))

      tl_hTp = tl_y*(c1+y*(c2+c3*y))+y*(tl_y*(c2+c3*y)+y*c3*tl_y)

      if(z.ne.0.d0) then

        yz = y*z
        tl_yz = tl_y*z+y*tl_z

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

!       fTp = fTp+z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))+y*(z*(d6+            &
!    &           z*(d7+z*(d8+(d9+d10*z)*z)))+y*(z*(d11+z*(d12+          &
!    &          z*(d13+z*(d14+d15*z))))+y*(z*(d16+z*(d17+z*(d18+        &
!    &          d19*z)))+y*(y*(d20*y*z+z*(d21+(d22+d23*z)*z))+z*(d24+   &
!    &              z*(d25+(d26+d27*z)*z))))))

        fac1=z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))
        fac2=z*(d6+z*(d7+z*(d8+(d9+d10*z)*z)))
        fac3=z*(d11+z*(d12+z*(d13+z*(d14+d15*z))))
        fac4=z*(d16+z*(d17+z*(d18+d19*z)))
        fac5=y*(d20*y*z+z*(d21+(d22+d23*z)*z))
        fac6=z*(d24+z*(d25+(d26+d27*z)*z))

        fTp = fTp+fac1+y*(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))

        tl_fac1=tl_z*(d1+z*(d2+z*(d3+z*(d4+d5*z))))+                    &
     &                     z*(tl_z*(d2+z*(d3+z*(d4+d5*z)))+             &
     &                           z*(tl_z*(d3+z*(d4+d5*z))+              &
     &                                 z*(tl_z*(d4+d5*z)+z*d5*tl_z)     &
     &                             )                                    &
     &                       )

        tl_fac2=tl_z*(d6+z*(d7+z*(d8+(d9+d10*z)*z)))+                   &
     &             z*(tl_z*(d7+z*(d8+(d9+d10*z)*z))+                    &
     &             z*(tl_z*(d8+(d9+d10*z)*z)+                           &
     &                   z*((d9+d10*z)*tl_z+d10*z*tl_z)) )

        tl_fac3=tl_z*(d11+z*(d12+z*(d13+z*(d14+d15*z))))+               &
     &             z*(tl_z*(d12+z*(d13+z*(d14+d15*z)))+                 &
     &                   z*(tl_z*(d13+z*(d14+d15*z))+                   &
     &                         z*(tl_z*(d14+d15*z)+z*d15*tl_z)          &
     &                     )                                            &
     &               )

        tl_fac4=tl_z*(d16+z*(d17+z*(d18+d19*z)))+                       &
     &             z*(tl_z*(d17+z*(d18+d19*z))+                         &
     &                   z*(tl_z*(d18+d19*z)+z*d18*tl_z)                &
     &               )

        tl_fac5=tl_y*(d20*y*z+z*(d21+(d22+d23*z)*z))+                   &
     &             y*(d20*tl_y*z+d20*y*tl_z+tl_z*(d21+(d22+d23*z)*z)+   &
     &                      z*((d22+d23*z)*tl_z+d23*z*tl_z)             &
     &               )

        tl_fac6=tl_z*(d24+z*(d25+(d26+d27*z)*z))+                       &
     &             z*(tl_z*(d25+(d26+d27*z)*z)+                         &
     &                   z*((d26+d27*z)*tl_z+d27*z*tl_z)                &
     &               )

        tl_fTp = tl_fTp+tl_fac1+                                        &
     &            tl_y*(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))+          &
     &               y*(tl_fac2+tl_y*(fac3+y*(fac4+y*(fac5+fac6)))+     &
     &                    y*(tl_fac3+tl_y*(fac4+y*(fac5+fac6))+         &
     &                y*(tl_fac4+tl_y*(fac5+fac6)+y*(tl_fac5+tl_fac6))) &
     &                   )
       

        e1=-729.116529735046d0
        e2=343.956902961561d0
        e3=-124.687671116248d0
        e4=31.656964386073d0
        e5=-7.04658803315449d0
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

!       gTp = gTp+z*(e1+z*(e2+z*(e3+(e4+e5*z)*z)))+y*(z*(e6+            &
!    &            z*(e7+z*(e8+z*(e9+e10*z))))+y*(y*z*(e11+z*(e12+       &
!    &            e13*z))+z*(e14+z*(e15+z*(e16+e17*z)))))

        fac1=z*(e1+z*(e2+z*(e3+(e4+e5*z)*z)))
        fac2=z*(e6+z*(e7+z*(e8+z*(e9+e10*z))))
        fac3=y*z*(e11+z*(e12+e13*z))
        fac4=z*(e14+z*(e15+z*(e16+e17*z)))

        gTp = gTp+fac1+y*(fac2+y*(fac3+fac4))

        tl_fac1=tl_z*(e1+z*(e2+z*(e3+(e4+e5*z)*z)))+                    &
     &             z*(tl_z*(e2+z*(e3+(e4+e5*z)*z))+                     &
     &                   z*(tl_z*(e3+(e4+e5*z)*z)+                      &
     &                         z*((e4+e5*z)*tl_z+e5*z*tl_z))            &
     &               )

        tl_fac2=tl_z*(e6+z*(e7+z*(e8+z*(e9+e10*z))))+                   &
     &             z*(tl_z*(e7+z*(e8+z*(e9+e10*z)))+                    &
     &                   z*(tl_z*(e8+z*(e9+e10*z))+                     &
     &                         z*(tl_z*(e9+e10*z)+z*e10*tl_z))          &
     &               )

        tl_fac3=(tl_y*z+y*tl_z)*(e11+z*(e12+e13*z))+                    &
     &           y*z*(tl_z*(e12+e13*z)+z*e13*tl_z)

        tl_fac4=tl_z*(e14+z*(e15+z*(e16+e17*z)))+                       &
     &             z*(tl_z*(e15+z*(e16+e17*z))+                         &
     &                   z*(tl_z*(e16+e17*z)+z*e17*tl_z))

        tl_gTp = tl_gTp+tl_fac1+tl_y*(fac2+y*(fac3+fac4))+              &
     &              y*(tl_fac2+tl_y*(fac3+fac4)+y*(tl_fac3+tl_fac4))

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

!       hTp =hTp+z*(f1+z*(f2+f3*z))+y*(y*(f4*z+f5*y*z)+z*(f6+(f7+f8*z)*z))

        fac1=z*(f1+z*(f2+f3*z))
        fac2=y*(f4*z+f5*y*z)
        fac3=z*(f6+(f7+f8*z)*z)

        hTp =hTp+fac1+y*(fac2+fac3)

        tl_fac1=tl_z*(f1+z*(f2+f3*z))+z*(tl_z*(f2+f3*z)+z*f3*tl_z)

        tl_fac2=tl_y*(f4*z+f5*y*z)+y*(f4*tl_z+f5*(y*tl_z+tl_y*z))

        tl_fac3=tl_z*(f6+(f7+f8*z)*z)+z*((f7+f8*z)*tl_z+f8*z*tl_z)

        tl_hTp =tl_hTp+tl_fac1+tl_y*(fac2+fac3)+y*(tl_fac2+tl_fac3)

      end if
       y = 2.5d-2*t
       tl_y = 2.5d-2*tl_t

       z = 1.0d-4*p
       tl_z = 1.0d-4*tl_p

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

        tl_fTp =tl_fTp-tl_y*(g2+y*(g3+y*(g4+y*(g5+(g6+g7*y)*y))))-      &
     &                    y*(tl_y*(g3+y*(g4+y*(g5+(g6+g7*y)*y)))+       &
     &                          y*(tl_y*(g4+y*(g5+(g6+g7*y)*y))+        &
     &                                y*(tl_y*(g5+(g6+g7*y)*y)+         &
     &                                y*((g6+g7*y)*tl_y+g7*y*tl_y)      &
     &                                  )                               &
     &                            )                                     &
     &                      )


        h1=-1858.920033948178d0
        h2=781.281858144429d0
        h3=-388.6250910633612d0
        h4=87.18719211065d0

!       gTp = gTp - y*(-1858.920033948178d0 + y*(781.281858144429d0 +   &
!    &       y*(-388.6250910633612d0 + 87.18719211065d0*y)))

        gTp=gTp- y*(h1+y*(h2+y*(h3+h4*y)))

        tl_gTp=tl_gTp-tl_y*(h1+y*(h2+y*(h3+h4*y)))-                     &
     &                   y*(tl_y*(h2+y*(h3+h4*y))+                      &
     &                         y*(tl_y*(h3+h4*y)+y*h4*tl_y))

        p1=317.440355256842d0
        p2=-202.5696441786141d0
        p3=67.5605099586024d0

!       hTp = hTp - y*(317.440355256842d0 + y*(-202.5696441786141d0 +   &
!    &                                     67.5605099586024d0*y))

        hTp=hTp-y*(p1+y*(p2+p3*y))

        tl_hTp=tl_hTp-tl_y*(p1+y*(p2+p3*y))-y*(tl_y*(p2+p3*y)+y*p3*tl_y)

        kTp = 22.6683558512829d0*(z1-z)

        tl_kTp = 22.6683558512829d0*(tl_z1-tl_z)

      if(z.ne.0.d0) then

        yz = y*z
        tl_yz = tl_y*z+y*tl_z

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

!        fTp=fTp-(z*(q1+z*(q2+z*(q3+z*(q4+q5*z))))+y*(z*(q6+z*(q7+      &
!    &            z*(q8+(q9+q10*z)*z)))+y*(z*(q11+z*(q12+z*(q13+z*(q14+ &
!    &            q15*z))))+y*(z*(q16+z*(q17+z*(q18+q19*z)))+           &
!    &            y*(y*(q20*y*z+z*(q21+(q22+q23*z)*z))+z*(q24+z*(q25+   &
!    &            (q26+q27*z)*z)))))))

         fac1=z*(q1+z*(q2+z*(q3+z*(q4+q5*z))))
         fac2=z*(q6+z*(q7+z*(q8+(q9+q10*z)*z)))
         fac3=z*(q11+z*(q12+z*(q13+z*(q14+q15*z))))
         fac4=z*(q16+z*(q17+z*(q18+q19*z)))
         fac5=y*(q20*y*z+z*(q21+(q22+q23*z)*z))
         fac6=z*(q24+z*(q25+(q26+q27*z)*z))

         fTp=fTp-(fac1+y*(fac2+y*(fac3+y*(fac4+y*(fac5+fac6)))))

         tl_fac1=tl_z*(q1+z*(q2+z*(q3+z*(q4+q5*z))))+                   &
     &              z*(tl_z*(q2+z*(q3+z*(q4+q5*z)))+                    &
     &                    z*(tl_z*(q3+z*(q4+q5*z))+                     &
     &                          z*(tl_z*(q4+q5*z)+tl_z*q5*z)            &
     &                      )                                           &
     &                )

         tl_fac2=tl_z*(q6+z*(q7+z*(q8+(q9+q10*z)*z)))+                  &
     &              z*(tl_z*(q7+z*(q8+(q9+q10*z)*z))+                   &
     &                    z*(tl_z*(q8+(q9+q10*z)*z)+                    &
     &                          z*((q9+q10*z)*tl_z+q10*z*tl_z)          &
     &                      )                                           &
     &                )

         tl_fac3=tl_z*(q11+z*(q12+z*(q13+z*(q14+q15*z))))+              &
     &              z*(tl_z*(q12+z*(q13+z*(q14+q15*z)))+                &
     &                    z*(tl_z*(q13+z*(q14+q15*z))+                  &
     &                          z*(tl_z*(q14+q15*z)+z*q15*tl_z)         &
     &                      )                                           &
     &                )

         tl_fac4=tl_z*(q16+z*(q17+z*(q18+q19*z)))+                      &
     &              z*(tl_z*(q17+z*(q18+q19*z))+                        &
     &                    z*(tl_z*(q18+q19*z)+z*q19*tl_z)   )

         tl_fac5=tl_y*(q20*y*z+z*(q21+(q22+q23*z)*z))+                  &
     &              y*(q20*(tl_y*z+y*tl_z)+                             &
     &                     tl_z*(q21+(q22+q23*z)*z)+                    &
     &                        z*((q22+q23*z)*tl_z+q23*tl_z*z) )

         tl_fac6=tl_z*(q24+z*(q25+(q26+q27*z)*z))+                      &
     &              z*(tl_z*(q25+(q26+q27*z)*z)+                        &
     &                    z*((q26+q27*z)*tl_z+q27*z*tl_z)     )

         tl_fTp=tl_fTp-                                                 &
     &             (tl_fac1+tl_y*(fac2+y*(fac3+y*(fac4+y*(fac5+fac6))))+&
     &                  y*(tl_fac2+tl_y*(fac3+y*(fac4+y*(fac5+fac6)))+  &
     &                  y*(tl_fac3+tl_y*(fac4+y*(fac5+fac6))+           &
     &               y*(tl_fac4+tl_y*(fac5+fac6)+y*(tl_fac5+tl_fac6)))  &
     &                               )                                  &
     &                 )

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

!        gTp=gTp-(z*(r1+z*(r2+z*(r3+(r4+r5*z)*z)))+y*(z*(r6+z*(r7+      &
!    &            z*(r8+z*(r9+r10*z))))+y*(y*z*(r11+z*(r12+r13*z))+     &
!    &            z*(r14+z*(r15+z*(r16+r17*z))))))
      
         fac1=z*(r1+z*(r2+z*(r3+(r4+r5*z)*z)))
         fac2=z*(r6+z*(r7+z*(r8+z*(r9+r10*z))))
         fac3=y*z*(r11+z*(r12+r13*z))
         fac4=z*(r14+z*(r15+z*(r16+r17*z)))

         gTp=gTp-(fac1+y*(fac2+y*(fac3+fac4)))

         tl_fac1=tl_z*(r1+z*(r2+z*(r3+(r4+r5*z)*z)))+                   &
     &              z*(tl_z*(r2+z*(r3+(r4+r5*z)*z))+                    &
     &                    z*(tl_z*(r3+(r4+r5*z)*z)+                     &
     &                          z*((r4+r5*z)*tl_z+r5*z*tl_z)            &
     &                      )                                           &
     &                )

         tl_fac2=tl_z*(r6+z*(r7+z*(r8+z*(r9+r10*z))))+                  &
     &              z*(tl_z*(r7+z*(r8+z*(r9+r10*z)))+                   &
     &                    z*(tl_z*(r8+z*(r9+r10*z))+                    &
     &                          z*(tl_z*(r9+r10*z)+z*r10*tl_z)          &
     &                      )                                           &
     &                )

         tl_fac3=(tl_y*z+y*tl_z)*(r11+z*(r12+r13*z))+                   &
     &            y*z*(tl_z*(r12+r13*z)+tl_z*r13*z)

         tl_fac4=tl_z*(r14+z*(r15+z*(r16+r17*z)))+                      &
     &              z*(tl_z*(r15+z*(r16+r17*z))+                        &
     &                    z*(tl_z*(r16+r17*z)+tl_z*r17*z)   )

         tl_gTp=tl_gTp-(tl_fac1+tl_y*(fac2+y*(fac3+fac4))+              &
     &                            y*(tl_fac2+tl_y*(fac3+fac4)+          &
     &                                          y*(tl_fac3+tl_fac4) ) )

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

!        hTp =hTp-(z*(s1+z*(s2+s3*z))+y*(y*(s4*z+s5*y*z)+               &
!    &                z*(s6+(s7+s8*z)*z)))

         fac1=z*(s1+z*(s2+s3*z))
         fac2=y*(s4*z+s5*y*z)
         fac3=z*(s6+(s7+s8*z)*z)

         hTp =hTp-(fac1+y*(fac2+fac3))

         tl_fac1=tl_z*(s1+z*(s2+s3*z))+
     &              z*(tl_z*(s2+s3*z)+tl_z*s3*z)

         tl_fac2=tl_y*(s4*z+s5*y*z)+y*(s4*tl_z+s5*(tl_y*z+y*tl_z))

         tl_fac3=tl_z*(s6+(s7+s8*z)*z)+z*((s7+s8*z)*tl_z+s8*z*tl_z)

         tl_hTp =tl_hTp-(tl_fac1+tl_y*(fac2+fac3)+y*(tl_fac2+tl_fac3))

      end if

      entropy_diff_F = 2.5d-2*(fTp+x2*(gTp+x*hTp+x2*kTp))

      tl_entropy_diff_F = 2.5d-2*(tl_fTp+tl_x2*(gTp+x*hTp+x2*kTp)+       &
     &              x2*(tl_gTp+tl_x*hTp+x*tl_hTp+tl_x2*kTp+x2*tl_kTp))

      RETURN
      END FUNCTION tl_entropy_diff_F

      REAL FUNCTION tl_de_dt_F(s,tl_s,t,tl_t,p,tl_p)
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

      x2 = 2.5d-2*s 
      tl_x2 = 2.5d-2*tl_s 

      x = dsqrt(x2) 
      tl_x=0.5d0*tl_x2/x

      y = 2.5d-2*t 
      tl_y = 2.5d-2*tl_t 

      z = 1.0d-4*p
      tl_z = 1.0d-4*tl_p

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

!     de_dt = 24715.571866078d0 + x2*(-1858.920033948178d0 +            &
!    &                x*(317.440355256842d0 + y*(-405.1392883572282d0 + &
!    &          202.6815298758072d0*y)) + y*(1562.563716288858d0 +      &
!    &          y*(-1165.8752731900836d0 + 348.7487684426d0*y))) +      &
!    &         y*(-4420.4472249096725d0 + y*(1778.231237203896d0 +      &
!    &            y*(-1160.5182516851419d0 + (569.531539542516d0 -      &
!    &                              128.13429152494615d0*y)*y)))

!     de_dt =a1+x2*(a2+x*(a3+y*(a4+a5*y))+y*(a6+y*(a7+a8*y)))+          &
!    &         y*(a9+y*(a10+y*(a11+(a12+a13*y)*y)))

      fac1=x*(a3+y*(a4+a5*y))
      tl_fac1=tl_x*(a3+y*(a4+a5*y))+x*(tl_y*(a4+a5*y)+y*a5*tl_y)

      fac2=y*(a6+y*(a7+a8*y))
      tl_fac2=tl_y*(a6+y*(a7+a8*y))+y*(tl_y*(a7+a8*y)+y*a8*tl_y)

      fac3=y*(a10+y*(a11+(a12+a13*y)*y))
      tl_fac3=tl_y*(a10+y*(a11+(a12+a13*y)*y))+                         &
     &           y*(tl_y*(a11+(a12+a13*y)*y)+                           &
     &                 y*((a12+a13*y)*tl_y+a13*y*tl_y) )

!     de_dt =a1+x2*(a2+fac1+fac2)+y*(a9+fac3)

      tl_de_dt =tl_x2*(a2+fac1+fac2)+x2*(tl_fac1+tl_fac2)+              &
     &           tl_y*(a9+fac3)+y*tl_fac3

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

!       de_dt = de_dt +                                                 &
!    &           z*(-2910.0729080936d0 + x2*(1721.528607567954d0 +      &
!    &        y*(-4165.4688847996085d0 + 3572.7449038462437d0*y) +      &
!    &            x*(-766.116132004952d0 + (2761.9195908075417d0 -      &
!    &         2814.78225133626d0*y)*y)) + y*(4035.04669887042d0 +      &
!    &            y*(-4775.621344883664d0 + y*(3892.3662123519d0 +      &
!    &      y*(-1905.341809925355d0 + 404.50541014508605d0*y)))) +      &
!    &          z*(1513.116771538718d0 + x2*(-674.819060538734d0 +      &
!    &              108.3834525034224d0*x + (1229.337851789418d0 -      &
!    &        896.713693665072d0*y)*y) + y*(-2996.162344914912d0 +      &
!    &            y*(3621.784567462512d0 + y*(-2410.4130980405d0 +      &
!    &          668.691951421377d0*y))) + z*(-546.959324647056d0 +      &
!    &             x2*(356.629112415276d0 - 51.2796974779828d0*x +      &
!    &         y*(-681.370187043564d0 + 437.84750280190565d0*y)) +      &
!    &         y*(1437.2719839264719d0 + y*(-1826.356460806092d0 +      &
!    &        (1105.446104680304d0 - 245.11816254543362d0*y)*y)) +      &
!    &             z*(111.1208127634436d0 + x2*(-88.4080716616d0 +      &
!    &           66.7696405958478d0*y) + y*(-292.8075111563232d0 +      &
!    &         (316.49805267936244d0 - 129.6381336154442d0*y)*y) +      &
!    &             (-8.68841343834394d0 + 15.84003094423364d0*x2 +      &
!    &                              9.978426372534301d0*y)*z))))

!       de_dt = de_dt +                                                 &
!    &           z*(b1+x2*(b2+y*(b3+b4*y)+x*(b5+(b6+b7*y)*y))+          &
!    &              y*(b8+y*(b9+y*(b10+y*(b11+b12*y))))+                &
!    &          z*(b13+x2*(b14+b15*x+(b16+b17*y)*y)+y*(b18+             &
!    &            y*(b19+y*(b20+b21*y)))+z*(b22+x2*(b23+b24*x+          &
!    &         y*(b25+b26*y))+y*(b27+y*(b28+(b29+b30*y)*y))+            &
!    &             z*(b31+x2*(b32+b33*y)+y*(b34+(b35+b36*y)*y)+         &
!    &             (b37+b38*x2+b39*y)*z))))

        fac1=x2*(b2+y*(b3+b4*y)+x*(b5+(b6+b7*y)*y))
        tl_fac1=tl_x2*(b2+y*(b3+b4*y)+x*(b5+(b6+b7*y)*y))+              &
     &             x2*(tl_y*(b3+b4*y)+y*b4*tl_y+tl_x*(b5+(b6+b7*y)*y)+  &
     &                                    x*((b6+b7*y)*tl_y+b7*y*tl_y))
 
        fac2=y*(b8+y*(b9+y*(b10+y*(b11+b12*y))))
        tl_fac2=tl_y*(b8+y*(b9+y*(b10+y*(b11+b12*y))))+                 &
     &             y*(tl_y*(b9+y*(b10+y*(b11+b12*y)))+                  &
     &                   y*(tl_y*(b10+y*(b11+b12*y))+                   &
     &                         y*(tl_y*(b11+b12*y)+y*b12*tl_y) ) )


        fac3=x2*(b14+b15*x+(b16+b17*y)*y)
        tl_fac3=tl_x2*(b14+b15*x+(b16+b17*y)*y)+                        &
     &             x2*(b15*tl_x+(b16+b17*y)*tl_y+b17*y*tl_y)

        fac4=y*(b18+y*(b19+y*(b20+b21*y)))
        tl_fac4=tl_y*(b18+y*(b19+y*(b20+b21*y)))+                       &
     &             y*(tl_y*(b19+y*(b20+b21*y))+                         &
     &                   y*(tl_y*(b20+b21*y)+y*b21*tl_y) )

        fac5=x2*(b23+b24*x+y*(b25+b26*y))
        tl_fac5=tl_x2*(b23+b24*x+y*(b25+b26*y))+                        &
     &             x2*(b24*tl_x+tl_y*(b25+b26*y)+y*b26*tl_y)

        fac6=y*(b27+y*(b28+(b29+b30*y)*y))
        tl_fac6=tl_y*(b27+y*(b28+(b29+b30*y)*y))+                       &
     &             y*(tl_y*(b28+(b29+b30*y)*y)+                         &
     &                   y*((b29+b30*y)*tl_y+b30*y*tl_y)  )

        fac7=x2*(b32+b33*y)
        tl_fac7=tl_x2*(b32+b33*y)+x2*b33*tl_y

        fac8=y*(b34+(b35+b36*y)*y)
        tl_fac8=tl_y*(b34+(b35+b36*y)*y)+y*((b35+b36*y)*tl_y+b36*y*tl_y)

        fac9=(b37+b38*x2+b39*y)*z
        tl_fac9=(b37+b38*x2+b39*y)*tl_z+(b38*tl_x2+b39*tl_y)*z

!       de_dt = de_dt +                                                 &
!    &       z*(b1+fac1+fac2+z*(b13+fac3+fac4+z*(b22+fac5+fac6+         &
!    &                     z*(b31+fac7+fac8+fac9))))
 
        tl_de_dt = tl_de_dt +                                           &
     &       tl_z*(b1+fac1+fac2+z*(b13+fac3+fac4+z*(b22+fac5+fac6+      &
     &                     z*(b31+fac7+fac8+fac9))))+
     &          z*(tl_fac1+tl_fac2+tl_z*(b13+fac3+fac4+z*(b22+fac5+fac6+&
     &                     z*(b31+fac7+fac8+fac9)))+                    &
     &                                z*(tl_fac3+tl_fac4+               &
     &                  tl_z*(b22+fac5+fac6+z*(b31+fac7+fac8+fac9))+    &
     &                   z*(tl_fac5+tl_fac6+tl_z*(b31+fac7+fac8+fac9)+  &
     &                          z*(tl_fac7+tl_fac8+tl_fac9)  )  )  )

      end if

!     de_dt_F = 6.25d-4*de_dt

      tl_de_dt_F = 6.25d-4*tl_de_dt

      RETURN
      END FUNCTION tl_de_dt_F

      REAL FUNCTION tl_t_from_ct(s, tl_s, ct, tl_ct, p, tl_p)
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
      tl_theta = tl_theta_from_ct(s,tl_s,ct,tl_ct)

!AMM I THINK THIS IS WRONG      t_from_ct = theta_from_t(s,theta,pr0,p)
!     t_from_ct = theta_from_t(s,theta,p,pr0)
!AMM I THINK THIS IS WRONG      tl_t_from_ct = tl_theta_from_t(s,tl_s,theta,tl_theta,pr0,p,tl_p,pr0)
      tl_t_from_ct = tl_theta_from_t(s,tl_s,theta,tl_theta,p,tl_p,pr0)

      RETURN
      END FUNCTION tl_t_from_ct

      REAL FUNCTION tl_theta_from_t(s,tl_s,t,tl_t,p,tl_p,pr)
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

      tl_pr=0.0d0

      a1=8.65483913395442d-6
      a2=1.41636299744881d-6
      a3=7.38286467135737d-9
      a4=-8.38241357039698d-6
      a5=2.83933368585534d-8
      a6=1.77803965218656d-8
      a7=1.71155619208233d-10

!     th0 = t+(p-pr)*( 8.65483913395442d-6   -                          &
!    &             s*  1.41636299744881d-6   -                          &
!    &        (p+pr)*  7.38286467135737d-9   +                          &
!    &             t*(-8.38241357039698d-6   +                          &
!    &             s*  2.83933368585534d-8   +                          &
!    &             t*  1.77803965218656d-8   +                          &
!    &          (p+pr)*  1.71155619208233d-10))

      th0 = t+(p-pr)*(a1-s*a2-(p+pr)*a3+t*(a4+s*a5+t*a6+(p+pr)*a7))

!     tl_th0 = tl_t+tl_p*( 8.65483913395442d-6   -                      &
!    &             s*  1.41636299744881d-6   -                          &
!    &        (p+pr)*  7.38286467135737d-9   +                          &
!    &             t*(-8.38241357039698d-6   +                          &
!    &             s*  2.83933368585534d-8   +                          &
!    &             t*  1.77803965218656d-8   +                          &
!    &          (p+pr)*  1.71155619208233d-10))+                        &
!    &          (p-pr)*(                                                &
!    &             -tl_s*  1.41636299744881d-6                          &
!    &             -tl_p*  7.38286467135737d-9 +                        &
!    &             tl_t*(-8.38241357039698d-6   +                       &
!    &             s*  2.83933368585534d-8   +                          &
!    &             t*  1.77803965218656d-8   +                          &
!    &          (p+pr)*  1.71155619208233d-10)                          &
!    &             t*(                                                  &
!    &             tl_s*  2.83933368585534d-8 +                         &
!    &             tl_t*  1.77803965218656d-8 +                         &
!    &          tl_p*  1.71155619208233d-10))

      tl_th0 = tl_t+tl_p*(a1-s*a2-(p+pr)*a3+t*(a4+s*a5+t*a6+(p+pr)*a7))+&
     &          (p-pr)*(-tl_s*a2-tl_p*a3+tl_t*(a4+s*a5+t*a6+(p+pr)*a7)+ &
     &                               t*(tl_s*a5+tl_t*a6+tl_p*a7))

      de_dt = 13.6d0
      tl_de_dt = 0.0d0
      nloops = 2                                      ! default
!    NOTE: nloops = 1 gives theta with a maximum error of 5.48x10^-06
!          nloops = 2 gives theta with a maximum error of 2.84x10^-14
!AMM      n = 1
!AMM      do while(n.le.nloops)
      do n=1,nloops

       dentropy = entropy_diff_F(s,t,p,th0,pr)
       tl_dentropy = tl_entropy_diff_F(s,tl_s,t,tl_t,p,tl_p,th0,tl_th0, &
     &                                     pr,tl_pr)

       theta = th0-dentropy/de_dt
       tl_theta = tl_th0-tl_dentropy/de_dt+                             &
     &                            tl_de_dt*dentropy/(de_dt*de_dt)

       theta = 0.5d0*(theta+th0)
       tl_theta = 0.5d0*(tl_theta+tl_th0)

       de_dt = de_dt_F(s,theta,pr)
       tl_de_dt = tl_de_dt_F(s,tl_s,theta,tl_theta,pr,tl_pr)

       theta = th0-dentropy/de_dt
       tl_theta = tl_th0-tl_dentropy/de_dt+                             &
     &                           tl_de_dt*dentropy/(de_dt*de_dt)

!AMM       n = n+1 
       th0 = theta
       tl_th0 = tl_theta

      end do
      theta_from_t = th0
      tl_theta_from_t = tl_th0

      RETURN
      END FUNCTION tl_theta_from_t

      REAL FUNCTION tl_fp_theta(s,tl_s,p,tl_p,sat)
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

      a1=2.5180516744541290d-03
      a2=-5.8545863698926184d-02
      a3=2.2979985780124325d-03
      a4=3.0086338218235500d-04
      a5=-7.0023530029351803d-04
      a6=8.4149607219833806d-09
      a7=1.1845857563107403d-11

      sqrts = dsqrt(s)
      tl_sqrts = 0.5d0*tl_s/sqrts

!     tf_num =                    2.5180516744541290d-03    +           &
!    &                        s*(-5.8545863698926184d-02    +           &
!    &                    sqrts*( 2.2979985780124325d-03    -           &
!    &                   sqrts *  3.0086338218235500d-04))  +           &
!    &                        p*(-7.0023530029351803d-04    +           &
!    &                        p*( 8.4149607219833806d-09    +           &
!    &                       s *  1.1845857563107403d-11));           

      tf_num =a1+s*(a2+sqrts*(a3-sqrts*a4))+p*(a5+p*(a6+s*a7))

      tl_tf_num=tl_s*(a2+sqrts*(a3-sqrts*a4))+                          &
     &             s*(tl_sqrts*(a3-sqrts*a4)-tl_sqrts*sqrts*a4)+        &
     &          tl_p*(a5+p*(a6+s*a7))+p*(tl_p*(a6+s*a7)+p*tl_s*a7)

      b1=1.0000000000000000d+00
      b2=-3.8493266309172074d-05
      b3=9.1686537446749641d-10
      b4=1.3632481944285909d-06

!     tf_den =                          1.0000000000000000d+00    +     &
!    &                              p*(-3.8493266309172074d-05    +     &
!    &                             p *  9.1686537446749641d-10)   +     &
!    &                      s*s*sqrts*  1.3632481944285909d-06

      tf_den =b1+p*(b2+p*b3)+s*s*sqrts*b4

      tl_tf_den =tl_p*(b2+p*b3)+p*tl_p*b3+2.0d0*tl_s*s*sqrts*b4+        &
     &                                    s*s*tl_sqrts*b4

      fp_theta = tf_num/tf_den

      tl_fp_theta = tl_tf_num/tf_den-tl_tf_den*fp_theta/tf_den

      if(sat.eq.'air-sat') then
      fp_theta = fp_theta          -2.5180516744541290d-03    +         &
     &                             s *  1.428571428571429d-05

      tl_fp_theta=tl_fp_theta+tl_s*1.428571428571429d-05
      elseif(sat.eq.'air-free') then
      continue
      else
      print *, '***  Error in fp_theta.f90: invalid third argument  ***'
      print *
      stop
      endif
      return

      END FUNCTION tl_fp_theta

# else

      REAL FUNCTION tl_dTemp(Sal,tl_Sal,Temp,tl_Temp,Pres,tl_Pres)
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
      real(r8), intent(in) :: tl_Sal, tl_Temp, tl_Pres
      real(r8)             :: s0,a0,a1,a2,a3,b0,b1,c0,c1,c2,c3
      real(r8)             :: d0,d1,e0,e1,e2,ds,tl_ds
      real(r8)             :: fac1, fac2, tl_fac1, tl_fac2

      data s0 /35.0D0/
      data a0,a1,a2,a3 /3.5803D-5, 8.5258D-6, -6.8360D-8, 6.6228D-10/
      data b0,b1       /1.8932D-6, -4.2393D-8/
      data c0,c1,c2,c3 /1.8741D-8, -6.7795D-10, 8.7330D-12, -5.4481D-14/
      data d0,d1       /-1.1351D-10, 2.7759D-12/
      data e0,e1,e2    /-4.6206D-13,  1.8676D-14, -2.1687D-16/

      ds = Sal-s0
      tl_ds = tl_Sal

!     dTemp = ( ( (e2*Temp + e1)*Temp + e0 )*Pres                       &
!    &      + ( (d1*Temp + d0)*ds                                       &
!    &      + ( (c3*Temp + c2)*Temp + c1 )*Temp + c0 ) )*Pres           &
!    &      + (b1*Temp + b0)*ds +  ( (a3*Temp + a2)*Temp + a1 )*Temp    &
!    &      + a0

      fac1=((e2*Temp+e1)*Temp+e0)*Pres
      tl_fac1=((e2*Temp+e1)*Temp+e0)*tl_Pres+((e2*Temp+e1)*tl_Temp+     &
     &                                         e2*Temp*tl_Temp)*Pres

      fac2=(d1*Temp + d0)*ds + ( (c3*Temp + c2)*Temp + c1 )*Temp + c0
      tl_fac2=d1*tl_Temp*ds+(d1*Temp+d0)*tl_ds+                         &
     &                  ((c3*Temp+c2)*tl_Temp+c3*Temp*tl_Temp)*Temp+    &
     &                  ((c3*Temp+c2)*Temp+c1)*tl_Temp

!     dTemp = (fac1+fac2)*Pres+(b1*Temp + b0)*ds+                       &
!    &             ( (a3*Temp + a2)*Temp + a1 )*Temp + a0

      tl_dTemp = (tl_fac1+tl_fac2)*Pres+(fac1+fac2)*tl_Pres+            &
     &           b1*tl_Temp*ds+(b1*Temp+b0)*tl_ds+                      &
     &             ((a3*Temp+a2)*Temp+a1)*tl_Temp+                      &
     &             ((a3*Temp+a2)*tl_Temp+a3*Temp*tl_Temp)*Temp



      RETURN
      END FUNCTION tl_dTemp

      REAL FUNCTION tl_thetaa(Sal,tl_Sal,Temp,tl_Temp,Pres,tl_Pres,   &
     &                             RPres,tl_RPres)
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
      real(r8), intent(in) ::  tl_Sal,tl_Temp,tl_Pres,tl_RPres
      real(r8)             ::  p,t,dp,dt,q,ct2,ct3,cq2a,cq2b,cq3a,cq3b
      real(r8)             ::  tl_p,tl_t,tl_dp,tl_dt,tl_q

      data ct2 ,ct3  /0.29289322 ,  1.707106781/
      data cq2a,cq2b /0.58578644 ,  0.121320344/
      data cq3a,cq3b /3.414213562, -4.121320344/

      p  = Pres
      tl_p  = tl_Pres

      t  = Temp
      tl_t  = tl_Temp

      dp = RPres-Pres
      tl_dp = tl_RPres-tl_Pres

      dt = dp*dTemp(Sal,t,p)
      tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)

      t  = t +0.5_r8*dt
      tl_t  = tl_t +0.5_r8*tl_dt

      q = dt
      tl_q = tl_dt

      p  = p +0.5_r8*dp
      tl_p  = tl_p +0.5_r8*tl_dp

      dt = dp*dTemp(Sal,t,p)
      tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)

      t  = t + ct2*(dt-q)
      tl_t  = tl_t + ct2*(tl_dt-tl_q)

      q  = cq2a*dt + cq2b*q
      tl_q  = cq2a*tl_dt + cq2b*tl_q

      dt = dp*dTemp(Sal,t,p)
      tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)

      t  = t + ct3*(dt-q)
      tl_t  = tl_t + ct3*(tl_dt-tl_q)

      q  = cq3a*dt + cq3b*q
      tl_q  = cq3a*tl_dt + cq3b*tl_q

      p  = RPres
      tl_p  = tl_RPres

      dt = dp*dTemp(Sal,t,p)
      tl_dt = tl_dp*dTemp(Sal,t,p)+dp*tl_dTemp(Sal,tl_Sal,t,tl_t,p,tl_p)

      thetaa = t + (dt-q-q)/6.0_r8
      tl_thetaa = tl_t + (tl_dt-tl_q-tl_q)/6.0_r8

      END FUNCTION tl_thetaa


! *********************************************************************
      SUBROUTINE tl_potit(Sal,tl_Sal,theta,tl_theta,Pres,tl_Pres,RPres, &
     &                                tl_RPres,Temp,tl_Temp,i,j)
! *********************************************************************
! Calculates from the salinity (sal, psu), potential temperature 
! (theta, degC) and reference pressure (pres, dbar) the in-situ 
! temperaure (Temp_insitu, degC) related to the in-situ pressure 
! (rfpres, dbar) with the help of an iterative method.

      USE mod_kinds

      integer, intent(in)   :: i, j
      real(r8), intent(in)  :: Sal, Pres,theta,RPres
      real(r8), intent(in)  :: tl_Sal, tl_Pres,tl_theta,tl_Rpres
      real(r8), intent(out) :: Temp
      real(r8), intent(out) :: tl_Temp

      integer               :: ind
      real(r8)              :: tpmd, theta1, thetad, epsi
      real(r8)              :: tl_thetad, tl_epsi

      data tpmd / 0.001_r8 /

      epsi = 0.0_r8
      tl_epsi = 0.0_r8
      do ind=1,100

      Temp   = theta+epsi
      tl_Temp   = tl_theta+tl_epsi


      thetad  = thetaa(Sal,Temp,Pres,RPres)-theta
      tl_thetad  = tl_thetaa(Sal,tl_Sal,Temp,tl_Temp,Pres,tl_Pres,      &
     &                                       RPres,tl_RPres)-tl_theta

!   Perform 100 iterations in TL routine always.
!     IF(abs(thetad).lt.tpmd) return
       epsi = epsi-thetad
       tl_epsi = tl_epsi-tl_thetad

      enddo
!     write(6,*) ' WARNING!',                                           &
!    & ' in-situ temperature calculation has not converged!', i,j
      RETURN
      END SUBROUTINE tl_potit
! *********************************************************************

# endif

#endif
      END MODULE tl_iceshelf_vbc_mod
