      SUBROUTINE tl_t3dmix2 (ng, tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2013 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This subroutine computes tangent linear horizontal harmonic mixing  !
!  of tracers along isopycnic surfaces.                                !
!                                                                      !
!  BASIC STATE variables needed: diff2, Hz, rho, t, z_r                !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef CLIMA_TS_MIX
      USE mod_clima
#endif
#ifdef DIAGNOSTICS_TS
!!    USE mod_diags
#endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iTLM, 26)
#endif
      CALL tl_t3dmix2_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#endif
     &                      GRID(ng) % om_v,                            &
     &                      GRID(ng) % on_u,                            &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % tl_z_r,                          &
     &                      MIXING(ng) % diff2,                         &
     &                      OCEAN(ng) % rho,                            &
     &                      OCEAN(ng) % tl_rho,                         &
#ifdef CLIMA_TS_MIX
     &                      CLIMA(ng) % tclm,                           &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                      DIAGS(ng) % DiaTwrk,                        &
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)
#ifdef PROFILE
      CALL wclock_off (ng, iTLM, 26)
#endif
      RETURN
      END SUBROUTINE tl_t3dmix2
!
!***********************************************************************
      SUBROUTINE tl_t3dmix2_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nnew,                           &
#ifdef MASKING
     &                            umask, vmask,                         &
#endif
     &                            om_v, on_u, pm, pn,                   &
     &                            Hz, tl_Hz,                            &
     &                            z_r, tl_z_r,                          &
     &                            diff2,                                &
     &                            rho, tl_rho,                          &
#ifdef CLIMA_TS_MIX
     &                            tclm,                                 &
#endif
#ifdef DIAGNOSTICS_TS
!!   &                            DiaTwrk,                              &
#endif
     &                            t, tl_t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: diff2(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# ifdef CLIMA_TS_MIX
      real(r8), intent(in) :: tclm(LBi:,LBj:,:,:)
# endif
      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_rho(LBi:,LBj:,:)
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
# endif

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: diff2(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: rho(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# ifdef CLIMA_TS_MIX
      real(r8), intent(in) :: tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng))
# endif
      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: tl_rho(LBi:UBi,LBj:UBj,N(ng))
# ifdef DIAGNOSTICS_TS
!!    real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
!!   &                                   NDT)
# endif
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k, k1, k2

      real(r8), parameter :: eps = 0.5_r8
      real(r8), parameter :: small = 1.0E-14_r8
      real(r8), parameter :: slope_max = 0.0001_r8
      real(r8), parameter :: strat_min = 0.1_r8

      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3, tl_cff4

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: FS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dRde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dRdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTdr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dTdx

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_FS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dRde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dRdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dTde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dTdr
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dTdx

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute horizontal harmonic diffusion along isopycnic surfaces.
!-----------------------------------------------------------------------
!
!  Compute horizontal and density gradients.  Notice the recursive
!  blocking sequence.  The vertical placement of the gradients is:
!
!        dTdx,dTde(:,:,k1) k     rho-points
!        dTdx,dTde(:,:,k2) k+1   rho-points
!          FS,dTdr(:,:,k1) k-1/2   W-points
!          FS,dTdr(:,:,k2) k+1/2   W-points
!
      T_LOOP : DO itrc=1,NT(ng)
        k2=1
        K_LOOP : DO k=0,N(ng)
          k1=k2
          k2=3-k1
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(pm(i,j)+pm(i-1,j))
#ifdef MASKING
                cff=cff*umask(i,j)
#endif
                dRdx(i,j,k2)=cff*(rho(i  ,j,k+1)-                       &
     &                            rho(i-1,j,k+1))
                tl_dRdx(i,j,k2)=cff*(tl_rho(i  ,j,k+1)-                 &
     &                               tl_rho(i-1,j,k+1))
#ifdef CLIMA_TS_MIX
                dTdx(i,j,k2)=cff*((t(i  ,j,k+1,nrhs,itrc)-              &
     &                             tclm(i  ,j,k+1,itrc))-               &
     &                            (t(i-1,j,k+1,nrhs,itrc)-              &
     &                             tclm(i-1,j,k+1,itrc)))
#else
                dTdx(i,j,k2)=cff*(t(i  ,j,k+1,nrhs,itrc)-               &
     &                            t(i-1,j,k+1,nrhs,itrc))
#endif
                tl_dTdx(i,j,k2)=cff*(tl_t(i  ,j,k+1,nrhs,itrc)-         &
     &                               tl_t(i-1,j,k+1,nrhs,itrc))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(pn(i,j)+pn(i,j-1))
#ifdef MASKING
                cff=cff*vmask(i,j)
#endif
                dRde(i,j,k2)=cff*(rho(i,j  ,k+1)-                       &
     &                            rho(i,j-1,k+1))
                tl_dRde(i,j,k2)=cff*(tl_rho(i,j  ,k+1)-                 &
     &                               tl_rho(i,j-1,k+1))
#ifdef CLIMA_TS_MIX
                dTde(i,j,k2)=cff*((t(i,j  ,k+1,nrhs,itrc)-              &
     &                             tclm(i,j  ,k+1,itrc))-               &
     &                            (t(i,j-1,k+1,nrhs,itrc)-              &
     &                             tclm(i,j-1,k+1,itrc)))
#else
                dTde(i,j,k2)=cff*(t(i,j  ,k+1,nrhs,itrc)-               &
     &                            t(i,j-1,k+1,nrhs,itrc))
#endif
                tl_dTde(i,j,k2)=cff*(tl_t(i,j  ,k+1,nrhs,itrc)-         &
     &                               tl_t(i,j-1,k+1,nrhs,itrc))
              END DO
            END DO
          END IF
          IF ((k.eq.0).or.(k.eq.N(ng))) THEN
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
                dTdr(i,j,k2)=0.0_r8
                tl_dTdr(i,j,k2)=0.0_r8
                FS(i,j,k2)=0.0_r8
                tl_FS(i,j,k2)=0.0_r8
              END DO
            END DO
          ELSE
            DO j=Jstr-1,Jend+1
              DO i=Istr-1,Iend+1
#if defined MAX_SLOPE
                cff1=SQRT(dRdx(i,j,k2)**2+dRdx(i+1,j,k2)**2+            &
     &                    dRdx(i,j,k1)**2+dRdx(i+1,j,k1)**2+            &
     &                    dRde(i,j,k2)**2+dRde(i,j+1,k2)**2+            &
     &                    dRde(i,j,k1)**2+dRde(i,j+1,k1)**2)
                IF (cff1.ne.0.0_r8) THEN
                  tl_cff1=(dRdx(i  ,j,k2)*tl_dRdx(i  ,j,k2)+            &
     &                     dRdx(i+1,j,k2)*tl_dRdx(i+1,j,k2)+            &
     &                     dRdx(i  ,j,k1)*tl_dRdx(i  ,j,k1)+            &
     &                     dRdx(i+1,j,k1)*tl_dRdx(i+1,j,k1)+            &
     &                     dRde(i,j  ,k2)*tl_dRde(i,j  ,k2)+            &
     &                     dRde(i,j+1,k2)*tl_dRde(i,j+1,k2)+            &
     &                     dRde(i,j  ,k1)*tl_dRde(i,j  ,k1)+            &
     &                     dRde(i,j+1,k1)*tl_dRde(i,j+1,k1))/cff1
                ELSE
                  tl_cff1=0.0_r8
                END IF
                cff2=0.25_r8*slope_max*                                 &
     &               (z_r(i,j,k+1)-z_r(i,j,k))*cff1
                tl_cff2=0.25_r8*slope_max*                              &
     &                  ((tl_z_r(i,j,k+1)-tl_z_r(i,j,k))*cff1+          &
     &                   (z_r(i,j,k+1)-z_r(i,j,k))*tl_cff1)
                cff3=MAX(rho(i,j,k)-rho(i,j,k+1),small)
                tl_cff3=(0.5_r8+SIGN(0.5_r8,rho(i,j,k)-rho(i,j,k+1)-    &
     &                                      small))*                    &
     &                  (tl_rho(i,j,k)-tl_rho(i,j,k+1))
                cff4=MAX(cff2,cff3)
                tl_cff4=(0.5_r8+SIGN(0.5_r8,cff2-cff3))*tl_cff2+        &
     &                  (0.5_r8-SIGN(0.5_r8,cff2-cff3))*tl_cff3
                cff=-1.0_r8/cff4
                tl_cff=cff*cff*tl_cff4
#elif defined MIN_STRAT
                cff1=MAX(rho(i,j,k)-rho(i,j,k+1),                       &
     &                   strat_min*(z_r(i,j,k+1)-z_r(i,j,k)))
                tl_cff1=(0.5_r8+SIGN(0.5_r8,                            &
     &                               rho(i,j,k)-rho(i,j,k+1)-           &
     &                               strat_min*(z_r(i,j,k+1)-           &
     &                                          z_r(i,j,k  ))))*        &
     &                  (tl_rho(i,j,k)-tl_rho(i,j,k+1))+                &
     &                  (0.5_r8-SIGN(0.5_r8,                            &
     &                               rho(i,j,k)-rho(i,j,k+1)-           &
     &                               strat_min*(z_r(i,j,k+1)-           &
     &                                          z_r(i,j,k  ))))*        &
     &                  (strat_min*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k  )))
                cff=-1.0_r8/cff1
                tl_cff=cff*cff*tl_cff1
#else
                cff1=MAX(rho(i,j,k)-rho(i,j,k+1),eps)
                tl_cff1=(0.5_r8+SIGN(0.5_r8,                            &
     &                               rho(i,j,k)-rho(i,j,k+1)-eps))*     &
     &                  (tl_rho(i,j,k)-tl_rho(i,j,k+1))
                cff=-1.0_r8/cff1
                tl_cff=cff*cff*tl_cff1
#endif
#ifdef CLIMA_TS_MIX
                dTdr(i,j,k2)=cff*((t(i,j,k+1,nrhs,itrc)-                &
     &                             tclm(i,j,k+1,itrc))-                 &
     &                            (t(i,j,k  ,nrhs,itrc)-                &
     &                             tclm(i,j,k  ,itrc)))
#else
                dTdr(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-                 &
     &                            t(i,j,k  ,nrhs,itrc))
#endif
#ifdef CLIMA_TS_MIX
                tl_dTdr(i,j,k2)=tl_cff*((t(i,j,k+1,nrhs,itrc)-          &
     &                                   tclm(i,j,k+1,itrc))-           &
     &                                  (t(i,j,k  ,nrhs,itrc)-          &
     &                                   tclm(i,j,k  ,itrc)))+          &
     &                          cff*(tl_t(i,j,k+1,nrhs,itrc)-           &
     &                               tl_t(i,j,k  ,nrhs,itrc))
#else
                tl_dTdr(i,j,k2)=tl_cff*(t(i,j,k+1,nrhs,itrc)-           &
     &                                  t(i,j,k  ,nrhs,itrc))+          &
     &                          cff*(tl_t(i,j,k+1,nrhs,itrc)-           &
     &                               tl_t(i,j,k  ,nrhs,itrc))
#endif
                FS(i,j,k2)=cff*(z_r(i,j,k+1)-z_r(i,j,k))
                tl_FS(i,j,k2)=tl_cff*(z_r(i,j,k+1)-z_r(i,j,k))+         &
     &                        cff*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))
              END DO
            END DO
          END IF
!
!  Compute components of the rotated tracer flux (T m4/s) along
!  isopycnic surfaces.
!
          IF (k.gt.0) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.25_r8*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*        &
     &              on_u(i,j)
!>              FX(i,j)=cff*                                            &
!>   &                  (Hz(i,j,k)+Hz(i-1,j,k))*                        &
!>   &                  (dTdx(i,j,k1)-                                  &
!>   &                   0.5_r8*(MAX(dRdx(i,j,k1),0.0_r8)*              &
!>   &                              (dTdr(i-1,j,k1)+                    &
!>   &                               dTdr(i  ,j,k2))+                   &
!>   &                           MIN(dRdx(i,j,k1),0.0_r8)*              &
!>   &                              (dTdr(i-1,j,k2)+                    &
!>   &                               dTdr(i,j,k1))))
!>
                tl_FX(i,j)=cff*                                         &
     &                     (((tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*             &
     &                       (dTdx(i,j,k1)-                             &
     &                        0.5_r8*(MAX(dRdx(i,j,k1),0.0_r8)*         &
     &                                   (dTdr(i-1,j,k1)+               &
     &                                    dTdr(i  ,j,k2))+              &
     &                                MIN(dRdx(i,j,k1),0.0_r8)*         &
     &                                   (dTdr(i-1,j,k2)+               &
     &                                    dTdr(i  ,j,k1)))))+           &
     &                      ((Hz(i,j,k)+Hz(i-1,j,k))*                   &
     &                       (tl_dTdx(i,j,k1)-                          &
     &                        0.5_r8*(MAX(dRdx(i,j,k1),0.0_r8)*         &
     &                                   (tl_dTdr(i-1,j,k1)+            &
     &                                    tl_dTdr(i  ,j,k2))+           &
     &                                MIN(dRdx(i,j,k1),0.0_r8)*         &
     &                                   (tl_dTdr(i-1,j,k2)+            &
     &                                    tl_dTdr(i  ,j,k1)))-          &
     &                        0.5_r8*((0.5_r8+                          &
     &                                 SIGN(0.5_r8, dRdx(i,j,k1)))*     &
     &                                tl_dRdx(i,j,k1)*                  &
     &                                (dTdr(i-1,j,k1)+dTdr(i,j,k2))+    &
     &                                (0.5_r8+                          &
     &                                 SIGN(0.5_r8,-dRdx(i,j,k1)))*     &
     &                                tl_dRdx(i,j,k1)*                  &
     &                                (dTdr(i-1,j,k2)+dTdr(i,j,k1))))))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.25_r8*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*        &
     &              om_v(i,j)
!>              FE(i,j)=cff*                                            &
!>   &                  (Hz(i,j,k)+Hz(i,j-1,k))*                        &
!>   &                  (dTde(i,j,k1)-                                  &
!>   &                   0.5_r8*(MAX(dRde(i,j,k1),0.0_r8)*              &
!>   &                              (dTdr(i,j-1,k1)+                    &
!>   &                               dTdr(i,j  ,k2))+                   &
!>   &                           MIN(dRde(i,j,k1),0.0_r8)*              &
!>   &                              (dTdr(i,j-1,k2)+                    &
!>   &                               dTdr(i,j  ,k1))))
!>
                tl_FE(i,j)=cff*                                         &
     &                     (((tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*             &
     &                       (dTde(i,j,k1)-                             &
     &                        0.5_r8*(MAX(dRde(i,j,k1),0.0_r8)*         &
     &                                   (dTdr(i,j-1,k1)+               &
     &                                    dTdr(i,j  ,k2))+              &
     &                                MIN(dRde(i,j,k1),0.0_r8)*         &
     &                                   (dTdr(i,j-1,k2)+               &
     &                                    dTdr(i,j  ,k1)))))+           &
     &                      ((Hz(i,j,k)+Hz(i,j-1,k))*                   &
     &                       (tl_dTde(i,j,k1)-                          &
     &                        0.5_r8*(MAX(dRde(i,j,k1),0.0_r8)*         &
     &                                   (tl_dTdr(i,j-1,k1)+            &
     &                                    tl_dTdr(i,j  ,k2))+           &
     &                                MIN(dRde(i,j,k1),0.0_r8)*         &
     &                                   (tl_dTdr(i,j-1,k2)+            &
     &                                    tl_dTdr(i,j  ,k1)))-          &
     &                        0.5_r8*((0.5_r8+                          &
     &                                 SIGN(0.5_r8, dRde(i,j,k1)))*     &
     &                                tl_dRde(i,j,k1)*                  &
     &                                (dTdr(i,j-1,k1)+dTdr(i,j,k2))+    &
     &                                (0.5_r8+                          &
     &                                 SIGN(0.5_r8,-dRde(i,j,k1)))*     &
     &                                tl_dRde(i,j,k1)*                  &
     &                                (dTdr(i,j-1,k2)+dTdr(i,j,k1))))))
              END DO
            END DO
            IF (k.lt.N(ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  cff1=MAX(dRdx(i  ,j,k1),0.0_r8)
                  cff2=MAX(dRdx(i+1,j,k2),0.0_r8)
                  cff3=MIN(dRdx(i  ,j,k2),0.0_r8)
                  cff4=MIN(dRdx(i+1,j,k1),0.0_r8)
                  tl_cff1=(0.5_r8+SIGN(0.5_r8, dRdx(i  ,j,k1)))*        &
     &                    tl_dRdx(i  ,j,k1)
                  tl_cff2=(0.5_r8+SIGN(0.5_r8, dRdx(i+1,j,k2)))*        &
     &                    tl_dRdx(i+1,j,k2)
                  tl_cff3=(0.5_r8+SIGN(0.5_r8,-dRdx(i  ,j,k2)))*        &
     &                    tl_dRdx(i  ,j,k2)
                  tl_cff4=(0.5_r8+SIGN(0.5_r8,-dRdx(i+1,j,k1)))*        &
     &                    tl_dRdx(i+1,j,k1)
                  cff=cff1*(cff1*dTdr(i,j,k2)-dTdx(i  ,j,k1))+          &
     &                cff2*(cff2*dTdr(i,j,k2)-dTdx(i+1,j,k2))+          &
     &                cff3*(cff3*dTdr(i,j,k2)-dTdx(i  ,j,k2))+          &
     &                cff4*(cff4*dTdr(i,j,k2)-dTdx(i+1,j,k1))
                  tl_cff=tl_cff1*(cff1*dTdr(i  ,j,k2)-                  &
     &                            dTdx(i  ,j,k1))+                      &
     &                   tl_cff2*(cff2*dTdr(i,j,k2)-                    &
     &                            dTdx(i+1,j,k2))+                      &
     &                   tl_cff3*(cff3*dTdr(i,j,k2)-                    &
     &                            dTdx(i  ,j,k2))+                      &
     &                   tl_cff4*(cff4*dTdr(i,j,k2)-                    &
     &                            dTdx(i+1,j,k1))+                      &
     &                   cff1*(tl_cff1*dTdr(i,j,k2)+                    &
     &                         cff1*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTdx(i  ,j,k1))+                      &
     &                   cff2*(tl_cff2*dTdr(i,j,k2)+                    &
     &                         cff2*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTdx(i+1,j,k2))+                      &
     &                   cff3*(tl_cff3*dTdr(i,j,k2)+                    &
     &                         cff3*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTdx(i  ,j,k2))+                      &
     &                   cff4*(tl_cff4*dTdr(i,j,k2)+                    &
     &                         cff4*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTdx(i+1,j,k1))
                  cff1=MAX(dRde(i,j  ,k1),0.0_r8)
                  cff2=MAX(dRde(i,j+1,k2),0.0_r8)
                  cff3=MIN(dRde(i,j  ,k2),0.0_r8)
                  cff4=MIN(dRde(i,j+1,k1),0.0_r8)
                  tl_cff1=(0.5_r8+SIGN(0.5_r8, dRde(i,j  ,k1)))*        &
     &                    tl_dRde(i,j  ,k1)
                  tl_cff2=(0.5_r8+SIGN(0.5_r8, dRde(i,j+1,k2)))*        &
     &                    tl_dRde(i,j+1,k2)
                  tl_cff3=(0.5_r8+SIGN(0.5_r8,-dRde(i,j  ,k2)))*        &
     &                    tl_dRde(i,j  ,k2)
                  tl_cff4=(0.5_r8+SIGN(0.5_r8,-dRde(i,j+1,k1)))*        &
     &                    tl_dRde(i,j+1,k1)
                  cff=cff+                                              &
     &                cff1*(cff1*dTdr(i,j,k2)-dTde(i,j  ,k1))+          &
     &                cff2*(cff2*dTdr(i,j,k2)-dTde(i,j+1,k2))+          &
     &                cff3*(cff3*dTdr(i,j,k2)-dTde(i,j  ,k2))+          &
     &                cff4*(cff4*dTdr(i,j,k2)-dTde(i,j+1,k1))
                  tl_cff=tl_cff+                                        &
     &                   tl_cff1*(cff1*dTdr(i,j,k2)-                    &
     &                            dTde(i,j  ,k1))+                      &
     &                   tl_cff2*(cff2*dTdr(i,j,k2)-                    &
     &                            dTde(i,j+1,k2))+                      &
     &                   tl_cff3*(cff3*dTdr(i,j,k2)-                    &
     &                            dTde(i,j  ,k2))+                      &
     &                   tl_cff4*(cff4*dTdr(i,j,k2)-                    &
     &                            dTde(i,j+1,k1))+                      &
     &                   cff1*(tl_cff1*dTdr(i,j,k2)+                    &
     &                         cff1*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTde(i,j  ,k1))+                      &
     &                   cff2*(tl_cff2*dTdr(i,j,k2)+                    &
     &                         cff2*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTde(i,j+1,k2))+                      &
     &                   cff3*(tl_cff3*dTdr(i,j,k2)+                    &
     &                         cff3*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTde(i,j  ,k2))+                      &
     &                   cff4*(tl_cff4*dTdr(i,j,k2)+                    &
     &                         cff4*tl_dTdr(i,j,k2)-                    &
     &                         tl_dTde(i,j+1,k1))
!>                FS(i,j,k2)=0.5_r8*cff*diff2(i,j,itrc)*FS(i,j,k2)
!>
                  tl_FS(i,j,k2)=0.5_r8*diff2(i,j,itrc)*                 &
     &                          (tl_cff*FS(i,j,k2)+                     &
     &                           cff*tl_FS(i,j,k2))
                END DO
              END DO
            END IF
!
!  Time-step harmonic, isopycnic diffusion term (m Tunits).
!
            DO j=Jstr,Jend
              DO i=Istr,Iend
!>              cff=dt(ng)*pm(i,j)*pn(i,j)*                             &
!>   &                     (FX(i+1,j)-FX(i,j)+                          &
!>   &                      FE(i,j+1)-FE(i,j))+                         &
!>   &              dt(ng)*(FS(i,j,k2)-FS(i,j,k1))
!>
                tl_cff=dt(ng)*pm(i,j)*pn(i,j)*                          &
     &                        (tl_FX(i+1,j)-tl_FX(i,j)+                 &
     &                         tl_FE(i,j+1)-tl_FE(i,j))+                &
     &                 dt(ng)*(tl_FS(i,j,k2)-tl_FS(i,j,k1))
!>              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff
!>
                tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff
#ifdef DIAGNOSTICS_TS
!!              DiaTwrk(i,j,k,itrc,iThdif)=cff
#endif
              END DO
            END DO
          END IF
        END DO K_LOOP
      END DO T_LOOP
      RETURN
      END SUBROUTINE tl_t3dmix2_tile
