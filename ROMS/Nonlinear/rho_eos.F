#include "cppdefs.h"
      MODULE rho_eos_mod
#ifdef SOLVE3D
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rho_eos
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rho_eos (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 14)
# endif
      CALL rho_eos_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef VAR_RHO_2D
     &                   GRID(ng) % Hz,                                 &
# endif
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   OCEAN(ng) % t,                                 &
# ifdef VAR_RHO_2D
     &                   COUPLING(ng) % rhoA,                           &
     &                   COUPLING(ng) % rhoS,                           &
# endif
# ifdef BV_FREQUENCY
     &                   MIXING(ng) % bvf,                              &
# endif
# if defined LMD_SKPP    || defined LMD_BKPP         || \
     defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                   MIXING(ng) % alpha,                            &
     &                   MIXING(ng) % beta,                             &
#  ifdef LMD_DDMIX
     &                   MIXING(ng) % alfaobeta,                        &
#  endif
# endif
# if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
     &                   OCEAN(ng) % pden,                              &
# endif
     &                   OCEAN(ng) % rho)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 14)
# endif

      RETURN
      END SUBROUTINE rho_eos

# ifdef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
#  ifdef MASKING
     &                         rmask,                                   &
#  endif
#  ifdef VAR_RHO_2D
     &                         Hz,                                      &
#  endif
     &                         z_r, z_w, t,                             &
#  ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
#  endif
#  ifdef BV_FREQUENCY
     &                         bvf,                                     &
#  endif
#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                         alpha, beta,                             &
#   ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#   endif
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
     &                         pden,                                    &
#  endif
     &                         rho)
!***********************************************************************
!
      USE mod_param
#  ifdef GSW_EOS
      USE mod_gswcoef
#  else
      USE mod_eoscoef
#  endif
      USE mod_scalars
#  ifdef SEDIMENT
      USE mod_sediment
#  endif
!
      USE exchange_2d_mod
      USE exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k

      real(r8) :: SedDen, Tp, Tpr10, Ts, Tt, sqrtTs
#  ifdef BV_FREQUENCY
      real(r8) :: bulk_dn, bulk_up, den_dn, den_up
#  endif
      real(r8) :: cff, cff1, cff2

      real(r8), dimension(0:9) :: C
#  ifdef EOS_TDERIVATIVE
      real(r8), dimension(0:9) :: dCdT(0:9)

      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: DbulkDT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dden1DT
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Scof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Tcof
      real(r8), dimension(IminS:ImaxS,N(ng)) :: wrk
#  endif
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk0
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bulk2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den
      real(r8), dimension(IminS:ImaxS,N(ng)) :: den1

      real (r8) :: pp, rz0, rz1, rz2, rz3, sau, ss, ctu, zu, deltas
      real (r8) :: sa, ct, p, r1

#  include "set_bounds.h"
!
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!=======================================================================
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
#ifdef GSW_EOS

! CHECK VALUES (SA=30g/kg, CT=10degC, p=1e3dbar):
!  rho = 1027.45140
!  a   = 0.179649406
!  b   = 0.765554495
!  r1  = 1.00447333
!  rdot= 1022.87574
!
! REFERENCES:
!  Roquet, F., Madec, G., McDougall, T. J., Barker, P. M., 2014:
!Accurate 
!   polynomial expressions for the density and specific volume of 
!   seawater using the TEOS-10 standard. Ocean Modelling.
!  McDougall, T. J., D. R. Jackett, D. G. Wright and R. Feistel, 2003: 
!   Accurate and computationally efficient algorithms for potential 
!   temperature and density of seawater.  Journal of Atmospheric and 
!   Oceanic Technology, 20, 730-741. 
!
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!
            ct=MAX(-5.0_r8,t(i,j,k,nrhs,itemp))
#  ifdef SALINITY
            sa=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
#  else
            sa=0.0_r8
#  endif
            p=-z_r(i,j,k)
! reduced variables
            sau = 40.0_r8*gsw_ups
            ctu = 40.0_r8
            zu = 1e4_r8
            deltas = 32.0_r8
            ss   = sqrt ( (sa+deltas)/sau )
            tt   = ct / ctu
            pp   = p / zu
            r1=(((((r15*pp+r14)*pp+r13)*pp+r12)*pp+r11)*pp+r10)*pp+     &
     &         1.0_r8
! stiffened density
            rz3=r013*tt + r103*ss + r003
            rz2=(r022*tt+r112*ss+r012)*tt+(r202*ss+r102)*ss+r002
            rz1=(((r041*tt+r131*ss+r031)*tt+                            &
     &          (r221*ss+r121)*ss+r021)*tt+                             &
     &          ((r311*ss+r211)*ss+r111)*ss+r011)*tt+                   &
     &          (((r401*ss+r301)*ss+r201)*ss+r101)*ss+r001

            den1(i,k)=(((((r060*tt+r150*ss+r050)*tt+                    &
     &          (r240*ss+r140)*ss+r040)*tt+                             &
     &          ((r330*ss+r230)*ss+r130)*ss+r030)*tt+                   &
     &          (((r420*ss+r320)*ss+r220)*ss+r120)*ss+r020)*tt+         &
     &        ((((r510*ss+r410)*ss+r310)*ss+r210)*ss+r110)*ss+r010)*tt+ &
     &    (((((r600*ss+r500)*ss+r400)*ss+r300)*ss+r200)*ss+r100)*ss+r000

            den(i,k)=(((rz3*pp+rz2)*pp+rz1)*pp+den1(i,k))-1000.0_r8

# if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!! INCOMPLETE AND WRONG ALLOCATIONS OF ALPHA, BETA. NEEDS ALPHAONBETA
! thermal expansion a
            alpha(i,k)=r1*(                                             &
               ((alp003*pp+                                             &
     &          alp012*tt + alp102*ss+alp002)*pp+                       &
     &          ((alp031*tt+                                            &
     &          alp121*ss+alp021)*tt+                                   &
     &          (alp211*ss+alp111)*ss+alp011)*tt+                       &
     &          ((alp301*ss+alp201)*ss+alp101)*ss+alp001)*pp+           &
     &          ((((alp050*tt+                                          &
     &          alp140*ss+alp040)*tt+                                   &
     &          (alp230*ss+alp130)*ss+alp030)*tt+                       &
     &          ((alp320*ss+alp220)*ss+alp120)*ss+alp020)*tt+           &
     &        (((alp410*ss+alp310)*ss+alp210)*ss+alp110)*ss+alp010)*tt+ &
     & ((((alp500*ss+alp400)*ss+alp300)*ss+alp200)*ss+alp100)*ss+alp000)

! haline contraction b
            beta(i,k)=((bet003*pp+                                      &
     &        bet012*tt + bet102*ss+bet002)*pp+                         &
     &        ((bet031*tt+                                              &
     &        bet121*ss+bet021)*tt+                                     &
     &        (bet211*ss+bet111)*ss+bet011)*tt+                         &
     &        ((bet301*ss+bet201)*ss+bet101)*ss+bet001)*pp+             &
     &        ((((bet050*tt+                                            &
     &        bet140*ss+bet040)*tt+                                     &
     &        (bet230*ss+bet130)*ss+bet030)*tt+                         &
     &        ((bet320*ss+bet220)*ss+bet120)*ss+bet020)*tt+             &
     &   (((bet410*ss+bet310)*ss+bet210)*ss+bet110)*ss+bet010)*tt+      &
     & ((((bet500*ss+bet400)*ss+bet300)*ss+bet200)*ss+bet100)*ss+bet000
 
             beta(i,k)=b(i,k)*r1/ss
#endif
          END DO
        END DO
#else
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!
            Tt=MAX(-5.0_r8,t(i,j,k,nrhs,itemp))
#  ifdef SALINITY
            Ts=MAX(0.0_r8,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
#  else
            Ts=0.0_r8
            sqrtTs=0.0_r8
#  endif
            Tp=z_r(i,j,k)
            Tpr10=0.1_r8*Tp

!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
#  ifdef EOS_TDERIVATIVE
!
            dCdT(0)=Q01+Tt*(2.0_r8*Q02+Tt*(3.0_r8*Q03+Tt*(4.0_r8*Q04+   &
     &                      Tt*5.0_r8*Q05)))
            dCdT(1)=U01+Tt*(2.0_r8*U02+Tt*(3.0_r8*U03+Tt*4.0_r8*U04))
            dCdT(2)=V01+Tt*2.0_r8*V02
#  endif
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)

#  ifdef EOS_TDERIVATIVE
!
!  Compute d(den1)/d(S) and d(den1)/d(T) derivatives used in the
!  computation of thermal expansion and saline contraction
!  coefficients.
!
            Dden1DS(i,k)=C(1)+1.5_r8*C(2)*sqrtTs+2.0_r8*W00*Ts
            Dden1DT(i,k)=dCdT(0)+Ts*(dCdT(1)+sqrtTs*dCdT(2))
#  endif
!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
#  ifdef EOS_TDERIVATIVE
!
            dCdT(3)=A01+Tt*(2.0_r8*A02+Tt*(3.0_r8*A03+Tt*4.0_r8*A04))
            dCdT(4)=B01+Tt*(2.0_r8*B02+Tt*3.0_r8*B03)
            dCdT(5)=D01+Tt*2.0_r8*D02
            dCdT(6)=E01+Tt*(2.0_r8*E02+Tt*3.0_r8*E03)
            dCdT(7)=F01+Tt*2.0_r8*F02
            dCdT(8)=G02+Tt*2.0_r8*G03
            dCdT(9)=H01+Tt*2.0_r8*H02
#  endif
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!  Compute d(bulk)/d(S) and d(bulk)/d(T) derivatives used
!  in the computation of thermal expansion and saline contraction
!  coefficients.
!
            DbulkDS(i,k)=C(4)+sqrtTs*1.5_r8*C(5)-                       &
     &                   Tp*(C(7)+sqrtTs*1.5_r8*G00-Tp*C(9))
            DbulkDT(i,k)=dCdT(3)+Ts*(dCdT(4)+sqrtTs*dCdT(5))-           &
     &                   Tp*(dCdT(6)+Ts*dCdT(7)-                        &
     &                       Tp*(dCdT(8)+Ts*dCdT(9)))
#  endif
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0_r8/(bulk(i,k)+Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
#  if defined SEDIMENT && defined SED_DENS
            SedDen=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               t(i,j,k,nrhs,itrc)*                                &
     &               (Srho(ised,ng)-den(i,k))*cff1
            END DO
            den(i,k)=den(i,k)+SedDen
#  endif
            den(i,k)=den(i,k)-1000.0_r8
#  ifdef MASKING
            den(i,k)=den(i,k)*rmask(i,j)
#  endif
          END DO
        END DO

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  (rhoS) used in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=den(i,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=den(i,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif

#  if defined BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points:
!
!                  bvf = - g/rho d(rho)/d(z).
!
!  The density anomaly difference is computed by lowering/rising the
!  water parcel above/below adiabatically at W-point depth "z_w".
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
            bulk_up=bulk0(i,k+1)-                                       &
     &              z_w(i,j,k)*(bulk1(i,k+1)-                           &
     &                          bulk2(i,k+1)*z_w(i,j,k))
            bulk_dn=bulk0(i,k  )-                                       &
     &              z_w(i,j,k)*(bulk1(i,k  )-                           &
     &                          bulk2(i,k  )*z_w(i,j,k))
            cff1=1.0_r8/(bulk_up+0.1_r8*z_w(i,j,k))
            cff2=1.0_r8/(bulk_dn+0.1_r8*z_w(i,j,k))
            den_up=cff1*(den1(i,k+1)*bulk_up)
            den_dn=cff2*(den1(i,k  )*bulk_dn)
            bvf(i,j,k)=-g*(den_up-den_dn)/                              &
     &                 (0.5_r8*(den_up+den_dn)*                         &
     &                  (z_r(i,j,k+1)-z_r(i,j,k)))
          END DO
        END DO
        DO i=IstrT,IendT
!!        bvf(i,j,0)=bvf(i,j,1)
!!        bvf(i,j,N(ng))=bvf(i,j,N(ng)-1)
          bvf(i,j,0)=0.0_r8
          bvf(i,j,N(ng))=0.0_r8
        END DO
#  endif

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
#   ifdef LMD_DDMIX
        DO k=1,N(ng)
#   else
        DO k=N(ng),N(ng)
#   endif
          DO i=IstrT,IendT
            Tpr10=0.1_r8*z_r(i,j,k)
!
!  Compute thermal expansion and saline contraction coefficients.
!
            cff=bulk(i,k)+Tpr10
            cff1=Tpr10*den1(i,k)
            cff2=bulk(i,k)*cff
            wrk(i,k)=(den(i,k)+1000.0_r8)*cff*cff
            Tcof(i,k)=-(DbulkDT(i,k)*cff1+                              &
     &                  Dden1DT(i,k)*cff2)
            Scof(i,k)= (DbulkDS(i,k)*cff1+                              &
     &                  Dden1DS(i,k)*cff2)
#   ifdef LMD_DDMIX
            alfaobeta(i,j,k)=Tcof(i,k)/Scof(i,k)
#   endif
          END DO
          IF (k.eq.N(ng)) THEN
            DO i=IstrT,IendT
              cff=1.0_r8/wrk(i,N(ng))
              alpha(i,j)=cff*Tcof(i,N(ng))
              beta (i,j)=cff*Scof(i,N(ng))
            END DO
          END IF
        END DO
#  endif
#endif
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) and potential
!  density anomaly (kg/m3 - 1000) referenced to the surface into global
!  arrays. Notice that this is done in a separate (i,k) DO-loops to
!  facilitate the adjoint.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=den(i,k)
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
            pden(i,j,k)=(den1(i,k)-1000.0_r8)
#   ifdef MASKING
            pden(i,j,k)=pden(i,j,k)*rmask(i,j)
#   endif
#  endif
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)

#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          pden)
#  endif

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
#   ifdef LMD_DDMIX
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          alfaobeta)
#   endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          beta)
#  endif

#  ifdef VAR_RHO_2D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
#  endif

#  ifdef BV_FREQUENCY
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          bvf)
#  endif
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rho)

#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pden)
#   endif

#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
#    ifdef LMD_DDMIX
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alfaobeta)
#    endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alpha, beta)
#   endif

#   ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoA, rhoS)
#   endif

#   ifdef BV_FREQUENCY
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bvf)
#   endif
#  endif

      RETURN
      END SUBROUTINE rho_eos_tile
# endif

# ifndef NONLIN_EOS
!
!***********************************************************************
      SUBROUTINE rho_eos_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
#  ifdef MASKING
     &                         rmask,                                   &
#  endif
#  ifdef VAR_RHO_2D
     &                         Hz,                                      &
#  endif
     &                         z_r, z_w, t,                             &
#  ifdef VAR_RHO_2D
     &                         rhoA, rhoS,                              &
#  endif
#  ifdef BV_FREQUENCY
     &                         bvf,                                     &
#  endif
#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
     &                         alpha, beta,                             &
#   ifdef LMD_DDMIX
     &                         alfaobeta,                               &
#   endif
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
     &                         pden,                                    &
#  endif
     &                         rho)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
#  ifdef SEDIMENT
      USE mod_sediment
#  endif
!
      USE exchange_2d_mod
      USE exchange_3d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:,LBj:)
      real(r8), intent(out) :: rhoS(LBi:,LBj:)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:,LBj:,0:)
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:,LBj:)
      real(r8), intent(out) :: beta(LBi:,LBj:)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:,LBj:,0:)
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: rho(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef VAR_RHO_2D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   ifdef VAR_RHO_2D
      real(r8), intent(out) :: rhoA(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: rhoS(LBi:UBi,LBj:UBj)
#   endif
#   ifdef BV_FREQUENCY
      real(r8), intent(out) :: bvf(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
      real(r8), intent(out) :: alpha(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: beta(LBi:UBi,LBj:UBj)
#    ifdef LMD_DDMIX
      real(r8), intent(out) :: alfaobeta(LBi:UBi,LBj:UBj,0:N(ng))
#    endif
#   endif
#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      real(r8), intent(out) :: pden(LBi:UBi,LBj:UBj,N(ng))
#   endif
      real(r8), intent(out) :: rho(LBi:UBi,LBj:UBj,N(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, ised, itrc, j, k
      real(r8) :: SedDen, cff, cff1, cff2

#  include "set_bounds.h"
!
!=======================================================================
!  Linear equation of state.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
!  equation of state.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO k=1,N(ng)
          DO i=IstrT,IendT
            rho(i,j,k)=R0(ng)-                                          &
     &                 R0(ng)*Tcoef(ng)*(t(i,j,k,nrhs,itemp)-T0(ng))
#  ifdef SALINITY
            rho(i,j,k)=rho(i,j,k)+                                      &
     &                 R0(ng)*Scoef(ng)*(t(i,j,k,nrhs,isalt)-S0(ng))
#  endif
#  if defined SEDIMENT && defined SED_DENS
            SedDen=0.0_r8
            DO ised=1,NST
              itrc=idsed(ised)
              cff1=1.0_r8/Srho(ised,ng)
              SedDen=SedDen+                                            &
     &               t(i,j,k,nrhs,itrc)*                                &
     &               (Srho(ised,ng)-rho(i,j,k))*cff1
            END DO
            rho(i,j,k)=rho(i,j,k)+SedDen
#  endif
            rho(i,j,k)=rho(i,j,k)-1000.0_r8
#  ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#  endif
#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
            pden(i,j,k)=rho(i,j,k)
#  endif
          END DO
        END DO

#  ifdef VAR_RHO_2D
!
!-----------------------------------------------------------------------
!  Compute vertical averaged density (rhoA) and density perturbation
!  used (rhoS) in barotropic pressure gradient.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          cff1=rho(i,j,N(ng))*Hz(i,j,N(ng))
          rhoS(i,j)=0.5_r8*cff1*Hz(i,j,N(ng))
          rhoA(i,j)=cff1
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrT,IendT
            cff1=rho(i,j,k)*Hz(i,j,k)
            rhoS(i,j)=rhoS(i,j)+Hz(i,j,k)*(rhoA(i,j)+0.5_r8*cff1)
            rhoA(i,j)=rhoA(i,j)+cff1
          END DO
        END DO
        cff2=1.0_r8/rho0
        DO i=IstrT,IendT
          cff1=1.0_r8/(z_w(i,j,N(ng))-z_w(i,j,0))
          rhoA(i,j)=cff2*cff1*rhoA(i,j)
          rhoS(i,j)=2.0_r8*cff1*cff1*cff2*rhoS(i,j)
        END DO
#  endif

#  ifdef BV_FREQUENCY
!
!-----------------------------------------------------------------------
!  Compute Brunt-Vaisala frequency (1/s2) at horizontal RHO-points
!  and vertical W-points.
!-----------------------------------------------------------------------
!
        DO k=1,N(ng)-1
          DO i=IstrT,IendT
            bvf(i,j,k)=-gorho0*(rho(i,j,k+1)-rho(i,j,k))/               &
     &                         (z_r(i,j,k+1)-z_r(i,j,k))
          END DO
        END DO
#  endif

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Compute thermal expansion (1/Celsius) and saline contraction
!  (1/PSU) coefficients.
!-----------------------------------------------------------------------
!
        DO i=IstrT,IendT
          alpha(i,j)=ABS(Tcoef(ng))
#   ifdef SALINITY
          beta(i,j)=ABS(Scoef(ng))
#   else
          beta(i,j)=0.0_r8
#   endif
        END DO
#   ifdef LMD_DDMIX
!
!  Compute ratio of thermal expansion and saline contraction
!  coefficients.
!
        IF (Scoef(ng).eq.0.0_r8) THEN
          cff=1.0_r8
        ELSE
          cff=1.0_r8/Scoef(ng)
        END IF
        DO k=1,N(ng)
          DO i=IstrT,IendT
            alfaobeta(i,j,k)=cff*Tcoef(ng)
          END DO
        END DO
#   endif
#  endif
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rho)

#  if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          pden)
#  endif

#  if defined LMD_SKPP    || defined LMD_BKPP         || \
      defined BULK_FLUXES || defined BALANCE_OPERATOR
#   ifdef LMD_DDMIX
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          alfaobeta)
#   endif
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          alpha)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          beta)
#  endif

#  ifdef VAR_RHO_2D
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoA)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rhoS)
#  endif

#  ifdef BV_FREQUENCY
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          bvf)
#  endif
      END IF

#  ifdef DISTRIBUTE
!
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rho)

#   if defined LMD_SKPP || defined LMD_BKPP || defined DIAGNOSTICS
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pden)
#   endif

#   if defined LMD_SKPP    || defined LMD_BKPP         || \
       defined BULK_FLUXES || defined BALANCE_OPERATOR
#    ifdef LMD_DDMIX
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alfaobeta)
#    endif
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    alpha, beta)
#   endif

#   ifdef VAR_RHO_2D
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rhoA, rhoS)
#   endif

#   ifdef BV_FREQUENCY
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bvf)
#   endif
#  endif

      RETURN
      END SUBROUTINE rho_eos_tile
# endif
#endif
      END MODULE rho_eos_mod
