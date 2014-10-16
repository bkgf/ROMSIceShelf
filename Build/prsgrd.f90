      MODULE prsgrd_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the baroclinic hydrostatic pressure gradient  !
!  term.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: prsgrd
      CONTAINS
      SUBROUTINE prsgrd (ng, tile)
!
!svn $Id$
!***********************************************************************
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine evaluates the nonlinear  baroclinic,  hydrostatic   !
!  pressure gradient term using a  nonconservative  Density-Jacobian   !
!  scheme,  based on  cubic polynomial fits for  "rho" and  "z_r" as   !
!  functions of nondimensional coordinates (XI,ETA,s), that is,  its   !
!  respective array indices. The  cubic polynomials  are monotonized   !
!  by using  harmonic mean instead of linear averages to interpolate   !
!  slopes. This scheme retains exact anti-symmetry:                    !
!                                                                      !
!        J(rho,z_r)=-J(z_r,rho).                                       !
!                                                                      !
!  If parameter OneFifth (below) is set to zero,  the scheme becomes   !
!  identical to standard Jacobian.                                     !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Shchepetkin A.F and J.C. McWilliams, 2003:  A method for          !
!      computing horizontal pressure gradient force in an ocean        !
!      model with non-aligned vertical coordinate, JGR, 108,           !
!      1-34.                                                           !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
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
      CALL wclock_on (ng, iNLM, 23)
      CALL prsgrd_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nrhs(ng),                                       &
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % z_r,                                 &
     &                  GRID(ng) % z_w,                                 &
     &                  GRID(ng) % zice,                                &
     &                  OCEAN(ng) % rho,                                &
     &                  OCEAN(ng) % ru,                                 &
     &                  OCEAN(ng) % rv)
      CALL wclock_off (ng, iNLM, 23)
      RETURN
      END SUBROUTINE prsgrd
!
!***********************************************************************
      SUBROUTINE prsgrd_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs,                                     &
     &                        om_v, on_u,                               &
     &                        Hz, z_r, z_w,                             &
     &                        zice,                                     &
     &                        rho,                                      &
     &                        ru, rv)
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
      integer, intent(in) :: nrhs
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: OneFifth = 0.2_r8
      real(r8), parameter :: OneTwelfth = 1.0_r8/12.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8), parameter :: drhodz = 0.00478_r8
      real(r8) :: GRho, GRho0,  HalfGRho
      real(r8) :: cff, cff1, cff2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: P
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aux
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZx
!
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
!-----------------------------------------------------------------------
!  Preliminary step (same for XI- and ETA-components:
!-----------------------------------------------------------------------
!
      GRho=g/rho0
      GRho0=1000.0_r8*GRho
      HalfGRho=0.5_r8*GRho
!
      DO j=JstrV-1,Jend
        DO k=1,N(ng)-1
          DO i=IstrU-1,Iend
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
          END DO
        END DO
        DO i=IstrU-1,Iend
          dR(i,N(ng))=dR(i,N(ng)-1)
          dZ(i,N(ng))=dZ(i,N(ng)-1)
          dR(i,0)=dR(i,1)
          dZ(i,0)=dZ(i,1)
        END DO
        DO k=N(ng),1,-1
          DO i=IstrU-1,Iend
            cff=2.0_r8*dR(i,k)*dR(i,k-1)
            IF (cff.gt.eps) THEN
              dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
            ELSE
              dR(i,k)=0.0_r8
            END IF
            dZ(i,k)=2.0_r8*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
          END DO
        END DO
        DO i=IstrU-1,Iend
          cff1=1.0_r8/(z_r(i,j,N(ng))-z_r(i,j,N(ng)-1))
          cff2=0.5_r8*(rho(i,j,N(ng))-rho(i,j,N(ng)-1))*                &
     &         (z_w(i,j,N(ng))-z_r(i,j,N(ng)))*cff1
          P(i,j,N(ng))=GRho0*(z_w(i,j,N(ng))-zice(i,j))-                &
     &                 GRho*(rho(i,j,N(ng))+0.5_r8*drhodz*zice(i,j))*   &
     &                 zice(i,j)+                                       &
     &                 GRho*(rho(i,j,N(ng))+cff2)*                      &
     &                 (z_w(i,j,N(ng))-z_r(i,j,N(ng)))
        END DO
        DO k=N(ng)-1,1,-1
          DO i=IstrU-1,Iend
            P(i,j,k)=P(i,j,k+1)+                                        &
     &               HalfGRho*((rho(i,j,k+1)+rho(i,j,k))*               &
     &                         (z_r(i,j,k+1)-z_r(i,j,k))-               &
     &                         OneFifth*                                &
     &                         ((dR(i,k+1)-dR(i,k))*                    &
     &                          (z_r(i,j,k+1)-z_r(i,j,k)-               &
     &                           OneTwelfth*                            &
     &                           (dZ(i,k+1)+dZ(i,k)))-                  &
     &                          (dZ(i,k+1)-dZ(i,k))*                    &
     &                          (rho(i,j,k+1)-rho(i,j,k)-               &
     &                           OneTwelfth*                            &
     &                           (dR(i,k+1)+dR(i,k)))))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute XI-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N(ng),1,-1
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend+1
            aux(i,j)=z_r(i,j,k)-z_r(i-1,j,k)
            FC(i,j)=rho(i,j,k)-rho(i-1,j,k)
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff=2.0_r8*aux(i,j)*aux(i+1,j)
            IF (cff.gt.eps) THEN
              cff1=1.0_r8/(aux(i,j)+aux(i+1,j))
              dZx(i,j)=cff*cff1
            ELSE
              dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i+1,j)
            IF (cff1.gt.eps) THEN
              cff2=1.0_r8/(FC(i,j)+FC(i+1,j))
              dRx(i,j)=cff1*cff2
            ELSE
              dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            ru(i,j,k,nrhs)=on_u(i,j)*0.5_r8*                            &
     &                     (Hz(i,j,k)+Hz(i-1,j,k))*                     &
     &                     (P(i-1,j,k)-P(i,j,k)-                        &
     &                      HalfGRho*                                   &
     &                      ((rho(i,j,k)+rho(i-1,j,k))*                 &
     &                       (z_r(i,j,k)-z_r(i-1,j,k))-                 &
     &                        OneFifth*                                 &
     &                        ((dRx(i,j)-dRx(i-1,j))*                   &
     &                         (z_r(i,j,k)-z_r(i-1,j,k)-                &
     &                          OneTwelfth*                             &
     &                          (dZx(i,j)+dZx(i-1,j)))-                 &
     &                         (dZx(i,j)-dZx(i-1,j))*                   &
     &                         (rho(i,j,k)-rho(i-1,j,k)-                &
     &                          OneTwelfth*                             &
     &                          (dRx(i,j)+dRx(i-1,j))))))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N(ng),1,-1
        DO j=JstrV-1,Jend+1
          DO i=Istr,Iend
            aux(i,j)=z_r(i,j,k)-z_r(i,j-1,k)
            FC(i,j)=rho(i,j,k)-rho(i,j-1,k)
          END DO
        END DO
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff=2.0_r8*aux(i,j)*aux(i,j+1)
            IF (cff.gt.eps) THEN
              cff1=1.0_r8/(aux(i,j)+aux(i,j+1))
              dZx(i,j)=cff*cff1
            ELSE
              dZx(i,j)=0.0_r8
            END IF
            cff1=2.0_r8*FC(i,j)*FC(i,j+1)
            IF (cff1.gt.eps) THEN
              cff2=1.0_r8/(FC(i,j)+FC(i,j+1))
              dRx(i,j)=cff1*cff2
            ELSE
              dRx(i,j)=0.0_r8
            END IF
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            rv(i,j,k,nrhs)=om_v(i,j)*0.5_r8*                            &
     &                     (Hz(i,j,k)+Hz(i,j-1,k))*                     &
     &                     (P(i,j-1,k)-P(i,j,k)-                        &
     &                      HalfGRho*                                   &
     &                      ((rho(i,j,k)+rho(i,j-1,k))*                 &
     &                       (z_r(i,j,k)-z_r(i,j-1,k))-                 &
     &                        OneFifth*                                 &
     &                        ((dRx(i,j)-dRx(i,j-1))*                   &
     &                         (z_r(i,j,k)-z_r(i,j-1,k)-                &
     &                          OneTwelfth*                             &
     &                          (dZx(i,j)+dZx(i,j-1)))-                 &
     &                         (dZx(i,j)-dZx(i,j-1))*                   &
     &                         (rho(i,j,k)-rho(i,j-1,k)-                &
     &                          OneTwelfth*                             &
     &                          (dRx(i,j)+dRx(i,j-1))))))
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
      END MODULE prsgrd_mod
