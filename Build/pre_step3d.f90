      MODULE pre_step3d_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine initialize computations for new time step of the    !
!  3D primitive variables.                                             !
!                                                                      !
!  Both n-1 and n-2 time-step contributions of the  Adams/Bashforth    !
!  scheme are added here to u and v at time index "nnew", since the    !
!  right-hand-side  arrays ru and rv at  n-2  will be overwriten in    !
!  subsequent calls to routines within the 3D engine.                  !
!                                                                      !
!  It also computes the time  "n"  vertical viscosity and diffusion    !
!  contributions of the Crank-Nicholson implicit scheme because the    !
!  thicknesses "Hz" will be overwriten at the end of the  2D engine    !
!  (barotropic mode) computations.                                     !
!                                                                      !
!  The actual time step will be carried out in routines "step3d_uv"    !
!  and "step3d_t".                                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: pre_step3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE pre_step3d (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
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
      CALL wclock_on (ng, iNLM, 22)
      CALL pre_step3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nrhs(ng), nstp(ng), nnew(ng),               &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % Huon,                            &
     &                      GRID(ng) % Hvom,                            &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % z_w,                             &
     &                      FORCES(ng) % btflx,                         &
     &                      FORCES(ng) % bustr,                         &
     &                      FORCES(ng) % bvstr,                         &
     &                      FORCES(ng) % stflx,                         &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr,                         &
     &                      MIXING(ng) % Akt,                           &
     &                      MIXING(ng) % Akv,                           &
     &                      OCEAN(ng) % W,                              &
     &                      OCEAN(ng) % ru,                             &
     &                      OCEAN(ng) % rv,                             &
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v)
      CALL wclock_off (ng, iNLM, 22)
      RETURN
      END SUBROUTINE pre_step3d
!
!***********************************************************************
      SUBROUTINE pre_step3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nrhs, nstp, nnew,                     &
     &                            rmask, umask, vmask,                  &
     &                            pm, pn,                               &
     &                            Hz, Huon, Hvom,                       &
     &                            z_r, z_w,                             &
     &                            btflx, bustr, bvstr,                  &
     &                            stflx, sustr, svstr,                  &
     &                            Akt, Akv,                             &
     &                            W,                                    &
     &                            ru, rv,                               &
     &                            t, u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: btflx(LBi:,LBj:,:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(in) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: rv(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
!
!  Local variable declarations.
!
      integer :: i, indx, is, itrc, j, k, ltrc
      real(r8), parameter :: Gamma = 1.0_r8/6.0_r8
      real(r8), parameter :: eps = 1.0E-16_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
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
!=======================================================================
!  Tracer equation(s).
!=======================================================================
!
!-----------------------------------------------------------------------
!  Compute intermediate tracer at n+1/2 time-step, t(i,j,k,3,itrc).
!-----------------------------------------------------------------------
!
!  Compute time rate of change of intermediate tracer due to
!  horizontal advection.
!
      T_LOOP1 :DO itrc=1,NT(ng)
        K_LOOP: DO k=1,N(ng)
!
!  Fourth-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,nstp,itrc)-                             &
     &                t(i-1,j,k,nstp,itrc)
              FX(i,j)=FX(i,j)*umask(i,j)
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,nstp,itrc)+                            &
     &                 t(i  ,j,k,nstp,itrc)-                            &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,nstp,itrc)-                             &
     &                t(i,j-1,k,nstp,itrc)
              FE(i,j)=FE(i,j)*vmask(i,j)
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,nstp,itrc)+                            &
     &                 t(i,j  ,k,nstp,itrc)-                            &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
            END DO
          END DO
!
!  Time-step horizontal advection (m Tunits).
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.5_r8*dt(ng)
            cff1=1.0_r8
            cff2=0.0_r8
          ELSE
            cff=(1.0_r8-Gamma)*dt(ng)
            cff1=0.5_r8+Gamma
            cff2=0.5_r8-Gamma
          END IF
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j,k,3,itrc)=Hz(i,j,k)*(cff1*t(i,j,k,nstp,itrc)+       &
     &                                   cff2*t(i,j,k,nnew,itrc))-      &
     &                        cff*pm(i,j)*pn(i,j)*                      &
     &                        (FX(i+1,j)-FX(i,j)+                       &
     &                         FE(i,j+1)-FE(i,j))
            END DO
          END DO
        END DO K_LOOP
      END DO T_LOOP1
!
!  Compute artificial continuity equation (same for all tracers) and
!  load it into private array DC (1/m). Notice pipelined J-loop.
!
      J_LOOP1 : DO j=Jstr,Jend
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff=0.5_r8*dt(ng)
        ELSE
          cff=(1.0_r8-Gamma)*dt(ng)
        END IF
        DO k=1,N(ng)
          DO i=Istr,Iend
            DC(i,k)=1.0_r8/(Hz(i,j,k)-                                  &
     &                      cff*pm(i,j)*pn(i,j)*                        &
     &                      (Huon(i+1,j,k)-Huon(i,j,k)+                 &
     &                       Hvom(i,j+1,k)-Hvom(i,j,k)+                 &
     &                      (W(i,j,k)-W(i,j,k-1))))
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Compute time rate of change of intermediate tracer due to vertical
!  advection.  Impose artificial continuity equation.
!-----------------------------------------------------------------------
!
        T_LOOP2: DO itrc=1,NT(ng)
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,nstp,itrc)+                      &
     &                       t(i,j,k+1,nstp,itrc))-                     &
     &                 cff3*(t(i,j,k-1,nstp,itrc)+                      &
     &                       t(i,j,k+2,nstp,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,nstp,itrc)+                           &
     &               cff2*t(i,j,2,nstp,itrc)-                           &
     &               cff3*t(i,j,3,nstp,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,nstp,itrc)+               &
     &                     cff2*t(i,j,N(ng)-1,nstp,itrc)-               &
     &                     cff3*t(i,j,N(ng)-2,nstp,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
!
! Time-step vertical advection of tracers (Tunits).
!
          IF (iic(ng).eq.ntfirst(ng)) THEN
            cff=0.5_r8*dt(ng)
          ELSE
            cff=(1.0_r8-Gamma)*dt(ng)
          END IF
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=cff*pm(i,j)*pn(i,j)
              t(i,j,k,3,itrc)=DC(i,k)*                                  &
     &                        (t(i,j,k,3,itrc)-                         &
     &                         cff1*(FC(i,k)-FC(i,k-1)))
            END DO
          END DO
        END DO T_LOOP2
      END DO J_LOOP1
!
!-----------------------------------------------------------------------
!  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
!-----------------------------------------------------------------------
!
!  Compute vertical diffusive fluxes "FC" of the tracer fields at
!  current time step n, and at horizontal RHO-points and vertical
!  W-points.  Notice that the vertical diffusion coefficients for
!  passive tracers is the same as that for salinity (ltrc=NAT).
!
      DO j=Jstr,Jend
        cff3=dt(ng)*(1.0_r8-lambda)
        DO itrc=1,NT(ng)
          ltrc=MIN(NAT,itrc)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FC(i,k)=cff3*cff*Akt(i,j,k,ltrc)*                         &
     &                (t(i,j,k+1,nstp,itrc)-                            &
     &                 t(i,j,k  ,nstp,itrc))
            END DO
          END DO
!
!  Apply bottom and surface tracer flux conditions.
!
          DO i=Istr,Iend
            FC(i,0)=dt(ng)*btflx(i,j,itrc)
            FC(i,N(ng))=dt(ng)*stflx(i,j,itrc)
          END DO
!
!  Compute new tracer field (m Tunits).
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=Hz(i,j,k)*t(i,j,k,nstp,itrc)
              cff2=FC(i,k)-FC(i,k-1)
              t(i,j,k,nnew,itrc)=cff1+cff2
            END DO
          END DO
        END DO
      END DO
!
!=======================================================================
!  3D momentum equation in the XI-direction.
!=======================================================================
!
!  Compute U-component viscous vertical momentum fluxes "FC" at
!  current time-step n, and at horizontal U-points and vertical
!  W-points.
!
      J_LOOP2: DO j=Jstr,Jend
        cff3=dt(ng)*(1.0_r8-lambda)
        DO k=1,N(ng)-1
          DO i=IstrU,Iend
            cff=1.0_r8/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-                    &
     &                  z_r(i,j,k  )-z_r(i-1,j,k  ))
            FC(i,k)=cff3*cff*(u(i,j,k+1,nstp)-u(i,j,k,nstp))*           &
     &              (Akv(i,j,k)+Akv(i-1,j,k))
          END DO
        END DO
!
!  Apply bottom and surface stresses, if so is prescribed.
!
        DO i=IstrU,Iend
          FC(i,0)=dt(ng)*bustr(i,j)
          FC(i,N(ng))=dt(ng)*sustr(i,j)
        END DO
!
!  Compute new U-momentum (m m/s).
!
        cff=dt(ng)*0.25_r8
        DO i=IstrU,Iend
          DC(i,0)=cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
        END DO
        indx=3-nrhs
        IF (iic(ng).eq.ntfirst(ng)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff1=u(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              cff2=FC(i,k)-FC(i,k-1)
              u(i,j,k,nnew)=cff1+cff2
            END DO
          END DO
        ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff1=u(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              cff2=FC(i,k)-FC(i,k-1)
              cff3=0.5_r8*DC(i,0)
              u(i,j,k,nnew)=cff1-                                       &
     &                      cff3*ru(i,j,k,indx)+                        &
     &                      cff2
            END DO
          END DO
        ELSE
          cff1= 5.0_r8/12.0_r8
          cff2=16.0_r8/12.0_r8
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff3=u(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              cff4=FC(i,k)-FC(i,k-1)
              u(i,j,k,nnew)=cff3+                                       &
     &                      DC(i,0)*(cff1*ru(i,j,k,nrhs)-               &
     &                               cff2*ru(i,j,k,indx))+              &
     &                      cff4
            END DO
          END DO
        END IF
!
!=======================================================================
!  3D momentum equation in the ETA-direction.
!=======================================================================
!
!  Compute V-component viscous vertical momentum fluxes "FC" at
!  current time-step n, and at horizontal V-points and vertical
!  W-points.
!
        IF (j.ge.JstrV) THEN
          cff3=dt(ng)*(1.0_r8-lambda)
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-                  &
     &                    z_r(i,j,k  )-z_r(i,j-1,k  ))
              FC(i,k)=cff3*cff*(v(i,j,k+1,nstp)-v(i,j,k,nstp))*         &
     &                (Akv(i,j,k)+Akv(i,j-1,k))
            END DO
          END DO
!
!  Apply bottom and surface stresses, if so is prescribed.
!
          DO i=Istr,Iend
            FC(i,0)=dt(ng)*bvstr(i,j)
            FC(i,N(ng))=dt(ng)*svstr(i,j)
          END DO
!
!  Compute new V-momentum (m m/s).
!
          cff=dt(ng)*0.25_r8
          DO i=Istr,Iend
            DC(i,0)=cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
          END DO
          IF (iic(ng).eq.ntfirst(ng)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=v(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                cff2=FC(i,k)-FC(i,k-1)
                v(i,j,k,nnew)=cff1+cff2
              END DO
            END DO
          ELSE IF (iic(ng).eq.(ntfirst(ng)+1)) THEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=v(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                cff2=FC(i,k)-FC(i,k-1)
                cff3=0.5_r8*DC(i,0)
                v(i,j,k,nnew)=cff1-                                     &
     &                        cff3*rv(i,j,k,indx)+                      &
     &                        cff2
              END DO
            END DO
          ELSE
            cff1= 5.0_r8/12.0_r8
            cff2=16.0_r8/12.0_r8
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff3=v(i,j,k,nstp)*0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                cff4=FC(i,k)-FC(i,k-1)
                v(i,j,k,nnew)=cff3+                                     &
     &                        DC(i,0)*(cff1*rv(i,j,k,nrhs)-             &
     &                                 cff2*rv(i,j,k,indx))+            &
     &                        cff4
              END DO
            END DO
          END IF
        END IF
      END DO J_LOOP2
!
!=======================================================================
!  Apply intermediate tracers lateral boundary conditions.
!=======================================================================
!
      DO itrc=1,NT(ng)
        CALL t3dbc_tile (ng, tile, itrc,                                &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, 3,                                       &
     &                   t)
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,3,itrc))
        END IF
      END DO
      RETURN
      END SUBROUTINE pre_step3d_tile
      END MODULE pre_step3d_mod
