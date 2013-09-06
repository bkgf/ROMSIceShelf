      MODULE u2dbc_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated U-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: u2dbc, u2dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE u2dbc (ng, tile, kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout
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
      CALL u2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 krhs(ng), kstp(ng), kout,                        &
     &                 OCEAN(ng) % ubar,                                &
     &                 OCEAN(ng) % vbar,                                &
     &                 OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE u2dbc
!
!***********************************************************************
      SUBROUTINE u2dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       krhs, kstp, kout,                          &
     &                       ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, kout
!
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, know
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8) :: Ce, Cx, Zx
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8) :: cff, cff1, cff2, cff3, dt2d, dUde, dUdt, dUdx, tau
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
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
      IF (iif(ng).eq.1) THEN
        know=krhs
        dt2d=dtfast(ng)
      ELSE IF (PREDICTOR_2D_STEP(ng)) THEN
        know=krhs
        dt2d=2.0_r8*dtfast(ng)
      ELSE
        know=kstp
        dt2d=dtfast(ng)
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isUbar,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Istr  ,j)=ubar(Istr  ,j  ,know)-                       &
     &                     ubar(Istr  ,j-1,know)
            grad(Istr+1,j)=ubar(Istr+1,j  ,know)-                       &
     &                     ubar(Istr+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dUdt=ubar(Istr+1,j,know)-ubar(Istr+1,j,kout)
              dUdx=ubar(Istr+1,j,kout)-ubar(Istr+2,j,kout)
              IF (LBC(iwest,isUbar,ng)%nudging) THEN
                IF ((dUdt*dUdx).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,iwest)
                ELSE
                  tau=M2obc_out(ng,iwest)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Istr+1,j  )+                              &
     &                   grad(Istr+1,j+1))).gt.0.0_r8) THEN
                dUde=grad(Istr+1,j  )
              ELSE
                dUde=grad(Istr+1,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=dUdt*dUdx
              Ce=0.0_r8
              ubar(Istr,j,kout)=(cff*ubar(Istr  ,j,know)+               &
     &                           Cx *ubar(Istr+1,j,kout)-               &
     &                           MAX(Ce,0.0_r8)*grad(Istr,j  )-         &
     &                           MIN(Ce,0.0_r8)*grad(Istr,j+1))/        &
     &                          (cff+Cx)
              IF (LBC(iwest,isUbar,ng)%nudging) THEN
                ubar(Istr,j,kout)=ubar(Istr,j,kout)+                    &
     &                            tau*(BOUNDARY(ng)%ubar_west(j)-       &
     &                                 ubar(Istr,j,know))
              END IF
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
            END IF
          END DO
!
!  Western edge, Flather boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%Flather) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              bry_val=BOUNDARY(ng)%ubar_west(j)
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)-                     &
     &                        ABS(GRID(ng)%zice(Istr-1,j))+             &
     &                        GRID(ng)%h(Istr  ,j)-                     &
     &                        ABS(GRID(ng)%zice(Istr  ,j))))
              Cx=SQRT(g*cff)
              ubar(Istr,j,kout)=bry_val-                                &
     &                          Cx*(0.5_r8*(zeta(Istr-1,j,know)+        &
     &                                      zeta(Istr  ,j,know))-       &
     &                              BOUNDARY(ng)%zeta_west(j))
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
            END IF
          END DO
!
!  Western edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (LBC(iwest,isUbar,ng)%Shchepetkin) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              bry_val=BOUNDARY(ng)%ubar_west(j)
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)+                 &
     &                            GRID(ng)%h(Istr  ,j)))
              cff1=SQRT(g*cff)
              Cx=dt2d*cff1*0.5_r8*(GRID(ng)%pm(Istr-1,j)+               &
     &                             GRID(ng)%pm(Istr  ,j))
              Zx=(0.5_r8+Cx)*zeta(Istr  ,j,know)+                       &
     &           (0.5_r8-Cx)*zeta(Istr-1,j,know)
              IF (Cx.gt.Co) THEN
                cff2=(1.0_r8-Co/Cx)**2
                cff3=zeta(Istr,j,kout)+                                 &
     &               Cx*zeta(Istr-1,j,know)-                            &
     &               (1.0_r8-Cx)*zeta(Istr,j,know)
                Zx=Zx+cff2*cff3
              END IF
              ubar(Istr,j,kout)=0.5_r8*                                 &
     &                          ((1.0_r8-Cx)*ubar(Istr,j,know)+         &
     &                           Cx*ubar(Istr+1,j,know)+                &
     &                           bry_val-                               &
     &                           cff1*(Zx-BOUNDARY(ng)%zeta_west(j)))
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar(Istr,j,kout)=BOUNDARY(ng)%ubar_west(j)
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar(Istr,j,kout)=ubar(Istr+1,j,kout)
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
            END IF
          END DO
!
!  Western edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%reduced) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              IF (LBC(iwest,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(Istr,j,know)-                          &
     &                      BOUNDARY(ng)%zeta_west(j))*                 &
     &                  0.5_r8*GRID(ng)%pm(Istr,j)
              ELSE
                bry_pgr=-g*(zeta(Istr  ,j,know)-                        &
     &                      zeta(Istr-1,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Istr-1,j)+                  &
     &                          GRID(ng)%pm(Istr  ,j))
              END IF
              bry_cor=0.0_r8
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Istr-1,j)-                     &
     &                        ABS(GRID(ng)%zice(Istr-1,j))+             &
     &                        zeta(Istr-1,j,kout)+                      &
     &                        GRID(ng)%h(Istr  ,j)-                     &
     &                        ABS(GRID(ng)%zice(Istr  ,j))+             &
     &                        zeta(Istr  ,j,kout)))
              bry_str=cff*(FORCES(ng)%sustr(Istr,j)-                    &
     &                     FORCES(ng)%bustr(Istr,j))
              ubar(Istr,j,kout)=ubar(Istr,j,know)+                      &
     &                          dt2d*(bry_pgr+                          &
     &                                bry_cor+                          &
     &                                bry_str)
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*                      &
     &                          GRID(ng)%umask(Istr,j)
            END IF
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isUbar,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              ubar(Istr,j,kout)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isUbar,ng)%radiation) THEN
          DO j=Jstr,Jend+1
            grad(Iend  ,j)=ubar(Iend  ,j  ,know)-                       &
     &                     ubar(Iend  ,j-1,know)
            grad(Iend+1,j)=ubar(Iend+1,j  ,know)-                       &
     &                     ubar(Iend+1,j-1,know)
          END DO
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dUdt=ubar(Iend,j,know)-ubar(Iend  ,j,kout)
              dUdx=ubar(Iend,j,kout)-ubar(Iend-1,j,kout)
              IF (LBC(ieast,isUbar,ng)%nudging) THEN
                IF ((dUdt*dUdx).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,ieast)
                ELSE
                  tau=M2obc_out(ng,ieast)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(Iend,j  )+                                &
     &                   grad(Iend,j+1))).gt.0.0_r8) THEN
                dUde=grad(Iend,j)
              ELSE
                dUde=grad(Iend,j+1)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=dUdt*dUdx
              Ce=0.0_r8
              ubar(Iend+1,j,kout)=(cff*ubar(Iend+1,j,know)+             &
     &                             Cx *ubar(Iend  ,j,kout)-             &
     &                             MAX(Ce,0.0_r8)*grad(Iend+1,j  )-     &
     &                             MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/    &
     &                            (cff+Cx)
              IF (LBC(ieast,isUbar,ng)%nudging) THEN
                ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%ubar_east(j)-     &
     &                                   ubar(Iend+1,j,know))
              END IF
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, Flather boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%Flather) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              bry_val=BOUNDARY(ng)%ubar_east(j)
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)-                     &
     &                        ABS(GRID(ng)%zice(Iend  ,j))+             &
     &                        GRID(ng)%h(Iend+1,j)-                     &
     &                        ABS(GRID(ng)%zice(Iend+1,j))))
              Cx=SQRT(g*cff)
              ubar(Iend+1,j,kout)=bry_val+                              &
     &                            Cx*(0.5_r8*(zeta(Iend  ,j,know)+      &
     &                                        zeta(Iend+1,j,know))-     &
     &                                BOUNDARY(ng)%zeta_east(j))
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (LBC(ieast,isUbar,ng)%Shchepetkin) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              bry_val=BOUNDARY(ng)%ubar_east(j)
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)+                 &
     &                            GRID(ng)%h(Iend+1,j)))
              cff1=SQRT(g*cff)
              Cx=dt2d*cff1*0.5_r8*(GRID(ng)%pm(Iend  ,j)+               &
     &                             GRID(ng)%pm(Iend+1,j))
              Zx=(0.5_r8+Cx)*zeta(Iend  ,j,know)+                       &
     &           (0.5_r8-Cx)*zeta(Iend+1,j,know)
              IF (Cx.gt.Co) THEN
                cff2=(1.0_r8-Co/Cx)**2
                cff3=zeta(Iend,j,kout)+                                 &
     &               Cx*zeta(Iend+1,j,know)-                            &
     &               (1.0_r8-Cx)*zeta(Iend,j,know)
                Zx=Zx+cff2*cff3
              END IF
              ubar(Iend+1,j,kout)=0.5_r8*                               &
     &                            ((1.0_r8-Cx)*ubar(Iend+1,j,know)+     &
     &                             Cx*ubar(Iend,j,know)+                &
     &                             bry_val+                             &
     &                             cff1*(Zx-BOUNDARY(ng)%zeta_east(j)))
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%clamped) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar(Iend+1,j,kout)=BOUNDARY(ng)%ubar_east(j)
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%gradient) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar(Iend+1,j,kout)=ubar(Iend,j,kout)
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%reduced) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              IF (LBC(ieast,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_east(j)-                  &
     &                      zeta(Iend,j,know))*                         &
     &                  0.5_r8*GRID(ng)%pm(Iend,j)
              ELSE
                bry_pgr=-g*(zeta(Iend+1,j,know)-                        &
     &                      zeta(Iend  ,j,know))*                       &
     &                  0.5_r8*(GRID(ng)%pm(Iend  ,j)+                  &
     &                          GRID(ng)%pm(Iend+1,j))
              END IF
              bry_cor=0.0_r8
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(Iend  ,j)-                     &
     &                        ABS(GRID(ng)%zice(Iend  ,j))+             &
     &                        zeta(Iend  ,j,kout)+                      &
     &                        GRID(ng)%h(Iend+1,j)-                     &
     &                        ABS(GRID(ng)%zice(Iend+1,j))+             &
     &                        zeta(Iend+1,j,kout)))
              bry_str=cff*(FORCES(ng)%sustr(Iend+1,j)-                  &
     &                     FORCES(ng)%bustr(Iend+1,j))
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,know)+                  &
     &                            dt2d*(bry_pgr+                        &
     &                                  bry_cor+                        &
     &                                  bry_str)
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%umask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isUbar,ng)%closed) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              ubar(Iend+1,j,kout)=0.0_r8
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isUbar,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jstr-1)=ubar(i+1,Jstr-1,know)-                       &
     &                     ubar(i  ,Jstr-1,know)
            grad(i,Jstr  )=ubar(i+1,Jstr  ,know)-                       &
     &                     ubar(i  ,Jstr  ,know)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dUdt=ubar(i,Jstr,know)-ubar(i,Jstr  ,kout)
              dUde=ubar(i,Jstr,kout)-ubar(i,Jstr+1,kout)
              IF (LBC(isouth,isUbar,ng)%nudging) THEN
                IF ((dUdt*dUde).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,isouth)
                ELSE
                  tau=M2obc_out(ng,isouth)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jstr)+                                &
     &                   grad(i  ,Jstr))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jstr)
              ELSE
                dUdx=grad(i  ,Jstr)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=0.0_r8
              Ce=dUdt*dUde
              ubar(i,Jstr-1,kout)=(cff*ubar(i,Jstr-1,know)+             &
     &                             Ce *ubar(i,Jstr  ,kout)-             &
     &                             MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/    &
     &                            (cff+Ce)
              IF (LBC(isouth,isUbar,ng)%nudging) THEN
                ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)+                &
     &                              tau*(BOUNDARY(ng)%ubar_south(i)-    &
     &                                   ubar(i,Jstr-1,know))
              END IF
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, Chapman boundary condition.
!
        ELSE IF (LBC(isouth,isUbar,ng)%Flather.or.                      &
     &           LBC(isouth,isUbar,ng)%reduced.or.                      &
     &           LBC(isouth,isUbar,ng)%Shchepetkin) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jstr)+                   &
     &                         GRID(ng)%pn(i  ,Jstr))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jstr)-                     &
     &                        ABS(GRID(ng)%zice(i-1,Jstr))+             &
     &                        GRID(ng)%h(i  ,Jstr)-                     &
     &                        ABS(GRID(ng)%zice(i  ,Jstr))))
              Ce=cff*cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              ubar(i,Jstr-1,kout)=cff2*(ubar(i,Jstr-1,know)+            &
     &                                  Ce*ubar(i,Jstr,kout))
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isUbar,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar(i,Jstr-1,kout)=BOUNDARY(ng)%ubar_south(i)
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUbar,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              ubar(i,Jstr-1,kout)=ubar(i,Jstr,kout)
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
            END IF
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isUbar,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%south(i)) THEN
              ubar(i,Jstr-1,kout)=gamma2(ng)*ubar(i,Jstr,kout)
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*                  &
     &                            GRID(ng)%umask(i,Jstr-1)
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isUbar,ng)%radiation) THEN
          DO i=IstrU-1,Iend
            grad(i,Jend  )=ubar(i+1,Jend  ,know)-                       &
     &                     ubar(i  ,Jend  ,know)
            grad(i,Jend+1)=ubar(i+1,Jend+1,know)-                       &
     &                     ubar(i  ,Jend+1,know)
          END DO
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dUdt=ubar(i,Jend,know)-ubar(i,Jend  ,kout)
              dUde=ubar(i,Jend,kout)-ubar(i,Jend-1,kout)
              IF (LBC(inorth,isUbar,ng)%nudging) THEN
                IF ((dUdt*dUde).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,inorth)
                ELSE
                  tau=M2obc_out(ng,inorth)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
              IF ((dUdt*(grad(i-1,Jend)+                                &
     &                   grad(i  ,Jend))).gt.0.0_r8) THEN
                dUdx=grad(i-1,Jend)
              ELSE
                dUdx=grad(i  ,Jend)
              END IF
              cff=MAX(dUdx*dUdx+dUde*dUde,eps)
              Cx=0.0_r8
              Ce=dUdt*dUde
              ubar(i,Jend+1,kout)=(cff*ubar(i,Jend+1,know)+             &
     &                             Ce *ubar(i,Jend  ,kout)-             &
     &                             MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/    &
     &                            (cff+Ce)
              IF (LBC(inorth,isUbar,ng)%nudging) THEN
                ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)+                &
     &                              tau*(BOUNDARY(ng)%ubar_north(i)-    &
     &                                   ubar(i,Jend+1,know))
              END IF
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, Chapman boundary condition.
!
        ELSE IF (LBC(inorth,isUbar,ng)%Flather.or.                      &
     &           LBC(inorth,isUbar,ng)%reduced.or.                      &
     &           LBC(inorth,isUbar,ng)%Shchepetkin) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pn(i-1,Jend)+                   &
     &                         GRID(ng)%pn(i  ,Jend))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(i-1,Jend)-                     &
     &                        ABS(GRID(ng)%zice(i-1,Jend))+             &
     &                        GRID(ng)%h(i  ,Jend)-                     &
     &                        ABS(GRID(ng)%zice(i  ,Jend))))
              Ce=cff*cff1
              cff2=1.0_r8/(1.0_r8+Ce)
              ubar(i,Jend+1,kout)=cff2*(ubar(i,Jend+1,know)+            &
     &                                  Ce*ubar(i,Jend,kout))
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isUbar,ng)%clamped) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar(i,Jend+1,kout)=BOUNDARY(ng)%ubar_north(i)
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUbar,ng)%gradient) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              ubar(i,Jend+1,kout)=ubar(i,Jend,kout)
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isUbar,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO i=Imin,Imax
            IF (LBC_apply(ng)%north(i)) THEN
              ubar(i,Jend+1,kout)=gamma2(ng)*ubar(i,Jend,kout)
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%umask(i,Jend+1)
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            ubar(Istr,Jstr-1,kout)=0.5_r8*(ubar(Istr+1,Jstr-1,kout)+    &
     &                                     ubar(Istr  ,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            ubar(Iend+1,Jstr-1,kout)=0.5_r8*(ubar(Iend  ,Jstr-1,kout)+  &
     &                                       ubar(Iend+1,Jstr  ,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            ubar(Istr,Jend+1,kout)=0.5_r8*(ubar(Istr  ,Jend  ,kout)+    &
     &                                     ubar(Istr+1,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            ubar(Iend+1,Jend+1,kout)=0.5_r8*(ubar(Iend+1,Jend  ,kout)+  &
     &                                       ubar(Iend  ,Jend+1,kout))
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Impose wetting and drying conditions.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,j,kout))*              &
     &                    GRID(ng)%umask_wet(Istr,j)
              cff=0.5_r8*GRID(ng)%umask_wet(Istr,j)*cff1+               &
     &            cff2*(1.0_r8-cff1)
              ubar(Istr,j,kout)=ubar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=Jstr,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,j,kout))*            &
     &                    GRID(ng)%umask_wet(Iend+1,j)
              cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              ubar(Iend+1,j,kout)=ubar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
      END IF
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=IstrU,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(i,Jstr-1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(i,Jstr-1,kout))*            &
     &                    GRID(ng)%umask_wet(i,Jstr-1)
              cff=0.5_r8*GRID(ng)%umask_wet(i,Jstr-1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              ubar(i,Jstr-1,kout)=ubar(i,Jstr-1,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff1=ABS(ABS(GRID(ng)%umask_wet(i,Jend+1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,ubar(i,Jend+1,kout))*            &
     &                    GRID(ng)%umask_wet(i,Jend+1)
              cff=0.5_r8*GRID(ng)%umask_wet(i,Jend+1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              ubar(i,Jend+1,kout)=ubar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
      END IF
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,Jstr-1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,Jstr-1,kout))*           &
     &                  GRID(ng)%umask_wet(Istr,Jstr-1)
            cff=0.5_r8*GRID(ng)%umask_wet(Istr,Jstr-1)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            ubar(Istr,Jstr-1,kout)=ubar(Istr,Jstr-1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,Jstr-1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,Jstr-1,kout))*         &
     &                  GRID(ng)%umask_wet(Iend+1,Jstr-1)
            cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,Jstr-1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            ubar(Iend+1,Jstr-1,kout)=ubar(Iend+1,Jstr-1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Istr,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Istr,Jend+1,kout))*           &
     &                  GRID(ng)%umask_wet(Istr,Jend+1)
            cff=0.5_r8*GRID(ng)%umask_wet(Istr,Jend+1)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            ubar(Istr,Jend+1,kout)=ubar(Istr,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%umask_wet(Iend+1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,ubar(Iend+1,Jend+1,kout))*         &
     &                  GRID(ng)%umask_wet(Iend+1,Jend+1)
            cff=0.5_r8*GRID(ng)%umask_wet(Iend+1,Jend+1)+cff1+          &
     &          cff2*(1.0_r8-cff1)
            ubar(Iend+1,Jend+1,kout)=ubar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE u2dbc_tile
      END MODULE u2dbc_mod
