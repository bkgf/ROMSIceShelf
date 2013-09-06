      MODULE v2dbc_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for vertically     !
!  integrated V-velocity.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: v2dbc, v2dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE v2dbc (ng, tile, kout)
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
      CALL v2dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 krhs(ng), kstp(ng), kout,                        &
     &                 OCEAN(ng) % ubar,                                &
     &                 OCEAN(ng) % vbar,                                &
     &                 OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE v2dbc
!
!***********************************************************************
      SUBROUTINE v2dbc_tile (ng, tile,                                  &
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
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Jmin, Jmax
      integer :: i, j, know
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8) :: Ce, Cx, Ze
      real(r8) :: bry_pgr, bry_cor, bry_str, bry_val
      real(r8):: cff, cff1, cff2, cff3, dt2d, dVde, dVdt, dVdx, tau
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
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isVbar,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jstr  )=vbar(i  ,Jstr  ,know)-                       &
     &                     vbar(i-1,Jstr  ,know)
            grad(i,Jstr+1)=vbar(i  ,Jstr+1,know)-                       &
     &                     vbar(i-1,Jstr+1,know)
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              dVdt=vbar(i,Jstr+1,know)-vbar(i,Jstr+1,kout)
              dVde=vbar(i,Jstr+1,kout)-vbar(i,Jstr+2,kout)
              IF (LBC(isouth,isVbar,ng)%nudging) THEN
                IF ((dVdt*dVde).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,isouth)
                ELSE
                  tau=M2obc_out(ng,isouth)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i  ,Jstr+1)+                              &
     &                   grad(i+1,Jstr+1))).gt.0.0_r8) THEN
                dVdx=grad(i  ,Jstr+1)
              ELSE
                dVdx=grad(i+1,Jstr+1)
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=0.0_r8
              Ce=dVdt*dVde
              vbar(i,Jstr,kout)=(cff*vbar(i,Jstr  ,know)+               &
     &                           Ce *vbar(i,Jstr+1,kout)-               &
     &                           MAX(Cx,0.0_r8)*grad(i  ,Jstr)-         &
     &                           MIN(Cx,0.0_r8)*grad(i+1,Jstr))/        &
     &                          (cff+Ce)
              IF (LBC(isouth,isVbar,ng)%nudging) THEN
                vbar(i,Jstr,kout)=vbar(i,Jstr,kout)+                    &
     &                            tau*(BOUNDARY(ng)%vbar_south(i)-      &
     &                                 vbar(i,Jstr,know))
              END IF
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
            END IF
          END DO
!
!  Southern edge, Flather boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              bry_val=BOUNDARY(ng)%vbar_south(i)
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)-                     &
     &                        ABS(GRID(ng)%zice(i,Jstr-1))+             &
     &                        GRID(ng)%h(i,Jstr  )-                     &
     &                        ABS(GRID(ng)%zice(i,Jstr  ))))
              Ce=SQRT(g*cff)
              vbar(i,Jstr,kout)=bry_val-                                &
     &                          Ce*(0.5_r8*(zeta(i,Jstr-1,know)+        &
     &                                      zeta(i,Jstr  ,know))-       &
     &                              BOUNDARY(ng)%zeta_south(i))
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
            END IF
          END DO
!
!  Southern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (LBC(isouth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              bry_val=BOUNDARY(ng)%vbar_south(i)
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)+                 &
     &                            GRID(ng)%h(i,Jstr  )))
              cff1=SQRT(g*cff)
              Ce=dt2d*cff1*0.5_r8*(GRID(ng)%pn(i,Jstr-1)+               &
     &                             GRID(ng)%pn(i,Jstr  ))
              Ze=(0.5_r8+Ce)*zeta(i,Jstr  ,know)+                       &
     &           (0.5_r8-Ce)*zeta(i,Jstr-1,know)
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                cff3=zeta(i,Jstr,kout)+                                 &
     &               Ce*zeta(i,Jstr-1,know)-                            &
     &               (1.0_r8-Ce)*zeta(i,Jstr,know)
                Ze=Ze+cff2*cff3
              END IF
              vbar(i,Jstr,kout)=0.5_r8*                                 &
     &                          ((1.0_r8-Ce)*vbar(i,Jstr,know)+         &
     &                           Ce*vbar(i,Jstr+1,know)+                &
     &                           bry_val-                               &
     &                           cff1*(Ze-BOUNDARY(ng)%zeta_south(i)))
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
            END IF
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              vbar(i,Jstr,kout)=BOUNDARY(ng)%vbar_south(i)
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
            END IF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              vbar(i,Jstr,kout)=vbar(i,Jstr+1,kout)
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
            END IF
          END DO
!
!  Southern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              IF (LBC(isouth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(zeta(i,Jstr,know)-                          &
     &                      BOUNDARY(ng)%zeta_south(i))*                &
     &                  0.5_r8*GRID(ng)%pn(i,Jstr)
              ELSE
                bry_pgr=-g*(zeta(i,Jstr  ,know)-                        &
     &                      zeta(i,Jstr-1,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jstr-1)+                  &
     &                          GRID(ng)%pn(i,Jstr  ))
              END IF
              bry_cor=0.0_r8
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jstr-1)-                     &
     &                        ABS(GRID(ng)%zice(i,Jstr-1))+             &
     &                        zeta(i,Jstr-1,kout)+                      &
     &                        GRID(ng)%h(i,Jstr  )-                     &
     &                        ABS(GRID(ng)%zice(i,Jstr  ))+             &
     &                        zeta(i,Jstr  ,kout)))
              bry_str=cff*(FORCES(ng)%svstr(i,Jstr)-                    &
     &                     FORCES(ng)%bvstr(i,Jstr))
              vbar(i,Jstr,kout)=vbar(i,Jstr,know)+                      &
     &                          dt2d*(bry_pgr+                          &
     &                                bry_cor+                          &
     &                                bry_str)
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*                      &
     &                          GRID(ng)%vmask(i,Jstr)
            END IF
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              vbar(i,Jstr,kout)=0.0_r8
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
        IF (LBC(inorth,isVbar,ng)%radiation) THEN
          DO i=Istr,Iend+1
            grad(i,Jend  )=vbar(i  ,Jend  ,know)-                       &
     &                     vbar(i-1,Jend  ,know)
            grad(i,Jend+1)=vbar(i  ,Jend+1,know)-                       &
     &                     vbar(i-1,Jend+1,know)
          END DO
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              dVdt=vbar(i,Jend,know)-vbar(i,Jend  ,kout)
              dVde=vbar(i,Jend,kout)-vbar(i,Jend-1,kout)
              IF (LBC(inorth,isVbar,ng)%nudging) THEN
                IF ((dVdt*dVde).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,inorth)
                ELSE
                  tau=M2obc_out(ng,inorth)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(i  ,Jend)+                                &
     &                   grad(i+1,Jend))).gt.0.0_r8) THEN
                dVdx=grad(i  ,Jend)
              ELSE
                dVdx=grad(i+1,Jend)
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=0.0_r8
              Ce=dVdt*dVde
              vbar(i,Jend+1,kout)=(cff*vbar(i,Jend+1,know)+             &
     &                             Ce *vbar(i,Jend  ,kout)-             &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-     &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/    &
     &                            (cff+Ce)
              IF (LBC(inorth,isVbar,ng)%nudging) THEN
                vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)+                &
     &                              tau*(BOUNDARY(ng)%vbar_north(i)-    &
     &                                   vbar(i,Jend+1,know))
              END IF
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, Flather boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%Flather) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              bry_val=BOUNDARY(ng)%vbar_north(i)
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )-                     &
     &                        ABS(GRID(ng)%zice(i,Jend  ))+             &
     &                        GRID(ng)%h(i,Jend+1)-                     &
     &                        ABS(GRID(ng)%zice(i,Jend+1))))
              Ce=SQRT(g*cff)
              vbar(i,Jend+1,kout)=bry_val+                              &
     &                            Ce*(0.5_r8*(zeta(i,Jend  ,know)+      &
     &                                        zeta(i,Jend+1,know))-     &
     &                                BOUNDARY(ng)%zeta_north(i))
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, Shchepetkin boundary condition (Maison et al., 2010).
!
        ELSE IF (LBC(inorth,isVbar,ng)%Shchepetkin) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              bry_val=BOUNDARY(ng)%vbar_north(i)
              cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )+                 &
     &                            GRID(ng)%h(i,Jend+1)))
              cff1=SQRT(g*cff)
              Ce=dt2d*cff1*0.5_r8*(GRID(ng)%pn(i,Jend  )+               &
     &                             GRID(ng)%pn(i,Jend+1))
              Ze=(0.5_r8+Ce)*zeta(i,Jend  ,know)+                       &
     &           (0.5_r8-Ce)*zeta(i,Jend+1,know)
              IF (Ce.gt.Co) THEN
                cff2=(1.0_r8-Co/Ce)**2
                cff3=zeta(i,Jend,kout)+                                 &
     &               Ce*zeta(i,Jend+1,know)-                            &
     &               (1.0_r8-Ce)*zeta(i,Jend,know)
                Ze=Ze+cff2*cff3
              END IF
              vbar(i,Jend+1,kout)=0.5_r8*                               &
     &                            ((1.0_r8-Ce)*vbar(i,Jend+1,know)+     &
     &                             Ce*vbar(i,Jend,know)+                &
     &                             bry_val+                             &
     &                             cff1*(Ze-BOUNDARY(ng)%zeta_north(i)))
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%clamped) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=BOUNDARY(ng)%vbar_north(i)
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%gradient) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=vbar(i,Jend,kout)
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, reduced-physics boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%reduced) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              IF (LBC(inorth,isFsur,ng)%acquire) THEN
                bry_pgr=-g*(BOUNDARY(ng)%zeta_north(i)-                 &
     &                      zeta(i,Jend,know))*                         &
     &                  0.5_r8*GRID(ng)%pn(i,Jend)
              ELSE
                bry_pgr=-g*(zeta(i,Jend+1,know)-                        &
     &                      zeta(i,Jend  ,know))*                       &
     &                  0.5_r8*(GRID(ng)%pn(i,Jend  )+                  &
     &                          GRID(ng)%pn(i,Jend+1))
              END IF
              bry_cor=0.0_r8
          cff=1.0_r8/(0.5_r8*(GRID(ng)%h(i,Jend  )-                     &
     &                        ABS(GRID(ng)%zice(i,Jend  ))+             &
     &                        zeta(i,Jend  ,kout)+                      &
     &                        GRID(ng)%h(i,Jend+1)-                     &
     &                        ABS(GRID(ng)%zice(i,Jend+1))+             &
     &                        zeta(i,Jend+1,kout)))
              bry_str=cff*(FORCES(ng)%svstr(i,Jend+1)-                  &
     &                     FORCES(ng)%bvstr(i,Jend+1))
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,know)+                  &
     &                            dt2d*(bry_pgr+                        &
     &                                  bry_cor+                        &
     &                                  bry_str)
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*                  &
     &                            GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isVbar,ng)%closed) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              vbar(i,Jend+1,kout)=0.0_r8
            END IF
          END DO
        END IF
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
        IF (LBC(iwest,isVbar,ng)%radiation) THEN
          DO j=JstrV-1,Jend
            grad(Istr-1,j)=vbar(Istr-1,j+1,know)-                       &
     &                     vbar(Istr-1,j  ,know)
            grad(Istr  ,j)=vbar(Istr  ,j+1,know)-                       &
     &                     vbar(Istr  ,j  ,know)
          END DO
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              dVdt=vbar(Istr,j,know)-vbar(Istr  ,j,kout)
              dVdx=vbar(Istr,j,kout)-vbar(Istr+1,j,kout)
              IF (LBC(iwest,isVbar,ng)%nudging) THEN
                IF ((dVdt*dVdx).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,iwest)
                ELSE
                  tau=M2obc_out(ng,iwest)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Istr,j-1)+                                &
     &                   grad(Istr,j  ))).gt.0.0_r8) THEN
                dVde=grad(Istr,j-1)
              ELSE
                dVde=grad(Istr,j  )
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=dVdt*dVdx
              Ce=0.0_r8
              vbar(Istr-1,j,kout)=(cff*vbar(Istr-1,j,know)+             &
     &                             Cx *vbar(Istr  ,j,kout)-             &
     &                             MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-     &
     &                             MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/    &
     &                            (cff+Cx)
              IF (LBC(iwest,isVbar,ng)%nudging) THEN
                vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%vbar_west(j)-     &
     &                                   vbar(Istr-1,j,know))
              END IF
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
            END IF
          END DO
!
!  Western edge, Chapman boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%Flather.or.                       &
     &           LBC(iwest,isVbar,ng)%reduced.or.                       &
     &           LBC(iwest,isVbar,ng)%Shchepetkin) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Istr,j-1)+                   &
     &                         GRID(ng)%pm(Istr,j  ))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Istr,j-1)-                     &
     &                        ABS(GRID(ng)%zice(Istr,j-1))+             &
     &                        GRID(ng)%h(Istr,j  )-                     &
     &                        ABS(GRID(ng)%zice(Istr,j  ))))
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              vbar(Istr-1,j,kout)=cff2*(vbar(Istr-1,j,know)+            &
     &                                  Cx*vbar(Istr,j,kout))
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
            END IF
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=BOUNDARY(ng)%vbar_west(j)
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
            END IF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=vbar(Istr,j,kout)
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
            END IF
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%west(j)) THEN
              vbar(Istr-1,j,kout)=gamma2(ng)*vbar(Istr,j,kout)
              vbar(Istr-1,j,kout)=vbar(Istr-1,j,kout)*                  &
     &                            GRID(ng)%vmask(Istr-1,j)
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
        IF (LBC(ieast,isVbar,ng)%radiation) THEN
          DO j=JstrV-1,Jend
            grad(Iend  ,j)=vbar(Iend  ,j+1,know)-                       &
     &                     vbar(Iend  ,j  ,know)
            grad(Iend+1,j)=vbar(Iend+1,j+1,know)-                       &
     &                     vbar(Iend+1,j  ,know)
          END DO
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              dVdt=vbar(Iend,j,know)-vbar(Iend  ,j,kout)
              dVdx=vbar(Iend,j,kout)-vbar(Iend-1,j,kout)
              IF (LBC(ieast,isVbar,ng)%nudging) THEN
                IF ((dVdt*dVdx).lt.0.0_r8) THEN
                  tau=M2obc_in(ng,ieast)
                ELSE
                  tau=M2obc_out(ng,ieast)
                END IF
                tau=tau*dt2d
              END IF
              IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
              IF ((dVdt*(grad(Iend,j-1)+                                &
     &                   grad(Iend,j  ))).gt.0.0_r8) THEN
                dVde=grad(Iend,j-1)
              ELSE
                dVde=grad(Iend,j  )
              END IF
              cff=MAX(dVdx*dVdx+dVde*dVde,eps)
              Cx=dVdt*dVdx
              Ce=0.0_r8
              vbar(Iend+1,j,kout)=(cff*vbar(Iend+1,j,know)+             &
     &                             Cx *vbar(Iend  ,j,kout)-             &
     &                             MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-     &
     &                             MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/    &
     &                            (cff+Cx)
              IF (LBC(ieast,isVbar,ng)%nudging) THEN
                vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)+                &
     &                              tau*(BOUNDARY(ng)%vbar_east(j)-     &
     &                                   vbar(Iend+1,j,know))
              END IF
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, Chapman boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%Flather.or.                       &
     &           LBC(ieast,isVbar,ng)%reduced.or.                       &
     &           LBC(ieast,isVbar,ng)%Shchepetkin) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff=dt2d*0.5_r8*(GRID(ng)%pm(Iend,j-1)+                   &
     &                         GRID(ng)%pm(Iend,j  ))
          cff1=SQRT(g*0.5_r8*(GRID(ng)%h(Iend,j-1)-                     &
     &                        ABS(GRID(ng)%zice(Iend,j-1))+             &
     &                        GRID(ng)%h(Iend,j  )-                     &
     &                        ABS(GRID(ng)%zice(Iend,j  ))))
              Cx=cff*cff1
              cff2=1.0_r8/(1.0_r8+Cx)
              vbar(Iend+1,j,kout)=cff2*(vbar(Iend+1,j,know)+            &
     &                                  Cx*vbar(Iend,j,kout))
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%clamped) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=BOUNDARY(ng)%vbar_east(j)
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isVbar,ng)%gradient) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=vbar(Iend,j,kout)
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
            END IF
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isVbar,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO j=Jmin,Jmax
            IF (LBC_apply(ng)%east(j)) THEN
              vbar(Iend+1,j,kout)=gamma2(ng)*vbar(Iend,j,kout)
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*                  &
     &                            GRID(ng)%vmask(Iend+1,j)
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
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            vbar(Istr-1,Jstr,kout)=0.5_r8*(vbar(Istr  ,Jstr  ,kout)+    &
     &                                     vbar(Istr-1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            vbar(Iend+1,Jstr,kout)=0.5_r8*(vbar(Iend  ,Jstr  ,kout)+    &
     &                                     vbar(Iend+1,Jstr+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            vbar(Istr-1,Jend+1,kout)=0.5_r8*(vbar(Istr-1,Jend  ,kout)+  &
     &                                       vbar(Istr  ,Jend+1,kout))
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            vbar(Iend+1,Jend+1,kout)=0.5_r8*(vbar(Iend+1,Jend  ,kout)+  &
     &                                      vbar(Iend  ,Jend+1,kout))
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
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%west(j)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,j,kout))*            &
     &                    GRID(ng)%vmask_wet(Istr-1,j)
              cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(Istr,j,kout)=vbar(Istr,j,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          DO j=JstrV,Jend
            IF (LBC_apply(ng)%east(j)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,j))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,j,kout))*            &
     &                    GRID(ng)%vmask_wet(Iend+1,j)
              cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,j)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(Iend+1,j,kout)=vbar(Iend+1,j,kout)*cff
            END IF
          END DO
        END IF
      END IF
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%south(i)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jstr))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jstr,kout))*              &
     &                    GRID(ng)%vmask_wet(i,Jstr)
              cff=0.5_r8*GRID(ng)%vmask_wet(i,Jstr)*cff1+               &
     &            cff2*(1.0_r8-cff1)
              vbar(i,Jstr,kout)=vbar(i,Jstr,kout)*cff
            END IF
          END DO
        END IF
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          DO i=Istr,Iend
            IF (LBC_apply(ng)%north(i)) THEN
              cff1=ABS(ABS(GRID(ng)%vmask_wet(i,Jend+1))-1.0_r8)
              cff2=0.5_r8+DSIGN(0.5_r8,vbar(i,Jend+1,kout))*            &
     &                    GRID(ng)%vmask_wet(i,Jend+1)
              cff=0.5_r8*GRID(ng)%vmask_wet(i,Jend+1)*cff1+             &
     &            cff2*(1.0_r8-cff1)
              vbar(i,Jend+1,kout)=vbar(i,Jend+1,kout)*cff
            END IF
          END DO
        END IF
      END IF
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr  )) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jstr))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jstr,kout))*           &
     &                  GRID(ng)%vmask_wet(Istr-1,Jstr)
            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jstr)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            vbar(Istr-1,Jstr,kout)=vbar(Istr-1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr  )) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jstr))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jstr,kout))*           &
     &                  GRID(ng)%vmask_wet(Iend+1,Jstr)
            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jstr)*cff1+            &
     &          cff2*(1.0_r8-cff1)
            vbar(Iend+1,Jstr,kout)=vbar(Iend+1,Jstr,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Istr-1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Istr-1,Jend+1,kout))*         &
     &                  GRID(ng)%vmask_wet(Istr-1,Jend+1)
            cff=0.5_r8*GRID(ng)%vmask_wet(Istr-1,Jend+1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            vbar(Istr-1,Jend+1,kout)=vbar(Istr-1,Jend+1,kout)*cff
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            cff1=ABS(ABS(GRID(ng)%vmask_wet(Iend+1,Jend+1))-1.0_r8)
            cff2=0.5_r8+DSIGN(0.5_r8,vbar(Iend+1,Jend+1,kout))*         &
     &                  GRID(ng)%vmask_wet(Iend+1,Jend+1)
            cff=0.5_r8*GRID(ng)%vmask_wet(Iend+1,Jend+1)*cff1+          &
     &          cff2*(1.0_r8-cff1)
            vbar(Iend+1,Jend+1,kout)=vbar(Iend+1,Jend+1,kout)*cff
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE v2dbc_tile
      END MODULE v2dbc_mod
