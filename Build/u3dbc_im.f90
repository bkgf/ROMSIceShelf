      MODULE u3dbc_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  U-velocity.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: u3dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE u3dbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
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
      CALL u3dbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 OCEAN(ng) % u)
      RETURN
      END SUBROUTINE u3dbc
!
!***********************************************************************
      SUBROUTINE u3dbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       u)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
!
!  Local variable declarations.
!
      integer :: Imin, Imax
      integer :: i, j, k
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dUde, dUdt, dUdx, tau
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
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Istr  ,j)=u(Istr  ,j  ,k,nstp)-                      &
     &                       u(Istr  ,j-1,k,nstp)
              grad(Istr+1,j)=u(Istr+1,j  ,k,nstp)-                      &
     &                       u(Istr+1,j-1,k,nstp)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dUdt=u(Istr+1,j,k,nstp)-u(Istr+1,j,k,nout)
                dUdx=u(Istr+1,j,k,nout)-u(Istr+2,j,k,nout)
                IF (LBC(iwest,isUvel,ng)%nudging) THEN
                  IF ((dUdt*dUdx).lt.0.0_r8) THEN
                    tau=M3obc_in(ng,iwest)
                  ELSE
                    tau=M3obc_out(ng,iwest)
                  END IF
                  tau=tau*dt(ng)
                END IF
                IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(Istr+1,j  )+                            &
     &                     grad(Istr+1,j+1))).gt.0.0_r8) THEN
                  dUde=grad(Istr+1,j  )
                ELSE
                  dUde=grad(Istr+1,j+1)
                END IF
                cff=MAX(dUdx*dUdx+dUde*dUde,eps)
                Cx=dUdt*dUdx
                Ce=0.0_r8
                u(Istr,j,k,nout)=(cff*u(Istr  ,j,k,nstp)+               &
     &                            Cx *u(Istr+1,j,k,nout)-               &
     &                            MAX(Ce,0.0_r8)*grad(Istr,j  )-        &
     &                            MIN(Ce,0.0_r8)*grad(Istr,j+1))/       &
     &                           (cff+Cx)
                IF (LBC(iwest,isUvel,ng)%nudging) THEN
                  u(Istr,j,k,nout)=u(Istr,j,k,nout)+                    &
     &                             tau*(BOUNDARY(ng)%u_west(j,k)-       &
     &                                  u(Istr,j,k,nstp))
                END IF
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                u(Istr,j,k,nout)=BOUNDARY(ng)%u_west(j,k)
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                u(Istr,j,k,nout)=u(Istr+1,j,k,nout)
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isUvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                u(Istr,j,k,nout)=0.0_r8
              END IF
            END DO
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
        IF (LBC(ieast,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Iend  ,j)=u(Iend  ,j  ,k,nstp)-                      &
     &                       u(Iend  ,j-1,k,nstp)
              grad(Iend+1,j)=u(Iend+1,j  ,k,nstp)-                      &
     &                       u(Iend+1,j-1,k,nstp)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dUdt=u(Iend,j,k,nstp)-u(Iend  ,j,k,nout)
                dUdx=u(Iend,j,k,nout)-u(Iend-1,j,k,nout)
                IF (LBC(ieast,isUvel,ng)%nudging) THEN
                  IF ((dUdt*dUdx).lt.0.0_r8) THEN
                    tau=M3obc_in(ng,ieast)
                  ELSE
                    tau=M3obc_out(ng,ieast)
                  END IF
                  tau=tau*dt(ng)
                END IF
                IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(Iend,j  )+                              &
     &                     grad(Iend,j+1))).gt.0.0_r8) THEN
                  dUde=grad(Iend,j  )
                ELSE
                  dUde=grad(Iend,j+1)
                END IF
                cff=MAX(dUdx*dUdx+dUde*dUde,eps)
                Cx=dUdt*dUdx
                Ce=0.0_r8
                u(Iend+1,j,k,nout)=(cff*u(Iend+1,j,k,nstp)+             &
     &                              Cx *u(Iend  ,j,k,nout)-             &
     &                              MAX(Ce,0.0_r8)*grad(Iend+1,j  )-    &
     &                              MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/   &
     &                             (cff+Cx)
                IF (LBC(ieast,isUvel,ng)%nudging) THEN
                  u(Iend+1,j,k,nout)=u(Iend+1,j,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%u_east(j,k)-     &
     &                                    u(Iend+1,j,k,nstp))
                END IF
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                u(Iend+1,j,k,nout)=BOUNDARY(ng)%u_east(j,k)
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                u(Iend+1,j,k,nout)=u(Iend,j,k,nout)
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isUvel,ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                u(Iend+1,j,k,nout)=0.0_r8
              END IF
            END DO
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
        IF (LBC(isouth,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
              grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-                      &
     &                       u(i  ,Jstr-1,k,nstp)
              grad(i,Jstr  )=u(i+1,Jstr  ,k,nstp)-                      &
     &                       u(i  ,Jstr  ,k,nstp)
            END DO
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dUdt=u(i,Jstr,k,nstp)-u(i,Jstr  ,k,nout)
                dUde=u(i,Jstr,k,nout)-u(i,Jstr+1,k,nout)
                IF (LBC(isouth,isUvel,ng)%nudging) THEN
                  IF ((dUdt*dUde).lt.0.0_r8) THEN
                    tau=M3obc_in(ng,isouth)
                  ELSE
                    tau=M3obc_out(ng,isouth)
                  END IF
                  tau=tau*dt(ng)
                END IF
                IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(i-1,Jstr)+                              &
     &                     grad(i  ,Jstr))).gt.0.0_r8) THEN
                  dUdx=grad(i-1,Jstr)
                ELSE
                  dUdx=grad(i  ,Jstr)
                END IF
                cff=MAX(dUdx*dUdx+dUde*dUde,eps)
                Cx=0.0_r8
                Ce=dUdt*dUde
                u(i,Jstr-1,k,nout)=(cff*u(i,Jstr-1,k,nstp)+             &
     &                              Ce *u(i,Jstr  ,k,nout)-             &
     &                              MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/   &
     &                             (cff+Ce)
                IF (LBC(isouth,isUvel,ng)%nudging) THEN
                  u(i,Jstr-1,k,nout)=u(i,Jstr-1,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%u_south(i,k)-    &
     &                                    u(i,Jstr-1,k,nstp))
                END IF
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                u(i,Jstr-1,k,nout)=BOUNDARY(ng)%u_south(i,k)
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                u(i,Jstr-1,k,nout)=u(i,Jstr,k,nout)
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(isouth,isUvel,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO k=1,N(ng)
            DO i=Imin,Imax
              IF (LBC_apply(ng)%south(i)) THEN
                u(i,Jstr-1,k,nout)=gamma2(ng)*u(i,Jstr,k,nout)
              END IF
            END DO
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
        IF (LBC(inorth,isUvel,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=IstrU-1,Iend
              grad(i,Jend  )=u(i+1,Jend  ,k,nstp)-                      &
     &                       u(i  ,Jend  ,k,nstp)
              grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-                      &
     &                       u(i  ,Jend+1,k,nstp)
            END DO
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dUdt=u(i,Jend,k,nstp)-u(i,Jend  ,k,nout)
                dUde=u(i,Jend,k,nout)-u(i,Jend-1,k,nout)
                IF (LBC(inorth,isUvel,ng)%nudging) THEN
                  IF ((dUdt*dUde).lt.0.0_r8) THEN
                    tau=M3obc_in(ng,inorth)
                  ELSE
                    tau=M3obc_out(ng,inorth)
                  END IF
                  tau=tau*dt(ng)
                END IF
                IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
                IF ((dUdt*(grad(i-1,Jend)+                              &
     &                     grad(i  ,Jend))).gt.0.0_r8) THEN
                  dUdx=grad(i-1,Jend)
                ELSE
                  dUdx=grad(i  ,Jend)
                END IF
                cff=MAX(dUdx*dUdx+dUde*dUde,eps)
                Cx=0.0_r8
                Ce=dUdt*dUde
                u(i,Jend+1,k,nout)=(cff*u(i,Jend+1,k,nstp)+             &
     &                              Ce *u(i,Jend  ,k,nout)-             &
     &                              MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-    &
     &                              MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/   &
     &                             (cff+Ce)
                IF (LBC(inorth,isUvel,ng)%nudging) THEN
                  u(i,Jend+1,k,nout)=u(i,Jend+1,k,nout)+                &
     &                               tau*(BOUNDARY(ng)%u_north(i,k)-    &
     &                                    u(i,Jend+1,k,nstp))
                END IF
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isUvel,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                u(i,Jend+1,k,nout)=BOUNDARY(ng)%u_north(i,k)
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isUvel,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=IstrU,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                u(i,Jend+1,k,nout)=u(i,Jend,k,nout)
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
        ELSE IF (LBC(inorth,isUvel,ng)%closed) THEN
          IF (EWperiodic(ng)) THEN
            Imin=IstrU
            Imax=Iend
          ELSE
            Imin=Istr
            Imax=IendR
          END IF
          DO k=1,N(ng)
            DO i=Imin,Imax
              IF (LBC_apply(ng)%north(i)) THEN
                u(i,Jend+1,k,nout)=gamma2(ng)*u(i,Jend,k,nout)
              END IF
            END DO
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
            DO k=1,N(ng)
              u(Istr,Jstr-1,k,nout)=0.5_r8*(u(Istr+1,Jstr-1,k,nout)+    &
     &                                      u(Istr  ,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=1,N(ng)
              u(Iend+1,Jstr-1,k,nout)=0.5_r8*(u(Iend  ,Jstr-1,k,nout)+  &
     &                                        u(Iend+1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr  ).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
              u(Istr,Jend+1,k,nout)=0.5_r8*(u(Istr  ,Jend  ,k,nout)+    &
     &                                      u(Istr+1,Jend+1,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
              u(Iend+1,Jend+1,k,nout)=0.5_r8*(u(Iend+1,Jend  ,k,nout)+  &
     &                                        u(Iend  ,Jend+1,k,nout))
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE u3dbc_tile
      END MODULE u3dbc_mod
