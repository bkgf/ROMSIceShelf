      MODULE t3dbc_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: t3dbc_tile
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE t3dbc (ng, tile, nout, itrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout, itrc
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
      CALL t3dbc_tile (ng, tile, itrc,                                  &
     &                 LBi, UBi, LBj, UBj, N(ng), NT(ng),               &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nstp(ng), nout,                                  &
     &                 OCEAN(ng)% t)
      RETURN
      END SUBROUTINE t3dbc
!
!***********************************************************************
      SUBROUTINE t3dbc_tile (ng, tile, itrc,                            &
     &                       LBi, UBi, LBj, UBj, UBk, UBt,              &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp, nout,                                &
     &                       t)
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
      integer, intent(in) :: ng, tile, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nout
!
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau
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
        IF (LBC(iwest,isTvar(itrc),ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Istr-1,j)=t(Istr-1,j  ,k,nstp,itrc)-                 &
     &                       t(Istr-1,j-1,k,nstp,itrc)
              grad(Istr  ,j)=t(Istr  ,j  ,k,nstp,itrc)-                 &
     &                       t(Istr  ,j-1,k,nstp,itrc)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dTdt=t(Istr,j,k,nstp,itrc)-t(Istr  ,j,k,nout,itrc)
                dTdx=t(Istr,j,k,nout,itrc)-t(Istr+1,j,k,nout,itrc)
                IF (LBC(iwest,isTvar(itrc),ng)%nudging) THEN
                  tau=Tobc_out(itrc,ng,iwest)
                  IF ((dTdt*dTdx).lt.0.0_r8) THEN
                    tau=Tobc_in(itrc,ng,iwest)
                  ELSE
                    tau=tau*dt(ng)
                  END IF
                END IF
                IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(Istr,j  )+                              &
     &                     grad(Istr,j+1))).gt.0.0_r8) THEN
                  dTde=grad(Istr,j  )
                ELSE
                  dTde=grad(Istr,j+1)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=dTdt*dTdx
                Ce=0.0_r8
                t(Istr-1,j,k,nout,itrc)=(cff*t(Istr-1,j,k,nstp,itrc)+   &
     &                                   Cx *t(Istr  ,j,k,nout,itrc)-   &
     &                                   MAX(Ce,0.0_r8)*                &
     &                                      grad(Istr-1,j  )-           &
     &                                   MIN(Ce,0.0_r8)*                &
     &                                      grad(Istr-1,j+1))/          &
     &                                  (cff+Cx)
                IF (LBC(iwest,isTvar(itrc),ng)%nudging) THEN
                  t(Istr-1,j,k,nout,itrc)=t(Istr-1,j,k,nout,itrc)+      &
     &                                    tau*                          &
     &                                    (BOUNDARY(ng)%t_west(j,k,     &
     &                                                         itrc)-   &
     &                                     t(Istr-1,j,k,nstp,itrc))
                END IF
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                t(Istr-1,j,k,nout,itrc)=BOUNDARY(ng)%t_west(j,k,itrc)
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                t(Istr-1,j,k,nout,itrc)=t(Istr,j,k,nout,itrc)
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
        IF (LBC(ieast,isTvar(itrc),ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend+1
              grad(Iend  ,j)=t(Iend  ,j  ,k,nstp,itrc)-                 &
     &                       t(Iend  ,j-1,k,nstp,itrc)
              grad(Iend+1,j)=t(Iend+1,j  ,k,nstp,itrc)-                 &
     &                       t(Iend+1,j-1,k,nstp,itrc)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dTdt=t(Iend,j,k,nstp,itrc)-t(Iend  ,j,k,nout,itrc)
                dTdx=t(Iend,j,k,nout,itrc)-t(Iend-1,j,k,nout,itrc)
                IF (LBC(ieast,isTvar(itrc),ng)%nudging) THEN
                  tau=Tobc_out(itrc,ng,ieast)
                  IF ((dTdt*dTdx).lt.0.0_r8) THEN
                    tau=Tobc_in(itrc,ng,ieast)
                  ELSE
                    tau=tau*dt(ng)
                  END IF
                END IF
                IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(Iend,j  )+                              &
     &                     grad(Iend,j+1))).gt.0.0_r8) THEN
                  dTde=grad(Iend,j  )
                ELSE
                  dTde=grad(Iend,j+1)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=dTdt*dTdx
                Ce=0.0_r8
                t(Iend+1,j,k,nout,itrc)=(cff*t(Iend+1,j,k,nstp,itrc)+   &
     &                                   Cx *t(Iend  ,j,k,nout,itrc)-   &
     &                                   MAX(Ce,0.0_r8)*                &
     &                                      grad(Iend+1,j  )-           &
     &                                   MIN(Ce,0.0_r8)*                &
     &                                      grad(Iend+1,j+1))/          &
     &                                  (cff+Cx)
                IF (LBC(ieast,isTvar(itrc),ng)%nudging) THEN
                  t(Iend+1,j,k,nout,itrc)=t(Iend+1,j,k,nout,itrc)+      &
     &                                    tau*                          &
     &                                    (BOUNDARY(ng)%t_east(j,k,     &
     &                                                         itrc)-   &
     &                                     t(Iend+1,j,k,nstp,itrc))
                END IF
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                t(Iend+1,j,k,nout,itrc)=BOUNDARY(ng)%t_east(j,k,itrc)
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                t(Iend+1,j,k,nout,itrc)=t(Iend,j,k,nout,itrc)
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
        IF (LBC(isouth,isTvar(itrc),ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=t(i  ,Jstr  ,k,nstp,itrc)-                 &
     &                       t(i-1,Jstr  ,k,nstp,itrc)
              grad(i,Jstr-1)=t(i  ,Jstr-1,k,nstp,itrc)-                 &
     &                       t(i-1,Jstr-1,k,nstp,itrc)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dTdt=t(i,Jstr,k,nstp,itrc)-t(i,Jstr  ,k,nout,itrc)
                dTde=t(i,Jstr,k,nout,itrc)-t(i,Jstr+1,k,nout,itrc)
                IF (LBC(isouth,isTvar(itrc),ng)%nudging) THEN
                  tau=Tobc_out(itrc,ng,isouth)
                  IF ((dTdt*dTde).lt.0.0_r8) THEN
                    tau=Tobc_in(itrc,ng,isouth)
                  ELSE
                    tau=tau*dt(ng)
                  END IF
                END IF
                IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(i  ,Jstr)+                              &
     &                     grad(i+1,Jstr))).gt.0.0_r8) THEN
                  dTdx=grad(i  ,Jstr)
                ELSE
                  dTdx=grad(i+1,Jstr)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=0.0_r8
                Ce=dTdt*dTde
                t(i,Jstr-1,k,nout,itrc)=(cff*t(i,Jstr-1,k,nstp,itrc)+   &
     &                                   Ce *t(i,Jstr  ,k,nout,itrc )-  &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      grad(i  ,Jstr-1)-           &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      grad(i+1,Jstr-1))/          &
     &                                  (cff+Ce)
                IF (LBC(isouth,isTvar(itrc),ng)%nudging) THEN
                  t(i,Jstr-1,k,nout,itrc)=t(i,Jstr-1,k,nout,itrc)+      &
     &                                    tau*                          &
     &                                    (BOUNDARY(ng)%t_south(i,k,    &
     &                                                          itrc)-  &
     &                                     t(i,Jstr-1,k,nstp,itrc))
                END IF
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                t(i,Jstr-1,k,nout,itrc)=BOUNDARY(ng)%t_south(i,k,itrc)
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                t(i,Jstr-1,k,nout,itrc)=t(i,Jstr,k,nout,itrc)
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
        IF (LBC(inorth,isTvar(itrc),ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=t(i  ,Jend  ,k,nstp,itrc)-                 &
     &                       t(i-1,Jend  ,k,nstp,itrc)
              grad(i,Jend+1)=t(i  ,Jend+1,k,nstp,itrc)-                 &
     &                       t(i-1,Jend+1,k,nstp,itrc)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dTdt=t(i,Jend,k,nstp,itrc)-t(i,Jend  ,k,nout,itrc)
                dTde=t(i,Jend,k,nout,itrc)-t(i,Jend-1,k,nout,itrc)
                IF (LBC(inorth,isTvar(itrc),ng)%nudging) THEN
                  tau=Tobc_out(itrc,ng,inorth)
                  IF ((dTdt*dTde).lt.0.0_r8) THEN
                    tau=Tobc_in(itrc,ng,inorth)
                  ELSE
                    tau=tau*dt(ng)
                  END IF
                END IF
                IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
                IF ((dTdt*(grad(i  ,Jend)+                              &
     &                     grad(i+1,Jend))).gt.0.0_r8) THEN
                  dTdx=grad(i  ,Jend)
                ELSE
                  dTdx=grad(i+1,Jend)
                END IF
                cff=MAX(dTdx*dTdx+dTde*dTde,eps)
                Cx=0.0_r8
                Ce=dTdt*dTde
                t(i,Jend+1,k,nout,itrc)=(cff*t(i,Jend+1,k,nstp,itrc)+   &
     &                                   Ce *t(i,Jend  ,k,nout,itrc)-   &
     &                                   MAX(Cx,0.0_r8)*                &
     &                                      grad(i  ,Jend+1)-           &
     &                                   MIN(Cx,0.0_r8)*                &
     &                                      grad(i+1,Jend+1))/          &
     &                                  (cff+Ce)
                IF (LBC(inorth,isTvar(itrc),ng)%nudging) THEN
                  t(i,Jend+1,k,nout,itrc)=t(i,Jend+1,k,nout,itrc)+      &
     &                                    tau*                          &
     &                                    (BOUNDARY(ng)%t_north(i,k,    &
     &                                                          itrc)-  &
     &                                     t(i,Jend+1,k,nstp,itrc))
                END IF
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isTvar(itrc),ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                t(i,Jend+1,k,nout,itrc)=BOUNDARY(ng)%t_north(i,k,itrc)
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isTvar(itrc),ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isTvar(itrc),ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                t(i,Jend+1,k,nout,itrc)=t(i,Jend,k,nout,itrc)
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
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=1,N(ng)
              t(Istr-1,Jstr-1,k,nout,itrc)=0.5_r8*                      &
     &                                     (t(Istr  ,Jstr-1,k,nout,     &
     &                                        itrc)+                    &
     &                                      t(Istr-1,Jstr  ,k,nout,     &
     &                                        itrc))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=1,N(ng)
              t(Iend+1,Jstr-1,k,nout,itrc)=0.5_r8*                      &
     &                                     (t(Iend  ,Jstr-1,k,nout,     &
     &                                        itrc)+                    &
     &                                      t(Iend+1,Jstr  ,k,nout,     &
     &                                        itrc))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=1,N(ng)
              t(Istr-1,Jend+1,k,nout,itrc)=0.5_r8*                      &
     &                                     (t(Istr-1,Jend  ,k,nout,     &
     &                                        itrc)+                    &
     &                                      t(Istr  ,Jend+1,k,nout,     &
     &                                        itrc))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=1,N(ng)
              t(Iend+1,Jend+1,k,nout,itrc)=0.5_r8*                      &
     &                                     (t(Iend+1,Jend  ,k,nout,     &
     &                                        itrc)+                    &
     &                                      t(Iend  ,Jend+1,k,nout,     &
     &                                        itrc))
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE t3dbc_tile
      END MODULE t3dbc_mod
