      MODULE exchange_3d_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  These routines apply periodic boundary conditions to generic        !
!  3D fields.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng                      Nested grid number.                      !
!     tile                    Domain partition.                        !
!     LBi                     I-dimension Lower bound.                 !
!     UBi                     I-dimension Upper bound.                 !
!     LBj                     J-dimension Lower bound.                 !
!     UBj                     J-dimension Upper bound.                 !
!     LBk                     J-dimension Lower bound.                 !
!     UBk                     J-dimension Upper bound.                 !
!     A                       3D field.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A                       Processed 3D field.                      !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!     exchange_p3d_tile       periodic conditions at PSI-points        !
!     exchange_r3d_tile       periodic conditions at RHO-points        !
!     exchange_u3d_tile       periodic conditions at U-points          !
!     exchange_v3d_tile       periodic conditions at V-points          !
!     exchange_w3d_tile       periodic conditions at W-points          !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
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
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        EW_exchange=.TRUE.
      ELSE
        EW_exchange=.FALSE.
      END IF
      IF (NSperiodic(ng)) THEN
        NS_exchange=.TRUE.
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=Jstr
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=Istr
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
              A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
              A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
              A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
              A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
              A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
              A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
              A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
              A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE exchange_p3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
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
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        EW_exchange=.TRUE.
      ELSE
        EW_exchange=.FALSE.
      END IF
      IF (NSperiodic(ng)) THEN
        NS_exchange=.TRUE.
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
              A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
              A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
              A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
              A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
              A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
              A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
              A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
              A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE exchange_r3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
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
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        EW_exchange=.TRUE.
      ELSE
        EW_exchange=.FALSE.
      END IF
      IF (NSperiodic(ng)) THEN
        NS_exchange=.TRUE.
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=Istr
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
              A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
              A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
              A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
              A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
              A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
              A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
              A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
              A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE exchange_u3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
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
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        EW_exchange=.TRUE.
      ELSE
        EW_exchange=.FALSE.
      END IF
      IF (NSperiodic(ng)) THEN
        NS_exchange=.TRUE.
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=Jstr
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
              A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
              A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
              A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
              A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
              A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
              A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
              A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
              A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE exchange_v3d_tile
!
!***********************************************************************
      SUBROUTINE exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              A)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      logical :: EW_exchange
      logical :: NS_exchange
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
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
!  Determine processing switches.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        EW_exchange=.TRUE.
      ELSE
        EW_exchange=.FALSE.
      END IF
      IF (NSperiodic(ng)) THEN
        NS_exchange=.TRUE.
      ELSE
        NS_exchange=.FALSE.
      END IF
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng)) THEN
        IF (NSperiodic(ng)) THEN
          Jmin=Jstr
          Jmax=Jend
        ELSE
          Jmin=JstrR
          Jmax=JendR
        END IF
!
        IF (EW_exchange) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(Lm(ng)+1,j,k)=A(1,j,k)
                A(Lm(ng)+2,j,k)=A(2,j,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO j=Jmin,Jmax
                  A(Lm(ng)+3,j,k)=A(3,j,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO j=Jmin,Jmax
                A(-2,j,k)=A(Lm(ng)-2,j,k)
                A(-1,j,k)=A(Lm(ng)-1,j,k)
                A( 0,j,k)=A(Lm(ng)  ,j,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  North-South periodic boundary conditions.
!-----------------------------------------------------------------------
!
      IF (NSperiodic(ng)) THEN
        IF (EWperiodic(ng)) THEN
          Imin=Istr
          Imax=Iend
        ELSE
          Imin=IstrR
          Imax=IendR
        END IF
!
        IF (NS_exchange) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,Mm(ng)+1,k)=A(i,1,k)
                A(i,Mm(ng)+2,k)=A(i,2,k)
              END DO
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                DO i=Imin,Imax
                  A(i,Mm(ng)+3,k)=A(i,3,k)
                END DO
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO k=LBk,UBk
              DO i=Imin,Imax
                A(i,-2,k)=A(i,Mm(ng)-2,k)
                A(i,-1,k)=A(i,Mm(ng)-1,k)
                A(i, 0,k)=A(i,Mm(ng)  ,k)
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).and.NSperiodic(ng)) THEN
        IF (EW_exchange.and.NS_exchange) THEN
          IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,Mm(ng)+1,k)=A(1,1,k)
              A(Lm(ng)+1,Mm(ng)+2,k)=A(1,2,k)
              A(Lm(ng)+2,Mm(ng)+1,k)=A(2,1,k)
              A(Lm(ng)+2,Mm(ng)+2,k)=A(2,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+1,Mm(ng)+3,k)=A(1,3,k)
                A(Lm(ng)+2,Mm(ng)+3,k)=A(2,3,k)
                A(Lm(ng)+3,Mm(ng)+1,k)=A(3,1,k)
                A(Lm(ng)+3,Mm(ng)+2,k)=A(3,2,k)
                A(Lm(ng)+3,Mm(ng)+3,k)=A(3,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,Mm(ng)+1,k)=A(Lm(ng)-2,1,k)
              A(-1,Mm(ng)+1,k)=A(Lm(ng)-1,1,k)
              A( 0,Mm(ng)+1,k)=A(Lm(ng)  ,1,k)
              A(-2,Mm(ng)+2,k)=A(Lm(ng)-2,2,k)
              A(-1,Mm(ng)+2,k)=A(Lm(ng)-1,2,k)
              A( 0,Mm(ng)+2,k)=A(Lm(ng)  ,2,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(-2,Mm(ng)+3,k)=A(Lm(ng)-2,3,k)
                A(-1,Mm(ng)+3,k)=A(Lm(ng)-1,3,k)
                A( 0,Mm(ng)+3,k)=A(Lm(ng)  ,3,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
            DO k=LBk,UBk
              A(Lm(ng)+1,-2,k)=A(1,Mm(ng)-2,k)
              A(Lm(ng)+1,-1,k)=A(1,Mm(ng)-1,k)
              A(Lm(ng)+1, 0,k)=A(1,Mm(ng)  ,k)
              A(Lm(ng)+2,-2,k)=A(2,Mm(ng)-2,k)
              A(Lm(ng)+2,-1,k)=A(2,Mm(ng)-1,k)
              A(Lm(ng)+2, 0,k)=A(2,Mm(ng)  ,k)
            END DO
            IF (NghostPoints.eq.3) THEN
              DO k=LBk,UBk
                A(Lm(ng)+3,-2,k)=A(3,Mm(ng)-2,k)
                A(Lm(ng)+3,-1,k)=A(3,Mm(ng)-1,k)
                A(Lm(ng)+3, 0,k)=A(3,Mm(ng)  ,k)
              END DO
            END IF
          END IF
          IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
            DO k=LBk,UBk
              A(-2,-2,k)=A(Lm(ng)-2,Mm(ng)-2,k)
              A(-2,-1,k)=A(Lm(ng)-2,Mm(ng)-1,k)
              A(-2, 0,k)=A(Lm(ng)-2,Mm(ng)  ,k)
              A(-1,-2,k)=A(Lm(ng)-1,Mm(ng)-2,k)
              A(-1,-1,k)=A(Lm(ng)-1,Mm(ng)-1,k)
              A(-1, 0,k)=A(Lm(ng)-1,Mm(ng)  ,k)
              A( 0,-2,k)=A(Lm(ng)  ,Mm(ng)-2,k)
              A( 0,-1,k)=A(Lm(ng)  ,Mm(ng)-1,k)
              A( 0, 0,k)=A(Lm(ng)  ,Mm(ng)  ,k)
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE exchange_w3d_tile
      END MODULE exchange_3d_mod
