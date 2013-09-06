      MODULE shapiro_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group         Kate Hedstrom   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This package contains shapiro filter routines for order 2 and       !
!  reduced order at the boundary and mask edges.                       !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    shapirp2d_tile       Shapiro filter for 2D fields.                !
!    shapirp3d_tile       Shapiro filter for 3D fields.                !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE shapiro2d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Amask,                                 &
     &                           A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      real(r8), intent(in) :: Amask(LBi:,LBj:)
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk2
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
!  Shapiro filter requested 2D field.
!-----------------------------------------------------------------------
!
!  This subroutine will apply a Shapiro filter of order 2 (defined
!  as twice the order in Shapiro (1970), with N even) to an array, A.
!  The order of the filter is reduced at the boundaries and at the
!  mask edges, if any.
!
!  Initialize filter in the Y-direction.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
          Awrk1(i,j)=0.25_r8*                                           &
     &               (A(i,j-1)*Amask(i,j-1)+                            &
     &                A(i,j+1)*Amask(i,j+1)-                            &
     &                2.0_r8*A(i,j)*Amask(i,j))*                        &
     &               Amask(i,j-1)*Amask(i,j+1)*Amask(i,j)
        END DO
      END DO
!
!  Add the changes to the field.
!
      DO j=Jstr,Jend
        DO i=Istr-1,Iend+1
          Awrk2(i,j)=A(i,j)+Awrk1(i,j)
        END DO
      END DO
!
!  Initialize filter in the X-direction.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Awrk1(i,j)=0.25_r8*                                           &
     &               (Awrk2(i-1,j)*Amask(i-1,j)+                        &
     &                Awrk2(i+1,j)*Amask(i+1,j)-                        &
     &                2.0_r8*Awrk2(i,j)*Amask(i,j))*                    &
     &               Amask(i-1,j)*Amask(i+1,j)*Amask(i,j)
        END DO
      END DO
!
!  Add changes to field.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          A(i,j)=Awrk2(i,j)+Awrk1(i,j)
        END DO
      END DO
      RETURN
      END SUBROUTINE shapiro2d_tile
!
!***********************************************************************
      SUBROUTINE shapiro3d_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Amask,                                 &
     &                           A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: Amask(LBi:,LBj:)
      real(r8), intent(inout) :: A(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk1
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Awrk2
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
!  Shapiro filter requested 3D field.
!-----------------------------------------------------------------------
!
!  This subroutine will apply a Shapiro filter of order 2 (defined
!  as twice the order in Shapiro (1970), with N even) to an array, A.
!  The order of the filter is reduced at the boundaries and at the
!  mask edges, if any.
!
!  Initialize filter in the Y-direction.
!
      DO k=LBk,UBk
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (A(i,j-1,k)*Amask(i,j-1)+                        &
     &                  A(i,j+1,k)*Amask(i,j+1)-                        &
     &                  2.0_r8*A(i,j,k)*Amask(i,j))*                    &
     &                 Amask(i,j-1)*Amask(i,j+1)*Amask(i,j)
          END DO
        END DO
!
!  Add the changes to the field.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            Awrk2(i,j)=A(i,j,k)+Awrk1(i,j)
          END DO
        END DO
!
!  Initialize filter in the X-direction.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Awrk1(i,j)=0.25_r8*                                         &
     &                 (Awrk2(i-1,j)*Amask(i-1,j)+                      &
     &                  Awrk2(i+1,j)*Amask(i+1,j)-                      &
     &                  2.0_r8*Awrk2(i,j)*Amask(i,j))*                  &
     &                 Amask(i-1,j)*Amask(i+1,j)*Amask(i,j)
          END DO
        END DO
!
!  Add changes to field.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            A(i,j,k)=Awrk2(i,j)+Awrk1(i,j)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE shapiro3d_tile
      END MODULE shapiro_mod
