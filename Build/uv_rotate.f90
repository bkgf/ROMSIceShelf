      MODULE uv_rotate_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!=======================================================================
!                                                                      !
!  These routines average momentum component to RHO-points and then    !
!  rotates from (XI,ETA) coordinates to geographical Eastward and      !
!  Northward directions.                                               !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: uv_rotate2d
      PUBLIC  :: uv_rotate3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE uv_rotate2d (ng, tile, add, Lboundary,                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        CosAngler, SinAngler,                     &
     &                        rmask_io,                                 &
     &                        Uinp, Vinp, Uout, Vout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
!
!  Imported variable declarations.
!
      logical, intent(in) :: add, Lboundary
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
      real(r8), intent(in) :: rmask_io(LBi:,LBj:)
      real(r8), intent(in) :: Uinp(LBi:,LBj:)
      real(r8), intent(in) :: Vinp(LBi:,LBj:)
      real(r8), intent(inout) :: Uout(LBi:,LBj:)
      real(r8), intent(inout) :: Vout(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: Urho, Vrho
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
!  Rotate 2D vector components to Eastward and Northward directions.
!-----------------------------------------------------------------------
!
      IF (add) THEN
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Urho=0.5_r8*(Uinp(i,j)+Uinp(i+1,j))
            Vrho=0.5_r8*(Vinp(i,j)+Vinp(i,j+1))
            Uout(i,j)=Uout(i,j)+                                        &
     &                Urho*CosAngler(i,j)-                              &
     &                Vrho*SinAngler(i,j)
            Vout(i,j)=Vout(i,j)+                                        &
     &                Vrho*CosAngler(i,j)+                              &
     &                Urho*SinAngler(i,j)
            Uout(i,j)=Uout(i,j)*rmask_io(i,j)
            Vout(i,j)=Vout(i,j)*rmask_io(i,j)
          END DO
        END DO
      ELSE
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Urho=0.5_r8*(Uinp(i,j)+Uinp(i+1,j))
            Vrho=0.5_r8*(Vinp(i,j)+Vinp(i,j+1))
            Uout(i,j)=Urho*CosAngler(i,j)-                              &
     &                Vrho*SinAngler(i,j)
            Vout(i,j)=Vrho*CosAngler(i,j)+                              &
     &                Urho*SinAngler(i,j)
            Uout(i,j)=Uout(i,j)*rmask_io(i,j)
            Vout(i,j)=Vout(i,j)*rmask_io(i,j)
          END DO
        END DO
      END IF
!
!  Exchange boundary data, if applicable.
!
      IF (Lboundary) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Uout)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Vout)
        END IF
      END IF
      RETURN
      END SUBROUTINE uv_rotate2d
!
!***********************************************************************
      SUBROUTINE uv_rotate3d (ng, tile, add, Lboundary,                 &
     &                        LBi, UBi, LBj, UBj, LBk, UBk,             &
     &                        CosAngler, SinAngler,                     &
     &                        rmask_io,                                 &
     &                        Uinp, Vinp, Uout, Vout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
!
!  Imported variable declarations.
!
      logical, intent(in) :: add, Lboundary
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
      real(r8), intent(in) :: rmask_io(LBi:,LBj:)
      real(r8), intent(in) :: Uinp(LBi:,LBj:,LBk:)
      real(r8), intent(in) :: Vinp(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: Uout(LBi:,LBj:,LBk:)
      real(r8), intent(inout) :: Vout(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: Urho, Vrho
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
!  Rotate 3D vector components to Eastward and Northward directions.
!-----------------------------------------------------------------------
!
      IF (add) THEN
        DO k=LBk,UBk
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Urho=0.5_r8*(Uinp(i,j,k)+Uinp(i+1,j,k))
              Vrho=0.5_r8*(Vinp(i,j,k)+Vinp(i,j+1,k))
              Uout(i,j,k)=Uout(i,j,k)+                                  &
     &                    Urho*CosAngler(i,j)-                          &
     &                    Vrho*SinAngler(i,j)
              Vout(i,j,k)=Vout(i,j,k)+                                  &
     &                    Vrho*CosAngler(i,j)+                          &
     &                    Urho*SinAngler(i,j)
              Uout(i,j,k)=Uout(i,j,k)*rmask_io(i,j)
              Vout(i,j,k)=Vout(i,j,k)*rmask_io(i,j)
            END DO
          END DO
        END DO
      ELSE
        DO k=LBk,UBk
          DO j=Jstr,Jend
            DO i=Istr,Iend
              Urho=0.5_r8*(Uinp(i,j,k)+Uinp(i+1,j,k))
              Vrho=0.5_r8*(Vinp(i,j,k)+Vinp(i,j+1,k))
              Uout(i,j,k)=Urho*CosAngler(i,j)-                          &
     &                    Vrho*SinAngler(i,j)
              Vout(i,j,k)=Vrho*CosAngler(i,j)+                          &
     &                    Urho*SinAngler(i,j)
              Uout(i,j,k)=Uout(i,j,k)*rmask_io(i,j)
              Vout(i,j,k)=Vout(i,j,k)*rmask_io(i,j)
            END DO
          END DO
        END DO
      END IF
!
!  Exchange boundary data.
!
      IF (Lboundary) THEN
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            Uout)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            Vout)
        END IF
      END IF
      RETURN
      END SUBROUTINE uv_rotate3d
      END MODULE uv_rotate_mod
