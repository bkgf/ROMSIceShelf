      MODULE state_dotprod_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the dot product between two model states:     !
!                                                                      !
!      DotProd(0:NstateVars) = < s1, s2 >                              !
!                                                                      !
!  where                                                               !
!                                                                      !
!      DotProd(0)           All state variable dot product             !
!      DotProd(isUvel)      3D U-momentum contribution                 !
!      DotProd(isVvel)      3D V-momentum contribution                 !
!      DotProd(isTvar(:))   Tracer-type variables contribution         !
!      DotProd(isFsur)      Free-surface contribution                  !
!                                                                      !
!=======================================================================
!
      implicit none
      PUBLIC  :: state_dotprod
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_dotprod (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          NstateVars, DotProd,                    &
     &                          rmask, umask, vmask,                    &
     &                          s1_t, s2_t,                             &
     &                          s1_u, s2_u,                             &
     &                          s1_v, s2_v,                             &
     &                          s1_zeta, s2_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: NstateVars
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: s1_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s2_t(LBi:,LBj:,:,:)
      real(r8), intent(in) :: s1_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_u(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s2_v(LBi:,LBj:,:)
      real(r8), intent(in) :: s1_zeta(LBi:,LBj:)
      real(r8), intent(in) :: s2_zeta(LBi:,LBj:)
!
      real(r8), intent(out), dimension(0:NstateVars) :: DotProd
!
!  Local variable declarations.
!
      integer :: NSUB, i, j, k
      integer :: ir, it
      real(r8) :: cff
      real(r8), dimension(0:NstateVars) :: my_DotProd
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
!  Compute dot product between S1 and S2 model state trajectories.
!-----------------------------------------------------------------------
!
      DO i=0,NstateVars
        my_DotProd(i)=0.0_r8
      END DO
!
!  Free-surface.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          cff=s1_zeta(i,j)*s2_zeta(i,j)
          cff=cff*rmask(i,j)
          my_DotProd(0)=my_DotProd(0)+cff
          my_DotProd(isFsur)=my_DotProd(isFsur)+cff
        END DO
      END DO
!
!  3D U-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            cff=s1_u(i,j,k)*s2_u(i,j,k)
            cff=cff*umask(i,j)
            my_DotProd(0)=my_DotProd(0)+cff
            my_DotProd(isUvel)=my_DotProd(isUvel)+cff
          END DO
        END DO
      END DO
!
!  3D V-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            cff=s1_v(i,j,k)*s2_v(i,j,k)
            cff=cff*vmask(i,j)
            my_DotProd(0)=my_DotProd(0)+cff
            my_DotProd(isVvel)=my_DotProd(isVvel)+cff
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO it=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              cff=s1_t(i,j,k,it)*s2_t(i,j,k,it)
              cff=cff*rmask(i,j)
              my_DotProd(0)=my_DotProd(0)+cff
              my_DotProd(isTvar(it))=my_DotProd(isTvar(it))+cff
            END DO
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Perform parallel global reduction operations.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
      IF (tile_count.eq.0) THEN
        DO i=0,NstateVars
          DotProd(i)=0.0_r8
        END DO
      END IF
      DO i=0,NstateVars
        DotProd(i)=DotProd(i)+my_DotProd(i)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
      END IF
      RETURN
      END SUBROUTINE state_dotprod
      END MODULE state_dotprod_mod
