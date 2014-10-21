      MODULE mod_boundary
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
!  zeta_west      Free-surface (m) western boundary conditions.        !
!  zetaG_west     Latest two-time snapshots of input free-surface (m)  !
!                   western boundary data.                             !
!  zeta_east      Free-surface (m) eastern boundary conditions.        !
!  zetaG_east     Latest two-time snapshots of input free-surface (m)  !
!                   eastern boundary data.                             !
!  zeta_south     Free-surface (m) southern boundary conditions.       !
!  zetaG_south    Latest two-time snapshots of input free-surface (m)  !
!                   southern boundary data.                            !
!  zeta_north     Free-surface (m) northern boundary conditions.       !
!  zetaG_north    Latest two-time snapshots of input free-surface (m)  !
!                   northern boundary data.                            !
!  ubar_west      2D u-momentum (m/s) western boundary conditions.     !
!  vbar_west      2D v-momentum (m/s) western boundary conditions.     !
!  ubarG_west     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vbarG_west     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) western boundary data.                       !
!  ubar_east      2D u-momentum (m/s) eastern boundary conditions.     !
!  vbar_east      2D v-momentum (m/s) eastern boundary conditions.     !
!  ubarG_east     Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vbarG_east     Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) eastern boundary data.                       !
!  ubar_south     2D u-momentum (m/s) southern boundary conditions.    !
!  vbar_south     2D v-momentum (m/s) southern boundary conditions.    !
!  ubarG_south    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  vbarG_south    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) southern boundary data.                      !
!  ubar_north     2D u-momentum (m/s) northern boundary conditions.    !
!  vbar_north     2D v-momentum (m/s) northern boundary conditions.    !
!  ubarG_north    Latest two-time snapshots of input 2D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  vbarG_north    Latest two-time snapshots of input 2D v-momentum     !
!                   (m/s) northern boundary data.                      !
!  u_west         3D u-momentum (m/s) western boundary conditions.     !
!  v_west         3D v-momentum (m/s) western boundary conditions.     !
!  uG_west        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) western boundary data.                       !
!  vG_west        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) western boundary data.                       !
!  u_east         3D u-momentum (m/s) eastern boundary conditions.     !
!  v_east         3D v-momentum (m/s) eastern boundary conditions.     !
!  uG_east        Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) eastern boundary data.                       !
!  vG_east        Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) eastern boundary data.                       !
!  u_south        3D u-momentum (m/s) southern boundary conditions.    !
!  v_south        3D v-momentum (m/s) southern boundary conditions.    !
!  uG_south       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) southern boundary data.                      !
!  vG_south       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) southern boundary data.                      !
!  u_north        3D u-momentum (m/s) northern boundary conditions.    !
!  v_north        3D v-momentum (m/s) northern boundary conditions.    !
!  uG_north       Latest two-time snapshots of input 3D u-momentum     !
!                   (m/s) northern boundary data.                      !
!  vG_north       Latest two-time snapshots of input 3D v-momentum     !
!                   (m/s) northern boundary data.                      !
!  t_west         Tracer (T units) western boundary conditions.        !
!  tG_west        Latest two-time snapshots of input tracer (Tunits)   !
!                   western boundary data.                             !
!  t_east         Tracer (T units) eastern boundary conditions.        !
!  tG_east        Latest two-time snapshots of input tracer (Tunits)   !
!                   eastern boundary data.                             !
!  t_south        Tracer (T units) southern boundary conditions.       !
!  tG_south       Latest two-time snapshots of input tracer (Tunits)   !
!                   southern boundary data.                            !
!  t_north        Tracer (T units) northern boundary conditions.       !
!  tG_north       Latest two-time snapshots of input tracer (Tunits)   !
!                   northern boundary data.                            !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
!
!-----------------------------------------------------------------------
!  Lateral boundary condition apply switches.
!-----------------------------------------------------------------------
!
!  The following switches are used to control which grid points are
!  processed by the lateral boundary conditions. These switches are
!  set to TRUE by default.  However in composite grids, the points
!  processed by nesting are set to FALSE to allow mixed boundary
!  conditions along the grid edges.
!
        TYPE T_APPLY
          logical, pointer :: west(:)
          logical, pointer :: east(:)
          logical, pointer :: south(:)
          logical, pointer :: north(:)
        END TYPE
        TYPE (T_APPLY), allocatable :: LBC_apply(:)
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions structure.
!-----------------------------------------------------------------------
!
        TYPE T_BOUNDARY
!
!  Nonlinear model state.
!
          real(r8), pointer :: zeta_west(:)
          real(r8), pointer :: zetaG_west(:,:)
          real(r8), pointer :: zeta_east(:)
          real(r8), pointer :: zetaG_east(:,:)
          real(r8), pointer :: zeta_south(:)
          real(r8), pointer :: zetaG_south(:,:)
          real(r8), pointer :: zeta_north(:)
          real(r8), pointer :: zetaG_north(:,:)
          real(r8), pointer :: ubar_west(:)
          real(r8), pointer :: vbar_west(:)
          real(r8), pointer :: ubarG_west(:,:)
          real(r8), pointer :: vbarG_west(:,:)
          real(r8), pointer :: ubar_east(:)
          real(r8), pointer :: vbar_east(:)
          real(r8), pointer :: ubarG_east(:,:)
          real(r8), pointer :: vbarG_east(:,:)
          real(r8), pointer :: ubar_south(:)
          real(r8), pointer :: vbar_south(:)
          real(r8), pointer :: ubarG_south(:,:)
          real(r8), pointer :: vbarG_south(:,:)
          real(r8), pointer :: ubar_north(:)
          real(r8), pointer :: vbar_north(:)
          real(r8), pointer :: ubarG_north(:,:)
          real(r8), pointer :: vbarG_north(:,:)
          real(r8), pointer :: u_west(:,:)
          real(r8), pointer :: v_west(:,:)
          real(r8), pointer :: uG_west(:,:,:)
          real(r8), pointer :: vG_west(:,:,:)
          real(r8), pointer :: u_east(:,:)
          real(r8), pointer :: v_east(:,:)
          real(r8), pointer :: uG_east(:,:,:)
          real(r8), pointer :: vG_east(:,:,:)
          real(r8), pointer :: u_south(:,:)
          real(r8), pointer :: v_south(:,:)
          real(r8), pointer :: uG_south(:,:,:)
          real(r8), pointer :: vG_south(:,:,:)
          real(r8), pointer :: u_north(:,:)
          real(r8), pointer :: v_north(:,:)
          real(r8), pointer :: uG_north(:,:,:)
          real(r8), pointer :: vG_north(:,:,:)
          real(r8), pointer :: t_west(:,:,:)
          real(r8), pointer :: tG_west(:,:,:,:)
          real(r8), pointer :: t_east(:,:,:)
          real(r8), pointer :: tG_east(:,:,:,:)
          real(r8), pointer :: t_south(:,:,:)
          real(r8), pointer :: tG_south(:,:,:,:)
          real(r8), pointer :: t_north(:,:,:)
          real(r8), pointer :: tG_north(:,:,:,:)
        END TYPE T_BOUNDARY
        TYPE (T_BOUNDARY), allocatable ::BOUNDARY(:)
      CONTAINS
      SUBROUTINE allocate_boundary (ng)
!
!=======================================================================
!                                                                      !
!  This routine initializes all variables in the module for all nested !
!  grids.  Currently, there is not parallel tiling in boundary arrays. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: my_tile
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Se dimension ranges. Notice that the boundary arrays are dimensioned
!  with the global dimensions of grid. That is, no tiling ranges in
!  distributed-memory. This is done to facilitate processing.
!
      my_tile=-1                           ! for global values
      LBi=BOUNDS(ng)%LBi(my_tile)
      UBi=BOUNDS(ng)%UBi(my_tile)
      LBj=BOUNDS(ng)%LBj(my_tile)
      UBj=BOUNDS(ng)%UBj(my_tile)
!
!  Allocate structures.
!
      IF (ng.eq.1) THEN
        allocate ( LBC_apply(Ngrids) )
        allocate ( BOUNDARY(Ngrids) )
      END IF
!
!  Lateral boundary conditions apply switches.  These switches need to
!  be initilized to TRUE here because 'initialize_boundary' is called
!  several times in adjoint-based application to clear state arrays.
!  These switches are part of the application grid and will be set to
!  FALSE elsewhere, if the boundary point is assigned by a nested grid.
!
      allocate ( LBC_apply(ng) % west(LBj:UBj) )
      LBC_apply(ng) % west = .TRUE.
      allocate ( LBC_apply(ng) % east(LBj:UBj) )
      LBC_apply(ng) % east = .TRUE.
      allocate ( LBC_apply(ng) % south(LBi:UBi) )
      LBC_apply(ng) % south = .TRUE.
      allocate ( LBC_apply(ng) % north(LBi:UBi) )
      LBC_apply(ng) % north = .TRUE.
!
!  Nonlinear model state.
!
      IF (LBC(iwest,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_west(LBj:UBj) )
        allocate ( BOUNDARY(ng) % zetaG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_east(LBj:UBj) )
        allocate ( BOUNDARY(ng) % zetaG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_south(LBi:UBi) )
        allocate ( BOUNDARY(ng) % zetaG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_north(LBi:UBi) )
        allocate ( BOUNDARY(ng) % zetaG_north(LBi:UBi,2) )
      END IF
!
      IF (LBC(iwest,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_west(LBj:UBj) )
        allocate ( BOUNDARY(ng) % ubarG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_east(LBj:UBj) )
        allocate ( BOUNDARY(ng) % ubarG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_south(LBi:UBi) )
        allocate ( BOUNDARY(ng) % ubarG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_north(LBi:UBi) )
        allocate ( BOUNDARY(ng) % ubarG_north(LBi:UBi,2) )
      END IF
!
      IF (LBC(iwest,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_west(LBj:UBj) )
        allocate ( BOUNDARY(ng) % vbarG_west(LBj:UBj,2) )
      END IF
      IF (LBC(ieast,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_east(LBj:UBj) )
        allocate ( BOUNDARY(ng) % vbarG_east(LBj:UBj,2) )
      END IF
      IF (LBC(isouth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_south(LBi:UBi) )
        allocate ( BOUNDARY(ng) % vbarG_south(LBi:UBi,2) )
      END IF
      IF (LBC(inorth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_north(LBi:UBi) )
        allocate ( BOUNDARY(ng) % vbarG_north(LBi:UBi,2) )
      END IF
!
      IF (LBC(iwest,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_west(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % uG_west(LBj:UBj,N(ng),2) )
      END IF
      IF (LBC(ieast,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_east(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % uG_east(LBj:UBj,N(ng),2) )
      END IF
      IF (LBC(isouth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_south(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % uG_south(LBi:UBi,N(ng),2) )
      END IF
      IF (LBC(inorth,isUvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % u_north(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % uG_north(LBi:UBi,N(ng),2) )
      END IF
!
      IF (LBC(iwest,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_west(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % vG_west(LBj:UBj,N(ng),2) )
      END IF
      IF (LBC(ieast,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_east(LBj:UBj,N(ng)) )
        allocate ( BOUNDARY(ng) % vG_east(LBj:UBj,N(ng),2) )
      END IF
      IF (LBC(isouth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_south(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % vG_south(LBi:UBi,N(ng),2) )
      END IF
      IF (LBC(inorth,isVvel,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % v_north(LBi:UBi,N(ng)) )
        allocate ( BOUNDARY(ng) % vG_north(LBi:UBi,N(ng),2) )
      END IF
!
      IF (ANY(LBC(iwest,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_west(LBj:UBj,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % tG_west(LBj:UBj,N(ng),2,NT(ng)) )
      END IF
      IF (ANY(LBC(ieast,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_east(LBj:UBj,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % tG_east(LBj:UBj,N(ng),2,NT(ng)) )
      END IF
      IF (ANY(LBC(isouth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_south(LBi:UBi,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % tG_south(LBi:UBi,N(ng),2,NT(ng)) )
      END IF
      IF (ANY(LBC(inorth,isTvar(:),ng)%acquire)) THEN
        allocate ( BOUNDARY(ng) % t_north(LBi:UBi,N(ng),NT(ng)) )
        allocate ( BOUNDARY(ng) % tG_north(LBi:UBi,N(ng),2,NT(ng)) )
      END IF
      RETURN
      END SUBROUTINE allocate_boundary
      SUBROUTINE initialize_boundary (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8
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
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_west = IniVal
          BOUNDARY(ng) % zetaG_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_east = IniVal
          BOUNDARY(ng) % zetaG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_south = IniVal
          BOUNDARY(ng) % zetaG_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_north = IniVal
          BOUNDARY(ng) % zetaG_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_west = IniVal
          BOUNDARY(ng) % ubarG_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_east = IniVal
          BOUNDARY(ng) % ubarG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_south = IniVal
          BOUNDARY(ng) % ubarG_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_north = IniVal
          BOUNDARY(ng) % ubarG_north = IniVal
          BOUNDARY(ng) % vbarG_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_west = IniVal
          BOUNDARY(ng) % vbarG_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_east = IniVal
          BOUNDARY(ng) % vbarG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_south = IniVal
          BOUNDARY(ng) % vbarG_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_north = IniVal
          BOUNDARY(ng) % vbarG_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_west = IniVal
          BOUNDARY(ng) % uG_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_east = IniVal
          BOUNDARY(ng) % uG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_south = IniVal
          BOUNDARY(ng) % uG_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUvel,ng)%acquire) THEN
          BOUNDARY(ng) % u_north = IniVal
          BOUNDARY(ng) % uG_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_west = IniVal
          BOUNDARY(ng) % vG_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_east = IniVal
          BOUNDARY(ng) % vG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_south = IniVal
          BOUNDARY(ng) % vG_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVvel,ng)%acquire) THEN
          BOUNDARY(ng) % v_north = IniVal
          BOUNDARY(ng) % vG_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      ANY(LBC(iwest,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_west = IniVal
          BOUNDARY(ng) % tG_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      ANY(LBC(ieast,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_east = IniVal
          BOUNDARY(ng) % tG_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      ANY(LBC(isouth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_south = IniVal
          BOUNDARY(ng) % tG_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      ANY(LBC(inorth,isTvar(:),ng)%acquire)) THEN
          BOUNDARY(ng) % t_north = IniVal
          BOUNDARY(ng) % tG_north = IniVal
        END IF
      END IF
      RETURN
      END SUBROUTINE initialize_boundary
      END MODULE mod_boundary
