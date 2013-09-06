      MODULE mod_ocean
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  2D Primitive Variables.                                             !
!                                                                      !
!  rubar        Right-hand-side of 2D U-momentum equation (m4/s2).     !
!  rvbar        Right-hand-side of 2D V-momentum equation (m4/s2).     !
!  rzeta        Right-hand-side of free surface equation (m3/s).       !
!  ubar         Vertically integrated U-momentum component (m/s).      !
!  vbar         Vertically integrated V-momentum component (m/s).      !
!  zeta         Free surface (m).                                      !
!                                                                      !
!  3D Primitive Variables.                                             !
!                                                                      !
!  pden         Potential Density anomaly (kg/m3).                     !
!  rho          Density anomaly (kg/m3).                               !
!  ru           Right-hand-side of 3D U-momentum equation (m4/s2).     !
!  rv           Right hand side of 3D V-momentum equation (m4/s2).     !
!  t            Tracer type variables (active and passive).            !
!  u            3D U-momentum component (m/s).                         !
!  v            3D V-momentum component (m/s).                         !
!  W            S-coordinate (omega*Hz/mn) vertical velocity (m3/s).   !
!                                                                      !
!  Biology Variables.                                                  !
!                                                                      !
!  pH           Surface concentration of hydrogen ions.                !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_OCEAN
!
!  Nonlinear model state.
!
          real(r8), pointer :: rubar(:,:,:)
          real(r8), pointer :: rvbar(:,:,:)
          real(r8), pointer :: rzeta(:,:,:)
          real(r8), pointer :: ubar(:,:,:)
          real(r8), pointer :: vbar(:,:,:)
          real(r8), pointer :: zeta(:,:,:)
          real(r8), pointer :: pden(:,:,:)
          real(r8), pointer :: rho(:,:,:)
          real(r8), pointer :: ru(:,:,:,:)
          real(r8), pointer :: rv(:,:,:,:)
          real(r8), pointer :: t(:,:,:,:,:)
          real(r8), pointer :: u(:,:,:,:)
          real(r8), pointer :: v(:,:,:,:)
          real(r8), pointer :: W(:,:,:)
          real(r8), pointer :: wvel(:,:,:)
        END TYPE T_OCEAN
        TYPE (T_OCEAN), allocatable :: OCEAN(:)
      CONTAINS
      SUBROUTINE allocate_ocean (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OCEAN(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( OCEAN(ng) % rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rzeta(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % zeta(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % rho(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % rv(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % W(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( OCEAN(ng) % wvel(LBi:UBi,LBj:UBj,0:N(ng)) )
      RETURN
      END SUBROUTINE allocate_ocean
      SUBROUTINE initialize_ocean (ng, tile, model)
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, rec
      integer :: itrc, k
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
!  Set array initialization range.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            OCEAN(ng) % rubar(i,j,1) = IniVal
            OCEAN(ng) % rubar(i,j,2) = IniVal
            OCEAN(ng) % rvbar(i,j,1) = IniVal
            OCEAN(ng) % rvbar(i,j,2) = IniVal
            OCEAN(ng) % rzeta(i,j,1) = IniVal
            OCEAN(ng) % rzeta(i,j,2) = IniVal
            OCEAN(ng) % ubar(i,j,1) = IniVal
            OCEAN(ng) % ubar(i,j,2) = IniVal
            OCEAN(ng) % ubar(i,j,3) = IniVal
            OCEAN(ng) % vbar(i,j,1) = IniVal
            OCEAN(ng) % vbar(i,j,2) = IniVal
            OCEAN(ng) % vbar(i,j,3) = IniVal
            OCEAN(ng) % zeta(i,j,1) = IniVal
            OCEAN(ng) % zeta(i,j,2) = IniVal
            OCEAN(ng) % zeta(i,j,3) = IniVal
          END DO
          DO k=1,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % pden(i,j,k) = IniVal
              OCEAN(ng) % rho(i,j,k) = IniVal
              OCEAN(ng) % u(i,j,k,1) = IniVal
              OCEAN(ng) % u(i,j,k,2) = IniVal
              OCEAN(ng) % v(i,j,k,1) = IniVal
              OCEAN(ng) % v(i,j,k,2) = IniVal
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Imin,Imax
              OCEAN(ng) % ru(i,j,k,1) = IniVal
              OCEAN(ng) % ru(i,j,k,2) = IniVal
              OCEAN(ng) % rv(i,j,k,1) = IniVal
              OCEAN(ng) % rv(i,j,k,2) = IniVal
              OCEAN(ng) % W(i,j,k) = IniVal
              OCEAN(ng) % wvel(i,j,k) = IniVal
            END DO
          END DO
          DO itrc=1,NT(ng)
            DO k=1,N(ng)
              DO i=Imin,Imax
                OCEAN(ng) % t(i,j,k,1,itrc) = IniVal
                OCEAN(ng) % t(i,j,k,2,itrc) = IniVal
                OCEAN(ng) % t(i,j,k,3,itrc) = IniVal
              END DO
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_ocean
      END MODULE mod_ocean
