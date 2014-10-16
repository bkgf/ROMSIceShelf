      MODULE mod_ncparam
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all the variables associated with input and    !
!  output  NetCDF  files.  The IO model is very generic and easy to    !
!  change or expand.  The NetCDF files can be in any language.  All    !
!  the IO information is managed using the following variables:        !
!                                                                      !
!  Vname      Input/output variables names and attributes:             !
!               Vname(1,*)  => field variable name.                    !
!               Vname(2,*)  => long-name attribute.                    !
!               Vname(3,*)  => units attribute.                        !
!               Vname(4,*)  => field type attribute.                   !
!               Vname(5,*)  => associated time variable name.          !
!  Tname      Input/output associated time variables names.            !
!                                                                      !
!  Cinfo      Input/output file names associated with each field       !
!                                                                      !
!  Linfo      Input/output fields logical information:                 !
!               Linfo(1,*)  => switch indicating grided data.          !
!               Linfo(2,*)  => switch indicating time cycling.         !
!               Linfo(3,*)  => switch indicating only one-time         !
!                              record available.                       !
!               Linfo(4,*)  => switch indication special record        !
!                              processing (like tides)                 !
!                                                                      !
!  Iinfo      Input/output fields integer information:                 !
!               Iinfo( 1,*) => variable grid type.                     !
!               Iinfo( 2,*) => field variable NetCDF ID.               !
!               Iinfo( 3,*) => associated time variable NetCDF ID.     !
!               Iinfo( 4,*) => number of time records.                 !
!               Iinfo( 5,*) => size of first  spatial dimension.       !
!               Iinfo( 6,*) => size of second spatial dimension.       !
!               Iinfo( 7,*) => size of third  spatial dimension.       !
!               Iinfo( 8,*) => rolling two-time levels index.          !
!               Iinfo( 9,*) => latest processed time record.           !
!               Iinfo(10,*) => number of field multi-files.            !
!                                                                      !
!  Finfo      Input/output field floating-point information:           !
!               Finfo( 1,*) => Starting time (days) of data.           !
!               Finfo( 2,*) => Ending time (days) of data.             !
!               Finfo( 3,*) => Data time lower bound (days) enclosing  !
!                                model starting time.                  !
!               Finfo( 4,*) => Data time upper bound (days) enclosing  !
!                                model starting time.                  !
!               Finfo( 5,*) => length (days) of time cycling.          !
!               Finfo( 6,*) => Scale to convert time to day units.     !
!               Finfo( 7,*) => latest monotonic time (sec).            !
!               Finfo( 8,*) => minimum value for current data.         !
!               Finfo( 9,*) => maximum value for current data.         !
!               Finfo(10,*) => value of scale_factor attribute if any. !
!  Fscale     Scale to convert input data to model units.              !
!  Fpoint     Latest two-time records of input point data.             !
!  Tintrp     Time (sec) of latest field snapshots used for            !
!               interpolation.                                         !
!  Vtime      Latest two-time values of processed input data.          !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
!  Maximum number of variables in a generic NetCDF file (MV) and
!  maximum number of variables in information arrays (NV).
!
        integer, parameter :: MV = 950
        integer, parameter :: NV = 950
!
!  Input/output grid-type variables.
!
        integer, parameter :: p2dvar = 1         ! 2D PSI-variable
        integer, parameter :: r2dvar = 2         ! 2D RHO-variable
        integer, parameter :: u2dvar = 3         ! 2D U-variable
        integer, parameter :: v2dvar = 4         ! 2D V-variable
        integer, parameter :: p3dvar = 5         ! 3D PSI-variable
        integer, parameter :: r3dvar = 6         ! 3D RHO-variable
        integer, parameter :: u3dvar = 7         ! 3D U-variable
        integer, parameter :: v3dvar = 8         ! 3D V-variable
        integer, parameter :: w3dvar = 9         ! 3D W-variable
        integer, parameter :: b3dvar = 10        ! 3D BED-sediment
!
!  Number of horizontal interior and boundary water points.
!
        integer, allocatable :: Nxyp(:)          ! PSI water points
        integer, allocatable :: Nxyr(:)          ! RHO water points
        integer, allocatable :: Nxyu(:)          ! U water points
        integer, allocatable :: Nxyv(:)          ! V water points
!
!  Number of horizontal interior water points only.
!
        integer, allocatable :: NwaterR(:)       ! RHO water points
        integer, allocatable :: NwaterU(:)       ! U water points
        integer, allocatable :: NwaterV(:)       ! V water points
!
!  Lower and upper bound ranges for RHO-type variables for processing
!  state vector and observations.
!
        integer, allocatable :: rILB(:)
        integer, allocatable :: rIUB(:)
        integer, allocatable :: rJLB(:)
        integer, allocatable :: rJUB(:)
!
        real(r8), allocatable :: rXmin(:)
        real(r8), allocatable :: rXmax(:)
        real(r8), allocatable :: rYmin(:)
        real(r8), allocatable :: rYmax(:)
!
!  Lower and upper bound ranges for U-type variables for processing
!  state vector and observations.
!
        integer, allocatable :: uILB(:)
        integer, allocatable :: uIUB(:)
        integer, allocatable :: uJLB(:)
        integer, allocatable :: uJUB(:)
!
        real(r8), allocatable :: uXmin(:)
        real(r8), allocatable :: uXmax(:)
        real(r8), allocatable :: uYmin(:)
        real(r8), allocatable :: uYmax(:)
!
!  Lower and upper bound ranges for V-type variables for processing
!  state vector and observations.
!
        integer, allocatable :: vILB(:)
        integer, allocatable :: vIUB(:)
        integer, allocatable :: vJLB(:)
        integer, allocatable :: vJUB(:)
!
        real(r8), allocatable :: vXmin(:)
        real(r8), allocatable :: vXmax(:)
        real(r8), allocatable :: vYmin(:)
        real(r8), allocatable :: vYmax(:)
!
!  Switches indicating which variables are written to output files.
!
        logical, allocatable :: Aout(:,:)     ! average file switches
        logical, allocatable :: Dout(:,:)     ! diagnostic file switches
        logical, allocatable :: Hout(:,:)     ! history file switches
        logical, allocatable :: Sout(:,:)     ! station file switches
!
!  Grid identification indices.
!
        integer  :: idXgrd = -1   ! XI-grid position
        integer  :: idYgrd = -2   ! ETA-grid position
        integer  :: idZgrd = -3   ! S-grid position
        integer  :: iddpth = -4   ! depth
        integer  :: idglon = -5   ! longitude
        integer  :: idglat = -6   ! latitude
!
!  Input/output identification indices.
!
        integer  :: idAice        ! fraction of cell covered by ice
        integer  :: idismr        ! iceshelf melt rate
        integer  :: iddraft       ! iceshelf draft
        integer  :: idbath        ! bathymetry
        integer  :: idCfra        ! cloud fraction
        integer  :: idCosW        ! COS(omega(k)*t)
        integer  :: idCos2        ! COS(omega(k)*t)*COS(omega(l)*t)
        integer  :: idDano        ! density anomaly
        integer  :: idDiff(2)     ! temperature and salinity diffusion
        integer  :: iddQdT        ! heat flux sensitivity to SST
        integer  :: idEmPf        ! E-P from bulk_flux.F
        integer  :: idevap        ! evaporation rate
        integer  :: idFsur        ! free-surface
        integer  :: idFsuD        ! detided free-surface
        integer  :: idFsuH        ! free-surface tide harmonics
        integer  :: idGhat(2)     ! KPP nonlocal transport
        integer  :: idHbbl        ! depth of bottom boundary layer
        integer  :: idHice        ! depth of ice cover
        integer  :: idHsbl        ! depth of surface boundary layer
        integer  :: idHsno        ! depth of snow cover
        integer  :: idKhor        ! convolution horizontal diffusion
        integer  :: idKver        ! convolution vertical diffusion
        integer  :: idLdwn        ! downwelling longwave radiation flux
        integer  :: idLhea        ! net latent heat flux
        integer  :: idLrad        ! net longwave radiation flux
        integer  :: idMadH        ! ADM interpolation weights
        integer  :: idMOMi        ! Initial model-observation misfit
        integer  :: idMOMf        ! final model-observation misfit
        integer  :: idMtke        ! turbulent kinetic energy
        integer  :: idMtls        ! turbulent length scale
        integer  :: idNLmi        ! initial NLM at observation locations
        integer  :: idNLmo        ! NLM at observations locations
        integer  :: idNobs        ! number of observations
        integer  :: idObsD        ! observations depth
        integer  :: idObsS        ! observations screening scale
        integer  :: idObsT        ! observations time
        integer  :: idObsX        ! observations X-grid location
        integer  :: idObsY        ! observations Y-grid location
        integer  :: idObsZ        ! observations Z-grid location
        integer  :: idOday        ! observations survey time
        integer  :: idOerr        ! observations error
        integer  :: idOtyp        ! observations type
        integer  :: idOval        ! observations value
        integer  :: idOvar        ! observations global variance
        integer  :: idOvel        ! omega vertical velocity
        integer  :: idQair        ! surface air humidity
        integer  :: idPair        ! surface air pressure
        integer  :: idPbar        ! streamfunction
        integer  :: idRdir        ! river runoff direction
        integer  :: idRepo        ! river runoff ETA-positions
        integer  :: idRflg        ! river runoff flag
        integer  :: idRtra        ! river runoff mass transport
        integer  :: idRuct        ! RHS of U-momentum coupling term
        integer  :: idRu2d        ! RHS of 2D U-momentum
        integer  :: idRu3d        ! RHS of total U-momentum
        integer  :: idRvct        ! RHS of V-momentum coupling term
        integer  :: idRv2d        ! RHS of 2D V-momentum
        integer  :: idRv3d        ! RHS of total V-momentum
        integer  :: idRxpo        ! river runoff XI-positions
        integer  :: idRvsh        ! river runoff transport profile
        integer  :: idRwet        ! wet/dry mask on RHO-points
        integer  :: idRzet        ! RHS of free-surface
        integer  :: idrain        ! rainfall rate
        integer  :: idragL        ! bottom linear drag coefficient
        integer  :: idragQ        ! bottom quadratic drag coefficient
        integer  :: idSdif        ! vertical S-diffusion coefficient
        integer  :: idSinW        ! SIN(omega(k)*t)
        integer  :: idSin2        ! SIN(omega(k)*t)*SIN(omega(l)*t)
        integer  :: idSrad        ! net shortwave radiation flux
        integer  :: idSSHc        ! SSH climatology
        integer  :: idSSHe        ! SSH error variance
        integer  :: idSSHo        ! SSH observations
        integer  :: idSSSc        ! SSS climatology
        integer  :: idSSTc        ! SST climatology
        integer  :: idSSTe        ! SST error variance
        integer  :: idSSTo        ! SST observations
        integer  :: idShea        ! net sensible heat flux
        integer  :: idSWCW        ! SIN(omega(k)*t)*COS(omega(l)*t)
        integer  :: idsfwf        ! surface freswater flux
        integer  :: idTLmo        ! TLM at observation locations
        integer  :: idTair        ! surface air temperature
        integer  :: idTdif        ! vertical T-diffusion coefficient
        integer  :: idTice        ! temperature of ice surface
        integer  :: idtime        ! ocean time
        integer  :: idTper        ! tidal period
        integer  :: idTvan        ! tidal current angle
        integer  :: idTvma        ! maximum tidal current
        integer  :: idTvmi        ! minimum tidal current
        integer  :: idTvph        ! tidal current phase
        integer  :: idTzam        ! tidal elevation amplitude
        integer  :: idTzph        ! tidal elevation phase
        integer  :: idu2dA        ! accumulated 2D U-velocity
        integer  :: idU2rs        ! 2D total U-radiation stress
        integer  :: idU3rs        ! 3D total U-radiation stress
        integer  :: idU2Sd        ! 2D U-Stokes drift velocity
        integer  :: idU3Sd        ! 3D U-Stokes drift velocity
        integer  :: idUads        ! 3D U-velocity adjoint sensitivity
        integer  :: idUair        ! surface U-wind
        integer  :: idUbar        ! 2D U-velocity
        integer  :: idUbas        ! 2D U-velocity adjoint sensitivity
        integer  :: idUbcl        ! 2D U-velocity climatology
        integer  :: idUbcs        ! bottom max U-momentum-wave stress
        integer  :: idUbed        ! bed load U-direction
        integer  :: idUbms        ! bottom U-momentum stress
        integer  :: idUbot        ! bed wave orbital U-velocity
        integer  :: idUbrs        ! bottom U-current stress
        integer  :: idUbtf        ! 2D U-velocity impulse forcing
        integer  :: idUbur        ! bottom U-velocity above bed
        integer  :: idUbws        ! bottom U-wave stress
        integer  :: idUclm        ! 3D U-velocity climatology
        integer  :: idUfx1        ! time averaged U-flux for 2D
        integer  :: idUfx2        ! time averaged U-flux for 3D
        integer  :: idUice        ! ice U-velocity
        integer  :: idUobs        ! 3D U-velocity observations
        integer  :: idUsms        ! surface U-momentum stress
        integer  :: idUsur        ! surface U-velocity observations
        integer  :: idUtlf        ! 3D U-velocity impulse forcing
        integer  :: idUVer        ! 3D velocity error variance
        integer  :: idUVse        ! surface velocity error variance
        integer  :: idUvel        ! 3D U-velocity
        integer  :: idUwet        ! wet/dry mask on U-points
        integer  :: idu2dD        ! detided 2D U-velocity
        integer  :: idu2dH        ! 2D U-velocity tide harmonics
        integer  :: idu2dE        ! 2D eastward velocity at RHO-points
        integer  :: idu3dD        ! detided 3D U-velocity
        integer  :: idu3dH        ! 3D U-velocity tide harmonics
        integer  :: idu3dE        ! 3D eastward velocity at RHO-points
        integer  :: idV2rs        ! 2D total V-radiation stress
        integer  :: idV3rs        ! 3D total V-radiation stress
        integer  :: idV2Sd        ! 2D U-Stokes drift velocity
        integer  :: idV3Sd        ! 3D U-Stokes drift velocity
        integer  :: idVads        ! 3D V-velocity adjoint sensitivity
        integer  :: idVair        ! surface V-wind
        integer  :: idVbar        ! 2D V-velocity
        integer  :: idVbas        ! 2D V-velocity adjoint sensitivity
        integer  :: idVbcl        ! 2D V-velocity climatology
        integer  :: idVbcs        ! bottom max V-current-wave stress
        integer  :: idVbed        ! bed load V-direction
        integer  :: idVbms        ! bottom V-momentum stress
        integer  :: idVbot        ! bed wave orbital V-velocity
        integer  :: idVbrs        ! bottom V-current stress
        integer  :: idVbtf        ! 2D V-velocity impulse forcing
        integer  :: idVbvr        ! bottom V-velocity above bed
        integer  :: idVbws        ! bottom V-wave stress
        integer  :: idVclm        ! 3D V-velocity climatology
        integer  :: idVfx1        ! 2D momentum time-averaged V-flux
        integer  :: idVfx2        ! 3D momentum time-averaged V-flux
        integer  :: idVice        ! ice V-velocity
        integer  :: idVmLS        ! vertical mixing length scale
        integer  :: idVmKK        ! Kinetic energy vertical mixing
        integer  :: idVmKP        ! Length scale vertical mixing
        integer  :: idVobs        ! 3D V-velocity observations
        integer  :: idVsms        ! surface V-momentum stress
        integer  :: idVsur        ! surface V-velocity observations
        integer  :: idVtlf        ! 3D V-velocity impulse forcing
        integer  :: idVvel        ! 3D V-velocity
        integer  :: idVvis        ! vertical viscosity coefficient
        integer  :: idVwet        ! wet/dry mask on V-points
        integer  :: idv2dD        ! detided 2D U-velocity
        integer  :: idv2dH        ! 2D U-velocity tide harmonics
        integer  :: idv2dN        ! 2D northward velocity at RHO-points
        integer  :: idv3dD        ! detided 3D U-velocity
        integer  :: idv3dH        ! 3D U-velocity tide harmonics
        integer  :: idv3dN        ! 3D northward velocity at RHO-points
        integer  :: idW2xx        ! 2D radiation stress, Sxx-component
        integer  :: idW2xy        ! 2D radiation stress, Sxy-component
        integer  :: idW2yy        ! 2D radiation stress, Syy-component
        integer  :: idW3xx        ! 3D radiation stress, Sxx-component
        integer  :: idW3xy        ! 3D radiation stress, Sxy-component
        integer  :: idW3yy        ! 3D radiation stress, Syy-component
        integer  :: idW3zx        ! 3D radiation stress, Szx-component
        integer  :: idW3zy        ! 3D radiation stress, Szy-component
        integer  :: idWamp        ! wind-induced wave amplitude
        integer  :: idWbrk        ! wind-induced wave breaking
        integer  :: idWdis        ! wind-induced wave dissipation
        integer  :: idWdir        ! wind-induced wave direction
        integer  :: idWlen        ! wind-induced wave length
        integer  :: idWptp        ! wind-induced surface wave period
        integer  :: idWpbt        ! wind-induced bottom wave period
        integer  :: idWorb        ! wind-induced wave orbital velocity
        integer  :: idWvel        ! true vertical velocity
        integer  :: idZoBL        ! bottom roughness length
        integer  :: idZads        ! Free-surface adjoint sensitivity
        integer  :: idZtlf        ! Free-surface impulse forcing
        integer  :: id2dPV        ! 2D potential vorticity
        integer  :: id2dRV        ! 2D relative vorticity
        integer  :: id3dPV        ! 3D potential vorticity
        integer  :: id3dRV        ! 3D relative vorticity
!
!  Input/output identification tracer indices.
!
        integer, allocatable :: idRtrc(:)    ! river runoff for tracers
        integer, allocatable :: idTads(:)    ! tracers adjoint sentivity
        integer, allocatable :: idTbot(:)    ! bottom flux for tracers
        integer, allocatable :: idTbry(:,:)  ! tracers boundary
        integer, allocatable :: idTclm(:)    ! tracers climatology
        integer, allocatable :: idTerr(:)    ! tracers error variance
        integer, allocatable :: idTobs(:)    ! tracers observations
        integer, allocatable :: idTsur(:)    ! surface flux for tracers
        integer, allocatable :: idTtlf(:)    ! tracers impulse forcing
!
!  Boundary conditions identification indices.
!
        integer  :: idU2bc(4)      ! 2D U-velocity boundary conditions
        integer  :: idU3bc(4)      ! 3D U-velocity boundary conditions
        integer  :: idV2bc(4)      ! 2D V-velocity boundary conditions
        integer  :: idV3bc(4)      ! 3D V-velocity boundary conditions
        integer  :: idZbry(4)      ! free-surface boundary conditions
!
!  Time-averaged quadratic terms IDs.
!
        integer  :: idU2av                    ! <ubar*ubar>
        integer  :: idV2av                    ! <vbar*vbar>
        integer  :: idZZav                    ! <zeta*zeta>
        integer  :: idHUav                    ! <Huon>
        integer  :: idHVav                    ! <Hvom>
        integer  :: idUUav                    ! <u*u>
        integer  :: idUVav                    ! <u*v>
        integer  :: idVVav                    ! <v*v>
        integer, allocatable :: iHUTav(:)     ! <Huon*t> for active tracers
        integer, allocatable :: iHVTav(:)     ! <Hvom*t> for active tracers
        integer, allocatable :: idTTav(:)     ! <t*t> for active tracers
        integer, allocatable :: idUTav(:)     ! <u*t> for active tracers
        integer, allocatable :: idVTav(:)     ! <v*t> for active tracers
!
!  State variables indices (order is important).
!
        integer  :: isFsur = 1                ! free-surface
        integer  :: isUbar = 2                ! 2D U-velocity
        integer  :: isVbar = 3                ! 2D V-velocity
        integer  :: isUvel = 4                ! 3D U-velocity
        integer  :: isVvel = 5                ! 3D V-velocity
        integer  :: isUstr                    ! surface u-stress
        integer  :: isVstr                    ! surface v-stress
        integer  :: isMtke                    ! turbulent kinetic energy
        integer, allocatable :: isTsur(:)     ! surface tracer flux
        integer, allocatable :: isTvar(:)     ! tracers
        integer, allocatable :: idBvar(:)     ! LBC variables indices
        integer, allocatable :: idSvar(:)     ! state vector indices
        integer, allocatable :: idSbry(:)     ! state boundaries indices
!
!  Generic state variables lateral boundary indices.
!
        integer  :: isBp2d                    ! 2D PSI-variables
        integer  :: isBr2d                    ! 2D RHO-variables
        integer  :: isBu2d                    ! 2D U-variables
        integer  :: isBv2d                    ! 2D V-variables
        integer  :: isBp3d                    ! 3D PSI-variables
        integer  :: isBr3d                    ! 3D RHO-variables
        integer  :: isBu3d                    ! 3D U-variables
        integer  :: isBv3d                    ! 3D V-variables
        integer  :: isBw3d                    ! 3D W-variables
!
!  Input forcing NetCDF files IDs.
!
        integer, allocatable :: ncFRCid(:,:)
!
!  Flags to create output files.
!
        integer, allocatable :: idefADJ(:)    ! adjoint file
        integer, allocatable :: idefAVG(:)    ! averages file
        integer, allocatable :: idefDIA(:)    ! diagnostics file
        integer, allocatable :: idefHIS(:)    ! history file
        integer, allocatable :: idefTLM(:)    ! tangent file
!
!  Output NetCDF variables IDs.
!
        integer, allocatable :: idTvar(:)     ! tracers variables
        integer, allocatable :: idTrcD(:)     ! detided tracer variables
        integer, allocatable :: idTrcH(:)     ! tracer variables hamonics
!
!  Input/Output information variables.
!
        logical, allocatable :: Linfo(:,:,:)
        integer, allocatable :: Iinfo(:,:,:)
        real(r8), allocatable :: Finfo(:,:,:)
        real(r8), allocatable :: Fpoint(:,:,:)
        real(r8), allocatable :: Fscale(:,:)
        real(r8), allocatable :: Tintrp(:,:,:)
        real(r8), allocatable :: Vtime(:,:,:)
        character (len=5  ) :: version = '3.7  '
        character (len=40 ) :: varnam(MV)
        character (len=44 ) :: date_str
        character (len=46 ) :: Tname(0:NV)
        character (len=100) :: Vname(5,0:NV)
        character (len=120) :: history
        character (len=256), allocatable :: Cinfo(:,:)
!
!  Source code root directory, cpp header file and directory, and
!  analytical expression directory.
!
        character (len=256) :: Rdir
        character (len=256) :: Hdir
        character (len=256) :: Hfile
        character (len=256) :: Adir
!
!  Analyical header file logical and names.
!
        logical :: Lanafile
        character (len=256), dimension(37) :: ANANAME
!
!  SVN revision and repository root URL.
!
        character (len=40 ) :: svn_rev
        character (len=256) :: svn_url
!
      CONTAINS
!
      SUBROUTINE allocate_ncparam
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate variables.
!-----------------------------------------------------------------------
!
      allocate ( Nxyp(Ngrids) )
      allocate ( Nxyr(Ngrids) )
      allocate ( Nxyu(Ngrids) )
      allocate ( Nxyv(Ngrids) )
      allocate ( NwaterR(Ngrids) )
      allocate ( NwaterU(Ngrids) )
      allocate ( NwaterV(Ngrids) )
      allocate ( rILB(Ngrids) )
      allocate ( rIUB(Ngrids) )
      allocate ( rJLB(Ngrids) )
      allocate ( rJUB(Ngrids) )
      allocate ( rXmin(Ngrids) )
      allocate ( rXmax(Ngrids) )
      allocate ( rYmin(Ngrids) )
      allocate ( rYmax(Ngrids) )
      allocate ( uILB(Ngrids) )
      allocate ( uIUB(Ngrids) )
      allocate ( uJLB(Ngrids) )
      allocate ( uJUB(Ngrids) )
      allocate ( uXmin(Ngrids) )
      allocate ( uXmax(Ngrids) )
      allocate ( uYmin(Ngrids) )
      allocate ( uYmax(Ngrids) )
      allocate ( vILB(Ngrids) )
      allocate ( vIUB(Ngrids) )
      allocate ( vJLB(Ngrids) )
      allocate ( vJUB(Ngrids) )
      allocate ( vXmin(Ngrids) )
      allocate ( vXmax(Ngrids) )
      allocate ( vYmin(Ngrids) )
      allocate ( vYmax(Ngrids) )
      allocate ( Aout(NV,Ngrids) )
      allocate ( Dout(NV,Ngrids) )
      allocate ( Hout(NV,Ngrids) )
      allocate ( Sout(NV,Ngrids) )
      allocate ( idRtrc(MT) )
      allocate ( idTads(MT) )
      allocate ( idTbot(MT) )
      allocate ( idTbry(4,MT) )
      allocate ( idTclm(MT) )
      allocate ( idTerr(MT) )
      allocate ( idTobs(MT) )
      allocate ( idTsur(MT) )
      allocate ( idTtlf(MT) )
      allocate ( iHUTav(MT) )
      allocate ( iHVTav(MT) )
      allocate ( idTTav(MT) )
      allocate ( idUTav(MT) )
      allocate ( idVTav(MT) )
      allocate ( isTsur(MT) )
      allocate ( isTvar(MT) )
      allocate ( idBvar(nLBCvar) )
      allocate ( idSvar(MAXVAL(NSV)+1) )
      allocate ( idSbry(MAXVAL(NSV)+1) )
      allocate ( ncFRCid(NV,Ngrids) )
      allocate ( idefADJ(Ngrids) )
      allocate ( idefAVG(Ngrids) )
      allocate ( idefDIA(Ngrids) )
      allocate ( idefHIS(Ngrids) )
      allocate ( idefTLM(Ngrids) )
      allocate ( idTvar(MT) )
      allocate ( idTrcD(NAT) )
      allocate ( idTrcH(NAT) )
      allocate ( Linfo(4,NV,Ngrids) )
      allocate ( Iinfo(10,NV,Ngrids) )
      allocate ( Finfo(10,NV,Ngrids) )
      allocate ( Fpoint(2,NV,Ngrids) )
      allocate ( Fscale(NV,Ngrids) )
      allocate ( Tintrp(2,NV,Ngrids) )
      allocate ( Vtime(MAX(2,MTC),NV,Ngrids) )
      allocate ( Cinfo(NV,Ngrids) )
      RETURN
      END SUBROUTINE allocate_ncparam
      SUBROUTINE initialize_ncparam
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_ncparam" for all nested grids.                                 !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Local variable declarations.
!
      logical :: load
      integer, parameter :: inp = 10
      integer :: Lvar, Ntiles, i, ic, ie, is, j, ng
      integer :: gtype, tile, varid
      real(r8), parameter :: spv = 0.0_r8
      real(r8) :: offset, scale
      character (len=120), dimension(7) :: Vinfo
!
!-----------------------------------------------------------------------
!  Initialize several variables.
!-----------------------------------------------------------------------
!
!  Initialize DOMAIN structure.
!
      DO ng=1,Ngrids
        DOMAIN(ng) % Eastern_Edge  = .FALSE.
        DOMAIN(ng) % Western_Edge  = .FALSE.
        DOMAIN(ng) % Northern_Edge = .FALSE.
        DOMAIN(ng) % Southern_Edge = .FALSE.
        DOMAIN(ng) % NorthEast_Corner = .FALSE.
        DOMAIN(ng) % NorthWest_Corner = .FALSE.
        DOMAIN(ng) % SouthEast_Corner = .FALSE.
        DOMAIN(ng) % SouthWest_Corner = .FALSE.
        DOMAIN(ng) % NorthEast_Test = .FALSE.
        DOMAIN(ng) % NorthWest_Test = .FALSE.
        DOMAIN(ng) % SouthEast_Test = .FALSE.
        DOMAIN(ng) % SouthWest_Test = .FALSE.
        DOMAIN(ng) % Xmin_psi = spv
        DOMAIN(ng) % Xmax_psi = spv
        DOMAIN(ng) % Ymin_psi = spv
        DOMAIN(ng) % Ymax_psi = spv
        DOMAIN(ng) % Xmin_rho = spv
        DOMAIN(ng) % Xmax_rho = spv
        DOMAIN(ng) % Ymin_rho = spv
        DOMAIN(ng) % Ymax_rho = spv
        DOMAIN(ng) % Xmin_u   = spv
        DOMAIN(ng) % Xmax_u   = spv
        DOMAIN(ng) % Ymin_u   = spv
        DOMAIN(ng) % Ymax_u   = spv
        DOMAIN(ng) % Xmin_v   = spv
        DOMAIN(ng) % Xmax_v   = spv
        DOMAIN(ng) % Ymin_v   = spv
        DOMAIN(ng) % Ymax_v   = spv
      END DO
!
!  Initialize NetCDF files creation flags.
!
      DO ng=1,Ngrids
        idefADJ(ng)=-1
        idefAVG(ng)=-1
        idefDIA(ng)=-1
        idefHIS(ng)=-1
        idefTLM(ng)=-1
      END DO
!
!  Analytical files switch and names.
!
      Lanafile=.TRUE.
      DO i=1,37
        DO j=1,LEN(ANANAME(1))
          ANANAME(i)(j:j)=' '
        END DO
      END DO
!
!  Set IDs for state some state variables.
!
      ic=5
      DO i=1,MT
        ic=ic+1
        isTvar(i)=ic
      END DO
!
!  Set generic lateral boundary indices for LBC structure.  Use the same
!  values of the state variables at the same C-grid location.  Generic
!  indices are used for testing periodicity.  The PSI-variables and
!  W-variables are assign the same value as the RHO-variables.
!
      isBp2d=isFsur                           ! 2D PSI-variables
      isBr2d=isFsur                           ! 2D RHO-variables
      isBu2d=isUbar                           ! 2D U-variables
      isBv2d=isVbar                           ! 2D V-variables
      isBp3d=isTvar(1)                        ! 3D PSI-variables
      isBr3d=isTvar(1)                        ! 3D RHO-variables
      isBu3d=isUvel                           ! 3D U-variables
      isBv3d=isVvel                           ! 3D V-variables
      isBw3d=isTvar(1)                        ! 3D W-variables
!
!  Initialize IO information variables.
!
      DO ng=1,Ngrids
        DO i=1,NV
          Linfo(1,i,ng)=.FALSE.
          Linfo(2,i,ng)=.FALSE.
          Linfo(3,i,ng)=.FALSE.
          Linfo(4,i,ng)=.FALSE.
          Aout(i,ng)=.FALSE.
          Dout(i,ng)=.FALSE.
          Hout(i,ng)=.FALSE.
          Sout(i,ng)=.FALSE.
          Iinfo(1,i,ng)=0
          Iinfo(2,i,ng)=-1
          Iinfo(3,i,ng)=-1
          Iinfo(4,i,ng)=0
          Iinfo(5,i,ng)=0
          Iinfo(6,i,ng)=0
          Iinfo(7,i,ng)=0
          Iinfo(8,i,ng)=0
          Iinfo(9,i,ng)=0
          Iinfo(10,i,ng)=0
          Finfo(1,i,ng)=0.0_r8
          Finfo(2,i,ng)=0.0_r8
          Finfo(3,i,ng)=0.0_r8
          Finfo(5,i,ng)=0.0_r8
          Finfo(6,i,ng)=0.0_r8
          Finfo(7,i,ng)=0.0_r8
          Finfo(10,i,ng)=1.0_r8
          Fscale(i,ng)=1.0_r8
          Fpoint(1,i,ng)=0.0_r8
          Fpoint(2,i,ng)=0.0_r8
          Tintrp(1,i,ng)=0.0_r8
          Tintrp(2,i,ng)=0.0_r8
          Vtime(1,i,ng)=0.0_r8
          Vtime(2,i,ng)=0.0_r8
          ncFRCid(i,ng)=-1
        END DO
      END DO
!
!  Set source code root directory, cpp header file and directory, and
!  analytical expression directory.
!
      Rdir="/home/bkgalton/data/hub/romsIceShelf"
      Hdir="/home/bkgalton/data/hub/romsIceShelf/ROMS/Include"
      Hfile="iceshelf_tides.h"
      Adir="/home/bkgalton/data/hub/romsIceShelf/ROMS/Functionals"
!
!-----------------------------------------------------------------------
!  Define names of variables for Input/Output NetCDF files.
!-----------------------------------------------------------------------
!
!  Open input variable information file.
!
      OPEN (inp, FILE=TRIM(varname), FORM='formatted', STATUS='old',    &
     &      ERR=10)
      GOTO 20
  10  IF (Master) WRITE(stdout,50) TRIM(varname)
      STOP
  20  CONTINUE
!
!  Read in variable information.  Ignore blank and comment [char(33)=!]
!  input lines.
!
      varid=0
      DO WHILE (.TRUE.)
        READ (inp,*,ERR=30,END=40) Vinfo(1)
        Lvar=LEN_TRIM(Vinfo(1))
!
!  Extract SVN Repository Root URL.
!
        IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).eq.CHAR(36))) THEN
          is=INDEX(Vinfo(1),'https')
          ie=INDEX(Vinfo(1),'/ROMS')-1
          IF (ie.gt.is) THEN
            svn_url=Vinfo(1)(is:ie)
          ELSE
            svn_url='https:://myroms.org/svn/src'
          END IF
          svn_rev="Unversioned directory"
!
!  Read in other variable information.
!
        ELSE IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).ne.CHAR(33))) THEN
          READ (inp,*,ERR=30) Vinfo(2)
          READ (inp,*,ERR=30) Vinfo(3)
          READ (inp,*,ERR=30) Vinfo(4)
          READ (inp,*,ERR=30) Vinfo(5)
          READ (inp,*,ERR=30) Vinfo(6)
          READ (inp,*,ERR=30) Vinfo(7)
          READ (inp,*,ERR=30) scale
!
!  Determine staggered C-grid variable.
!
          SELECT CASE (TRIM(ADJUSTL(Vinfo(7))))
            CASE ('p2dvar')
              gtype=p2dvar
            CASE ('r2dvar')
              gtype=r2dvar
            CASE ('u2dvar')
              gtype=u2dvar
            CASE ('v2dvar')
              gtype=v2dvar
            CASE ('p3dvar')
              gtype=p3dvar
            CASE ('r3dvar')
              gtype=r3dvar
            CASE ('u3dvar')
              gtype=u3dvar
            CASE ('v3dvar')
              gtype=v3dvar
            CASE ('w3dvar')
              gtype=w3dvar
            CASE ('b3dvar')
              gtype=b3dvar
            CASE DEFAULT
              gtype=0
          END SELECT
!
!  Assign identification indices.
!
          load=.TRUE.
          varid=varid+1
          SELECT CASE (TRIM(ADJUSTL(Vinfo(6))))
            CASE ('idtime')
              idtime=varid
            CASE ('idismr')
              idismr=varid
            CASE ('iddraft')
              iddraft=varid
            CASE ('idbath')
              idbath=varid
            CASE ('idFsur')
              idFsur=varid
            CASE ('idRzet')
              idRzet=varid
            CASE ('idUbar')
              idUbar=varid
            CASE ('idu2dE')
              idu2dE=varid
            CASE ('idRu2d')
              idRu2d=varid
            CASE ('idVbar')
              idVbar=varid
            CASE ('idv2dN')
              idv2dN=varid
            CASE ('idRv2d')
              idRv2d=varid
            CASE ('idUvel')
              idUvel=varid
            CASE ('idu3dE')
              idu3dE=varid
            CASE ('idRu3d')
              idRu3d=varid
            CASE ('idVvel')
              idVvel=varid
            CASE ('idv3dN')
              idv3dN=varid
            CASE ('idRv3d')
              idRv3d=varid
            CASE ('idWvel')
              idWvel=varid
            CASE ('idOvel')
              idOvel=varid
            CASE ('idDano')
              idDano=varid
            CASE ('idTvar(itemp)')
              idTvar(itemp)=varid
            CASE ('idTvar(isalt)')
              idTvar(isalt)=varid
            CASE ('idUsms')
              idUsms=varid
            CASE ('idVsms')
              idVsms=varid
            CASE ('idUbms')
              idUbms=varid
            CASE ('idVbms')
              idVbms=varid
            CASE ('idUbws')
              idUbws=varid
            CASE ('idUbcs')
              idUbcs=varid
            CASE ('idVbws')
              idVbws=varid
            CASE ('idVbcs')
              idVbcs=varid
            CASE ('idTsur(itemp)')
              idTsur(itemp)=varid
            CASE ('iddQdT')
              iddQdT=varid
            CASE ('idsfwf')
              idsfwf=varid
            CASE ('idTsur(isalt)')
              idTsur(isalt)=varid
            CASE ('idTbot(itemp)')
              idTbot(itemp)=varid
            CASE ('idTbot(isalt)')
              idTbot(isalt)=varid
            CASE ('idGhat(itemp)')
              idGhat(itemp)=varid
            CASE ('idGhat(isalt)')
              idGhat(isalt)=varid
            CASE ('idMtke')
              idMtke=varid
            CASE ('idMtls')
              idMtls=varid
            CASE ('idVvis')
              idVvis=varid
            CASE ('idTdif')
              idTdif=varid
              idDiff(itemp)=idTdif
            CASE ('idSdif')
              idSdif=varid
              idDiff(isalt)=idSdif
            CASE ('idVmLS')
              idVmLS=varid
            CASE ('idVmKK')
              idVmKK=varid
            CASE ('idVmKP')
              idVmKP=varid
            CASE ('idZbry(iwest)')
              idZbry(iwest)=varid
            CASE ('idZbry(ieast)')
              idZbry(ieast)=varid
            CASE ('idZbry(isouth)')
              idZbry(isouth)=varid
            CASE ('idZbry(inorth)')
              idZbry(inorth)=varid
            CASE ('idU2bc(iwest)')
              idU2bc(iwest)=varid
            CASE ('idU2bc(ieast)')
              idU2bc(ieast)=varid
            CASE ('idU2bc(isouth)')
              idU2bc(isouth)=varid
            CASE ('idU2bc(inorth)')
              idU2bc(inorth)=varid
            CASE ('idV2bc(iwest)')
              idV2bc(iwest)=varid
            CASE ('idV2bc(ieast)')
              idV2bc(ieast)=varid
            CASE ('idV2bc(isouth)')
              idV2bc(isouth)=varid
            CASE ('idV2bc(inorth)')
              idV2bc(inorth)=varid
            CASE ('idU3bc(iwest)')
              idU3bc(iwest)=varid
            CASE ('idU3bc(ieast)')
              idU3bc(ieast)=varid
            CASE ('idU3bc(isouth)')
              idU3bc(isouth)=varid
            CASE ('idU3bc(inorth)')
              idU3bc(inorth)=varid
            CASE ('idV3bc(iwest)')
              idV3bc(iwest)=varid
            CASE ('idV3bc(ieast)')
              idV3bc(ieast)=varid
            CASE ('idV3bc(isouth)')
              idV3bc(isouth)=varid
            CASE ('idV3bc(inorth)')
              idV3bc(inorth)=varid
            CASE ('idTbry(iwest,itemp)')
              idTbry(iwest,itemp)=varid
            CASE ('idTbry(ieast,itemp)')
              idTbry(ieast,itemp)=varid
            CASE ('idTbry(isouth,itemp)')
              idTbry(isouth,itemp)=varid
            CASE ('idTbry(inorth,itemp)')
              idTbry(inorth,itemp)=varid
            CASE ('idTbry(iwest,isalt)')
              idTbry(iwest,isalt)=varid
            CASE ('idTbry(ieast,isalt)')
              idTbry(ieast,isalt)=varid
            CASE ('idTbry(isouth,isalt)')
              idTbry(isouth,isalt)=varid
            CASE ('idTbry(inorth,isalt)')
              idTbry(inorth,isalt)=varid
            CASE ('idRwet')
              idRwet=varid
            CASE ('idUwet')
              idUwet=varid
            CASE ('idVwet')
              idVwet=varid
            CASE ('idPair')
              idPair=varid
            CASE ('idTair')
              idTair=varid
            CASE ('idQair')
              idQair=varid
            CASE ('idCfra')
              idCfra=varid
            CASE ('idSrad')
              idSrad=varid
            CASE ('idLdwn')
              idLdwn=varid
            CASE ('idLrad')
              idLrad=varid
            CASE ('idLhea')
              idLhea=varid
            CASE ('idShea')
              idShea=varid
            CASE ('idrain')
              idrain=varid
            CASE ('idEmPf')
              idEmPf=varid
            CASE ('idevap')
              idevap=varid
            CASE ('idUair')
              idUair=varid
            CASE ('idVair')
              idVair=varid
            CASE ('idWamp')
              idWamp=varid
            CASE ('idWbrk')
              idWbrk=varid
            CASE ('idWdis')
              idWdis=varid
            CASE ('idWdir')
              idWdir=varid
            CASE ('idWlen')
              idWlen=varid
            CASE ('idWptp')
              idWptp=varid
            CASE ('idWpbt')
              idWpbt=varid
            CASE ('idWorb')
              idWorb=varid
            CASE ('idW2xx')
              idW2xx=varid
            CASE ('idW2xy')
              idW2xy=varid
            CASE ('idW2yy')
              idW2yy=varid
            CASE ('idW3xx')
              idW3xx=varid
            CASE ('idW3xy')
              idW3xy=varid
            CASE ('idW3yy')
              idW3yy=varid
            CASE ('idW3zx')
              idW3zx=varid
            CASE ('idW3zy')
              idW3zy=varid
            CASE ('idU2rs')
              idU2rs=varid
            CASE ('idV2rs')
              idV2rs=varid
            CASE ('idU2Sd')
              idU2Sd=varid
            CASE ('idV2Sd')
              idV2Sd=varid
            CASE ('idU3rs')
              idU3rs=varid
            CASE ('idV3rs')
              idV3rs=varid
            CASE ('idU3Sd')
              idU3Sd=varid
            CASE ('idV3Sd')
              idV3Sd=varid
            CASE ('idTper')
              idTper=varid
            CASE ('idTzam')
              idTzam=varid
            CASE ('idTzph')
              idTzph=varid
            CASE ('idTvph')
              idTvph=varid
            CASE ('idTvan')
              idTvan=varid
            CASE ('idTvma')
              idTvma=varid
            CASE ('idTvmi')
              idTvmi=varid
            CASE ('idRxpo')
              idRxpo=varid
            CASE ('idRepo')
              idRepo=varid
            CASE ('idRdir')
              idRdir=varid
            CASE ('idRvsh')
              idRvsh=varid
            CASE ('idRtra')
              idRtra=varid
            CASE ('idRflg')
              idRflg=varid
            CASE ('idRtrc(itemp)')
              idRtrc(itemp)=varid
            CASE ('idRtrc(isalt)')
              idRtrc(isalt)=varid
            CASE ('idHsbl')
              idHsbl=varid
            CASE ('idHbbl')
              idHbbl=varid
            CASE ('idUbot')
              idUbot=varid
            CASE ('idVbot')
              idVbot=varid
            CASE ('idUbur')
              idUbur=varid
            CASE ('idVbvr')
              idVbvr=varid
            CASE ('idUbrs')
              idUbrs=varid
            CASE ('idVbrs')
              idVbrs=varid
            CASE ('idSSHc')
              idSSHc=varid
            CASE ('idUbcl')
              idUbcl=varid
            CASE ('idVbcl')
              idVbcl=varid
            CASE ('idUclm')
              idUclm=varid
            CASE ('idVclm')
              idVclm=varid
            CASE ('idSSSc')
              idSSSc=varid
            CASE ('idSSTc')
              idSSTc=varid
            CASE ('idSSHo')
              idSSHo=varid
            CASE ('idSSHe')
              idSSHe=varid
            CASE ('idUobs')
              idUobs=varid
            CASE ('idVobs')
              idVobs=varid
            CASE ('idUVer')
              idUVer=varid
            CASE ('idUsur')
              idUsur=varid
            CASE ('idVsur')
              idVsur=varid
            CASE ('idUVse')
              idUVse=varid
            CASE ('idSSTo')
              idSSTo=varid
            CASE ('idSSTe')
              idSSTe=varid
            CASE ('idTobs(itemp)')
              idTobs(itemp)=varid
            CASE ('idTerr(itemp)')
              idTerr(itemp)=varid
            CASE ('idTobs(isalt)')
              idTobs(isalt)=varid
            CASE ('idTerr(isalt)')
              idTerr(isalt)=varid
            CASE ('idU2av')
              idU2av=varid
            CASE ('idV2av')
              idV2av=varid
            CASE ('idZZav')
              idZZav=varid
            CASE ('idTTav(itrc)')
              load=.TRUE.
            CASE ('iHUTav(itrc)')
              load=.TRUE.
            CASE ('iHVTav(itrc)')
              load=.TRUE.
            CASE ('idUTav(itrc)')
              load=.TRUE.
            CASE ('idVTav(itrc)')
              load=.TRUE.
            CASE ('idHUav')
              idHUav=varid
            CASE ('idHVav')
              idHVav=varid
            CASE ('idUUav')
              idUUav=varid
            CASE ('idUVav')
              idUVav=varid
            CASE ('idVVav')
              idVVav=varid
            CASE ('idTTav(isalt)')
              idTTav(isalt)=varid
            CASE ('iHUTav(isalt)')
              iHUTav(isalt)=varid
            CASE ('iHVTav(isalt)')
              iHVTav(isalt)=varid
            CASE ('idUTav(isalt)')
              idUTav(isalt)=varid
            CASE ('idVTav(isalt)')
              idVTav(isalt)=varid
            CASE ('id2dPV')
              id2dPV=varid
            CASE ('id2dRV')
              id2dRV=varid
            CASE ('id3dPV')
              id3dPV=varid
            CASE ('id3dRV')
              id3dRV=varid
            CASE DEFAULT
              load=.FALSE.
          END SELECT
!
!  Load variable data into information arrays.
!
          IF (load) THEN
            load=.FALSE.
            IF (varid.gt.MV) THEN
              WRITE (stdout,60) MV, varid
              STOP
            END IF
            DO i=1,5
              Vname(i,varid)=TRIM(ADJUSTL(Vinfo(i)))
            END DO
            DO ng=1,Ngrids
              Iinfo(1,varid,ng)=gtype
              Fscale(varid,ng)=scale
            END DO
          ELSE
            varid=varid-1
          END IF
        END IF
      END DO
      GOTO 40
  30  WRITE (stdout,80) TRIM(ADJUSTL(Vinfo(1)))
      STOP
  40  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Set passive tracers surface flux metadata. The variable name is the
!  same as the basic tracer but with the _sflux suffix.
!-----------------------------------------------------------------------
!
      DO i=NAT+1,MT
        varid=varid+1
        IF (varid.gt.MV) THEN
          WRITE (stdout,60) MV, varid
          STOP
        END IF
        idTsur(i)=varid
        DO ng=1,Ngrids
          Fscale(varid,ng)=1.0_r8
          Iinfo(1,varid,ng)=r2dvar
        END DO
        WRITE (Vname(1,varid),'(a,a)')                                  &
     &        TRIM(ADJUSTL(Vname(1,idTvar(i)))), '_sflux'
        WRITE (Vname(2,varid),'(a,a)')                                  &
     &        TRIM(ADJUSTL(Vname(2,idTvar(i)))), ', surface flux'
        WRITE (Vname(3,varid),'(a,1x,a)')                               &
     &        TRIM(ADJUSTL(Vname(3,idTvar(i)))), 'meter second-1'
        WRITE (Vname(4,varid),'(a,a)')                                  &
     &        TRIM(Vname(1,varid)), ', scalar, series'
        WRITE (Vname(5,varid),'(a)')                                    &
     &        TRIM(ADJUSTL(Vname(5,idTvar(i))))
      END DO
!
!-----------------------------------------------------------------------
!  Set model lateral boundary variables index.
!-----------------------------------------------------------------------
!
      idBvar(isFsur)=idFsur
      idBvar(isUbar)=idUbar
      idBvar(isVbar)=idVbar
      ic=3
      idBvar(isUvel)=idUvel
      idBvar(isVvel)=idVvel
      ic=ic+2
      DO i=1,MT
        idBvar(isTvar(i))=idTvar(i)
        ic=ic+1
      END DO
!
!-----------------------------------------------------------------------
!  Set model state variables indices.
!-----------------------------------------------------------------------
!
      idSvar(isFsur)=idFsur
      idSvar(isUbar)=idUbar
      idSvar(isVbar)=idVbar
      ic=3
      idSvar(isUvel)=idUvel
      idSvar(isVvel)=idVvel
      ic=ic+2
      DO i=1,MT
        idSvar(isTvar(i))=idTvar(i)
        ic=ic+1
      END DO
!
  50  FORMAT (/,' MOD_NCPARAM - Unable to open variable information',   &
     &        ' file: ',/,15x,a,/,15x,'Default file is located in',     &
     &        ' source directory.')
  60  FORMAT (/,' MOD_NCPARAM - too small dimension ',                  &
     &        'parameter, MV = ',2i5,/,15x,                             &
     &        'change file  mod_ncparam.F  and recompile.')
  70  FORMAT (/,' MOD_NCPARM - Cannot load information for ',           &
     &        'variable: ',a,/,15x,'Need CASE construct for: ',a)
  80  FORMAT (/,' MOD_NCPARM - Error while reading information ',       &
     &        'for variable: ',a)
  90  FORMAT (a,i2.2)
 100  FORMAT (a,a,i2.2)
 110  FORMAT (a)
 120  FORMAT (a,a)
      RETURN
      END SUBROUTINE initialize_ncparam
      END MODULE mod_ncparam
