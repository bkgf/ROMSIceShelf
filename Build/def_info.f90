      SUBROUTINE def_info (ng, model, ncid, ncname, DimIDs)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer).                           !
!     model    Calling model identifier (integer).                     !
!     ncid     NetCDF file ID (integer).                               !
!     ncname   NetCDF file name (character).                           !
!     DimIDs   NetCDF dimensions IDs (integer vector):                 !
!                DimIDs( 1) => XI-dimension at RHO-points.             !
!                DimIDs( 2) => XI-dimension at U-points.               !
!                DimIDs( 3) => XI-dimension at V-points.               !
!                DimIDs( 4) => XI-dimension at PSI-points.             !
!                DimIDs( 5) => ETA-dimension at RHO-points.            !
!                DimIDs( 6) => ETA-dimension at U-points.              !
!                DimIDs( 7) => ETA-dimension at V-points.              !
!                DimIDs( 8) => ETA-dimension at PSI-points.            !
!                DimIDs( 9) => S-dimension at RHO-points.              !
!                DimIDs(10) => S-dimension at W-points.                !
!                DimIDs(11) => Number of tracers dimension.            !
!                DimIDs(12) => Unlimited time record dimension.        !
!                DimIDs(13) => Number of stations dimension.           !
!                DimIDs(14) => Boundary dimension.                     !
!                DimIDs(15) => Number of floats dimension.             !
!                DimIDs(16) => Number sediment bed layers dimension.   !
!                DimIDs(17) => Dimension 2D water RHO-points.          !
!                DimIDs(18) => Dimension 2D water U-points.            !
!                DimIDs(19) => Dimension 2D water V-points.            !
!                DimIDs(20) => Dimension 3D water RHO-points.          !
!                DimIDs(21) => Dimension 3D water U-points.            !
!                DimIDs(23) => Dimension 3D water W-points.            !
!                DimIDs(24) => Dimension sediment bed water points.    !
!                DimIDs(25) => Number of EcoSim phytoplankton groups.  !
!                DimIDs(26) => Number of EcoSim bateria groups.        !
!                DimIDs(27) => Number of EcoSim DOM groups.            !
!                DimIDs(28) => Number of EcoSim fecal groups.          !
!                DimIDs(29) => Number of state variables.              !
!                DimIDs(30) => Number of 3D variables time levels (2). !
!                DimIDs(31) => Number of 2D variables time levels (3). !
!                DimIDs(32) => Number of sediment tracers.             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings
!
      USE def_var_mod, ONLY : def_var
      USE strings_mod, ONLY : join_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      integer, intent(in) :: DimIDs(32)
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25
      integer :: brydim, i, ie, is, j, lstr, varid
      integer :: srdim, stadim, status, swdim, trcdim, usrdim
      integer :: iceshelfdim
      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: def_dim
      real(r8) :: Aval(6)
      character (len=11 ) :: frcatt
      character (len=50 ) :: tiling
      character (len=80 ) :: type
      character (len=120) :: Vinfo(Natt)
      character (len=1024) :: string
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
      iceshelfdim=DimIDs(5)
      stadim=DimIDs(13)
      brydim=DimIDs(14)
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(ncid.ne.GST(ng)%ncid)) THEN
        status=def_dim(ng, model, ncid, ncname, 'Nuser', 25, usrdim)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Define history global attribute.
!
        IF (LEN_TRIM(date_str).gt.0) THEN
          WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',       &
     &                                      TRIM( version),             &
     &                                      TRIM(date_str)
        ELSE
          WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version',              &
     &                               TRIM(version)
        END IF
!
!  Set tile decomposition global attribute.
!
        WRITE (tiling,10) NtileI(ng), NtileJ(ng)
!
!  Define file name global attribute.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'file',                &
     &                        TRIM(ncname))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define NetCDF format type.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'format',              &
     &                        'netCDF-3 64bit offset file')
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'format', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
        type='CF-1.4'
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'Conventions',         &
     &                        TRIM(type))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'Conventions', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file type global attribute.
!
        IF (ncid.eq.ADM(ng)%ncid) THEN
          type='ROMS/TOMS adjoint history file'
        ELSE IF (ncid.eq.AVG(ng)%ncid) THEN
          type='ROMS/TOMS nonlinear model averages file'
        ELSE IF (ncid.eq.DIA(ng)%ncid) THEN
          type='ROMS/TOMS diagnostics file'
        ELSE IF (ncid.eq.FLT(ng)%ncid) THEN
          type='ROMS/TOMS floats file'
        ELSE IF (ncid.eq.ERR(ng)%ncid) THEN
          type='ROMS/TOMS posterior analysis error covariance matrix'
        ELSE IF (ncid.eq.GST(ng)%ncid) THEN
          type='ROMS/TOMS GST check pointing restart file'
        ELSE IF (ncid.eq.HSS(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Hessian eigenvectors file'
        ELSE IF (ncid.eq.HIS(ng)%ncid) THEN
          type='ROMS/TOMS history file'
        ELSE IF (ncid.eq.ITL(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear model initial file'
        ELSE IF (ncid.eq.LCZ(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Lanczos vectors file'
        ELSE IF (ncid.eq.NRM(1,ng)%ncid) THEN
          type='ROMS/TOMS initial conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(2,ng)%ncid) THEN
          type='ROMS/TOMS model error covariance norm file'
        ELSE IF (ncid.eq.NRM(3,ng)%ncid) THEN
         type='ROMS/TOMS boundary conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(4,ng)%ncid) THEN
          type='ROMS/TOMS surface forcing error covariance norm file'
        ELSE IF (ncid.eq.RST(ng)%ncid) THEN
          type='ROMS/TOMS restart file'
        ELSE IF (ncid.eq.STA(ng)%ncid) THEN
          type='ROMS/TOMS station file'
        ELSE IF (ncid.eq.TLF(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear impulse forcing file'
        ELSE IF (ncid.eq.TLM(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear history file'
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'type',                &
     &                        TRIM(type))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define other global attributes to NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'title',               &
     &                        TRIM(title))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'title', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'rst_file',            &
     &                        TRIM(RST(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'rst_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).gt.0) THEN
              status=nf90_put_att(ncid, nf90_global, 'his_base',        &
     &                            TRIM(HIS(ng)%base))
            ELSE
              status=nf90_put_att(ncid, nf90_global, 'his_file',        &
     &                            TRIM(HIS(ng)%name))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) 'his_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (ObcData(ng)) THEN
          IF (exit_flag.eq.NoError) THEN
            CALL join_string (BRY(ng)%files, BRY(ng)%Nfiles,            &
     &                        string, lstr)
            status=nf90_put_att(ncid, nf90_global, 'bry_file',          &
     &                          string(1:lstr))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) 'bry_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  NLM Lateral boundary conditions.
!
        IF (exit_flag.eq.NoError) THEN
          CALL lbc_putatt (ng, ncid, ncname, 'NLM_LBC', LBC, status)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'NLM_LBC', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  SVN repository information.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_url',             &
     &                        TRIM(svn_url))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'svn_url', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_rev',             &
     &                        TRIM(svn_rev))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'svn_rev', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'code_dir',            &
     &                        TRIM(Rdir))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'code_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_dir',          &
     &                        TRIM(Hdir))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'header_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_file',         &
     &                        TRIM(Hfile))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'header_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Attributes describing platform and compiler
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'os',                  &
     &                        TRIM(my_os))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'os', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'cpu',                 &
     &                        TRIM(my_cpu))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'cpu', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_system',     &
     &                        TRIM(my_fort))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'compiler_system',            &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_command',    &
     &                        TRIM(my_fc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'compiler_command',           &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_flags',      &
     &                        TRIM(my_fflags))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'compiler_flags', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Tiling and history attributes.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'tiling',              &
     &                        TRIM(tiling))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'tiling', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'history',             &
     &                        TRIM(history))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'history', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Analytical header files used.
!
        IF (exit_flag.eq.NoError) THEN
          CALL join_string (ANANAME, 37, string, lstr)
          IF (lstr.gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'ana_file',          &
     &                          string(1:lstr))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) 'ana_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Activated CPP options.
!
        IF (exit_flag.eq.NoError) THEN
          lstr=LEN_TRIM(Coptions)-1
          status=nf90_put_att(ncid, nf90_global, 'CPP_options',         &
     &                        TRIM(Coptions(1:lstr)))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'CPP_options', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      IF (INT(time_ref).eq.-2) THEN
        Vinfo( 3)='days since 1968-05-23 00:00:00 GMT'
        Vinfo( 4)='gregorian'
      ELSE IF (INT(time_ref).eq.-1) THEN
        Vinfo( 3)='days since 0001-01-01 00:00:00'
        Vinfo( 4)='360_day'
      ELSE IF (INT(time_ref).eq.0) THEN
        Vinfo( 3)='days since 0001-01-01 00:00:00'
        Vinfo( 4)='julian'
      ELSE IF (time_ref.gt.0.0_r8) THEN
        WRITE (Vinfo( 3),'(a,1x,a)') 'days since', TRIM(r_text)
      END IF
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Horizontal mixing coefficients.
!
      Vinfo( 1)='nl_tnu2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='nl_visc2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo ,ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        Vinfo( 1)='FSobc_in'
        Vinfo( 2)='free-surface inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='FSobc_out'
        Vinfo( 2)='free-surface outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='M2obc_in'
        Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='M2obc_out'
        Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='Tobc_in'
        Vinfo( 2)='tracers inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='Tobc_out'
        Vinfo( 2)='tracers outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='M3obc_in'
        Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='M3obc_out'
        Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Various parameters.
!
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian. Writing characters in
!  parallel I/O is extremely inefficient.  It is better to write
!  this as an integer switch: 0=Cartesian, 1=spherical.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 9)='Cartesian'
      Vinfo(10)='spherical'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate parameters.
!
      Vinfo( 1)='Vtransform'
      Vinfo( 2)='vertical terrain-following transformation equation'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Vstretching'
      Vinfo( 2)='vertical terrain-following stretching function'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_rho, scalar'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_w, scalar'
      Vinfo(21)='ocean_s_coordinate'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_r, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_w, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(14)='user, scalar'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/usrdim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        Vinfo( 1)='h'
        Vinfo( 2)='bathymetry at RHO-points'
        Vinfo( 3)='meter'
        Vinfo(14)='bath, scalar'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Coriolis Parameter.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='f'
          Vinfo( 2)='Coriolis parameter at RHO-points'
          Vinfo( 3)='second-1'
          Vinfo(14)='coriolis, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='pm'
          Vinfo( 2)='curvilinear coordinate metric in XI'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pm, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo( 1)='pn'
          Vinfo( 2)='curvilinear coordinate metric in ETA'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pn, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_rho'
          Vinfo( 2)='longitude of RHO-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_rho, scalar'
          Vinfo(21)='longitude'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='lat_rho'
          Vinfo( 2)='latitude of RHO-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_rho, scalar'
          Vinfo(21)='latitude'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo,  ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_rho'
          Vinfo( 2)='x-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_rho, scalar'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='y_rho'
          Vinfo( 2)='y-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_rho, scalar'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_u'
          Vinfo( 2)='longitude of U-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_u, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='lat_u'
          Vinfo( 2)='latitude of U-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_u, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_u'
          Vinfo( 2)='x-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_u, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='y_u'
          Vinfo( 2)='y-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_u, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_v'
          Vinfo( 2)='longitude of V-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_v, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='lat_v'
          Vinfo( 2)='latitude of V-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_v, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_v'
          Vinfo( 2)='x-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_v, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='y_v'
          Vinfo( 2)='y-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_v, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_psi'
          Vinfo( 2)='longitude of PSI-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_psi, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='lat_psi'
          Vinfo( 2)='latitude of PSI-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_psi, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_psi'
          Vinfo( 2)='x-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_psi, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Vinfo( 1)='y_psi'
          Vinfo( 2)='y-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_psi, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Angle between XI-axis and EAST at RHO-points.
!
        Vinfo( 1)='angle'
        Vinfo( 2)='angle between XI-axis and EAST'
        Vinfo( 3)='radians'
        Vinfo(14)='angle, scalar'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo,  ncname)
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
  10  FORMAT (i3.3,'x',i3.3)
  20  FORMAT (/,' DEF_INFO - error while creating global attribute: ',  &
     &        a,/,12x,a)
  30  FORMAT ('frc_file_',i2.2)
      RETURN
      END SUBROUTINE def_info
