      SUBROUTINE regrid (ng, model, ncname, ncid,                       &
     &                   ncvname, ncvarid, gtype, iflag,                &
     &                   Nx, Ny, Finp, Amin, Amax,                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Imin, Imax, Jmin, Jmax,                        &
     &                   Xout, Yout, Fout)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates gridded data, Finp, to model locations    !
!  Xout and Yout.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
!     ncvname    NetCDF variable name (string)                         !
!     ncvarid    NetCDF variable ID (integer)                          !
!     gtype      C-grid type (integer)                                 !
!     iflag      Interpolation flag (integer, 0: linear, 1: cubic)     !
!     Nx         X-dimension size for gridded data, Finp (integer)     !
!     Ny         Y-dimension size for gridded data, Finp (integer)     !
!     Finp       Gridded data to interpolate from (real)               !
!     Amin       Gridded data minimum value (integer)                  !
!     Amax       Gridded data maximum value (integer)                  !
!     LBi        Fout I-dimension Lower bound (integer)                !
!     UBi        Fout I-dimension Upper bound (integer)                !
!     LBj        Fout J-dimension Lower bound (integer)                !
!     UBj        Fout J-dimension Upper bound (integer)                !
!     Imin       Fout starting data I-index (integer)                  !
!     Imax       Fout ending   data I-index (integer)                  !
!     Jmin       Fout starting data J-index (integer)                  !
!     Jmax       Fout ending   data J-index (integer)                  !
!     Xout       X-locations (longitude) to interpolate (real)         !
!     Yout       Y-locations (latitude)  to interpolate (real)         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated field (real)                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE interpolate_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, gtype, iflag
      integer, intent(in) :: Nx, Ny
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
!
      real(r8), intent(inout) :: Amin, Amax
!
      real(r8), intent(inout) :: Finp(Nx,Ny)
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname
!
!  Local variable declarations
!
      logical :: rectangular
      integer :: i, j
      integer :: Istr, Iend, Jstr, Jend
      real(r8), parameter :: IJspv = 0.0_r8
      real(r8) :: my_min, my_max, Xmin, Xmax, Ymin, Ymax
      real(r8), dimension(Nx,Ny) :: angle
      real(r8), dimension(Nx,Ny) :: Xinp
      real(r8), dimension(Nx,Ny) :: Yinp
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Iout
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Jout
!
!-----------------------------------------------------------------------
!  Get input variable coordinates.
!-----------------------------------------------------------------------
!
      CALL get_varcoords (ng, model, ncname, ncid,                      &
     &                    ncvname, ncvarid, Nx, Ny,                     &
     &                    Xmin, Xmax, Xinp, Ymin, Ymax, Yinp,           &
     &                    rectangular)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set input gridded data rotation angle.
!
      DO i=1,Nx
        DO j=1,Ny
          angle(i,j)=0.0_r8
        END DO
      END DO
!
!  Initialize local fractional coordinates arrays to avoid
!  deframentation.
!
      Iout=0.0_r8
      Jout=0.0_r8
!
!-----------------------------------------------------------------------
!  Check if gridded data contains model grid.
!-----------------------------------------------------------------------
!
      IF ((LonMin(ng).lt.Xmin).or.                                      &
     &    (LonMax(ng).gt.Xmax).or.                                      &
     &    (LatMin(ng).lt.Ymin).or.                                      &
     &    (LatMax(ng).gt.Ymax)) THEN
        IF (Master) THEN
          WRITE (stdout,10) Xmin, Xmax, Ymin, Ymax,                     &
     &                      LonMin(ng), LonMax(ng),                     &
     &                      LatMin(ng), LatMax(ng)
 10       FORMAT (/, ' REGRID - input gridded data does not contain',   &
     &               ' model grid:', /,                                 &
     &            /,10x,'Gridded:  LonMin = ',f9.4,' LonMax = ',f9.4,   &
     &            /,10x,'          LatMin = ',f9.4,' LatMax = ',f9.4,   &
     &            /,10x,'Model:    LonMin = ',f9.4,' LonMax = ',f9.4,   &
     &            /,10x,'          LatMin = ',f9.4,' LatMax = ',f9.4)
        END IF
        exit_flag=4
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Interpolate (bilinear or bicubic) to requested positions.
!-----------------------------------------------------------------------
!
!  Set tile starting and ending indices.
!
      Istr=Imin
      Iend=Imax
      Jstr=Jmin
      Jend=Jmax
!
!  Find fractional indices (Iout,Jout) of the grid cells in Finp
!  containing positions to intepolate.
!
      CALL hindices (ng, 1, Nx, 1, Ny, 1, Nx, 1, Ny,                    &
     &               angle, Xinp, Yinp,                                 &
     &               LBi, UBi, LBj, UBj,                                &
     &               Istr, Iend, Jstr, Jend,                            &
     &               Xout, Yout,                                        &
     &               Iout, Jout,                                        &
     &               IJspv, rectangular)
      IF (iflag.eq.linear) THEN
        CALL linterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Finp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, Xout, Yout,                         &
     &                  Fout, my_min, my_max)
      ELSE IF (iflag.eq.cubic) THEN
        CALL cinterp2d (ng, 1, Nx, 1, Ny,                               &
     &                  Xinp, Yinp, Finp,                               &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Iout, Jout, Xout, Yout,                         &
     &                  Fout, my_min, my_max)
      END IF
!
!  Compute global interpolated field minimum and maximum values.
!  Notice that gridded data values are overwritten.
!
      Amin=my_min
      Amax=my_max
      RETURN
      END SUBROUTINE regrid
