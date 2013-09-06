      MODULE nf_fwrite2d_bry_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 2D boundary array !
!  into an output NetCDF file.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF output file name (string)                    !
!     ncid         NetCDF file ID (integer)                            !
!     ncvname      NetCDF variable name (string)                       !
!     ncvarid      NetCDF variable ID (integer)                        !
!     tindex       NetCDF time record index to write (integer)         !
!     gtype        Grid type (integer)                                 !
!     LBij         IJ-dimension Lower bound (integer)                  !
!     UBij         IJ-dimension Upper bound (integer)                  !
!     Nrec         Number of boundary records (integer)                !
!     Ascl         Factor to scale field before writing (real)         !
!     A            Boundary field to write out (real)                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status       Error flag (integer).                               !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      FUNCTION nf_fwrite2d_bry (ng, model, ncname, ncid,                &
     &                          ncvname, ncvarid,                       &
     &                          tindex, gtype,                          &
     &                          LBij, UBij, Nrec,                       &
     &                          Ascl, A)  RESULT(status)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: LBij, UBij, Nrec
      real(r8), intent(in) :: Ascl
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname
      real(r8), intent(in) :: A(LBij:,:,:)
!
!  Local variable declarations.
!
      logical, dimension(4) :: bounded
      integer :: bc, i, ib, ic, ir, j, rc, tile
      integer :: IorJ, Imin, Imax, Jmin, Jmax, Npts
      integer :: Istr, Iend, Jstr, Jend
      integer, dimension(4) :: start, total
      integer :: status
      real(r8), parameter :: Aspv = 0.0_r8
      real(r8), dimension((UBij-LBij+1)*4*Nrec) :: Awrk
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
      tile=-1
      Imin=LBij
      Imax=UBij
      Jmin=LBij
      Jmax=UBij
      IorJ=IOBOUNDS(ng)%IorJ
      Npts=IorJ*4*Nrec
!
!  Get tile bounds.
!
      Istr=BOUNDS(ng)%Istr (tile)
      Iend=BOUNDS(ng)%Iend (tile)
      Jstr=BOUNDS(ng)%Jstr (tile)
      Jend=BOUNDS(ng)%Jend (tile)
!
!  Set switch to process boundary data by their associated tiles.
!
      bounded(iwest )=DOMAIN(ng)%Western_Edge(tile)
      bounded(ieast )=DOMAIN(ng)%Eastern_Edge(tile)
      bounded(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      bounded(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set NetCDF dimension counters for processing requested field.
!
      start(1)=1
      total(1)=IorJ
      start(2)=1
      total(2)=4
      start(3)=1
      total(3)=Nrec
      start(4)=tindex
      total(4)=1
!
!-----------------------------------------------------------------------
!  Pack and scale output data.
!-----------------------------------------------------------------------
!
      Awrk=Aspv
      DO ir=1,Nrec
        rc=(ir-1)*IorJ*4
        DO ib=1,4
          IF (bounded(ib)) THEN
            bc=(ib-1)*IorJ+rc
            IF ((ib.eq.iwest).or.(ib.eq.ieast)) THEN
              DO j=Jmin,Jmax
                ic=1+(j-LBij)+bc
                Awrk(ic)=A(j,ib,ir)*Ascl
              END DO
            ELSE IF ((ib.eq.isouth).or.(ib.eq.inorth)) THEN
              DO i=Imin,Imax
                ic=1+(i-LBij)+bc
                Awrk(ic)=A(i,ib,ir)*Ascl
              END DO
            END IF
          END IF
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
      IF (OutThread) THEN
        status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
      END IF
      RETURN
      END FUNCTION nf_fwrite2d_bry
      END MODULE nf_fwrite2d_bry_mod
