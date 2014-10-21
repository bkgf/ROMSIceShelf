      SUBROUTINE checkadj
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency with all available algorithms.                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ic = 0
      integer :: ifound
      character (len=40) :: string
!
!-----------------------------------------------------------------------
!  Report issues with various C-preprocessing options.
!-----------------------------------------------------------------------
!
      string=uppercase('ts_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF
      string=uppercase('ts_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF
      string=uppercase('uv_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF
      string=uppercase('uv_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string),                      &
     &                               'stability problems, WARNING'
      END IF
!
!-----------------------------------------------------------------------
!  Set execution error flag to stop execution.
!-----------------------------------------------------------------------
!
      IF (ic.gt.0) THEN
        exit_flag=5
      END IF
!
 10   FORMAT (/,' CHECKADJ - use caution when activating: ', a,/,12x,   &
     &        'REASON: ',a,'.')
      RETURN
      END SUBROUTINE checkadj
