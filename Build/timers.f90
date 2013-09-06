      SUBROUTINE wclock_on (ng, model, region)
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine turns on wall clock to meassure the elapsed time in    !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: iregion, MyModel
      integer :: my_getpid
      integer :: my_threadnum
      real(r8), dimension(2) :: wtime
      real(r8) :: my_wtime
!
!-----------------------------------------------------------------------
! Initialize timing for all threads.
!-----------------------------------------------------------------------
!
!  Insure that MyModel is not zero.  Notice that zero value is used to
!  indicate restart of the nonlinear model.
!
      MyModel=MAX(1,model)
      Cstr(region,MyModel,ng)=my_wtime(wtime)
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.0)) THEN
        DO iregion=1,Nregion
          Cend(iregion,MyModel,ng)=0.0_r8
          Csum(iregion,MyModel,ng)=0.0_r8
        END DO
        proc(1,MyModel,ng)=1
        proc(0,MyModel,ng)=my_getpid()
        WRITE (stdout,10) ' Thread #', MyThread,                        &
     &                    ' (pid=',proc(0,MyModel,ng),') is active.'
 10     FORMAT (a,i3,a,i8,a)
        thread_count=thread_count+1
        IF (thread_count.eq.numthreads) thread_count=0
      END IF
      RETURN
      END SUBROUTINE wclock_on
      SUBROUTINE wclock_off (ng, model, region)
!
!=======================================================================
!                                                                      !
!  This routine turns off wall clock to meassure the elapsed time in   !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: imodel, iregion, MyModel
      integer :: my_threadnum
      real(r8) :: percent, sumcpu, sumper, total
      real(r8), dimension(2) :: wtime
      real(r8) :: my_wtime
      character (len=14), dimension(4) :: label
!
!-----------------------------------------------------------------------
!  Compute elapsed wall time for all threads.
!-----------------------------------------------------------------------
!
!  Insure that MyModel is not zero.  Notice that zero value is used to
!  indicate restart of the nonlinear model.
!
      MyModel=MAX(1,model)
      IF (region.ne.0) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
      END IF
!
!  Report elapsed wall time.
!
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.1)) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
        DO imodel=1,4
          proc(1,imodel,ng)=0
        END DO
!
!  Report total elapsed time (seconds) for each CPU.
!
        WRITE (stdout,10) ' Thread #', MyThread, ' CPU:',               &
     &                    Cend(region,MyModel,ng)
 10     FORMAT (a,i3,a,f12.3)
!
! Report elapsed time profile for each region of the code.
!
        thread_count=thread_count+1
        DO imodel=1,4
          Csum(region,imodel,ng)=Csum(region,imodel,ng)+                &
     &                           Cend(region,imodel,ng)
          DO iregion=1,Nregion
            Csum(iregion,imodel,ng)=Csum(iregion,imodel,ng)+            &
     &                              Cend(iregion,imodel,ng)
          END DO
        END DO
        DO imodel=1,4
          IF (imodel.ne.MyModel) THEN
            DO iregion=1,Nregion
              Csum(region,imodel,ng)=Csum( region,imodel,ng)+           &
     &                               Csum(iregion,imodel,ng)
            END DO
          END IF
        END DO
        IF (thread_count.eq.numthreads) THEN
          thread_count=0
          IF (Master) THEN
            total=0.0_r8
            DO imodel=1,4
              total=total+Csum(region,imodel,ng)
            END DO
            WRITE (stdout,20) ' Total:', total
 20         FORMAT (a,8x,f14.3)
          END IF
!
!  Report profiling times.
!
          label(iNLM)='Nonlinear     '
          label(iTLM)='Tangent linear'
          label(iRPM)='Representer   '
          label(iADM)='Adjoint       '
          DO imodel=1,4
            IF (Master.and.(Csum(region,imodel,ng).gt.0.0_r8)) THEN
              WRITE (stdout,30) TRIM(label(imodel)),                    &
     &                          ' model elapsed time profile:'
 30           FORMAT (/,1x,a,a,/)
            END IF
            sumcpu=0.0_r8
            sumper=0.0_r8
            DO iregion=1,Mregion-1
              IF (Csum(iregion,imodel,ng).gt.0.0_r8) THEN
                percent=100.0_r8*Csum(iregion, imodel,ng)/              &
     &                           Csum( region,MyModel,ng)
                IF (Master) WRITE (stdout,40) Pregion(iregion),         &
     &                                        Csum(iregion,imodel,ng),  &
     &                                        percent
                sumcpu=sumcpu+Csum(iregion,imodel,ng)
                sumper=sumper+percent
              END IF
            END DO
 40         FORMAT (2x,a,t53,f14.3,2x,'(',f7.4,' %)')
            IF (Master.and.(Csum(region,imodel,ng).gt.0.0_r8)) THEN
              WRITE (stdout,50) sumcpu, sumper
 50           FORMAT (t47,'Total:',f14.3,2x,f8.4)
            END IF
          END DO
          IF (Master) WRITE (stdout,60) Csum(region,MyModel,ng)
  60      FORMAT (/,' All percentages are with respect to total time =',&
     &            5x,f12.3)
        END IF
      END IF
      RETURN
      END SUBROUTINE wclock_off
