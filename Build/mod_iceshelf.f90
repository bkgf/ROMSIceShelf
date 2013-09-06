      MODULE mod_iceshelf
!
!svn $Id$
!================================================== Hernan G. Arango ===
!    Copyright (c) 2002-2012 The ROMS/TOMS Group                       !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!============================================ Benjmain K. Galton-Fenzi==
!                                                                      !
!  Parameters for ice shelf model:                                     !
!  =============================                                       !
! COMMENT: List here		                                       !
!=======================================================================
!
      USE mod_param
      USE mod_kinds
!
      implicit none
      real(r8), parameter :: a = -0.057_r8
      real(r8), parameter :: b = 0.0939_r8
      real(r8), parameter :: c = 7.61e-4
      real(r8), parameter :: gamma = 0.0001_r8
      real(r8), parameter :: refSalt = 34.4_r8
      real(r8), parameter :: L=3.34e5_r8
      real(r8) :: temp_f
      real(r8), parameter :: trelax = 7.0_r8 * 86400.0_r8 ! 7 days
      real(r8), parameter :: saltMax = 35.0_r8
      real(r8), parameter :: saltMin = 34.0_r8
      real(r8), parameter :: sRateInc = 0.0085_r8
      real(r8), parameter :: sRateDec = 0.0283333_r8
      real(r8), parameter :: Tmax = 0.0_r8
      real(r8), parameter :: Tmin = -1.9_r8
      real(r8) :: tyear, sfcTemp, sfcSalt
      real(r8), parameter :: eps = 1.0E-14_r8
      END MODULE mod_iceshelf
