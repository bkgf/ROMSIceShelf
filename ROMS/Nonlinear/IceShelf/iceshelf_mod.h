#include "cppdefs.h"
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

#  if defined ICESHELF_2EQN_VBC
      real(r8), parameter :: gamma = 0.0001_r8
      real(r8), parameter :: refSalt = 34.4_r8
      real(r8), parameter :: L=3.34e5_r8
      real(r8) :: temp_f
#  elif defined ICESHELF_3EQN_VBC
      real(r8), parameter :: Pr = 13.8_r8
      real(r8), parameter :: Sc = 2432.2_r8
      real(r8), parameter :: Cd = 5.0e-3_r8
      real(r8), parameter :: visc = 1.95e-6_r8
      real(r8), parameter :: L = 3.33e5_r8
      real(r8), parameter :: small = 8.0e-6_r8 ! lower limit for diffusion alone
      real(r8), parameter :: dt_i = 1.54e-6
      real(r8), parameter :: cp_w = 3947.0_r8
      real(r8), parameter :: rho_i = 920.0_r8
      real(r8), parameter :: Ti = -20.0_r8
      real(r8), parameter :: Si = 0.0_r8
#  endif
#  ifdef ANA_SEAICE
      real(r8), parameter :: trelax = 3.0_r8 * 86400.0_r8 ! 3 days
      real(r8), parameter :: saltMax = 34.5_r8
      real(r8), parameter :: saltMin = 34.5_r8
      real(r8), parameter :: sRateInc = 0.0085_r8
      real(r8), parameter :: sRateDec = 0.0283333_r8
      real(r8), parameter :: Tmax = -1.85_r8
      real(r8), parameter :: Tmin = -1.85_r8
      real(r8) :: tyear, sfcTemp, sfcSalt
#   else
      real(r8), parameter :: trelax = 10.0_r8 * 86400.0_r8 !1 day
      real(r8), parameter :: sfcTemp = -1.85_r8
      real(r8), parameter :: sfcSalt = 34.5_r8
#   endif
      real(r8), parameter :: eps = 1.0E-14_r8
