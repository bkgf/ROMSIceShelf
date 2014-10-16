#include "cppdefs.h"
#if defined NONLINEAR && (defined SSH_TIDES && defined NODAL_TIDES)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:26:07:2007      *
! *                                                                            *
! *   REVISION    :  See revision 1, below                 BGF:16:10:2014      *
! *                                                                            *
! *   SOURCE      :  tpd_roms.f                                                *
! *   ROUTINE NAME:  tpd_roms                                                  *
! *   TYPE        :  MAIN                                                      *
! *                                                                            *
! *   FUNCTION    :  Tide Prediction Driver (TPD) subroutines for ROMS         *
! *                  based on tpd1604.f                                        *
! *   REVISION 1  : Implementation and conversion to be f90 compliant          *    
! *                                                                            *
! ******************************************************************************
! *                                                                            *
! *                                SOFTWARE LICENSING                          *
! *                                                                            *
! *                  Copyright (C) 2002 John Robert Hunter                     *
! *                                                                            *
! *                  This program is free software; you can redistribute       *
! *                  it and/or modify it under the terms of the GNU General    *
! *                  Public License as published by the Free Software          *
! *                  Foundation; either Version 2 of the license, or (at       *
! *                  your option) any later version.                           *
! *                                                                            *
! *                  This program is distributed in the hope that it will      *
! *                  be useful, but without any warranty; without even the     *
! *                  implied warranty of merchantability or fitness for a      *
! *                  particular purpose. See the GNU General Public License    *
! *                  for more details.                                         *
! *                                                                            *
! *                  A copy of the GNU General Public License is available     *
! *                  at http://www.gnu.org/copyleft/gpl.html or by writing     *
! *                  to the Free Software Foundation, Inc., 59 Temple Place    *
! *                  - Suite 330, Boston, MA 02111-1307, USA.                  *
! *                                                                            *
! ******************************************************************************
!   
      real*8 function tide_pred(iye,mon,idy,ihr,min,sec,tz_t,           &
     &                          htid,gtid,                              &
     &                          aconin,                                 &
     &                          conmax)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:22:01:1997      *
! *                                                                            *
! *   REVISION    :  Mod.to allow comments (prefix "#") in                     *
! *                  input file                            JRH:17:04:1998      *
! *                  Addition of time zones                JRH:17:04:1998      *
! *                  Interpolation of j and v based on GMT JRH:17:04:1998      *
! *                  tz_t changed to real*8                JRH:10:09:1999      *
! *                  Mod. to calc. of interp. weights      JRH:10:09:1999      *
! *                  Mod. to ignore "*" in harmonic                            *
! *                  constant file                         JRH:13:09:1999      *
! *                  Addition of many "save"s              JRH:31:07:2000      *
! *                  !heck Z0 occurs once and only once    JRH:31:07:2000      *
! *                  Mod. to allow any year and both "98"                      *
! *                  and "1998" type formats               JRH:06:09:2001      *
! *                  Modified for use with astronomical                        *
! *                  argument routine get_ast_arg          JRH:23:09:2004      *
! *                  Minor change to comment               JRH:09:12:2004      *
! *                  Number of constituents (including Z0)                     *
! *                  increased to 116                      JRH:09:12:2004      *
! *                  Unnecessary setting of ios removed    JRH:09:12:2004      *
! *                                                                            *
! *   SOUR!E      :  tpd1528.f                                                 *
! *   ROUTINE NAME:  tide_pred                                                 *
! *   TYPE        :  real*8 function                                           *
! *                                                                            *
! *   FUN!TION    :  General tide prediction subroutine                        *
! *                                                                            *
! *                  iye ........ year (1970-2030)                             *
! *                  mon ........ month (1-12)                                 *
! *                  idy ........ day of month (1-31)                          *
! *                  ihr ........ hour of day (0-24)                           *
! *                  min ........ minute of hour (0-60)                        *
! *                  sec ........ second of minute (0-60)                      *
! *                  tz_t ....... time zone of time (hours, positive eastwards)*
! *                  htid ....... input array of "h" harmonic constant         *
! *                  gtid ....... input array of "g" harmonic constant         *
! *                  aconin ..... input array of constituent names             *
! *                  conmax ..... max. number of constituent                   *
! *                                                                            *
! *                  NOTE that this subroutine uses input device number 7.     *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      integer conmax
!
      real*8 htid(conmax),gtid(conmax)
!
      real*8 tz_t
      real*4 sec
      integer iye,mon,idy,ihr,min
!
      character*6 aconin(conmax)
!
!     Local variables:
!
      integer pcon,ptry
      parameter(pcon=116,ptry=3)   ! >/= Maximum no. of constituents used (including Z0)
             ! No. of tries at different const. names
!
      real*8 h(pcon),g(pcon)
      real*8 j1(pcon),v1(pcon),vpv(pcon)
      real*8 j2(pcon),v2(pcon),v0(pcon)
      real*8 sigma(pcon)
!
      real*8 j_int,vpv_int
      real*8 j_day,j_day_1,j_day_2,tim
      real*8 tz_g
      real*8 d2r
      real*8 dum1,dum2,year
      real*8 w1,w2
!
      integer iye_ast_arg_old
      integer i,i_z0,n_z0,ncon
      integer ios,itry
      integer len,lenchar
      integer ifail_ast,ifail
!
      logical first,found,match,l_tz_g
!
      character*6 acon(pcon)
      character*6 atry(ptry,2)
!
      character*40 buf
      character*80 formt
!
      save iye_ast_arg_old,first
      save d2r,tz_g,l_tz_g,ncon,acon,h,g,i_z0
      save sigma,j1,v1,vpv,j2,v2,v0
!
      data iye_ast_arg_old/0/
!
      data atry(1,1),atry(1,2)/'sig1  ','sigma1'/
      data atry(2,1),atry(2,2)/'the1  ','theta1'/
      data atry(3,1),atry(3,2)/'lam2  ','lamda2'/
!
      data first/.true./
      data l_tz_g/.false./
!
      if(first) d2r=datan(1.d0)*4.d0/180.d0                 ! Degrees to radians
!
!     !heck year (based on email of David Blackman, 27/8/2004, but
!     NOTE THAT GREAT BRITAIN DID NOT ADOPT GREGORIAN !ALENDAR UNTIL 1752):
!
      if(iye.lt.1600) then
        stop
      endif
!
!     Set harmonic constants:
!
      h(1)=0.d0
      g(1)=0.d0
      acon(1)='z0    '
      i_z0=1
!
      ncon=conmax
!
      do i=2,ncon+1
        h(i)=htid(i-1)
        g(i)=gtid(i-1)
        acon(i)=aconin(i-1)
      end do
!
!     Calculate astronomical arguments, if necessary:
!
      if(iye.ne.iye_ast_arg_old) then
!
        iye_ast_arg_old=iye                                  ! Remember old year
!
        do i=1,ncon
!
          if(acon(i).eq.'z0    ') then                           ! Trap mean, Z0
            i_z0=i
          else
            itry=0
            found=.false.
            do while(.not.found)
              call get_ast_arg(iye,acon(i),sigma(i),j1(i),v1(i),vpv(i), &
     &                         ifail_ast)
              found=(ifail_ast.eq.0)                            ! .true. if O.K.
              match=(ifail_ast.eq.0)              ! If not O.K., try other names
!
!             write(6,5) acon(i)                                    ! Diagnostic
!   5         format('Trying constituent ',a6)                      ! Diagnostic
!
              do while(.not.match.and.itry.lt.ptry)
                itry=itry+1
                if(acon(i).eq.atry(itry,1)) then
                  acon(i)=atry(itry,2)
                  match=.true.
                else if(acon(i).eq.atry(itry,2)) then
                  acon(i)=atry(itry,1)
                  match=.true.
                endif
!               write(8,3) itry,i,acon(i)                           ! Diagnostic
!   3           format('itry,i,acon(i): ',2i5,1x,a6)                ! Diagnostic
              end do
!
              if(.not.found.and..not.match) then
                stop
              endif
            end do
!
            sigma(i)=sigma(i)/3600.d0            ! !onvert to degrees per second
            v0(i)=vpv(i)-v1(i)
!
          endif
!
        end do
!
        do i=1,ncon
!
          if(acon(i).ne.'z0    ') then                           ! Trap mean, Z0
            call get_ast_arg(iye+1,acon(i),dum1,j2(i),v2(i),dum2,       &
     &                       ifail_ast)
            if(ifail_ast.ne.0) then
              stop
            endif
!
          endif
!
        end do
!
      endif
!
!     Assume tz_g is the same as tz_t if not given in input file:
!
      if(.not.l_tz_g) tz_g=tz_t
!
!     Convert to time from start of year:
!
      call julday(iye,mon,idy,ihr,min,sec,j_day  ,ifail)
      if(ifail.ne.0) then
        stop
      endif
!
      call julday(iye,1  ,1  ,0  ,0  ,0. ,j_day_1,ifail)
      if(ifail.ne.0) then
        stop
      endif
!
      tim=(j_day-j_day_1)*86400.d0
!
!     Find length of current year:
!
      call julday(iye+1,1  ,1  ,0  ,0  ,0. ,j_day_2,ifail)
      if(ifail.ne.0) then
        stop
      endif
!
      year=(j_day_2-j_day_1)*86400.d0
!
!     Linearly interpolate nodal corrections through year,
!     and calculate tide:
!
      tide_pred=h(i_z0)                                               ! Mean, Z0
      w2=(tim-tz_t*3600.d0)/year                                   ! Time in GMT
      w1=1.d0-w2
      do i=1,ncon
        if(i.ne.i_z0) then
          j_int=j1(i)*w1+j2(i)*w2
          vpv_int=v0(i)+v1(i)*w1+v2(i)*w2
          tide_pred=tide_pred                                           &
     &      +h(i)*j_int*dcos((sigma(i)*(tim+(tz_g-tz_t)*3600.d0)        &
     &      +vpv_int-g(i))*d2r)
        endif
      end do
!
      return
      end
!
      subroutine changecase(ain,aout,n,itype)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:03:12:1990      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:1990      *
! *                                                                            *
! *   SOURCE      :  tpd1450.f                                                 *
! *   ROUTINE NAME:  changecase                                                *
! *   TYPE        :  subroutine                                                *
! *                                                                            *
! *   FUNCTION    :  Changes case of character variable                        *
! *                                                                            *
! *                  AIN ..... input variable                                  *
! *                  AOUT .... output variable                                 *
! *                  N ....... number of characters in AIN and AOUT            *
! *                  ITYPE ... 0 for change to lower case                      *
! *                            1 for change to upper case                      *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      integer n,itype
      character*(*) ain,aout
!
      integer idum,i
!
      if(itype.eq.0) then                                 ! Change to lower case
        do i=1,n
          idum=ichar(ain(i:i))
          if(idum.ge.65.and.idum.le.90) then
            idum=idum+32
            aout(i:i)=char(idum)
          else
            aout(i:i)=ain(i:i)
          endif
        end do
        return
      else if(itype.eq.1) then                            ! Change to upper case
        do i=1,n
          idum=ichar(ain(i:i))
          if(idum.ge.97.and.idum.le.122) then
            idum=idum-32
            aout(i:i)=char(idum)
          else
            aout(i:i)=ain(i:i)
          endif
        end do
        return
      else                                 ! ITYPE out of range ..... do nothing
        do i=1,n
          aout(i:i)=ain(i:i)
        end do
        return
      endif
      end
!
      subroutine error_handler(ifail,error_point)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  MODELLING                         REF:JRH:06:01:1995      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:1995      *
! *                                                                            *
! *   SOURCE      :  tpd1450.f                                                 *
! *   ROUTINE NAME:  error_handler                                             *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Handles error (ifail.ne.0)                                *
! *                                                                            *
! ******************************************************************************
      integer ifail,error_point
!
      if(ifail.ne.0) then
        write(6,1) ifail,error_point
    1   format(/' ifail returned as ',i5,' at error point ',i5,         &
     &          ' ..... program terminated'/)
        stop
      endif
      return
      end
!
      integer function idint2(x)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:23:04:1991      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:1991      *
! *                                                                            *
! *   SOURCE      :  tpd1450.f                                                 *
! *   ROUTINE NAME:  idint2                                                    *
! *   TYPE        :  integer function                                        *
! *                                                                            *
! *   FUNCTION    :  As INT, but trunctates towards -(infinity)                *
!                    This is D.P. version of INT2                              *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      real*8 x
      idint2=idint(x)
      if(dble(idint2).eq.x) return
      if(x.lt.0.d0) idint2=idint2-1
      return
      end 
!
      subroutine julday(iye,mon,idy,ihr,min,sec,julian,ifail)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:23:05:1991      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:1991      *
! *                                                                            *
! *   SOURCE      :  tpd1450.f                                                 *
! *   ROUTINE NAME:  julday                                                    *
! *   TYPE        :  subroutine                                                *
! *                                                                            *
! *   FUNCTION    :  Converts (year, month, day, hour, minute, second) to      *
! *                  Julian Day (based on Numerical Recipes)                   *
! *                                                                            *
! *                  IYE ........ Year                                         *
! *                  MON ........ Month                                        *
! *                  IDY ........ Day                                          *
! *                  IHR ........ Hour                                         *
! *                  MIN ........ Minute                                       *
! *                  SEC ........ Second                                       *
! *                  JULIAN ..... Julian Day (non-negative, but may be         *
! *                               non-integer)                                 *
! *                  IFAIL ...... 0 for successful execution, otherwise 1      *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      real*8 julian
      real*4 sec
      integer iye,mon,idy,ihr,min,ifail
      integer iyyy,jy,jm,igreg,ja,ijul
      integer idint2
      parameter (igreg=15+31*(10+12*1582))
!                              ..... Gregorian Calendar was adopted 15 Oct. 1582
      if(iye.eq.0.or.iye.lt.-4713) then                  ! There is no year zero
                                      ! Julian Day must be non-neagtive
        ifail=1
        return
      endif
      if(iye.lt.0) then
        iyyy=iye+1
      else
        iyyy=iye
      endif
      if(mon.gt.2) then
        jy=iyyy
        jm=mon+1
      else
        jy=iyyy-1
        jm=mon+13
      endif
      ijul=idint2(365.25d0*dble(jy))+idint2(30.6001d0*dble(jm))         &
     &       +idy+1720995
      if(idy+31*(mon+12*iyyy).ge.igreg) then
!                                    ..... Test for change to Gregorian !alendar
        ja=idint(0.01d0*dble(jy))
        ijul=ijul+2-ja+idint(0.25d0*dble(ja))
      endif
      julian=dble(ijul)                                                 &
     &         +dble(ihr)/24.d0+dble(min)/1440.d0+dble(sec)/86400.d0    &
     &         -0.5d0                         ! !orrection from midnight to noon
      if(julian.lt.0.d0) then                  ! Julian Day must be non-negative
        ifail=1
        return
      endif
      ifail=0
      return
      end
!
      integer function lenchar(c)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:09:01:1986      *
! *                                                                            *
! *   REVISION    :  Variable declarations                 JRH:30:09:1996      *
! *                                                                            *
! *   SOURCE      :  tpd1450.f                                                 *
! *   ROUTINE NAME:  lenchar                                                   *
! *   TYPE        :  integer function                                        *
! *                                                                            *
! *   FUNCTION    :  Returns length of CHARACTER  variable (defined by         *
! *                  removing blank characters from right-hand side).          *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      character*(*) c
!
      integer itot,i
!
      itot=len(c)
      do i=itot,1,-1
        if(c(i:i).ne.' ') go to 1
      end do
      lenchar=0                                           ! String is all blanks
      return
    1 lenchar=i
      return
      end
!
      subroutine caldat(iye,mon,idy,ihr,min,sec,julian)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:23:05:1991      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:1991      *
! *                                                                            *
! *   SOURCE      :  FORLIB.FVS                                                *
! *   ROUTINE NAME:  CALDAT                                                    *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Converts Julian Day to (year, month, day, hour, minute,   *
! *                  second) (based on Numerical Recipes)                      *
! *                                                                            *
! *                  IYE ........ Year                                         *
! *                  MON ........ Month                                        *
! *                  IDY ........ Day                                          *
! *                  IHR ........ Hour                                         *
! *                  MIN ........ Minute                                       *
! *                  SEC ........ Second                                       *
! *                  JULIAN ..... Julian Day (non-negative, but may be         *
! *                               non-integer)                                 *
! *                  IFAIL ...... 0 for successful execution, otherwise 1      *
! *                                                                            *
! ******************************************************************************
      real*8 julian
      real*4 sec
      integer iye,mon,idy,ihr,min
      integer igreg,jalpha,ja,jb,jc,jd,je,ijul
      parameter (igreg=2299161)               ! !ross-over to Gregorian !alendar
      if(julian.lt.0.d0) then                  ! Negative Julian Day not allowed
        stop
      endif
      ijul=idint(julian)                                   ! Integral Julian Day
      sec=sngl((julian-dble(ijul))*86400.d0)! Seconds from beginning of Jul. Day
      if(sec.ge.43200.) then                              ! In next calendar day
        ijul=ijul+1
        sec=sec-43200.                            ! Adjust from noon to midnight
      else                                                ! In same calendar day
        sec=sec+43200.                            ! Adjust from noon to midnight
      endif
      if(sec.ge.86400.) then              ! Final check to prevent time 24:00:00
        ijul=ijul+1
        sec=sec-86400.
      endif
      min=int(sec/60.)                  ! Integral minutes from beginning of day
      sec=sec-float(min*60)                   ! Seconds from beginning of minute
      ihr=min/60                          ! Integral hours from beginning of day
      min=min-ihr*60                   ! Integral minutes from beginning of hour
      if(ijul.ge.igreg)then                  ! !orrection for Gregorian !alendar
        jalpha=idint((dble(ijul-1867216)-0.25d0)/36524.25d0)
        ja=ijul+1+jalpha-idint(0.25d0*dble(jalpha))
      else                                                       ! No correction
        ja=ijul
      endif
      jb=ja+1524
      jc=idint(6680.d0+(dble(jb-2439870)-122.1d0)/365.25d0)
      jd=365*jc+idint(0.25d0*dble(jc))
      je=idint(dble(jb-jd)/30.6001d0)
      idy=jb-jd-idint(30.6001d0*dble(je))
      mon=je-1
      if(mon.gt.12)mon=mon-12
      iye=jc-4715
      if(mon.gt.2)iye=iye-1
      if(iye.le.0)iye=iye-1
      return
      end
!
!*******************************************************************************
!
!     The following is modified software provided by POL on 27 Aug 2004, with
!     the following email:
!
! From dlb@pol.ac.uk Mon Aug 30 15:25:04 2004
! Date: Fri, 27 Aug 2004 14:28:24 +0100
! From: David Blackman <dlb@pol.ac.uk>
! To: john.hunter@utas.edu.au
! Subject: Re: A favour
! 
! John,
! 
! Below is our code used in POL predictions; it follows Dr Catrwrights IHR
! LXII paper "Tidal Prediction and Modern Time Scales". The parameter
! delta T corresponds to the values in figure 1. Although you can predict
! back to 1600 the correct delta T is only included for year 1900 onwards.
! I do also have a version with slighly different values of delta T that
! go back to 1620 from Ian Vassie. These give differences in the 4th
! decimal place of V so are not really significant.
! 
! I hope you can understand the code!!
! 
! Regards,
! 
! David
!
!     For tests of these routines against previously provided astronomical
!     arguments, see: /home/johunter/newtides/make_pol_astr_arg
!
!*******************************************************************************
!
      subroutine get_ast_arg(iye,acon,sigma,j,va,vpv,ifail)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:21:09:2004      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:2004      *
! *                                                                            *
! *   SOURCE      :  tpd1604.f                                                 *
! *   ROUTINE NAME:  get_ast_arg                                               *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Returns angular frequency and astronomical arguments for  *
! *                  year iye                                                  *
! *                                                                            *
! *                  NOTE:                                                     *
! *                                                                            *
! *                  1. Number of astronomical arguments is 115, one more than *
! *                     in my previous software, with the addition of 2MN2.    *
! *                     All arrays are therefore dimensioned to 115 in this    *
! *                     set.                                                   *
! *                                                                            *
! *                  2. "v" has been renamed "va" as the argument of this      *
! *                     subroutine, so as not to conflict with POL's use of    *
! *                     "v".                                                   *
! *                                                                            *
! *                  3. All constituent names are lower case, and in my        *
! *                     standard form (i.e. alternative names have to be dealt *
! *                     with elsewhere).                                       *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      real*8 sigma,j,va,vpv
!
      integer iye,ifail
!
      character*6 acon
!     
      real*8 si(115),f(115),u(115),v(115)
!
      real*8 s,p,h,en,p1
!
      integer i,icon,iye_old
!
      character*6 cnam(115)
!
      save iye_old,cnam,si,f,u,v
!
      data iye_old/0/
!
      if(iye.ne.iye_old) then
!
        call inicm1(si)
        call inicm2(cnam)
!
        iye_old=iye
!
!     Calculate s, p, h, en and p1:
!
        call sphen(iye,1,s,p,h,en,p1)
!
!     Calculate u and f:
!
        call ufset(p,en,u,f)
!
!     Calculate v:
!
        call vset(s,p,p1,h,v)
!
      endif
!
!     Find constituent:
!
      i=1
      icon=0
!
      do while(i.le.115.and.icon.eq.0)
!
        if(acon.eq.cnam(i)) icon=i
        i=i+1
      end do
!
      if(icon.eq.0) then
        ifail=1
        return
      endif          
!
      sigma=si(icon)
      j=f(icon)
      va=u(icon)
      if(va.gt.180.d0) va=va-360.d0
      vpv=u(icon)+v(icon)
      if(vpv.lt.0.d0) vpv=vpv+360.d0
      if(vpv.gt.360.d0) vpv=vpv-360.d0
!
      ifail=0
      return
!
      end
!
      subroutine inicm1(si)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:21:09:2004      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:2004      *
! *                                                                            *
! *   SOURCE      :  tpd1604.f                                                 *
! *   ROUTINE NAME:  inicm1                                                    *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Returns angular frequencies, si (degrees/hour).           *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      real*8 si(115)
!
      si(1)   =  0.410686d-01
      si(2)   =  0.821373d-01
      si(3)   =  0.5443747d+00
      si(4)   =  0.10158958d+01
      si(5)   =  0.10980331d+01
      si(6)   =  0.128542862d+02
      si(7)   =  0.129271398d+02
      si(8)   =  0.133986609d+02
      si(9)   =  0.134715145d+02
      si(10)  =  0.139430356d+02
      si(11)  =  0.140251729d+02
      si(12)  =  0.144920521d+02
      si(13)  =  0.145695476d+02
      si(14)  =  0.149178647d+02
      si(15)  =  0.149589314d+02
      si(16)  =  0.150000000d+02
      si(17)  =  0.150410686d+02
      si(18)  =  0.150821353d+02
      si(19)  =  0.151232059d+02
      si(20)  =  0.155125897d+02
      si(21)  =  0.155854433d+02
      si(22)  =  0.160569644d+02
      si(23)  =  0.161391017d+02
      si(24)  =  0.273416965d+02
      si(25)  =  0.274238337d+02
      si(26)  =  0.278953548d+02
      si(27)  =  0.279682084d+02
      si(28)  =  0.284397295d+02
      si(29)  =  0.285125831d+02
      si(30)  =  0.289019669d+02
      si(31)  =  0.289841042d+02
      si(32)  =  0.290662415d+02
      si(33)  =  0.294556253d+02
      si(34)  =  0.295284789d+02
      si(35)  =  0.299589333d+02
      si(36)  =  0.300000000d+02
      si(37)  =  0.300410667d+02
      si(38)  =  0.300821373d+02
      si(39)  =  0.305443747d+02
      si(40)  =  0.306265120d+02
      si(41)  =  0.310158958d+02
      si(42)  =  0.429271398d+02
      si(43)  =  0.434761563d+02
      si(44)  =  0.439430356d+02
      si(45)  =  0.440251729d+02
      si(46)  =  0.450410686d+02
      si(47)  =  0.574238337d+02
      si(48)  =  0.579682084d+02
      si(49)  =  0.584397295d+02
      si(50)  =  0.589841042d+02
      si(51)  =  0.590662415d+02
      si(52)  =  0.600000000d+02
      si(53)  =  0.600821373d+02
      si(54)  =  0.864079380d+02
      si(55)  =  0.869523127d+02
      si(56)  =  0.874238337d+02
      si(57)  =  0.879682084d+02
      si(58)  =  0.880503457d+02
      si(59)  =  0.889841042d+02
      si(60)  =  0.890662415d+02
      si(61)  =  0.264079379d+02
      si(62)  =  0.268701754d+02
      si(63)  =  0.269523127d+02
      si(64)  =  0.275059710d+02
      si(65)  =  0.283575922d+02
      si(66)  =  0.299178627d+02
      si(67)  =  0.310887494d+02
      si(68)  =  0.423827651d+02
      si(69)  =  0.430092770d+02
      si(70)  =  0.445695475d+02
      si(71)  =  0.568701754d+02
      si(72)  =  0.569523127d+02
      si(73)  =  0.578860711d+02
      si(74)  =  0.719112441d+02
      si(75)  =  0.724602605d+02
      si(76)  =  0.730092771d+02
      si(77)  =  0.848476674d+02
      si(78)  =  0.853920422d+02
      si(79)  =  0.858542795d+02
      si(80)  =  0.859364168d+02
      si(81)  =  0.863258006d+02
      si(82)  =  0.864807915d+02
      si(83)  =  0.868701754d+02
      si(84)  =  0.874966873d+02
      si(85)  =  0.885125832d+02
      si(86)  =  0.885947204d+02
      si(87)  =  0.1148476674d+03
      si(88)  =  0.1153920422d+03
      si(89)  =  0.1159364168d+03
      si(90)  =  0.1164079379d+03
      si(91)  =  0.1169523127d+03
      si(92)  =  0.1170344500d+03
      si(93)  =  0.1175059710d+03
      si(94)  =  0.1179682084d+03
      si(95)  =  0.1180503457d+03
      si(96)  =  0.1459364168d+03
      si(97)  =  0.1469523127d+03
      si(98)  =  0.1743761463d+03
      si(99)  =  0.1749205210d+03
      si(100) =  0.1759364168d+03
      si(101) =  0.274966873d+02
      si(102) =  0.278860711d+02
      si(103) =  0.289430356d+02
      si(104) =  0.290251728d+02
      si(105) =  0.304715211d+02
      si(106) =  0.310980331d+02
      si(107) =  0.564079379d+02
      si(108) =  0.574966873d+02
      si(109) =  0.585125830d+02
      si(110) =  0.595284789d+02
      si(111) =  0.283986609d+02
      si(112) =  0.284807981d+02
      si(113) =  0.729271398d+02
      si(114) =  0.740251728d+02
      si(115) =  0.295284789d+02
!
      return
!
      end
!
      subroutine inicm2(cnam)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:21:09:2004      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:2004      *
! *                                                                            *
! *   SOURCE      :  tpd1604.f                                                 *
! *   ROUTINE NAME:  inicm2                                                    *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Returns constituent names, cnam (6 character).            *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
      character*6 cnam(115)
!
      cnam(1) = 'sa    '
      cnam(2) = 'ssa   '
      cnam(3) = 'mm    '
      cnam(4) = 'msf   '
      cnam(5) = 'mf    '
      cnam(6) = '2q1   '
      cnam(7) = 'sigma1'
      cnam(8) = 'q1    '
      cnam(9) = 'rho1  '
      cnam(10) = 'o1    '
      cnam(11) = 'mp1   '
      cnam(12) = 'm1    '
      cnam(13) = 'chi1  '
      cnam(14) = 'pi1   '
      cnam(15) = 'p1    '
      cnam(16) = 's1    '
      cnam(17) = 'k1    '
      cnam(18) = 'psi1  '
      cnam(19) = 'phi1  '
      cnam(20) = 'theta1'
      cnam(21) = 'j1    '
      cnam(22) = 'so1   '
      cnam(23) = 'oo1   '
      cnam(24) = 'oq2   '
      cnam(25) = 'mns2  '
      cnam(26) = '2n2   '
      cnam(27) = 'mu2   '
      cnam(28) = 'n2    '
      cnam(29) = 'nu2   '
      cnam(30) = 'op2   '
      cnam(31) = 'm2    '
      cnam(32) = 'mks2  '
      cnam(33) = 'lamda2'
      cnam(34) = 'l2    '
      cnam(35) = 't2    '
      cnam(36) = 's2    '
      cnam(37) = 'r2    '
      cnam(38) = 'k2    '
      cnam(39) = 'msn2  '
      cnam(40) = 'kj2   '
      cnam(41) = '2sm2  '
      cnam(42) = 'mo3   '
      cnam(43) = 'm3    '
      cnam(44) = 'so3   '
      cnam(45) = 'mk3   '
      cnam(46) = 'sk3   '
      cnam(47) = 'mn4   '
      cnam(48) = 'm4    '
      cnam(49) = 'sn4   '
      cnam(50) = 'ms4   '
      cnam(51) = 'mk4   '
      cnam(52) = 's4    '
      cnam(53) = 'sk4   '
      cnam(54) = '2mn6  '
      cnam(55) = 'm6    '
      cnam(56) = 'msn6  '
      cnam(57) = '2ms6  '
      cnam(58) = '2mk6  '
      cnam(59) = '2sm6  '
      cnam(60) = 'msk6  '
      cnam(61) = '2mn2s2'
      cnam(62) = '3m(sk2'
      cnam(63) = '3m2s2 '
      cnam(64) = 'mnk2s2'
      cnam(65) = 'snk2  '
      cnam(66) = '2sk2  '
      cnam(67) = '2ms2n2'
      cnam(68) = 'mq3   '
      cnam(69) = '2mp3  '
      cnam(70) = '2mq3  '
      cnam(71) = '3mk4  '
      cnam(72) = '3ms4  '
      cnam(73) = '2msk4 '
      cnam(74) = '3mk5  '
      cnam(75) = 'm5    '
      cnam(76) = '3mo5  '
      cnam(77) = '2mn)s6'
      cnam(78) = '3mns6 '
      cnam(79) = '4mk6  '
      cnam(80) = '4ms6  '
      cnam(81) = '2msnk6'
      cnam(82) = '2mv6  '
      cnam(83) = '3msk6 '
      cnam(84) = '4mn6  '
      cnam(85) = '3msn6 '
      cnam(86) = 'mkl6  '
      cnam(87) = '2(mn)8'
      cnam(88) = '3mn8  '
      cnam(89) = 'm8    '
      cnam(90) = '2msn8 '
      cnam(91) = '3ms8  '
      cnam(92) = '3mk8  '
      cnam(93) = 'msnk8 '
      cnam(94) = '2(ms)8'
      cnam(95) = '2msk8 '
      cnam(96) = '4ms10 '
      cnam(97) = '3m2s10'
      cnam(98) = '4msn12'
      cnam(99) = '5ms12 '
      cnam(100) = '4m2s12'
      cnam(101) = 'mvs2  '
      cnam(102) = '2mk2  '
      cnam(103) = 'ma2   '
      cnam(104) = 'mb2   '
      cnam(105) = 'msv2  '
      cnam(106) = 'skm2  '
      cnam(107) = '2mns4 '
      cnam(108) = 'mv4   '
      cnam(109) = '3mn4  '
      cnam(110) = '2msn4 '
      cnam(111) = 'na2   '
      cnam(112) = 'nb2   '
      cnam(113) = 'mso5  '
      cnam(114) = 'msk5  '
      cnam(115) = '2mn2  '
!
      return
!
      end
!
      subroutine sphen (year,vd,s,p,h,en,p1)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:21:09:2004      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:2004      *
! *                                                                            *
! *   SOURCE      :  tpd1604.f                                                 *
! *   ROUTINE NAME:  sphen                                                     *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Returns s, p, h, en and p1.                               *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
!     Checked for 2000 compliance by dlb on 23.09.1997
!
      real*8 s,p,h,en,p1
!
      integer year,vd
!
      real*8 cycle,t,td,delt(84),delta,deltat
!
      integer yr,ilc,icent,it,iday,iyd,ild,ipos
!
      data delt /-5.04d0,-3.90d0,-2.87d0,-0.58d0,0.71d0,1.80d0,         &
     &  3.08d0, 4.63d0, 5.86d0, 7.21d0, 8.58d0,10.50d0,12.10d0,         &
     & 12.49d0,14.41d0,15.59d0,15.81d0,17.52d0,19.01d0,18.39d0,         &
     & 19.55d0,20.36d0,21.01d0,21.81d0,21.76d0,22.35d0,22.68d0,         &
     & 22.94d0,22.93d0,22.69d0,22.94d0,23.20d0,23.31d0,23.63d0,         &
     & 23.47d0,23.68d0,23.62d0,23.53d0,23.59d0,23.99d0,23.80d0,         &
     & 24.20d0,24.99d0,24.97d0,25.72d0,26.21d0,26.37d0,26.89d0,         &
     & 27.68d0,28.13d0,28.94d0,29.42d0,29.66d0,30.29d0,30.96d0,         &
     & 31.09d0,31.59d0,31.52d0,31.92d0,32.45d0,32.91d0,33.39d0,         &
     & 33.80d0,34.23d0,34.73d0,35.40d0,36.14d0,36.99d0,37.87d0,         &
     & 38.75d0,39.70d0,40.70d0,41.68d0,42.82d0,43.96d0,45.00d0,         &
     & 45.98d0,47.00d0,48.03d0,49.10d0,50.10d0,50.97d0,51.81d0,         &
     & 52.57d0/
!
      cycle = 360.0d0
      ilc = 0
      icent = year/100
      yr = year - icent*100
      t = icent - 20
!
!     For the following equations
!     time origin is fixed at 00 hr of jan 1st,2000.
!     See notes by cartwright
!
      it = icent - 20
!
      if (it .lt. 0) then
        iday = it/4 -it
      else
        iday = (it+3)/4 - it
      endif
!
!     t is in julian century
!     Correction in gegorian calendar where only century year divisible
!     by 4 is leap year.
!
      td = 0.0d0
!
      if (yr .ne. 0) then
        iyd = 365*yr
        ild = (yr-1)/4
        if((icent - (icent/4)*4) .eq. 0) ilc = 1
        td = iyd + ild + ilc
      endif
!
      td = td + iday + vd - 1.5d0
      t = t + (td/36525.0d0)
!
      deltat=0.0d0
      ipos=year-1899
!
      if (ipos .gt. 0) then
!
        if (ipos .lt. 83) then
          delta = (delt(ipos+1)+delt(ipos))/2.0d0
        else
          delta= (65.0d0-50.5d0)/20.0d0*(year-1980)+50.5d0
        endif
!
        deltat = delta * 1.0d-6
!
      endif
!
      s  = 218.3165d0 + 481267.8813d0*t - 0.0016d0*t*t + 152.0d0*deltat
      h  = 280.4661d0 + 36000.7698d0*t  + 0.0003d0*t*t + 11.0d0*deltat
      p  =  83.3535d0 + 4069.0139d0*t   - 0.0103d0*t*t + deltat
      en = 125.0445d0 - 1934.1363d0*t   + 0.0021d0*t*t - deltat
      p1 = 282.9384d0 + 1.7195d0*t      + 0.0005d0*t*t
!
      s = dmod (s,cycle)
      if ( s .lt. 0.0d0) s = s+cycle
!
      h = dmod (h,cycle)
      if (h .lt. 0.0d0) h = h+cycle
!
      p = dmod (p,cycle)
      if (p .lt. 0.0d0) p = p+cycle
!
      en = dmod (en,cycle)
      if (en .lt. 0.0d0) en = en + cycle
!
      p1 = dmod (p1,cycle)
      if (p1 .lt. 0.0d0) p1 = p1 + cycle
!
      return
!
      end
!
      subroutine ufset (p,cn,u,f)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:21:09:2004      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:2004      *
! *                                                                            *
! *   SOURCE      :  tpd1604.f                                                 *
! *   ROUTINE NAME:  ufset                                                     *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Returns u and f.                                          *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
!     Version 92.01  24th january
!
!     Checked for 2000 compliance by dlb on 23.09.1997
!
      real*8 u(115),f(115)
!
      real*8 p,cn
!
      real*8 rad,deg,pi,nw,pw,w1,w2,w3,w4,w5,w6,a1,a2,a3,a4,x,y
!
      integer k
!
!     t is zero as compared to tifa.
!
      rad = 6.28318530717959d0/360.0d0
      deg = 360.0d0/6.28318530717959d0
      pi  = 6.28318530717959d0/2.0d0
!
      pw = p*rad
      nw = cn*rad
!
      w1 = dcos(nw)
      w2 = dcos(2.0d0*nw)
      w3 = dcos(3.0d0*nw)
      w4 = dsin(nw)
      w5 = dsin(2.0d0*nw)
      w6 = dsin(3.0d0*nw)
!
      a1 = pw-nw
      a2 = 2.0d0*pw
      a3 = a2-nw
      a4 = a2-2.0d0*nw
!
!     u's are computed in radians
!
      u(3)   =  0.0d0
      f(3)   =  1.0d0    -0.1300d0*w1 +0.0013d0*w2
!
      u(5)   =           -0.4143d0*w4 +0.0468d0*w5 -0.0066d0*w6
      f(5)   =  1.0429d0 +0.4135d0*w1 -0.004d0*w2
!
      u(10)  =            0.1885d0*w4 -0.0234d0*w5 +0.0033d0*w6
      f(10)  =  1.0089d0 +0.1871d0*w1 -0.0147d0*w2 +0.0014d0*w3
!
      x = 2.0d0*dcos(pw)+0.4d0*dcos(a1)
      y = dsin(pw)+0.2d0*dsin(a1)
!
      u(12)  = datan (y/x)
      if (x .lt. 0.0d0) u(12) = u(12)+pi
      f(12)  = dsqrt(x*x+y*y)
!
      u(17)  =           -0.1546d0*w4 +0.0119d0*w5 -0.0012d0*w6
      f(17)  =  1.0060d0 +0.1150d0*w1 -0.0088d0*w2 +0.0006d0*w3
!
      u(21)  =           -0.2258d0*w4 +0.0234d0*w5 -0.0033d0*w6
      f(21)  =  1.0129d0 +0.1676d0*w1 -0.0170d0*w2 +0.0016d0*w3
!
      f(23)  =  1.1027d0 +0.6504d0*w1 +0.0317d0*w2 -0.0014d0*w3
      u(23)  =           -0.6402d0*w4 +0.0702d0*w5 -0.0099d0*w6
!
      u(31)  =           -0.0374d0*w4
      f(31)  =  1.0004d0 -0.0373d0*w1 +0.0002d0*w2
!
      x  = 1.0d0-0.2505d0*dcos(a2)-0.1102d0*dcos(a3)                    &
     &     -0.0156d0*dcos(a4)-0.037d0*w1                   
      y  = -0.2505d0*dsin(a2)-0.1102d0*dsin(a3)                         &
     &     -0.0156d0*dsin(a4)-0.037d0*w4
!
      u(34)  = datan (y/x)
      if (x .lt. 0.0d0) u(34) = u(34)+pi
      f(34)  = dsqrt(x*x+y*y)
!
      u(38)  =           -0.3096d0*w4 +0.0119d0*w5 -0.0007d0*w6
      f(38)  =  1.0241d0 +0.2863d0*w1 +0.0083d0*w2 -0.0015d0*w3
!
      u(1)   =  0.0d0
      u(2)   =  0.0d0
      u(4)   = -u(31)
      u(6)   =  u(10)
      u(7)   =  u(10)
      u(8)   =  u(10)
      u(9)   =  u(10)
      u(11)  =  u(31)
      u(13)  =  u(21)
      u(14)  =  0.0d0
      u(15)  =  0.0d0
      u(16)  =  0.0d0
      u(18)  =  0.0d0
      u(19)  =  0.0d0
      u(20)  =  u(21)
      u(22)  = -u(10)
      u(24)  =  2.0d0*u(10)
      u(25)  =  2.0d0*u(31)
      u(26)  =  u(31)
      u(27)  =  u(31)
      u(28)  =  u(31)
      u(29)  =  u(31)
      u(30)  =  u(10)
      u(32)  =  u(31)+u(38)
      u(33)  =  u(31)
      u(35)  =  0.0d0
      u(36)  =  0.0d0
      u(37)  =  0.0d0
      u(39)  =  0.0d0
      u(40)  =  u(17)+u(21)
      u(41)  =  u(4)
      u(42)  =  u(31)+u(10)
      u(43)  =  u(31)*1.5
      u(44)  =  u(10)
      u(45)  =  u(31)+u(17)
      u(46)  =  u(17)
      u(47)  =  u(25)
      u(48)  =  u(25)
      u(49)  =  u(31)
      u(50)  =  u(31)
      u(51)  =  u(32)
      u(52)  =  0.0d0
      u(53)  =  u(38)
      u(54)  =  u(25)+u(31)
      u(55)  =  u(54)
      u(56)  =  u(25)
      u(57)  =  u(25)
      u(58)  =  u(25)+u(38)
      u(59)  =  u(31)
      u(60)  =  u(32)
      u(61)  =  0.0d0
      u(62)  =  u(54)-u(38)
      u(63)  =  u(54)
      u(64)  =  u(58)
      u(65)  =  u(31)-u(38)
      u(66)  = -u(38)
      u(67)  =  0.0d0
      u(68)  =  u(42)
      u(69)  =  u(25)
      u(70)  =  u(25)-u(10)
      u(71)  =  u(54)-u(38)
      u(72)  =  u(54)
      u(73)  =  u(25)-u(38)
      u(74)  =  u(54)-u(17)
      u(75)  =  2.5d0*u(31)
      u(76)  =  u(54)-u(10)
      u(77)  =  2.0d0*u(25)
      u(78)  =  u(77)
      u(79)  =  u(77)-u(38)
      u(80)  =  u(77)
      u(81)  =  u(71)
      u(82)  =  u(54)
      u(83)  =  u(71)
      u(84)  =  u(54)
      u(85)  =  u(25)
      u(86)  =  u(51)+u(34)
      u(87)  =  u(77)
      u(88)  =  u(77)
      u(89)  =  u(77)
      u(90)  =  u(54)
      u(91)  =  u(54)
      u(92)  =  u(54)+u(38)
      u(93)  =  u(58)
      u(94)  =  u(25)
      u(95)  =  u(58)
      u(96)  =  u(77)
      u(97)  =  u(54)
      u(98)  =  5.0d0*u(31)
      u(99)  =  u(98)
      u(100) =  u(77)
      u(101) =  u(25)
      u(102) =  u(73)
!
!     u(103),u(104) are changed according to
!     Dr cartwright's notes of nov 15,1977
!
      u(103) =  0.0d0
      u(104) =  0.0d0
      u(105) =  0.0d0
      u(106) = -u(65)
      u(107) =  u(54)
      u(108) =  u(25)
      u(109) =  0.0d0
      u(110) =  u(31)
      u(111) =  0.0d0
      u(112) =  0.0d0
      u(113) =  u(42)
      u(114) =  u(45)
      u(115) =  u(31)
!
!     Convert into degrees
!
      do k=1,115
        u(k)=mod(u(k)*deg,360.0d0)
        if (u(k) .lt. 0.0d0) u(k)=u(k)+360.0d0
      end do
!
      f(1)   =  1.0d0
      f(2)   =  1.0d0
      f(4)   =  f(31)
      f(6)   =  f(10)
      f(7)   =  f(10)
      f(8)   =  f(10)
      f(9)   =  f(10)
      f(11)  =  f(31)
      f(13)  =  f(21)
      f(14)  =  1.0d0
      f(15)  =  1.0d0
      f(16)  =  1.0d0
      f(18)  =  1.0d0
      f(19)  =  1.0d0
      f(20)  =  f(21)
      f(22)  =  f(10)
      f(24)  =  f(10)**2
      f(25)  =  f(31)**2
      f(26)  =  f(31)
      f(27)  =  f(31)
      f(28)  =  f(31)
      f(29)  =  f(31)
      f(30)  =  f(10)
      f(32)  =  f(31)*f(38)
      f(33)  =  f(31)
      f(35)  =  1.0d0
      f(36)  =  1.0d0
      f(37)  =  1.0d0
      f(39)  =  f(25)
      f(40)  =  f(17)*f(21)
      f(41)  =  f(31)
      f(42)  =  f(31)*f(10)
      f(43)  =  f(31)**1.5d0
      f(44)  =  f(10)
      f(45)  =  f(31)*f(17)
      f(46)  =  f(17)
      f(47)  =  f(25)
      f(48)  =  f(25)
      f(49)  =  f(31)
      f(50)  =  f(31)
      f(51)  =  f(32)
      f(52)  =  1.0d0
      f(53)  =  f(38)
      f(54)  =  f(25)*f(31)
      f(55)  =  f(54)
      f(56)  =  f(25)
      f(57)  =  f(25)
      f(58)  =  f(25)*f(38)
      f(59)  =  f(31)
      f(60)  =  f(32)
      f(61)  =  1.0d0
      f(62)  =  f(54)*f(38)
      f(63)  =  f(54)
      f(64)  =  f(58)
      f(65)  =  f(32)
      f(66)  =  f(38)
      f(67)  =  f(25)*f(25)
      f(68)  =  f(42)
      f(69)  =  f(25)
      f(70)  =  f(25)*f(10)
      f(71)  =  f(54)*f(38)
      f(72)  =  f(54)
      f(73)  =  f(58)
      f(74)  =  f(54)*f(17)
      f(75)  =  0.5d0*(f(25)+f(54))
      f(76)  =  f(54)*f(10)
      f(77)  =  f(67)
      f(78)  =  f(67)
      f(79)  =  f(67)*f(38)
      f(80)  =  f(67)
      f(81)  =  f(71)
      f(82)  =  f(54)
      f(83)  =  f(71)
      f(84)  =  f(54)*f(25)
      f(85)  =  f(67)
      f(86)  =  f(51)*f(34)
      f(87)  =  f(67)
      f(88)  =  f(67)
      f(89)  =  f(67)
      f(90)  =  f(54)
      f(91)  =  f(54)
      f(92)  =  f(71)
      f(93)  =  f(58)
      f(94)  =  f(25)
      f(95)  =  f(58)
      f(96)  =  f(67)
      f(97)  =  f(54)
      f(98)  =  f(84)
      f(99)  =  f(84)
      f(100) =  f(67)
      f(101) =  f(25)
      f(102) =  f(58)
!
!     f(103),f(104) are changed according to
!     Dr cartwright's notes of nov 15,1977
!
      f(103) =  1.0d0
      f(104) =  1.0d0
      f(105) =  1.0d0
      f(106) =  f(32)
      f(107) =  f(54)
      f(108) =  f(25)
      f(109) =  1.0d0
      f(110) =  f(54)
      f(111) =  1.0d0
      f(112) =  1.0d0
      f(113) =  f(42)
      f(114) =  f(45)
      f(115) =  f(54)
!
      return
!
      end
!
      subroutine vset (s,p,p1,h,v)
! ******************************************************************************
! *                                                                            *
! *                            FORTRAN SOURCE CODE                             *
! *                                                                            *
! *   PROGRAM SET :  UTILITY                           REF:JRH:21:09:2004      *
! *                                                                            *
! *   REVISION    :  -------------                         JRH:--:--:2004      *
! *                                                                            *
! *   SOURCE      :  tpd1604.f                                                 *
! *   ROUTINE NAME:  vset                                                      *
! *   TYPE        :  SUBROUTINE                                                *
! *                                                                            *
! *   FUNCTION    :  Returns v.                                                *
! *                                                                            *
! ******************************************************************************
!
      implicit none
!
!     Version 92.01   24nd january
!
!     Checked for 2000 compliance by dlb on 23.09.1997
!
      real*8 v(115)
!
      real*8 s,p,p1,h
!
      real*8 h2,h3,h4,p2,s2,s3,s4
!
      integer k
!
      h2=h+h
      h3=h2+h
      h4=h3+h
!
      s2=s+s
      s3=s2+s
      s4=s3+s
!
      p2=p+p
!
!     v's computed in degrees.
!
      v(1)   =  h
      v(2)   =  h2
      v(3)   =  s-p
      v(4)   =  s2-h2
      v(5)   =  s2
      v(6)   =  h-s4+p2+270.0d0
      v(7)   =  h3-s4+270.0d0
      v(8)   =  h-s3+p+270.0d0
      v(9)   =  h3-s3-p+270.0d0
      v(10)  =  h-s2+270.0d0
      v(11)  =  h3-s2+90.0d0
      v(12)  =  h-s+90.0d0
      v(13)  =  h3-s-p+90.0d0
      v(14)  =  p1-h2+270.0d0
      v(15)  =  270.0d0-h
      v(16)  =  180.0d0
      v(17)  =  h+90.0d0
      v(18)  =  h2-p1+90.0d0
      v(19)  =  h3+90.0d0
      v(20)  =  s-h+p+90.0d0
      v(21)  =  s+h-p+90.0d0
      v(23)  =  s2+h+90.0d0
      v(26)  =  h2-s4+p2
      v(27)  =  h4-s4
      v(28)  =  h2-s3+p
      v(29)  =  h4-s3-p
      v(31)  =  h2-s2
      v(32)  =  h4-s2
      v(33)  =  p-s+180.0d0
      v(34)  =  h2-s-p+180.0d0
      v(35)  =  p1-h
      v(36)  =  0.0d0
      v(37)  =  h-p1+180.0d0
      v(38)  =  h2
!
      v(22)  = -v(10)
      v(24)  =  v(10)+v(8)
      v(25)  =  v(31)+v(28)
      v(30)  =  v(10)+v(15)
      v(39)  =  v(31)-v(28)
      v(40)  =  v(17)+v(21)
      v(41)  = -v(31)
      v(42)  =  v(31)+v(10)
      v(43)  =  h3-s3+180.0d0
      v(44)  =  v(10)
      v(45)  =  v(31)+v(17)
      v(46)  =  v(17)
      v(47)  =  v(25)
      v(48)  =  v(31)+v(31)
      v(49)  =  v(28)
      v(50)  =  v(31)
      v(51)  =  v(31)+v(38)
      v(52)  =  0.0d0
      v(53)  =  v(38)
      v(54)  =  v(48)+v(28)
      v(55)  =  v(48)+v(31)
      v(56)  =  v(47)
      v(57)  =  v(48)
      v(58)  =  v(48)+v(38)
      v(59)  =  v(31)
      v(60)  =  v(51)
      v(61)  =  v(54)
      v(62)  =  v(55)-v(38)
      v(63)  =  v(55)
      v(64)  =  v(25)+v(38)
      v(65)  =  v(28)-v(38)
      v(66)  =  -v(38)
      v(67)  =  v(48)-v(28)-v(28)
      v(68)  =  v(31)+v(8)
      v(69)  =  v(48)-v(15)
      v(70)  =  v(48)-v(8)
      v(71)  =  v(62)
      v(72)  =  v(55)
      v(73)  =  v(48)-v(38)
      v(74)  =  v(55)-v(17)
      v(75)  =  v(31)+v(43)
      v(76)  =  v(55)-v(10)
      v(77)  =  v(54)+v(28)
      v(78)  =  v(55)+v(28)
!
      v(89)  =  v(48)+v(48)
!
      v(79)  =  v(89)-v(38)
      v(80)  =  v(89)
      v(81)  =  v(54)-v(38)
      v(82)  =  v(48)+v(29)
      v(83)  =  v(62)
      v(84)  =  v(89)-v(28)
      v(85)  =  v(55)-v(28)
      v(86)  =  v(51)+v(34)
      v(87)  =  v(77)
      v(88)  =  v(78)
!
      v(90)  =  v(54)
      v(91)  =  v(55)
      v(92)  =  v(55)+v(38)
      v(93)  =  v(64)
      v(94)  =  v(48)
      v(95)  =  v(58)
      v(96)  =  v(89)
      v(97)  =  v(55)
      v(98)  =  v(89)+v(28)
      v(99)  =  v(89)+v(31)
      v(100) =  v(89)
      v(101) =  v(31)+v(29)
      v(102) =  v(73)
!
!     v(103),v(104) are changed according to Dr cartwrights
!     notes of 15 nov,1977
!
      v(103) =  v(31)-h
      v(104) =  v(31)+h
      v(105) =  v(31)-v(29)
      v(106) =  v(38)-v(31)
      v(107) =  v(54)
      v(108) =  v(101)
      v(109) =  v(85)
      v(110) =  v(34)
      v(111) =  v(28)+v(35)
      v(112) =  v(28)-v(35)
      v(113) =  v(42)-270.0d0
      v(114) =  v(45)-90.0d0
!     v(115) =  v(48)-v(28)
      v(115) =  v(34)
!
      do k=1,115
        v(k)=dmod(v(k),360.0d0)
        if (v(k) .lt. 0.0d0) v(k)=v(k)+360.0d0
      end do
!
      return
!
      end
!
#endif
