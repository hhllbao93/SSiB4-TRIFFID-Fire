!
!This module is designed mainly for model time processing.
!  zzq, sept 28,2013
!
MODULE module_time
!                                                                       
      implicit none

contains
!====================================================================
      SUBROUTINE RADC2_time(YEAR,month,DAY,IDAY,TIME,DTT)
!
!=====================================================================
!
! ** SOLAR ZENITH ANGLE COMPUTATION; DOWNCOMING RADIATION AT BOTTOM.
!
      integer::year,month,iday,ndom,nDAYSPY
      real:: time,day,dtt
!
      call daynumofyear(year,nDAYSPY)
!
! ** JULIAN DAY AND TIME UPDATE; SKIP ON 1ST TIME STEP (INITIALIZED)
!
      TIME = TIME + DTT / 3600.
      IF ( TIME .GE. 23.999) THEN
           TIME = 0.0
           iday =iday +1
           DAY = DAY + 1
           IF ( DAY .GT. nDAYSPY ) then
            YEAR = YEAR + 1
            DAY = DAY - nDAYSPY
           ENDIF
       ENDIF
!
!       *** check for end of month
!
      call daynumofmonth(year,month,ndom)
      if(iday.gt.ndom) then
          month=month+1
          if(month .gt.12) month=1
          iday=1
      endif
!
      end subroutine
!
!====================================================================
      subroutine timeforward(dhour,nyy,nmm,ndd,dayhr,ierr)
!
!=====================================================================
!
      real dhour,dayhr
      integer nyy,nmm,ndd,ndom
      integer ierr
      integer ntotdays
      real tothrs
!
      ierr = 0
      if(dayhr.gt.24.or.ndd.gt.31.or.nmm.gt.12) then
         ierr=1
      endif
!
      tothrs=dhour+dayhr
      ntotdays =int(tothrs/24) 
      dayhr = tothrs - ntotdays*24
      ndd  =ndd + ntotdays
      call daynumofmonth(nyy,nmm,ndom)
!
      do while(ndd.gt.ndom)
        ndd = ndd - ndom
        nmm = nmm+1
        if(nmm.gt.12) then
            nmm=1
            nyy = nyy+1
        endif
        call daynumofmonth(nyy,nmm,ndom)
      enddo
!
      end subroutine
!=========================================================================

      subroutine jdaytodate(nyy,jday,nmm,ndd,ierr)
!
!=========================================================================
!     converts Julian day to month and date:
      integer nyy,jday
      integer nmm,ndd,ierr
      integer ndom

      ierr = 0
      if(jday.eq.0) then
         ierr =1
         nmm=0
         ndd=0
         return
      endif
!
      ndd = jday
      call daynumofmonth(nyy,nmm,ndom)
!
      do while(ndd.gt.ndom)
        ndd = ndd - ndom
        nmm = nmm+1
        if(nmm.gt.12) then
            nmm=1
            nyy = nyy+1
        endif
        call daynumofmonth(nyy,nmm,ndom)
      enddo

      end subroutine
!=========================================================================

      subroutine datetojday(nyy,nmm,ndd,jday,ierr)
!
!=========================================================================
!     converts date to Julian day:
      integer nyy,nmm,ndd
      integer jday,ierr
      integer im,ndom

      ierr = 0
         
      if(nmm.gt.12.or.ndd.gt.31) then
        ierr=1
      endif
    
      jday =0
      do im=1,nmm-1
       call daynumofmonth(nyy,im,ndom)
       jday =jday +ndom
      enddo
      jday = jday +ndd
      
      end subroutine
!=========================================================================
!
      subroutine daynumoffeb(nyy,ndom)
!
!=========================================================================
!
      integer nyy
      integer ndom
      integer i1,i2,i3

      i1=mod(nyy,4)
      i2=mod(nyy,100)
      i3=mod(nyy,400)

      if((i1.eq.0).and.(i2.ne.0.or.i3.eq.0)) then
        ndom=29
      else
        ndom=28
      endif
!
      end subroutine
!=========================================================================
!
      subroutine daynumofmonth(nyy,nmm,ndom)
!
!=========================================================================
      integer nyy,nmm
      integer ndom
      integer,dimension(12):: nmth     
      data nmth/31,28,31,30,31,30,31,31,30,31,30,31/

      if(nmm.lt.0.or.nmm.gt.12) then
         ndom = 0
         return
       endif

       ndom = nmth(nmm)
       if(nmm.eq.2) then
           call daynumoffeb(nyy,ndom)
       endif

      end subroutine
!=========================================================================
!
      subroutine daynumofyear(nyy,ndoy)
!
!=========================================================================
!
      integer nyy
      integer ndoy
      integer ndoffeb
!
      call daynumoffeb(nyy,ndoffeb)
      ndoy=337+ndoffeb
!
      end subroutine
!=========================================================================
!
      function nstepofyear(nyy,dtt) result(f_result)
!
!=========================================================================
      integer nyy ! year 
      real dtt    ! integration step in seconds
      integer nstpoy,f_result

      if(mod(nyy,4).eq.0) then
           nstpoy =366*24*3600/dtt
      else
           nstpoy =365*24*3600/dtt
      endif
!
      f_result = nstpoy
!
      end function
END MODULE
