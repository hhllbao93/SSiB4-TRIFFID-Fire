MODULE module_radc2
      use ssib4_module_comsconst
      use module_time
!
      implicit none
contains
!=======================================================================
!
      SUBROUTINE RADC2_cosz(loctmod,COSZ,IYEAR,DAY1,TIME1,DTT,ZLON1,ZLAT)
!
!=====================================================================
!
! ** SOLAR ZENITH ANGLE COMPUTATION; DOWNCOMING RADIATION AT BOTTOM.
!
!=====================================================================
       integer:: loctmod,nDAYSPY
       real:: year,day,time,zlon
       real:: COSZ,DAY1,TIME1,DTT,ZLON1,ZLAT
       real:: DECMAX,SOLS,SEASON,DEC,RFD,SS,COSHR,ANGLE
       real:: SIND,COSD,HAC,H,DAWN,DUSK,SR
       integer:: IYEAR
!
! ** SOLAR ZENITH ANGLE COMPUTATION; DOWNCOMING RADIATION AT BOTTOM.
! ** SOLAR DECLINATION CALCULATION
!
!zzq 10-30-2008, calculate local time
!
!------------------------------------------------------------
! loctmod =1 for that the input time is GMT,
! otherwise the input is the local time.
!
      year=IYEAR
      day=day1
      time =time1
      zlon=zlon1
      if(zlon.lt.0) zlon=360+zlon
!
      if(loctmod.eq.1) then
         time=time1+ZLON/15.
         call daynumofyear(iyear,nDAYSPY)
!
         IF ( TIME .GE. 23.999) THEN
           TIME = TIME-24.
           DAY = DAY + 1
           IF ( DAY .GT. nDAYSPY ) then
            YEAR = YEAR + 1
            DAY = DAY - nDAYSPY
           ENDIF
         ENDIF
      endif
!------------------------------------------------------------
!
! ** SOLAR DECLINATION CALCULATION
!
      DECMAX = PIE * ( 23.5 / 180.)
      SOLS   = ( 4141./24. ) + MOD( int(YEAR)+3, 4) * 0.25
!
      SEASON = ( DAY - SOLS ) / 365.2
      DEC    = DECMAX * COS ( 2. * PIE * SEASON )
!
      RFD  = PIE / 180.
      SIND = SIN( DEC )
      COSD = COS( DEC )
      HAC  = -TAN( ZLAT * RFD )*TAN( DEC )
      HAC  = AMIN1(HAC,1.0)
      HAC  = AMAX1(HAC,-1.0)
!
! **  H IS THE HALF-DAY LENGTH (IN RADIANS)
!
      H   = ACOS(HAC)
      DAWN= -H
      DUSK= +H
      SR  = 12.-(H/(15.*RFD))
      SS  = 12.+(H/(15.*RFD))
!
!     WRITE(6,53)SR,SS
!  53 FORMAT(1X,' SUNRISE IS AT ',F8.2,' HOURS AND SUNSET IS AT ',F8.2,
!    &' HOURS .')
!
! ** CALCULATION OF SOLAR ANGLE (SUNANG)
!
      COSHR = COS( - PIE + (TIME + 0.5*DTT/3600.) / 24. * 2. * PIE )
!
      ANGLE = SIN( ZLAT*RFD ) * SIND + COS ( ZLAT*RFD ) * COSD * COSHR
      ANGLE = ( ANGLE + ABS(ANGLE))/2.
      COSZ= ANGLE + 0.01745 * ( 1. - ANGLE )
!
      COSZ = AMAX1(0.01,COSZ)
!
      END SUBROUTINE
!=========================================================================
!
      function nstep_index(nyy,nmm,ndd,tdayhr,dtt) result(f_result)
!
!=========================================================================
!
      integer:: nyy,nmm,ndd
      real:: tdayhr,dtt
      integer:: jday,ierr,f_result
!
      call datetojday(nyy,nmm,ndd,jday,ierr)
      f_result=int(((jday-1)*24+tdayhr)*3600.0/dtt)+1
!
      end function
!=========================================================================
!
      function nstep_index_max(nyy,nstart,endhour,thour,dtt)  &
           result(f_result)
!
!=========================================================================
!
      integer:: nyy,nstart
      real:: endhour,thour,dtt
      integer ::nyrstepmax,n,f_result

      nyrstepmax =nstepofyear(nyy,dtt)
      n  = nstart + (endhour-thour)*3600./dtt
      if(n.gt.nyrstepmax) n  = nyrstepmax
      f_result = n
!
      end function

!=======================================================================
!
      SUBROUTINE RADC2(COSZ,YEAR,DAY,TIME,DTT,ZLAT)
!
!=====================================================================
!
! ** SOLAR ZENITH ANGLE COMPUTATION; DOWNCOMING RADIATION AT BOTTOM.
!
      integer:: year,iter
      real:: dayspy,ss,time,dtt,day
      real:: DECMAX,SOLS,SEASON,DEC,RFD,SIND,COSD,  &
           HAC,DAWN,DUSK,SR,COSHR,ANGLE,            &
           ZLAT,H,COSZ

!
      DAYSPY = 365.
      IF ( MOD( int(YEAR), 4 ) .EQ. 0 ) DAYSPY = 366.
!zzq      DAYSPY = 366.
!zzq      IF ( MOD( YEAR, 4 ) .EQ. 0 ) DAYSPY = 367.
!
! ** JULIAN DAY AND TIME UPDATE; SKIP ON 1ST TIME STEP (INITIALIZED)
!
!zzq oct 09,2008
      iter = 0 
!
      IF(ITER .EQ. 1)GO TO 10
      TIME = TIME + DTT / 3600.
      IF ( TIME .GE. 23.999) TIME = 0.0
      DAY = DAY +  DTT / 86400.
!
   10 CONTINUE
!
      IF ( DAY .GT. DAYSPY ) YEAR = YEAR + 1
      IF ( DAY .GT. DAYSPY ) DAY = DAY - DAYSPY
!
! ** SOLAR DECLINATION CALCULATION
!
      DECMAX = PIE * ( 23.5 / 180.)
      SOLS   = ( 4141./24. ) + MOD( int(YEAR)+3, 4) * 0.25
!
      SEASON = ( DAY - SOLS ) / 365.2
      DEC    = DECMAX * COS ( 2. * PIE * SEASON )
!
      RFD  = PIE / 180.
      SIND = SIN( DEC )
      COSD = COS( DEC )
      HAC  = -TAN( ZLAT * RFD )*TAN( DEC )
      HAC  = AMIN1(HAC,1.0)
      HAC  = AMAX1(HAC,-1.0)
!
! **  H IS THE HALF-DAY LENGTH (IN RADIANS)
!
      H   = ACOS(HAC)
      DAWN= -H
      DUSK= +H
      SR  = 12.-(H/(15.*RFD))
      SS  = 12.+(H/(15.*RFD))
!
!     WRITE(6,53)SR,SS
!  53 FORMAT(1X,' SUNRISE IS AT ',F8.2,' HOURS AND SUNSET IS AT ',F8.2,
!    &' HOURS .')
!
! ** CALCULATION OF SOLAR ANGLE (SUNANG)
!
      COSHR = COS( - PIE + (TIME + 0.5*DTT/3600.) / 24. * 2. * PIE )
!
      ANGLE = SIN( ZLAT*RFD ) * SIND + COS ( ZLAT*RFD ) * COSD * COSHR
      ANGLE = ( ANGLE + ABS(ANGLE))/2.
      COSZ= ANGLE + 0.01745 * ( 1. - ANGLE )
!
      COSZ = AMAX1(0.01,COSZ)
!
      END SUBROUTINE
!
!=====================================================================
!
        SUBROUTINE VEGMASK(SSIBMASK,SLMSK,NMAP,GWWW)
!
!=====================================================================
!       READ IN AND CHECK THE VEGETATION MAP       
        INTEGER,PARAMETER::LONB=360,LATB=180
        REAL,DIMENSION(LONB,LATB):: SSIBMASK,SLMSK
        REAL,DIMENSION(LONB,LATB,3)::GWWW
        INTEGER::NMAP,I,J
        REAL:: ZERO

        READ(NMAP) SSIBMASK
        ZERO=0.
        DO 121 J=1,LATB
        DO 121 I=1,LONB
        SSIBMASK(I,J)=NINT(SSIBMASK(I,J))
        IF(SLMSK(I,J).EQ.ZERO.AND.SSIBMASK(I,J).NE.ZERO) THEN
          WRITE(7,*) ' GCM MASK IS OCEAN AND SSIB MASK LAND AT I,J ',I,J
          WRITE(7,*) ' GCM MASK IS USED '
          WRITE(7,*) ' SUSPECT NMC LAND SEA MASK IS INCORRECT '
          SSIBMASK(I,J) = 0.0
        ENDIF
        IF(SLMSK(I,J).EQ.1.0.AND.SSIBMASK(I,J).EQ.ZERO) THEN
          WRITE(7,*) ' GCM MASK IS LAND AND SSIB MASK OCEAN AT I,J ',I,J
          WRITE(7,*) ' GROUNDCOVER PLANTED (VEG TYPE 7) '
          WRITE(7,*) ' SUSPECT GCM LAND SEA MASK IS INCORRECT '
          IF ((J.LE.2).OR.(J.GE.(LATB-1)).OR.(I.LE.2).OR.   &
              (I.GE.(LONB-1))) THEN
!      THE FOLLOWING LINE CAN BE MODIFIED  ACCORDING TO THE LATITUDE
!             IF (J.LE.12) THEN
!                SSIBMASK(I,J) = 10
!             ELSE
                 SSIBMASK(I,J)=7
                 GWWW(I,J,1)=0.5
                 GWWW(I,J,2)=0.5
                 GWWW(I,J,3)=0.5
!             END IF
              GO TO 2006
          END IF
           IF (SSIBMASK(I+1,J).GE.1) THEN
               SSIBMASK(I,J)=SSIBMASK(I+1,J)
                 GWWW(I,J,1)=SSIBMASK(I+1,J)
                 GWWW(I,J,2)=SSIBMASK(I+1,J)
                 GWWW(I,J,3)=SSIBMASK(I+1,J)
            ELSEIF (SSIBMASK(I-1,J).GE.1) THEN
               SSIBMASK(I,J)=SSIBMASK(I-1,J)
                 GWWW(I,J,1)=SSIBMASK(I-1,J)
                 GWWW(I,J,2)=SSIBMASK(I-1,J)
                 GWWW(I,J,3)=SSIBMASK(I-1,J)
            ELSEIF (SSIBMASK(I,J-1).GE.1) THEN
               SSIBMASK(I,J)=SSIBMASK(I,J+1)
                 GWWW(I,J,1)=SSIBMASK(I,J+1)
                 GWWW(I,J,2)=SSIBMASK(I,J+1)
                 GWWW(I,J,3)=SSIBMASK(I,J+1)
            ELSEIF (SSIBMASK(I,J+1).GE.1) THEN
               SSIBMASK(I,J)=SSIBMASK(I,J-1)
                 GWWW(I,J,1)=SSIBMASK(I,J-1)
                 GWWW(I,J,2)=SSIBMASK(I,J-1)
                 GWWW(I,J,3)=SSIBMASK(I,J-1)
            ELSEIF (SSIBMASK(I,J+2).GE.1) THEN
               SSIBMASK(I,J)=SSIBMASK(I,J+2)
                 GWWW(I,J,1)=SSIBMASK(I,J+2)
                 GWWW(I,J,2)=SSIBMASK(I,J+2)
                 GWWW(I,J,3)=SSIBMASK(I,J+2)
            ELSEIF (SSIBMASK(I,J-2).GE.1) THEN
               SSIBMASK(I,J)=SSIBMASK(I,J-2)
                 GWWW(I,J,1)=SSIBMASK(I,J-2)
                 GWWW(I,J,2)=SSIBMASK(I,J-2)
                 GWWW(I,J,3)=SSIBMASK(I,J-2)
            ELSE
               SSIBMASK(I,J)=7
                 GWWW(I,J,1)=0.5
                 GWWW(I,J,2)=0.5
                 GWWW(I,J,3)=0.5
            ENDIF
        ENDIF
 2006        CONTINUE
!       THE FOLLOWING LINE CHECKS THE LNAD ICE
        IF(SLMSK(I,J).EQ.2.0.AND.SSIBMASK(I,J).NE.ZERO) SSIBMASK(I,J)=13.0
  121   CONTINUE
!
! ..............................................................
      END SUBROUTINE
END MODULE
