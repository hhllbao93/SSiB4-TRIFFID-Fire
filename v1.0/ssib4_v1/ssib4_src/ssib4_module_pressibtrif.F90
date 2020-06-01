MODULE ssib4_module_pressibtrif      
      use ssib4_module_comsveg
      use ssib4_module_trifparms
      use ssib4_module_comsveg
      use ssib4_module_comssib4sfc
      use ssib4_module_vegcnvt
!
      IMPLICIT NONE
!
contains
       subroutine  PROG_expset0(month,xlon,xlat,ITYPE,FRAC,LAI,cs)
!
!-----------------------------------------------------------------------
       type(ssibdataios)::datafs
       integer:: MONTH,ITYPE,N
       real:: xlon,xlat
       real:: VEGCOVER,CS
       real,dimension(NTYPES)::FRAC
       real,dimension(6)::LAI,FRACT0
       data FRACT0 /1.0, 0.9 ,0.9 ,0.8,0.4,0.5/
!
         cs=10.0
          VEGCOVER=0.0
         do N=1,6
         if (N.eq.ITYPE) then
          FRAC(N)=FRACT0(ITYPE)
          else
          FRAC(N)=0.01
          end if
          VEGCOVER=VEGCOVER+FRAC(N)
        end do
        FRAC(7)=1-VEGCOVER

        DO N=1,6
          LAI(N)=3.0
        ENDDO

      END SUBROUTINE
!-----------------------------------------------------------------------
!
      subroutine  PROG_expset1(datafs,xlon,xlat,month,ITYPE,FRAC,LAI,CS)
!
!-----------------------------------------------------------------------
       type(ssibdataios)::datafs
       integer:: month,ITYPE,N
       real:: xlon,xlat
       real::VEGCOVER,CS
       real,dimension(NTYPES)::FRAC
!yliu09May2017 addtype       real,dimension(6)::LAI,FRACT0
       real,dimension(NTYPET)::LAI,FRACT0
!yliu09May2017 addtype       data FRACT0 /1.0, 0.9 ,0.9 ,0.8,0.4,0.5/
       data FRACT0 /1.0, 0.9 ,0.9 ,0.8,0.4,0.5,0.9/
!
         cs=10.0
         VEGCOVER=0.0
!yliu09May2017 addtype         DO N=1,6
         DO N=1,NTYPET
!zzq          FRAC(N)=0.14
!zzq          LAI(N)=4.0
!zzq 04-24-2008  both of the following can not be zero for TRIF
!
            FRAC(N)=1.0E-6
            LAI(N)=1.0E-6
!zzq            LAI(N)=LAI_MIN(N)
          if(ITYPE.eq.N) THEN
             FRAC(N) = VCOVER0(ITYPE,MONTH,1)
             LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
!zzq 2008-04-14
             if(itype.eq.5) then
               FRAC(n)=0.2
               LAI(n)=1.5
             endif
             if(itype.eq.3) LAI(n)=1.5
             if(itype.eq.6) LAI(n)=1.5
!yliu DBL 
             if(itype.eq.7) LAI(n)=5.0
          ENDIF
!
!zzq 2009May19
!zzq 2012Sep06          FRAC(N)=1/6.
          VEGCOVER=VEGCOVER+FRAC(N)
        ENDDO
        FRAC(8)=1-VEGCOVER
!
       call changeVT8(datafs,xlon,xlat,MONTH,ITYPE,FRAC,LAI)
!
      END SUBROUTINE
!-----------------------------------------------------------------------
!
      subroutine  PROG_expset2(datafs,xlon,xlat,month,ITYPE,FRAC,LAI,CS)
!
!-----------------------------------------------------------------------
       type(ssibdataios)::datafs
       integer:: month,ITYPE,N
       real:: xlon,xlat
       real::VEGCOVER,CS
       real,dimension(NTYPES)::FRAC
       real,dimension(NTYPET)::LAI,FRACT0
       data FRACT0 /1.0, 0.9 ,0.9 ,0.8,0.4,0.5,0.9/
!
         cs=10.0
         VEGCOVER=0.0
         DO N=1,NTYPET
            LAI(N)=1.0E-6
          if(ITYPE.eq.N) THEN
             LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
             if(itype.eq.5) then
               LAI(n)=1.5
             endif
             if(itype.eq.3) LAI(n)=1.5
             if(itype.eq.6) LAI(n)=1.5
!yliu DBL 
             if(itype.eq.7) LAI(n)=5.0
          ENDIF
!
!zzq 2009May19
!zzq 2012Sep06          FRAC(N)=1/6.
          VEGCOVER=VEGCOVER+FRAC(N)
        ENDDO
        FRAC(8)=1-VEGCOVER
!
      END SUBROUTINE
!=======================================================================
!=======================================================================
!
       subroutine init_ssib(datafs,MONTH,tm,em,tdd,prec)
!=======================================================================
!
      type(ssibdataios):: datafs
      integer:: N,MONTH,DAY     !Huilin add DAY Aug. 2019
      real:: tm,em,tdd,prec     !Huilin add prec Aug. 2019
!
       DO N=1,NTYPES
          datafs % td0(N) =tdd
          datafs % tg0(N) =tm
          datafs % tc0(N) =tm
          datafs % tdm(N) =tdd
          datafs % tgm(N) =tm
          datafs % tcm(N) =tm
          datafs % w01(N) =0.5
          datafs % w02(N) =0.5
          datafs % w03(N) =0.5
          datafs % wm1(N) =  datafs % w01(N)
          datafs % wm2(N) =  datafs % w02(N)
          datafs % wm3(N) =  datafs % w03(N)
          datafs % cap01(N) =0.0
          datafs % cap02(N) =0.0
          datafs % capm1(N) =0.0
          datafs % capm2(N) =0.0
          datafs % xta(N)=tm
          datafs % xea(N)=em
!yliu09May2017          datafs % xZFW(N) =1.0E-6
          datafs % xZFW(N) =0.8
          datafs % xCO2RA(N)=RBC0(N,MONTH)
          datafs % xCO2RB(N)=RBC0(N,MONTH)
          datafs % xfiltr(N)=0.0
       ENDDO
!yliu add for ssib4 ini only due to no data in the South Pole
      if (datafs%xvtype.eq.NTYPES) then
      DO N=1,NNPFT
         datafs%xveg_fracn(N) = 0.0
         datafs%xlai_trif(N)  = 0.0
         datafs%xHT_trif(N)   = 0.0
      ENDDO
      endif
      DO DAY=1,10
         datafs % xprec10(n) = prec
      ENDDO
      DO DAY=1,60
         datafs % xprec60(n) = prec
      ENDDO
      datafs % xprecd        = prec
!
      END SUBROUTINE
!=======================================================================
!
       subroutine init_trif(datafs,xlon,xlat,DTT,MONTH,IDAY,ITYPE)
!
!=======================================================================
!
      type(ssibdataios)::datafs
      real:: Cs,VEG_FRAC,FRAC_VS,WOOD,LEAF     !Huilin add LEAF initialization
      real,dimension(NPFT)  :: HT,LAI,SAI,GREEN_TRIF,WOOD1,LEAF1,ROOT1   !Huilin add LEAF1 ROOT1 WOOD1 initialization
      real,dimension(NTYPES):: FRAC_TRIF
      real:: xlon,xlat

      real:: DTT
      INTEGER:: ITYPE,MONTH,IDAY,N
      integer:: KTRIFDAY,ITRIFSTEP
!
!         CALL PROG_expset1(datafs,xlon,xlat,month,ITYPE,FRAC_TRIF,LAI,cs)
      do N=1,NTYPES
         FRAC_TRIF(N) = datafs % xVEG_FRACn(N)
      enddo
         CALL PROG_expset2(datafs,xlon,xlat,month,ITYPE,FRAC_TRIF,LAI,cs)

         KTRIFDAY = IDAY 
         ITRIFSTEP = 1
         datafs % xKTRIFDAY = KTRIFDAY
         datafs % xITRIFSTEP = ITRIFSTEP
         datafs % xCo2CS=CS
         datafs % xRESP_S_DR = 0.0
         datafs % xFRAC_VS = 1.0
!Huilin add fire output initialization
         datafs % xLEAF_LOSS_DR = 0.0
         datafs % xWOOD_LOSS_DR = 0.0
         datafs % xROOT_LOSS_DR = 0.0
!
!-----------------------------------------------------------------------
! Diagnose a Stem Area Index (SAI) using the balanced (Maximum) LAI
! from TRIFFID, and assuming that SAI is 0.5 of this value (consistent
! with BATS parameter table)
!-----------------------------------------------------------------------
!
       DO N=1,NPFT
          WOOD = A_WL(N)*(LAI(N)**B_WL(N))
!Huilin add leaf wood1 leaf1
          LEAF = SIGL(N)*LAI(N) 
          WOOD1(N) = WOOD
          LEAF1(N) = LEAF
          ROOT1(N) = LEAF
          HT(N) = WOOD / (A_WS(N) * ETA_SL(N))            &
                  * (A_WL(N)/WOOD)**(1.0/B_WL(N))
          SAI(N) = DSAI_DLAI(N) *                         &
                  ((A_WS(N)*ETA_SL(N)*HT(N)/A_WL(N))      &
                         **(1.0/(B_WL(N)-1)))
          green_trif(N)=1-SAI(N)/(LAI(N)+SAI(N))
       ENDDO
!
       DO N=1,NPFT
          datafs % xGPP_DR(N)= 0.0
          datafs % xNPP_DR(N)= 0.0
          datafs % xHT_TRIF(N) =HT(N)
          datafs % xLAI_TRIF(N)=LAI(N)
          datafs % xSAI_TRIF(N)=SAI(N)
          datafs % xGREEN_TRIF(N)=GREEN_TRIF(N)
          datafs % xRESP_W_DR(N)= 0.0
          datafs % xG_LEAF_DAY(N)= 0.0
          datafs % xG_LEAF_DR(N)= 0.0
          datafs % xPHEN(N)=0.0
          datafs % xfsmc(N)=0.5
          datafs % xfsmc_ave(N)=0.0
!yliu 16Jun2016
          datafs % xGREEN_SSIB(N)=GREEN_TRIF(N)
!Huilin add fire output initialization
          datafs % xleaf1(N)= LEAF1(N)
          datafs % xwood1(N)= WOOD1(N)
          datafs % xroot1(N)= ROOT1(N)
       ENDDO
!
      VEG_FRAC=0.0
      do n=1,NNPFT
          datafs % xVEG_FRACn(N)=FRAC_TRIF(N)
          VEG_FRAC=VEG_FRAC+FRAC_TRIF(N)
       enddo
        datafs % xVEG_FRAC = VEG_FRAC

!
      END SUBROUTINE
!-----------------------------------------------------------------------
      subroutine trif_prescribed(DDTT)
!-----------------------------------------------------------------------
     INTEGER:: ITYPE,ity
     REAL,parameter::                 &
       DAY_YEAR = 365.0,              & ! Number of days in a year (days).
       SEC_DAY  = 86400.,             & ! Number of seconds in a day (s).
       SEC_YEAR = DAY_YEAR*SEC_DAY      ! Number of seconds in a year (s).

    REAL C1,C11,DDTT
!-----------------------------------------------------------------------
! Defaults for run
!-----------------------------------------------------------------------
      TRIF_ON=.TRUE.
      PHEN_ON=.TRUE.
      STEP_DAY = SEC_DAY/DDTT
      VEG_EQUIL=.FALSE.
      DAY_TRIF=10
!zzq      VEG_EQUIL=.TRUE.
!zzq      DAY_TRIF=DAY_YEAR
      DAY_PHEN=1
      ANETC_DEF_KGC=2.0

!------------------------------------------------------------------------
! Set default ANETC and RDC to 20% of this (???)
!------------------------------------------------------------------------
      ANETC_DEF_KGC=ANETC_DEF_KGC/SEC_YEAR/12.0E-3
!zzq,05-12-2008   V_SAT=poros0(ITYPE)

!-----------------------------------------------------------------------
! Initialise defaults for deriving variables,
! prognostic variables and accumulations
!-----------------------------------------------------------------------
!zzq,05-09-2008    C1=-EXP(PH20(ITYPE,1))
!zzq,05-09-2008    C11=C1/PHSAT0(ITYPE)
!zzq,05-09-2008    v_wilt=EXP(-ALOG(C11)/BEE0(ITYPE))
!
       do ity=1,NTYPES
         C1=-EXP(PH20(ity,1))
         C11=C1/PHSAT0(ity)
         v_wilt0(ity)=EXP(-ALOG(C11)/BEE0(ity))
         V_SAT0(ity)=poros0(ity)
       enddo
!
!      print*,ITYPE,'wilt=',V_wilt0
!
      FRACA=0.0
      G_ANTH=0.0
!
      END SUBROUTINE
!=======================================================================
!
      SUBROUTINE VARINIT(TC,TGS,TD,WWW,CAPAC,YEAR,MONTH,DAY,iday,time,  &
                 ZLAT,zlon,                                             &
                 TC0,TGS0,TD0,WWW0,CAPAC0,TCM,TGSM,TDM,WWWM,CAPACM,     &
                 ITYPE, ZWINDN,DTT,TAN,itrifstep,ktrifday)
!                                                         12 AUGUST 2000
!=======================================================================
!     SET INITIAL CONDITIONS

!-----------------------------------------------------------------------
! Time parameters
!-----------------------------------------------------------------------
      integer:: itrifstep,ktrifday
      integer:: iday
      integer:: NPFT1,ITYPE
      real:: time
      real:: tc,tgs,td,zlat,zlon
      real:: month,day,year,dtt

      real,dimension(2)::CAPAC
      real,dimension(3)::WWW
      real,dimension(2,7)::CAPAC0,CAPACM
      real,dimension(3,7)::WWW0,WWWM
      real,dimension(7)::ZWINDN,TAN
      real,DIMENSION(7):: TC0,TGS0,TD0,TCM,TGSM,TDM

       do NPFT1=1,7
      CAPAC0(1,NPFT1)=0.
      CAPAC0(2,NPFT1)=0.
      WWW0(1,NPFT1)= 0.5
      WWW0(2,NPFT1)= 0.5
      WWW0(3,NPFT1)= 0.5
      TC0(NPFT1)=293.
      TGS0(NPFT1)=293.
      TD0(NPFT1)=285.
      CAPACM(1,NPFT1)=0.
      CAPACM(2,NPFT1)=0.
      WWWM(1,NPFT1)= 0.5
      WWWM(2,NPFT1)= 0.5
      WWWM(3,NPFT1)= 0.5
      TCM(NPFT1)=293.
      TGSM(NPFT1)=293.
      TDM(NPFT1)=285.
      TAN(NPFT1)=293.0
      end do
      ZWINDN(1)=64.0
      ZWINDN(2)=30.0
      ZWINDN(3)=2.0
      ZWINDN(4)=6.0
      ZWINDN(5)=6.0
      ZWINDN(6)=6.0
      ZWINDN(7)=2.0
      CAPAC(1)=0.
      CAPAC(2)=0.
      WWW(1)= 0.5
      WWW(2)= 0.5
      WWW(3)= 0.5
      TC=293.
      TGS=293.
      TD=285.
!============================================================
!
!     THE FOLLOWING LINES ARE NOT NEEDED FOR GCM
!
!=============================================================
!
      ZLAT=56.0    ! North
      ZLON=-98.30  ! West

      MONTH=11
      DAY=335.
      iday=1
      time=0.0
      YEAR=2002.
      ITYPE=2
      DTT=1800.

      itrifstep=1
      ktrifday=335
!
      END SUBROUTINE
END MODULE
