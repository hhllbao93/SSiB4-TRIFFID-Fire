!
! This code is transformed from F77, zzq, Sept26,2013
!
MODULE ssib4_module_ssibsfcbio
       use ssib4_module_comsconst
       use ssib4_module_ssibsub
       use ssib4_module_ssibco2
!
contains
      SUBROUTINE SFCBIO(                                &
!in ......................................................    
         DDTT,SUNANG, RADFRAC,RADN,RLWDOWN,             &
         PPL,PPC,TM,UMM1,VMM1,QM,EM,RHOAIR,             &
         PSURF,SIG,ZWIND,co2ca,L_VEG,ZLAT,MONTH,        &
! vegetation parameters
         TRAN,REF,GREEN,VCOVER,CHIL,                    &
         RSTPAR,VM0,TOPT,TL,TU,fc3,DEFAC,PH1,PH2,       &
         ZLT,Z0,XDD,Z2,Z1,RDC,RBC,ROOTD,SOREF,BEE,      &
         PHSAT,SATCO,SLOPE,POROS,ZDEPTH,                &
!inout.....................................................
         TD01,TG01,TC01,W01,CAPAC01,                    &
         TDM1,TGM1,TCM1,WM1,CAPACM1,                    &
         ra,crb,zfw1,filtr,ea,ta,                       &
!out.......................................................
         TASAVE,TDSAVE,RASAVE,                          &
         TD,TGS,TC,WWW,CAPAC,                           &
         ETMASS,EVAPSOIL,                               &
         EVAPWC,EVAPDC,EVAPSN,EVAPGX,                   &
         GHTFLX,xhlflx,xhsflx,                          &
         USTAR,DRAG,DRAGU,DRAGV,                        &
         TGEFF,BEDO,SALB,RADT,ROFF,SOILM,               &
         anet,gb,gs,                                    &
         CRd,rstfac,rsoil,shf,chf,                      &
         fsdown,fldown,fsup,flup,                       &
         umom,vmom,                                     &
         hlflx,hsflx,ulwsf1,evap,q2m,                   &
         plantr,cm,ch,fm,fh,crib,ZLT1,                  &
         xwcl,xwel,xwsl,xwcs,xwes,xwss,                 &
!yliu09May2017 PAR         ft1,ft2,hs)
         pilphr,                                        &   !Huilin for hr output  Feb 7th 2018
         ft1,ft2,hs,FPARBD,KPARBD)
!
!   THERE ARE A MAIN SSIB PROGRAM AND 12 SUBROUTINE FILES, WHICH ARE:        
!              CROPS
!              INTERC
!              NEWTON
!              RADAB
!              RASIT5
!              ROOT1
!              SFILT
!              STOMA1
!              STRES1
!              TEMRS1
!              UPDAT1
!              VEGOUT
!                                      YONGKANG XUE                     
!-----------------------------------------------------------------------
!                              INPUT  
!     DDTT:      TIME INTERVAL
!     SUNANG:   SOLAR ZENITH ANGLE
!     SWDOWN:   SHORT WAVE DOWN(W/M*M);
!     RADFRAC:  SHORT WAVE COMPONENTS (visible and near IR; direct and diffuse)
!     RLWDOWN:   LONG WAVE DOWN(W/M*M); 
!     PPL, PPC: LARGE SCALE AND CONVECTIVE PRECIPITATIONS AT THE TIME STEP (mm)
!     TM:       TEMPERETURE AT LOWEST MODEL LAYER (K)
!     UMM,VMM:  ZONAL AND MERIDIONAL WINDS AT LOWEST MODEL LAYER (m/S)
!     QM:       WATER VAPOR AT LOWEST MODEL LAYER;
!     PSURF:    SURFACE PRESSURE (pa)
!     ZWIND:    HEIGHT (m) OF LOWEST MODEL LAYER
!     ITYPE:    VEGETATION TYPE
!     ZLAT:     LATITUDE, SOUTH POLE IS -90 DEGREE AND NORTH POLE IS 90 DEGREE
!     MONTH:    MONTH
!     DAY:      CALENDER DATE
!     IYEAR:    YEAR
!                             OUTPUT
!     ETMASS:   EVAPORATION (mm/step)
!     ELATEN:   LATENT HEAT FLUX (w/m*m) 
!     EVAPSOIL,EVAPWC,EVAPDC,EVAPSN: LATENT HEAT FROM (SOIL, INTERCEPTION,
!               TRANSPIRATION, AND SNOW SURFACE)                  
!     HFLUX:    SENSIBLE HEAT FLUX(w/m*m)
!     GHTFLX:   GROUND HEAT FLUX(w/m*m) = CHF+SHF
!     USTAR:    FRICTION VELOCITY (m/s)
!     DRAG:     MOMENTUM FLUX (kg/m/s**2)
!     DRAGU:    U COMPONENT OF MOMENTUM FLUX (kg/m/s**2)
!     DRAGV:    V COMPONENT OF MOMENTUM FLUX (kg/m/s**2)
!     TGEFF:    RADIATIVE TEMPERATURE (K)
!     BEDO:     TOTAL ALBEDO
!     SALB:     ALBEDO FOR 4 COMPONENTS
!     RADT:     NET RADIATION AT CANOPY AND GROUND LEVELS
!     TGS:      SOIL SURFACE TEMPERATURE (K)
!     TC:       CANOPY TEMPERATURE (K)
!     TD:       DEEP SOIL TEMPERATURE (K)
!     TA:       TEMPERATURE AT CANOPY AIR SPACE (K)
!     CAPAC:    INTERCEPTION AT CANOPY (1) and SNOW DEPTH (2)
!     WWW:      SOIL MOISTURE
!     SOILM:    TOTAL SOIL WATER CONTENT
!     ROFF:     RUN OFF
!     
!----------------------------------------------------------------------
!=======================================================================
!  
!      include 'ssibcom.h'
!zzq      INTEGER NNPFT,NPFT 
!zzq      PARAMETER(NNPFT=7,NPFT=6)
!zzq      common /Anparam/fW 
!zzq      common /RSLT/wc,we,ws,al,gs,wc1,we1,ws1,all
!      common /Tfac/fT1,fT2,Respd
!zzq       common /Tfac/Respd
      REAL WWW(3), CAPAC(2), SATCAP(2)
      REAL W01(3),WM1(3),CAPAC01(2),CAPACM1(2)
      REAL TC01,TG01,TD01,TCM1,TGM1,TDM1
      REAL TRAN(2,3,2), REF(2,3,2),SOREF(3)
      REAL GREEN(2), VCOVER(2), ZLT(2), CHIL(2)
      REAL RSTPAR(2,3), TOPT(2), TL(2), TU(2), DEFAC(2)
      REAL PH1(2), PH2(2), RST(2), RSTFAC(2,4)
      REAL ROOTD(2), ZDEPTH(3), ROOTP(3),PHSOIL(3)
      REAL RADN(3,2), RADT(2), PAR(2), PD(2)
      REAL RADFRAC(2,2)
      REAL SALB(2,2), ALBEDO(2,3,2), EXTK(2,3,2)
      REAL YMATT(3),YMATQ(3)
!yliu09May2017 PAR
      REAL FPARBD, KPARBD
!
      real ra, crb,zfw1,TA,rhoair,bps
      real psurf
      real ZLAT
      integer MONTH      ! Huilin add month for ect Oct 2020
      LOGICAL L_VEG
      real pilphr                                             !Huilin for hr output  Feb 7th 2018
! 
!     The final albedo=original albedo+XADJ
      XADJ=0.
!       XADJ=-0.1
!     CTLPA controls stomatal resistance;  
!     Final stomatal resistance=ctlpa * stomatal resistance
      CTLPA=1.
!     NROOT controls root distribution. nroot=1: root uniformly distributes
!           in the soil layer; 
!     If NROOT not =1, root distribution is controled by rootp. 
      NROOT=1
!     if the numbers of iterations are larger than 3, the iteration stops
!zzq dec2013      ITRUNK=3
!     INTG=2
      INTG=1
      DTT =DDTT*FLOAT(INTG)
      TPREC=PPL
!
      UM=SQRT(UMM1**2+VMM1**2)
!zzq080321
      UM=AMAX1(UM,0.1)
!
!zzq... dec 25, 2007
!
!zzq_20080301
      SWDOWN = RADN(1,1)+RADN(1,2)+RADN(2,1)+RADN(2,2) 
!zzq.......
!zzq
      SWDOWN = AMAX1(SWDOWN,0.1)
      RADN(3,2)=rlwdown
      RADN(3,1)=0.
!khs
!khs....for CO2 routines
      RST(2)  =1.e4
!
      TD         =TDM1
      TGS        =TGM1 
      TC         =TCM1 
      CAPAC(1)   =CAPACM1(1)
      CAPAC(2)   =CAPACM1(2)
      WWW(1)     =WM1(1)    
      WWW(2)     =WM1(2)    
      WWW(3)     =WM1(3)    
!
!zzq      GASP = 287.05
!zzq      CPSY=CPAIR/(HLAT*0.622)
!zzq dec2013      AKAPPA = GASR/CPAIR
!zzq dec2013      SIGKI  =1.0E0 / EXP(AKAPPA*LOG(PMR/PSURF))
!zzq dec2013      BPS    =SIGKI
      BPS    =1.0E0/ SIG**AKAPPA
!
!cc   RHOAIR = 1.225
!                                                                       
!song 9.2.2011 CALL RADAB(TRAN,REF,GREEN,VCOVER,CHIL,ZLT,Z2,Z1,SOREF,
      CALL RADAB(TRAN,REF,GREEN,VCOVER,CHIL,ZLT,Z2,Z1,SOREF,  &
           WWW(1),POROS,                                      &
           TC,TGS,SATCAP,EXTK,CLOSS,GLOSS,THERMK,P1F,P2F,     &
           RADT,PAR,PD,SALB,ALBEDO,TGEFF,SUNANG,XADJ,CAPAC,   &
           RADN,BEDO,ZLWUP,RADFRAC,SWDOWN,SCOV2,              &
!yliu09May2017 PAR           fsdown,fldown,fsup,flup)
           fsdown,fldown,fsup,flup,FPARBD,KPARBD)

      CALL ROOT1(PHSAT,BEE,WWW,PHSOIL)
!          CRIB  = 1.0/(SQRT(U2)/RBC+ZLT(1)*.004)
!
!zzqa      if (L_VEG) then
!zzqa        if (crib.eq.0) then
!zzqa           CALL STOMA1(GREEN,VCOVER,CHIL,ZLT,PAR,PD,
!zzqa     &       EXTK,SUNANG,RST,RSTPAR, CTLPA,zfw1)
!zzqa            goto 87
!zzqa         end if

       CALL GMUPI(SUNANG,FC3,CHIL,ZLT,VCOVER,EXTK,PD,GREEN,  &
           GMU1,GMU2,BPIDR3,BPIDR4,BPIDF3,BPIDF4)

!khs..for CO2 routines
!Huilin    EM100 = EM*100.
           EM100 = EA*100.

      CALL ANGS2(PAR(1),TA,TC,EM100,CRB,Co2Ca,ZLT,VCOVER,PSURF,  &
           VM0,FC3,zFW1,PD(1),SUNANG,GMU1,GMU2,BPIDR3,BPIDR4,    &
           BPIDF3,BPIDF4,EXTK(1,1,1),TOPT(1),TU(1),TL(1),ANET,   &
           RST(1),ES,GB,GS,XWCL,XWEL,XWSL,XWCS,XWES,XWSS,        &
           CRd,ft1,ft2,hs)
!
!      Computes CO2 flux [umol m-2 s-1] from soil:
!zzq         Fcs = 6.46e-7*exp((Td-273.2-24.3)/10.*log(2.3)) ! Malhi et al. 1999
!      Computes CO2 flux from atmosphere to land surface:
!zzq        Fca = Fcs - Anet       ! [umol m-2 s-1]
!     Computes CO2 concentration [Pa] within canopy air space:
!zzq         Cacn=Co2Ca-(Anet-Fcs)*PSURF*1.4*ra/(44.6*273.2/Ta*PSURF/1.013e5)
!khs..additional output
!zzq      cs  = cacn - 1.4*psurf*anet/gb
!zzq      ci  = cs   - 1.6*psur*anet/gs
!zzqa      else
!zzqa           CALL STOMA1(GREEN,VCOVER,CHIL,ZLT,PAR,PD,
!zzqa     &       EXTK,SUNANG,RST,RSTPAR, CTLPA,zfw1)
!zzqa      endif
   87   continue
!zzq....For calculation of Fcs outside this subroutine
      TASAVE = TA
      TDSAVE = TD     
      RASAVE = RA
!zzq
!yliu09May2017      RST(1) = AMIN1(RST(1),1000.)
      RST(1) = AMIN1(RST(1),10000.)
!                                                                       
! ** WATER BALANCE CHECK                                                
      TOTWB = WWW(1) * POROS * ZDEPTH(1)      &                         
            + WWW(2) * POROS * ZDEPTH(2)      &                         
            + WWW(3) * POROS * ZDEPTH(3)      &                         
            + CAPAC(1) + CAPAC(2)                                       
!     
      CALL INTERC(DTT ,VCOVER,ZLT,TM,TC,TGS,CAPAC,WWW,PPC,PPL,ROFF,     & 
            ZDEPTH,POROS,CCX,CG,SATCO,SATCAP,SPWET,EXTK,RNOFFS,FILTR,   &
            SMELT)
!                                                                       
!zzq, TGS,TD,TC,TA,RST and RADT need to be intialized.

      CALL TEMRS1(DTT,TC,TGS,TD,TA,TM,QM,EM,PSURF,ZLAT,MONTH,WWW,CAPAC,SATCAP,DTC, &
            DTG,RA,RST,ZDEPTH,BEE,PHSAT,POROS,XDD,Z0,RDC,RBC,VCOVER,    &
            Z2,ZLT,DEFAC,TU,TL,TOPT,RSTFAC,NROOT,ROOTD,PHSOIL,ROOTP,    &
            PH1,PH2,ECT,ECI,EGT,EGI,EGS,HC,HG,EC,EG,EA,RADT,CHF,SHF,    &
            ALBEDO,ZLWUP,THERMK,RHOAIR,ZWIND,UM,USTAR,DRAG,CCX,CG,      &
!song 11.22.2011
            rd,                                                         &
            pilphr,                                                     & !Huilin for hr output  Feb 7th 2018
            BPS,XXT,XXQ,XAKMS,rcc,rib,CU,CT,flup,zfw1,crb,              &
            rst1,rst2,rst3,rst4,rsoil)
!                                                                       
      CALL UPDAT1(DTT ,TC,TGS,TD,CAPAC,DTC,DTG,ECT,ECI,EGT,EGI,         &
         EGS,EG,HC,HG,HFLUX,ETMASS,FILTR,SOILDIF,SOILDRA,ROFF,          &
         RNOFFB,RNOFFS,NROOT,ROOTD,ROOTP,POROS,BEE,SATCO,SLOPE,         &
         PHSAT,ZDEPTH,WWW,CCX,CG,CHF,SHF,SMELT)
!
      CALL SFILT(TD ,TGS ,TC ,WWW ,CAPAC ,                              &
                 TD01,TG01,TC01,W01,CAPAC01,                            &
                 TDM1,TGM1,TCM1,WM1,CAPACM1,INTG)
!                                                                       
      ENDWB = WWW(1) * POROS * ZDEPTH(1)                                &
            + WWW(2) * POROS * ZDEPTH(2)                                &
            + WWW(3) * POROS * ZDEPTH(3)                                &
            + CAPAC(1) + CAPAC(2) - TPREC/1000. + ETMASS/1000. + ROFF 
      ERROR = TOTWB - ENDWB                                             
      IF(ABS(ERROR) .GT. 0.0001) THEN
        WRITE(7,*) 'WARNING WATER BALANCE AFTER UPDATE ',error
      ENDIF
!                                                                       
      CBAL = RADT(1) - CHF - (ECT+HC+ECI)/DTT                           
      GBAL = RADT(2) - SHF - (EGT+EGI+HG+EGS)/DTT                       
      ZLHS = RADT(1) + RADT(2) - CHF - SHF
!zzq      ZRHS = HFLUX + (ECT + ECI + EGT + EGI + EGS)/DTT                  
      ZRHS = HFLUX/DTT + (ECT + ECI + EGT + EGI + EGS)/DTT             
!                                                                       
      IF(ABS (ZLHS - ZRHS) .GT. 1.) THEN                                
          ectw=ect/dtt
          eciw=eci/dtt
          egtw=egt/dtt
          egiw=egi/dtt
          egsw=egs/dtt                                              
      write(7,*)'warning energe balance' 
      write(7,*) 'swdn=',swdown,' radn(3,2)=',radn(3,2),' zlwup=',zlwup
      write(7,*) 'hflux=',hflux,' chf=',chf,' shf=',shf,' ectw=',ectw
      write(7,*) 'eciw=',eciw, ' egtw=',egtw, ' egiw=',egiw
      write(7,*) 'egsw=',egsw, ' zlhs=',zlhs, ' zrhs=',zrhs
      ENDIF
!
!	*** THE FOLLOWING LINES PROVIDE OUTPUT INFORMATION
!
        SOILM=(WWW(1) * POROS * ZDEPTH(1))   &
            + (WWW(2) * POROS * ZDEPTH(2))   &
            +  WWW(3) * POROS * ZDEPTH(3) 
!       UPDATING TM and QM
!zzq 080831       XXT=YMATT(1)+XXT/YMATT(3)
!zzq 080831       XXQ=YMATQ(1)+XXQ/YMATT(3)
!
!zzq 080831        FAC   =YMATT(3)*RHOAIR * DTT
!zzq        DTMDT =(YMATT(2)+HFLUX/FAC/(CPAIR*BPS))/XXT
!zzq 080831        DTMDT =(YMATT(2)+HFLUX/DTT/FAC/(CPAIR*BPS))/XXT
!zzq 080831        DQMDT =(YMATQ(2)+ETMASS/FAC           )/XXQ
!zzq 080831        DTM   =DTMDT   *DTT
!zzq 080831        DQM   =DQMDT   *DTT
!        TM1=TM
!        TM=TM+DTM
!        TM=TM1
!        QM   =QM+DQM
!rr
!zzq... taken from KHS
!
      UMOM=RHOAIR*CU*USTAR*UMM1
      VMOM=RHOAIR*CU*USTAR*VMM1
      HLFLX= ETMASS/RHOAIR/DTT
      HSFLX= HFLUX/CPAIR/RHOAIR/DTT
      ULWSF1=TGEFF*TGEFF*TGEFF*TGEFF*STEFAN
      Q2M=0.622*EA/(0.01*PSURF-EA)
      EVAP=ETMASS*HLAT
      PLANTR=RCC
      SHELEG=CAPAC(2)*1000.
      CM=(USTAR*USTAR)/(UM*UM)
      CH=1/(UM*RA)
      CRIB=RIB
      ZLT1=ZLT(1)
!
!fds Calculates turbulence terms (04/2012)
      FM=VKC/CU
      FH=VKC/CT
!zzq
!
        EVAPSOIL=EGS /DTT 
        EVAPWC=ECI /DTT 
        EVAPDC=ECT /DTT 
        EVAPSN=EGI /DTT
        EVAPGX=EGT /DTT
        ELATEN=EVAPSOIL+EVAPWC+EVAPDC+EVAPSN+EVAPGX
        GHTFLX=CHF+SHF
!KHS=================================================================
      xhsflx=(hc+hg)/dtt
      xhlflx=(ect+eci+egt+egi+egs)/dtt
!KHS=================================================================
        xrst1=rst(1)
        xrst2=rst(2)
        DRAGU=DRAG*UMM1/UM
        DRAGV=DRAG*VMM1/UM
  100   FORMAT (i8,1x,14(f7.2,1x))
!
       END SUBROUTINE
END MODULE 
