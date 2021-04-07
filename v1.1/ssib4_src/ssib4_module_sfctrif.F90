MODULE  ssib4_module_sfctrif       
      use ssib4_module_comsconst
      use ssib4_module_comssib4sfc
      use ssib4_module_comsveg
      use ssib4_module_trifparms
      use ssib4_module_weighting
      use ssib4_module_updatepara
      use ssib4_module_ssibsub
      use ssib4_module_sitriffid
      use ssib4_module_ssibsfcbio
      use ssib4_module_fireseason
      use ssib4_module_fireimpact
!
      IMPLICIT NONE

contains
      subroutine sfctrif(                                             &
! for ssib4
         dataios,                                                     &
         ISTRT,MONTH,JDAY,TIMESTEP,                                   &
         SUNANG,radn11,radn12,radn21,radn22,RLWDOWN,                  &
         PPL,PPC,TM,UMM1,VMM1,QM,EM,RHOAIR,KREC,                      &
         PSURF,SIG,ITYPE,co2ca,ZLAT,ZWIND,                            &
         LIGHT,PREC60,PREC10,PRECD,DP,AGR,GDP,PEATF,                  &!Huilin add fire variables 
!out.......................................................
         TD1,TG1,TC1,TA,                                              &
         WWW1,WWW2,WWW3,CAPAC1,CAPAC2,                                &
         SALB11,SALB12,SALB21,SALB22,                                 &
         radt1,radt2,                                                 &
         rstfac1,rstfac2,rstfac3,rstfac4,                             &
         ra1,rb1,rd1,rsoil,shf,chf,                                   &
         ETMASS,EVAPSOIL,                                             &
         EVAPWC,EVAPDC,EVAPSN,EVAPGX,                                 &
         GHTFLX,xhlflx,xhsflx,                                        &
         USTAR,DRAG,DRAGU,DRAGV,                                      &
         TGEFF,BEDO,roff,soilm,                                       &
         fsdown,fldown,fsup,flup,                                     &
         umom,vmom,                                                   &
         hlflx,hsflx,ulwsf1,evap,q2m,rh,pilphr,                       & ! Huilin Feb 7 2018
         plantr,canopy,sheleg,cm,ch,fm,fh,rib,ZLT1,                   &
         aveanet,avecrdn,fcs,fca,cacn,cs,ci,gb,gs,                    &
         xwcl,xwel,xwsl,xwcs,xwes,xwss,                               &
         aveft1,aveft2,avehs,                                         &
!TRIF output to GCM
         tvcover,tzlt,green_frac,                                     &
         GPP,NEP,NPP,RESP_P,RESP_S,                                   &
         LEAF1,ROOT1,WOOD1,DCVEG1,LIT_C,LIT_C_T,BAF_AVG60,            &
         AGB_AVG,FEO_AVG,FES_AVG,FE1_AVG,FE2_AVG,FD_AVG,FCLI,FSAT,    & 
         LH,IGN,FM_AVG,FT_AVG,FB_AVG,FM_FIRE,                         &
         BAF_AVG,BAF_PEAT,NFIRE_AVG,BURN_AVG,BURN_PEAT,BURN_FIRE,     & ! Huilin add fire output
         EM_CO2,EM_CO,EM_CH4,EM_NHMC,EM_H2,EM_NOX,EM_N2O,             &
         EM_PM25,EM_TPM,EM_TC,EM_OC,EM_BC,                            &
!yliu09May2017
         RESP_P_M,RESP_P_G,ANETC,RDC_TRIF)
!
!----------------------------------------------------------------------
      type(ssibdataios)::dataios
      real,parameter:: SEC_DAY = 86400.
      INTEGER:: MONTH, JDAY,ITYPE, NPFT1,IL,ISTRT,DAY,KREC
      real:: psurf,sig, zlat,zwind,light,prec60,prec10,precd,dp,agr,gdp,peatf            !Huilin add fire variables 
      real::                                                         &
        SUNANG, RLWDOWN,PPL,PPC,TM,UMM1,VMM1,QM, EM,                 &
        rhoair,radn11,radn12,radn21,radn22,radnsum
      real:: BAF_AVG60, RBB
      real,dimension(480)  :: BAF60

      real,dimension(2)   :: CAPAC,RADT
      real,dimension(3)   :: WWW
      real,dimension(2,2) :: SALB,RADFRAC
      real,dimension(2,4) :: RSTFAC
      real,dimension(3,2) :: RADN
!
      real :: xhlflx,xhsflx
      real :: fcs,fca,cs,ci,gb,gs,co2ca
      real :: fsdown,fldown,fsup,flup
      real :: umom,vmom,hlflx,hsflx,ulwsf1,evap,q2m
      real :: plantr,canopy,sheleg,cm,ch,fm,fh,rib,ZLT1
      real :: xwcl,xwel,xwsl,xwcs,xwes,xwss
      real :: capac1,capac2
      real :: www1,www2,www3,rh,pilphr                                     !Huilin Huang output pilphr Feb 7 2018
      real :: rstfac1,rstfac2,rstfac3,rstfac4
      real :: ra1,rb1,rd1,rsoil,shf,chf
      real :: TASAVE,TDSAVE,RASAVE
      real :: TG1,TC1,TD1,TGEFF
      real :: etmass,EVAPSOIL,EVAPWC
      real :: EVAPDC,EVAPSN,EVAPGX,GHTFLX,USTAR
      real :: DRAG,DRAGU,DRAGV,RADT1,RADT2,bedo,roff,soilm
      real :: salb11,salb12,salb21,salb22
      real :: aveft1,aveft2,avehs
      real :: aveanet,anetdgb,anetdgs,avecrdn,cacn,vfrac
!
      real,dimension(NNPFT) ::  xhlflxn,xhsflxn,gbn,gsn
      real,dimension(NNPFT) ::                              &
         fsdownn,fldownn,fsupn,flupn,                       &
         umomn,vmomn,hlflxn,hsflxn,ulwsf1n,evapn,q2mn,      &
         plantrn,cmn,chn,fmn,fhn,ribn,ZLT1n,                &
         xwcln,xweln,xwsln,xwcsn,xwesn,xwssn
      real,dimension(NNPFT):: WWW1N,WWW2N,WWW3N
      real,dimension(NNPFT):: CAPAC1N,CAPAC2N,TGEFFN
      real,dimension(NNPFT):: TGSN,TCN,TDN,TASAVN,TDSAVN,RASAVN
      real,dimension(NNPFT):: RSTFAC11N,RSTFAC12N,RSTFAC13N,RSTFAC14N
      real,dimension(NNPFT):: rsoiln,shfn,chfn
      real,dimension(NNPFT):: ETMASSN,EVAPSOILN,EVAPWCN,EVAPGXN,EVAPDCN,EVAPSNN
      real,dimension(NNPFT):: GHTFLXN,ANETN,USTARN,DRAGN,DRAGUN
      real,dimension(NNPFT):: DRAGVN,RADT1N,RADT2N,ROFFN,SOILMN,BEDON
      real,dimension(NNPFT):: salb11n,salb12N,salb21N,salb22N,PILPHRn   !Huilin output HR Feb 7 2018
      real,dimension(NNPFT):: ft1n,ft2n,hsn,TAN
      real,dimension(NNPFT):: ran,rbn,crdn

!
!-----------------------------------------------------------------------
! Inputs from SSiB
!-----------------------------------------------------------------------
      REAL,DIMENSION(NPFT)::    &
           ANETC,               & !IN Net canopy photosynthesis
           RDC_TRIF,            & !IN Canopy dark respiration(mol CO2/m2/s).
           FSMC,                & !IN Moisture availability factor.
           TSTAR                  !IN Surface temperature on veg tiles (K).
      REAL,DIMENSION(NPFT)::    &
           FPARn,KPARn
!-----------------------------------------------------------------------
! Prognostic variables
!-----------------------------------------------------------------------
      REAL,DIMENSION(NNPFT) ::  &
           HT_TRIF1,            & !WORK Canopy height (m)
           LAI1,                & !WORK Leaf area index.
           GREEN_TRIF1,         & !WORK GREEN FRACTION
           SAI1,                & !WORK Stem area index.
!zzq2012sept24
           GREEN_SSIB1,         & !WORK GREEN FRACTION
           LAID1                  !WORK dead leaf area index.
!
      real,dimension(NTYPES)::  &
           FRAC_TRIF              !WORK Vegetation coverage.

      REAL::                    & 
       TIMESTEP                   !IN Model TIMESTEP (s).
!
      real,dimension(NPFT)::    &
!zzq2008May12       STHU          !IN Unfrozen soil moisture as a fraction
       STHU                       !IN Unfrozen soil moisture as a fraction(mol CO2/m2/s).
      REAL::                    &
       T_TRIF                     !IN Air temperature (K).
      REAL,DIMENSION(NPFT) ::   &
!zzq2008May12                   & TS1, !IN Sub-surface temperature (K).
         TS1,                   & !IN Sub-surface temperature (K).
         GPP,                   & !OUT Gross Primary Productivity (kg C/m2/s).
         NPP,                   & !OUT Net Primary Productivity (kg C/m2/s).
         RESP_P,                & !OUT Plant respiration rate (kg C/m2/s).
         LEAF1,                 & !WORK Leaf area index.
         ROOT1,                 & !WORK Stem area index.
         WOOD1,                 & !WORK Stem area index.
         DCVEG1,                & !WORK Leaf area index.
         lit_c,                 & !WORK Net canopy photosynthesis (mol CO2/m2/s).
         PHEN,                  & !WORK Leaf area index.
         V_WILT,                & !WORK Volumetric soil moisture concentration
                                  ! below which stomata close (m3 H2O/m3 soil).
         v_sat                    ! volumetric soil moisture concentration
                                  !       at saturation (m3 h2o/m3 soil).
!yliu09May2017
      REAL,DIMENSION(NPFT)::    &
        RESP_P_M,RESP_P_G
      real::                    &
        NEP,                    & !OUT Net Ecosystem Productivity (kg C/m2/s).
        RESP_S,                 & !OUT Soil respiration rate (kg C/m2/s).
        lit_c_t                  !WORK Unfrozen soil moisture as a fraction of saturation.

      INTEGER::                 & 
         ITRIFSTEP,KTRIFDAY,N,NXTYPE  !Loop counters

      LOGICAL::                 &
         L_TRIF,                & !WORK .T. if vegetation to be updated on current TIMESTEP.
         L_PHEN,                & !WORK .T. if phenology to be updated  on current TIMESTEP.
         L_VEG,                 & !WORK .T. if vegetation coverage.
         L_FIRE                   !WORK .T. if fire is included Huilin Nov 2018 

      real:: TLAI,TSAI
      real:: FB_AVG,AGB_AVG,FES_AVG,FEO_AVG,FE1_AVG,FE2_AVG,FD_AVG
      real:: LH,IGN,FM_AVG,FT_AVG,FM_FIRE !Huilin May 11 for debug
      real:: FCLI,FSAT,BAF_PEAT                        !Huilin add for peatfire
      real:: BAF_AVG,NFIRE_AVG,BURN_AVG,BURN_PEAT      !Huilin add for fire
      real:: EM_CO2,EM_CO,EM_CH4,EM_NHMC,EM_H2,EM_NOX,EM_N2O   !Huilin add emission for gas and aerosols Nov 2018
      real:: EM_PM25,EM_TPM,EM_TC,EM_OC,EM_BC
      real:: RESP_S_DR,VEG_FRAC,FRAC_VS
      real:: SAI_MEAN,VEGCOVER
      real,dimension(2):: tvcover,zlt,tzlt,                     &
          green_frac,green,vcover,zvcover
      real,dimension(NPFT)::                                    &
           FSMCx,FSMC_ave,GPP_DR,NPP_DR,HT_TRIF,                &
           LAI,SAI,GREEN_TRIF,RESP_W_DR,G_LEAF_DAY,G_LEAF_DR

      REAL,DIMENSION(NPFT)::                                    &  
           BURN_FIRE,LEAF_LOSS_DR,WOOD_LOSS_DR,ROOT_LOSS_DR

!zzq2012sept24
      real,dimension(NPFT) :: LAID,GREEN_SSIB
      real,dimension(NNPFT):: veg_fracn
!
      real,dimension(2,3,2):: TRAN,REF
      real,dimension(2,3)  :: RSTPAR
      real,dimension(2)    :: TOPT,CHIL,TL,TU,DEFAC,PH1,PH2,ROOTD
      real,dimension(2)    :: CAPAC01,CAPACM1
      real,dimension(3)    :: SOREF,ZDEPTH,W01,WM1
      real poros,vm0,fc3,z0,xdd,z2,z1,     &
          rdc,rbc,bee,phsat,satco,slope
      real,dimension(NNPFT):: FPARBDn,KPARBDn
      real::                               &
          TD01,TG01,TC01,TDM1,TGM1,TCM1,   &
          RA,RB,zfw1,filtr,ea,TA
      real:: xneglect,FTIME_PHEN
!
      xneglect =1.0e-4
!-------------------------------------
!rr
! **  Read in vegetation parameters
!
      RADN(3,2) =RLWDOWN
      RADN(3,1) =0.
      SUNANG=AMAX1(SUNANG,0.01746)
!zzq Aug29,2008
      RADNSUM = RADN11 + RADN12 + RADN21 + RADN22
      if(sunang.ge.0.01746.and.RADNSUM.ne.0) then
         RADN(1,1)=RADN11
         RADN(1,2)=RADN12
         RADN(2,1)=RADN21
         RADN(2,2)=RADN22
         RADFRAC(1,1) = RADN(1,1) /RADNSUM
         RADFRAC(1,2) = RADN(1,2) /RADNSUM
         RADFRAC(2,1) = RADN(2,1) /RADNSUM
         RADFRAC(2,2) = RADN(2,2) /RADNSUM
      else
         RADN(1,1) = 0.
         RADN(1,2) = 0.
         RADN(2,1) = 0.
         RADN(2,2) = 0.
         RADFRAC(1,1) = 0.25
         RADFRAC(1,2) = 0.25
         RADFRAC(2,1) = 0.25
         RADFRAC(2,2) = 0.25
      endif

!--------------------------------------------------------

!zzq 2009Feb23  CALL VEGOUT1(POROS,ZDEPTH,ITYPE)
!
       KTRIFDAY = dataios % xKTRIFDAY
       ITRIFSTEP= dataios % xITRIFSTEP
       CS       = dataios % xCo2CS
       RESP_S_DR= dataios % xRESP_S_DR
       VEG_FRAC = dataios % xVEG_FRAC
       FRAC_VS  = dataios % xFRAC_VS

       DO N=1,NPFT
         GPP_DR(N)  = dataios % xGPP_DR(N)
         NPP_DR(N)  = dataios % xNPP_DR(N)
         HT_TRIF(N) = dataios % xHT_TRIF(N)
         LAI(N)     = dataios % xLAI_TRIF(N)
         SAI(N)     = dataios % xSAI_TRIF(N)
         GREEN_TRIF(N)= dataios % xGREEN_TRIF(N)
!zzq2012sept24
         LAID(N)    = dataios % xLAID_TRIF(N)
         GREEN_SSIB(N)= dataios % xGREEN_SSIB(N)
!
         RESP_W_DR(N) = dataios % xRESP_W_DR(N)
         G_LEAF_DAY(N)= dataios % xG_LEAF_DAY(N)
         G_LEAF_DR(N) = dataios % xG_LEAF_DR(N)
         PHEN(N)      = dataios % xPHEN(N)
         FSMC(N)      = dataios % xFSMC(N)
         FSMC_ave(N)  = dataios % xFSMC_ave(N)
!Huilin add fore fire output
         LEAF1(N)     = dataios % xLEAF1(N)
         WOOD1(N)     = dataios % xWOOD1(N)
         ROOT1(N)     = dataios % xROOT1(N)
         LEAF_LOSS_DR(N) = dataios % xLEAF_LOSS_DR(N)
         WOOD_LOSS_DR(N) = dataios % xWOOD_LOSS_DR(N)
         ROOT_LOSS_DR(N) = dataios % xROOT_LOSS_DR(N)
       ENDDO
!
       DO N=1,NNPFT
          FRAC_TRIF(N) =  dataios % xVEG_FRACN(N)
       ENDDO
!
       DO N=1,NPFT
         HT_TRIF1(N)=HT_TRIF(N)
         LAI1(N)=LAI(N)
         SAI1(N)=SAI(N)
         GREEN_TRIF1(N)=GREEN_TRIF(N)
!zzq2012sept24
         LAID1(N)=LAID(N)
         GREEN_SSIB1(N)=GREEN_SSIB(N)
       ENDDO

       HT_TRIF1(NNPFT)=0.1
       LAI1(NNPFT)=0.0001
       SAI1(NNPFT)=0.1
       GREEN_TRIF1(NNPFT)=0.0001
!
!zzq2012sept24
!
       LAID1(NNPFT)=0.1
       GREEN_SSIB1(NNPFT)=0.0001

!==========================================================================
!zzq2012sep26 update soil parameters
!zzq2012sep26
      CALL VEGOUT1(POROS,ZDEPTH,ITYPE)
!
      DO npft1=1,NNPFT
!no feedback
!zzq 2009-5-15          if(ITYPE.eq.Npft1) THEN
!zzq 2009-5-15             veg_fracn(Npft1) = 1.0
!zzq 2009-5-15          ELSE
!zzq 2009-5-15             veg_fracn(Npft1) = 0
!zzq 2009-5-15          ENDIF
!feedback
       veg_fracn(NPFT1) = dataios %  xVEG_FRACN(NPFT1)
      enddo

      BAF_AVG60         = SUM(BAF60)

!
      call update_soil(NNPFT,veg_fracn,month,itype,istrt)
!
!==========================================================================
!
      DO 1001 NPFT1=1,NNPFT
!
!zzq    if(veg_fracn(npft1).lt.xneglect) goto 1001

       L_VEG = .FALSE.
       IF(NPFT1.LE.NPFT) L_VEG=.TRUE.

       if(NPFT1.eq.1) NXTYPE=1
       if(NPFT1.eq.2) NXTYPE=2
       if(NPFT1.eq.3) NXTYPE=3
       if(NPFT1.eq.4) NXTYPE=4
       if(NPFT1.eq.5) NXTYPE=5
       if(NPFT1.eq.6) NXTYPE=6
       if(NPFT1.eq.7) NXTYPE=7
!yliu09May2017 addtype
       if(NPFT1.eq.8) NXTYPE=8

!zzq 2009Feb23
!zzq 2012Sep26     CALL VEGOUT1(POROS,ZDEPTH,ITYPE)
      call VEGOUT2(TRAN,REF,GREEN,VCOVER,CHIL,        &
           RSTPAR,VM0,TOPT,TL,TU,FC3,DEFAC,PH1,PH2,   &
           ZLT,Z0,XDD,Z2,Z1,RDC,RBC,ROOTD,            &
           SOREF,BEE, PHSAT,SATCO,SLOPE,              &
           MONTH,NXTYPE)
!
       ZVCOVER(1)=VCOVER(1)   !1.0
       ZVCOVER(2)=VCOVER(2)   !1.0
!
!zzq feedback from TRIF
!zzq 2012SEPT19 In response to FSales
!zzq 2012SEPT19 if npft1==NPFT, don't update zvcover and so forth. 
!
       if(npft1.LE.NPFT) then
         call updatepara(                              &
!zzq2012sept24           NXTYPE,HT_TRIF1(npft1),GREEN_TRIF1(npft1), &
!zzq2012sept24           LAI1(npft1),SAI1(npft1),                   &
            NXTYPE,HT_TRIF1(npft1),GREEN_SSIB1(npft1), & 
            LAI1(npft1),LAID1(npft1),                  &
            Z0,XDD,RBC,RDC,Z2,Z1,ZLT(1),GREEN(1),ZVCOVER(1))
       endif
!
        TDM1= dataios % TDM(npft1)
        TGM1= dataios % TGM(npft1)
        TCM1= dataios % TCM(npft1)
        TD01= dataios % TD0(npft1)
        TG01= dataios % TG0(npft1)
        TC01= dataios % TC0(npft1)
!
        CAPAC01(1)= dataios % CAP01(npft1)
        CAPAC01(2)= dataios % CAP02(npft1)
        CAPACM1(1)= dataios % CAPM1(npft1)
        CAPACM1(2)= dataios % CAPM2(npft1)
!
        w01(1)= dataios % W01(npft1)
        w01(2)= dataios % W02(npft1)
        w01(3)= dataios % W03(npft1)
        wM1(1)= dataios % WM1(npft1)
        wM1(2)= dataios % WM2(npft1)
        wM1(3)= dataios % WM3(npft1)

        zFW1= dataios % xZFW(npft1)
        RA  = dataios % xCO2RA(npft1)
        RB  = dataios % xCO2RB(npft1)
        EA  = dataios % xea(npft1)
        filtr      = dataios % xfiltr(npft1)
        TAN(npft1) = dataios % xTA(npft1)
!
       call SFCBIO(                                                    &
!in ......................................................            
         TIMESTEP,SUNANG, RADFRAC,RADN,RLWDOWN,                        &
         PPL,PPC,TM,UMM1,VMM1,QM,EM,RHOAIR,                            &
         PSURF,SIG,ZWIND,co2ca,L_VEG,ZLAT,MONTH,                       &
! vegetation parameters
         TRAN,REF,GREEN,ZVCOVER,CHIL,                                  &
         RSTPAR,VM0,TOPT,TL,TU,fc3,DEFAC,PH1,PH2,                      &
         ZLT,Z0,XDD,Z2,Z1,RDC,RBC,ROOTD,SOREF,BEE,                     &
         PHSAT,SATCO,SLOPE,POROS,ZDEPTH,                               &
!inout.....................................................
         TD01,TG01,TC01,W01,CAPAC01,                                   &
         TDM1,TGM1,TCM1,WM1,CAPACM1,                                   &
         ra,rb,zfw1,filtr,ea,tan(npft1),                               &
!out.......................................................
         TASAVN(npft1),  TDSAVN(npft1),  RASAVN(npft1),                &
         Tdn(npft1),     Tgsn(npft1),    Tcn(npft1),                   &
         WWW,    CAPAC,                                                &
         ETMASSN(npft1),                                               &
         EVAPSOILn(npft1), EVAPWCn(npft1), EVAPDCn(npft1),             &
         EVAPSNn(npft1),   EVAPGXn(npft1), GHTFLXn(npft1),             &
         Xhlflxn(NPFT1),   Xhsflxn(NPFT1), USTARn(npft1),              &
         Dragn(npft1),     Dragun(npft1),  Dragvn(npft1),              &
         Tgeffn(npft1),    Bedon(npft1),                               &
         SALB,  RADT,                                                  &
         Roffn(npft1),   soilmn(npft1), Anetn(npft1),                  &
         Gbn(npft1),     Gsn(npft1),    Crdn(npft1),                   &
         RSTFAC,                                                       &
         rsoiln(npft1),  shfn(npft1),    chfn(npft1),                  &
         Fsdownn(NPFT1), Fldownn(NPFT1), Fsupn(NPFT1),                 &
         Flupn(NPFT1),   Umomn(NPFT1),   Vmomn(NPFT1),                 &
         Hlflxn(NPFT1),  Hsflxn(NPFT1),  Ulwsf1n(NPFT1),               &
         Evapn(NPFT1),   Q2mn(NPFT1),    Plantrn(NPFT1),               &
         Cmn(NPFT1),     Chn(NPFT1),                                   &
         Fmn(NPFT1),     Fhn(NPFT1),                                   &
         RIbn(NPFT1),                                                  &
         Zlt1n(NPFT1),   Xwcln(NPFT1),   Xweln(NPFT1),                 &
         Xwsln(NPFT1),                                                 &
         Xwcsn(NPFT1),   Xwesn(NPFT1),   Xwssn(NPFT1),                 &
         PILPHRn(npft1),                                               &   !Huilin for hr output  Feb 7th 2018   
         Ft1n(npft1),    Ft2n(npft1),    hsn(npft1),                   &
         FPARBDn(npft1), KPARBDn(npft1)  )
!
         if (NPFT1.le.NPFT)then
             ANETC(NPFT1)=anetn(npft1)
             RDC_TRIF(NPFT1)=CRdn(npft1)
             TSTAR(NPFT1)=TCn(npft1)
             FSMCx(NPFT1) =RSTFAC(1,2)
         end if
!
         dataios % xZFW(npft1)  = zFW1
         dataios % xCO2RA(npft1)= RA
         dataios % xCO2RB(npft1)= RB
         dataios % xTA(npft1)  = TAN(npft1)
         dataios % xea(npft1)   = EA
         dataios % xfiltr(npft1)= filtr

         dataios % TDM(npft1) = TDM1
         dataios % TGM(npft1) = TGM1
         dataios % TCM(npft1) = TCM1
         dataios % TD0(npft1) = TD01
         dataios % TG0(npft1) = TG01
         dataios % TC0(npft1) = TC01
!
         dataios % CAP01(npft1) = CAPAC01(1)
         dataios % CAP02(npft1) = CAPAC01(2)
         dataios % CAPM1(npft1) = CAPACM1(1)
         dataios % CAPM2(npft1) = CAPACM1(2)
!
         dataios % W01(npft1) = W01(1)
         dataios % W02(npft1) = W01(2)
         dataios % W03(npft1) = W01(3)
         dataios % WM1(npft1) = WM1(1)
         dataios % WM2(npft1) = WM1(2)
         dataios % WM3(npft1) = WM1(3)
!
         ran(npft1)=ra
         rbn(npft1)=rb
         CAPAC1N(NPFT1)=CAPAC(1)
         CAPAC2N(NPFT1)=CAPAC(2)
!zzq20120907         www1n(NPFT1)=www(1)
!zzq20120907         www2n(NPFT1)=www(2)
!zzq20120907         www3n(NPFT1)=www(3)
         www1n(NPFT1)=www(1)*poros
         www2n(NPFT1)=www(2)*poros
         www3n(NPFT1)=www(3)*poros
         SALB11N(NPFT1)=SALB(1,1)
         SALB12N(NPFT1)=SALB(1,2)
         SALB21N(NPFT1)=SALB(2,1)
         SALB22N(NPFT1)=SALB(2,2)
         RSTFAC11N(NPFT1)=RSTFAC(1,1)
         RSTFAC12N(NPFT1)=RSTFAC(1,2)
         RSTFAC13N(NPFT1)=RSTFAC(1,3)
         RSTFAC14N(NPFT1)=RSTFAC(1,4)
         RADT1N(NPFT1)=RADT(1)
         RADT2N(NPFT1)=RADT(2)
 1001    continue
!
         RA1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,ran)
         RB1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,rbn)
         RD1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,crdn)
         TA     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TAN)
         TC1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TCN)
         TG1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TGSN)
         TD1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TDN)
         TASAVE = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TASAVN)
         TDSAVE = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TDSAVN)
         RASAVE = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,RASAVN)
         CAPAC1 = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,CAPAC1N)
         CAPAC2 = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,CAPAC2N)
         WWW1  = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,WWW1N)
         WWW2  = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,WWW2N)
         WWW3  = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,WWW3N)
         ETMASS    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,ETMASSN)
         EVAPSOIL  = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,EVAPSOILN)
         EVAPWC    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,EVAPWCN)
         EVAPDC    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,EVAPDCN)
         EVAPSN    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,EVAPSNN)
         EVAPGX    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,EVAPGXN)
         GHTFLX    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,GHTFLXN)
         xhlflx    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xhlflxN)
         xhsflx    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xhsflxN)
         USTAR     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,USTARN)
         DRAG      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,DRAGN)
         DRAGU     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,DRAGUN)
         DRAGV     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,DRAGVN)
         RADT1     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,RADT1N)
         RADT2     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,RADT2N)
         TGEFF     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,TGEFFN)
         bedo      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,bedon)
         roff      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,ROFFN)
         soilm     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,SOILMN)
         fsdown    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,fsdownN)
         fldown    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,fldownN)
         fsup      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,fsupN)
         flup      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,flupN)
         umom      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,umomN)
         vmom      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,vmomN)
         hlflx     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,hlflxN)
         hsflx     = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,hsflxN)
         ulwsf1    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,ulwsf1N)
         evap      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,evapN)
         q2m       = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,q2mN)
         PILPHR    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,PILPHRn)          !Huilin Feb 7 2018 
!
!  --- ... canopy water and snow water equivalent
          canopy = capac1 * 1000.0
          sheleg = capac2 * 1000.0
!
         cm        = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,cmN)
         ch        = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,chN)
         fm        = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,fmN)
         fh        = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,fhN)
         rib       = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,ribN)
         ZLT1      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,ZLT1N)
         xwcl      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xwclN)
         xwel      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xwelN)
         xwsl      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xwslN)
         xwcs      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xwcsN)
         xwes      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xwesN)
         xwss      = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,xwssN)
         salb11    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,salb11N)
         salb12    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,salb12N)
         salb21    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,salb21N)
         salb22    = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,salb22N)
         gb        = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,gbN)
         gs        = summ(1,NNPFT,NNPFT,veg_fracn,xneglect,gsN)
         aveanet   = summ(1,NPFT,NNPFT,veg_fracn,xneglect,anetN)
         avecrdn   = summ(1,NPFT,NNPFT,veg_fracn,xneglect,crdN)
         aveft1    = summ(1,NPFT,NNPFT,veg_fracn,xneglect,ft1N)
         aveft2    = summ(1,NPFT,NNPFT,veg_fracn,xneglect,ft2N)
         avehs     = summ(1,NPFT,NNPFT,veg_fracn,xneglect,hsN)
         RSTFAC1   = summ(1,NPFT,NNPFT,veg_fracn,xneglect,RSTFAC11N)
         RSTFAC2   = summ(1,NPFT,NNPFT,veg_fracn,xneglect,RSTFAC12N)
         RSTFAC3   = summ(1,NPFT,NNPFT,veg_fracn,xneglect,RSTFAC13N)
         RSTFAC4   = summ(1,NPFT,NNPFT,veg_fracn,xneglect,RSTFAC14N)
         rsoil     = summ(1,NPFT,NNPFT,veg_fracn,xneglect,rsoiln)
         shf       = summ(1,NPFT,NNPFT,veg_fracn,xneglect,shfn)
         chf       = summ(1,NPFT,NNPFT,veg_fracn,xneglect,chfn)
         plantr    = summ(1,NPFT,NNPFT,veg_fracn,xneglect,plantrn)
!

!yliu09May2017 -- bedo
       bedo = 0.
       vfrac = 0.
       do npft1 = 1,NNPFT
          if (bedon(NPFT1).ne.-999.) then
          bedo = bedo+FRAC_TRIF(NPFT1)*bedon(NPFT1)
          vfrac = vfrac+FRAC_TRIF(NPFT1)
          endif
       enddo
       if (vfrac.ne.0.) then
          bedo = bedo/vfrac
       else
          bedo = -999.
       endif
!!!--
       VFRAC = 0.
       TLAI =0.
       TSAI = 0.
       aveanet = 0.0
       avecrdn = 0.0
!
       DO npft1=1,NPFT
!zzq2009Feb13     ANETC(NPFT1)=FRAC_TRIF(npft1)*ANETC(npft1)
!zzq2009Feb13     RDC_TRIF(NPFT1)=FRAC_TRIF(npft1)*RDC_TRIF(npft1)
          TLAI = TLAI + FRAC_TRIF(NPFT1)*(LAI(NPFT1) + SAI(NPFT1))
          TSAI = TSAI + FRAC_TRIF(NPFT1)*SAI(NPFT1)
          aveanet = aveanet + FRAC_TRIF(npft1)*ANETC(npft1)
          avecrdn = avecrdn + FRAC_TRIF(npft1)*RDC_TRIF(Npft1)
          VFRAC = VFRAC +FRAC_TRIF(npft1)
        end do
!
        IF(VFRAC.NE.0) THEN
           TLAI = TLAI/VFRAC
           TSAI = TSAI/VFRAC
        ENDIF
!
!------------------------------------------------------------------------
! Determine whether phenology or vegetation dynamics need to be called
!------------------------------------------------------------------------

        L_TRIF=.FALSE.
        L_PHEN=.FALSE.

        IF (TRIF_ON .AND. (ITRIFSTEP .EQ. STEP_DAY)      &
            .AND. (MOD(KTRIFDAY,DAY_TRIF) .EQ. 0.0)) L_TRIF=.TRUE.

        IF (PHEN_ON .AND. (ITRIFSTEP .EQ. STEP_DAY)      &
            .AND. (MOD(KTRIFDAY,DAY_PHEN) .EQ. 0.0)) L_PHEN=.TRUE.
!
!zzq2008May12          T_TRIF = TA
!zzq2008May12          STHU = WWW2x
!zzq2008May12          TS1 = TD1
!
!zzq2009Feb19
        L_FIRE  = .True.
!----------------------------------------------------------------------- 
!                 Huilin add fire model 
!----------------------------------------------------------------------- 
        IF(L_FIRE) THEN
!----------------------------------------------------------------------- 
!                 Fire occurrence 
!----------------------------------------------------------------------- 

          CALL FIRESEASON(                                             &
!IN ........................................................        
            TIMESTEP,ZLAT,VFRAC,DP,LIGHT,AGR,GDP,PEATF,                &
            ZWIND,TM,EM,TG1,WWW1N,WWW2N,PREC60,                        & 
            FEO_AVG,FES_AVG,FE1_AVG,FE2_AVG,FD_AVG,                    &
            LH,IGN,FM_AVG,FT_AVG,FM_FIRE,                              &
            FRAC_TRIF,LEAF1,PHEN,WOOD1,                                &  
            LEAF_LOSS_DR,WOOD_LOSS_DR,                                 &  
!OUT........................................................        
            RH,FB_AVG,AGB_AVG,FCLI,FSAT,                               &
            BAF_AVG,BAF_PEAT,NFIRE_AVG                                 &
          )

!----------------------------------------------------------------------- 
!                 Fire impact 
!----------------------------------------------------------------------- 

          CALL FIREIMP(                                                &
!IN ....  ..................................................            
            TIMESTEP,ZLAT,DAY_TRIF,L_TRIF,                             &
            FRAC_TRIF,BAF_AVG,BAF_PEAT,                                &
            LEAF1,PHEN,WOOD1,ROOT1,CS,                                 &
!OUT....  ..................................................                     
            BURN_AVG,BURN_PEAT,BURN_FIRE,                              &
            LEAF_LOSS_DR,WOOD_LOSS_DR,ROOT_LOSS_DR,                    &
            EM_CO2,EM_CO,EM_CH4,EM_NHMC,EM_H2,EM_NOX,EM_N2O,           &
            EM_PM25,EM_TPM,EM_TC,EM_OC,EM_BC                           &
          )

        ELSE

          FB_AVG        = 0.0
          AGB_AVG       = 0.0 
          BAF_AVG       = 0.0 
          BAF_PEAT      = 0.0 
          NFIRE_AVG     = 0.0
          BURN_AVG      = 0.0
          BURN_PEAT     = 0.0
          BURN_FIRE     = 0.0
          LEAF_LOSS_DR  = 0.0
          WOOD_LOSS_DR  = 0.0
          ROOT_LOSS_DR  = 0.0  
          EM_CO2        = 0.0
          EM_CO         = 0.0
          EM_CH4        = 0.0
          EM_NHMC       = 0.0
          EM_H2         = 0.0
          EM_NOX        = 0.0
          EM_N2O        = 0.0          
          EM_PM25       = 0.0
          EM_TPM        = 0.0
          EM_TC         = 0.0
          EM_OC         = 0.0
          EM_BC         = 0.0         

        ENDIF
        DO DAY=1,479
           BAF60(DAY)   = dataios % xBAF60(DAY+1)
        END DO
        BAF60(480)      = BAF_AVG*3600*3
!----------------------------------------------------------------------- 
!                 Fire impact end 
!----------------------------------------------------------------------- 

        FTIME_PHEN=TIMESTEP/REAL(SEC_DAY*DAY_PHEN)
        DO npft1=1,NPFT
          fsmc_ave(npft1)= fsmc_ave(npft1) + fsmcx(npft1)*FTIME_PHEN 
        ENDDO
!
        IF(L_PHEN) then
            DO npft1=1,NPFT
               fsmc(npft1)= fsmc_ave(npft1)
               fsmc_ave(npft1)=0.0
            ENDDO
        endif
!
        DO npft1=1,NPFT
           STHU(npft1) = WWW2N(npft1) 
           TS1(npft1) = TDn(npft1)
           V_WILT(npft1) = V_WILT0(npft1)
           V_SAT(npft1) = V_SAT0(npft1)
!yliu09May2017 PAR --
           FPARn(npft1) = FPARBDn(npft1)
           KPARn(npft1) = KPARBDn(npft1)
!!!--
        ENDDO
!
!            print*,'anetc',anetc
!            print*,'RDC_TRIF',RDC_TRIF
!            print*,'FSMC',FSMC
!            print*,'sthu',sthu
!            print*,'T_TIRF',t_trif
!            print*,'Ts1',ts1
!            print*,'Tstar',tstar
!
          CALL SITRIFFID (     &                                 
               TIMESTEP,       & !IN (1) Model TIMESTEP
               L_PHEN,         & !IN (1) Logical
               DAY_PHEN,       & !IN (1) Time
               L_TRIF,         & !IN (1) Logical
               DAY_TRIF,       & !IN (1) Time
               VEG_EQUIL     , & !IN (1) Logical
               ANETC,          & !IN (NPFT) Anet fraction    ,input from ssib2
               RDC_TRIF,       & !IN (NPFT) RDC fraction     ,input from ssib2
               FSMC,           & !IN (NPFT) RSTFAC(1,2)      ,input from ssib2
!zzq2008May12            STHU,     &  !IN (1) www(2) average      ,input from ssib2 
!zzq2008May12            T_TRIF,   &  !IN (1) TA average          ,input from ssib2 
!zzq2008May12            TS1,      &  !IN (1) TD average          ,input from ssib2
               STHU,           & !IN (npft1) www(2)          ,input from ssib2 
               TS1,            & !IN (NPFT1) TD              ,input from ssib2
               TSTAR,          & !IN (NPFT) TC               ,input from ssib2
!zzq2008May12            V_SAT,     & !IN (1) SPOROS prescribed
!zzq2008May12            V_WILT,    & !IN (1) wilting point prescribed
               V_SAT,          & !IN (npft1) SPOROS prescribed
               V_WILT,         & !IN (npft1) wilting point prescribed
!yliu09May2017 PAR --             
               FPARn,          & !IN FPAR
               KPARn,          & !IN KPAR
!!!--
               FRACA,          & !IN (1) Agriculture factor prescribed
               G_ANTH,         & !IN (1) Anthropogenic factor prescribed
               AGR,            & !IN Huilin read agriculture fraction
               CS,             & !INOUT(1)     diagnosed
               FRAC_TRIF,      & !INOUT(NTYPES)diagnose
               HT_TRIF,        & !INOUT(NPFT)  diagnose
               LAI,            & !INOUT(NPFT)  diagnose
               G_LEAF_DAY,     & !INOUT(NPFT)  diagnose
               G_LEAF_DR,      & !INOUT(NPFT)  diagnose
               GPP_DR,         & !INOUT(NPFT)  accumulate& initialize
               NPP_DR,         & !INOUT(NPFT)  accumulate& initialize
               RESP_S_DR,      & !INOUT(1)     accumulate& initialize
               RESP_W_DR,      & !INOUT(NPFT)  accumulate& initialize
               VEG_FRAC,       & !INOUT(1)     accumulate& initialize
               FRAC_VS,        & !INOUT(1)     accumulate& initialize
               LEAF_LOSS_DR,   & !INOUT(NPFT) leaf loss due to fire Huilin
               WOOD_LOSS_DR,   & !INOUT(NPFT) wood loss due to fire Huilin
               ROOT_LOSS_DR,   & !INOUT(NPFT) rood loss due to fire Huilin   
               GPP,            & !OUT (NPFT)
               NEP,            & !OUT (1)
               NPP,            & !OUT (NPFT)
               RESP_P,         & !OUT (NPFT)
               RESP_S,         & !OUT (1)
               LEAF1,          & !INOUT (NPFT) Huilin change from OUT to INOUT to build carbon pool
               ROOT1,          & !INOUT (NPFT) Huilin change from OUT to INOUT to build carbon pool
               WOOD1,          & !INOUT (NPFT) Huilin change from OUT to INOUT to build carbon pool
               DCVEG1,         & !OUT (NPFT)
               LIT_C,          & !OUT (NPFT)
               LIT_C_T,        & !OUT (1)
               phen,           & !OUT (NPFT)
!zzq2012sept24
               laid,           & !OUT (NPFT), dead LAI
               green_ssib,     & !OUT (NPFT), green fraction
!yliu09May2017 --
               RESP_P_M,       &
               RESP_P_G        &
             )
!
!zzq fcs will not be calculated using empirical formula.
!zzq FCS = 6.46E-7*EXP((TDsave-273.2-24.3)/10.*LOG(2.3))
       Fcs=resp_s*1e3/12
       Fca=Fcs-aveanet
!zzq
       Cacn=Co2Ca-(aveanet-Fcs)*PSURF*   &
               1.4*rasave/(44.6*273.2/Tasave*PSURF/1.013e5)
!
       anetdgb = 0.0
       anetdgs = 0.0
       DO NPFT1 =1, NNPFT
           if(frac_trif(npft1).gt.0.0) then
            anetdgb =anetdgb +           &
                 frac_trif(npft1)*anetn(npft1)/gbn(npft1)
            anetdgs =anetdgs +           &
                 frac_trif(npft1)*anetn(npft1)/gsn(npft1)
          endif
       ENDDO
!
!zzq cs is calculated from TRIF
!zzq  cs  = cacn - 1.4*psurf*anetdgb
      ci  = cs   - 1.6*psurf*anetdgs
!
! Diagnose gridbox mean LAI and vegetation fraction for SSiB.
!------------------------------------------------------------------------
        IF (L_TRIF.OR.L_PHEN) THEN
!-----------------------------------------------------------------------
! Diagnose a Stem Area Index (SAI) using the balanced (Maximum) LAI
! from TRIFFID, and assuming that SAI is 0.5 of this value (consistent
! with BATS parameter table)
!-----------------------------------------------------------------------
            SAI_MEAN=0.0
            VEGCOVER=0.0
            DO N=1,NPFT
              SAI(N) = DSAI_DLAI(N) *                      &
                 ((A_WS(N)*ETA_SL(N)*HT_TRIF(N)/A_WL(N))   &
                         **(1.0/(B_WL(N)-1)))
              green_trif(N)=1-SAI(N)/(LAI(N)+SAI(N))
              SAI_MEAN=SAI_MEAN+SAI(N)*FRAC_TRIF(N)
              VEGCOVER=VEGCOVER+FRAC_TRIF(N)
            ENDDO

            FRAC_TRIF(NNPFT)=1-VEGCOVER
!----------------------------------------------------------------------
! Diagnose total LAI, as sum of green LAI and SAI
!----------------------------------------------------------------------
          tzlt(1)=SAI_MEAN
          tzlt(2)=0.0001
          tvcover(1)=0.0000
          tvcover(2)=0.0001
!
          DO N=1,NPFT
            tvcover(1)=tvcover(1)+FRAC_TRIF(N)
            tzlt(1)=tzlt(1)+LAI(N)*FRAC_TRIF(N)
          ENDDO
!          write(4,10)tzlt(1)
!10        format(f6.3)
          green_frac(1)=1-SAI_MEAN/tzlt(1)
          green_frac(2)=0.0001
        END IF
!
        ITRIFSTEP=ITRIFSTEP+1
!
!zzq this condition is modified according to the orignial, because
! the TIMESTEP is determined and the integrating number for one day 
! is also known. 
!
        if(ITRIFSTEP.GT.(24*3600/TIMESTEP)) then
          KTRIFDAY=KTRIFDAY+1
          ITRIFSTEP=1
        endif
!
       dataios % xKTRIFDAY =KTRIFDAY
       dataios % xITRIFSTEP=ITRIFSTEP
       dataios % xCo2CS    =CS
       dataios % xRESP_S_DR=RESP_S_DR
       dataios % xVEG_FRAC =VEG_FRAC
       dataios % xFRAC_VS  =FRAC_VS
       dataios % xBAF60    = BAF60
!

       DO N=1,NPFT
         dataios % xGPP_DR(N)    = GPP_DR(N)
         dataios % xNPP_DR(N)    = NPP_DR(N)
         dataios % xRESP_W_DR(N) = RESP_W_DR(N)
         dataios % xG_LEAF_DAY(N)= G_LEAF_DAY(N)
         dataios % xG_LEAF_DR(N) = G_LEAF_DR(N)
         dataios % xHT_TRIF(N)   = HT_TRIF(N)
         dataios % xLAI_TRIF(N)  = LAI(N)
         dataios % xSAI_TRIF(N)  = SAI(N)
         dataios % xGREEN_TRIF(N)= GREEN_TRIF(N)
!zzq2012sept24
         dataios % xLAID_TRIF(N)  = LAID(N)
         dataios % xGREEN_SSIB(N) = GREEN_SSIB(N)
!
         dataios % xPHEN(N)      = PHEN(N)
         dataios % xFSMC(N)      = FSMC(N)
         dataios % xFSMC_ave(N)  = FSMC_ave(N)
!Huilin add for fire output
         dataios % xLEAF1(N)     = LEAF1(N)    
         dataios % xWOOD1(N)     = WOOD1(N)    
         dataios % xROOT1(N)     = ROOT1(N)    
         dataios % xLEAF_LOSS_DR(N) = LEAF_LOSS_DR(N) 
         dataios % xWOOD_LOSS_DR(N) = WOOD_LOSS_DR(N)
         dataios % xROOT_LOSS_DR(N) = ROOT_LOSS_DR(N)
       ENDDO
!
       do n=1,NNPFT
          dataios % xVEG_FRACN(N) =FRAC_TRIF(N)
       enddo
      END SUBROUTINE
END MODULE
