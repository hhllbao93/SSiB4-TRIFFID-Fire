MODULE ssib4_module_sfcssib2       
      use ssib4_module_comsconst
      use ssib4_module_comssib4sfc
      use ssib4_module_weighting
      use ssib4_module_ssibsub
      use ssib4_module_ssibsfcbio
!
     IMPLICIT NONE

contains
      subroutine sfcssib2(                                         &
! for ssib2
         dataios,                                                  &
         ISTRT,MONTH,JDAY,DDTT,                                    &
         SUNANG,radn11,radn12,radn21,radn22,RLWDOWN,               &
         PPL1,PPC1,TM,UMM1,VMM1,QM,EM,RHOAIR,                      &
         PSURF,SIG,ITYPE,co2ca,ZLAT,ZWIND,                         &
!out.......................................................
         TD1,TG1,TC1,TA,                                           &
         WWW1,WWW2,WWW3,CAPAC1,CAPAC2,                             &
         SALB11,SALB12,SALB21,SALB22,                              &
         radt1,radt2,                                              &
         rstfac1,rstfac2,rstfac3,rstfac4,                          &
         ra1,rb1,rd1,rsoil,shf,chf,                                &
         ETMASS,EVAPSOIL,                                          &
         EVAPWC,EVAPDC,EVAPSN,EVAPGX,                              &
         GHTFLX,xhlflx,xhsflx,                                     &
         USTAR,DRAG,DRAGU,DRAGV,                                   &
         TGEFF,BEDO,roff,soilm,                                    &
         fsdown,fldown,fsup,flup,                                  &
         umom,vmom,                                                &
         hlflx,hsflx,ulwsf1,evap,q2m,rh,pilphr,                    &!Huilin Huang Feb 7 2018
         plantr,canopy,sheleg,cm,ch,fm,fh,rib,ZLT1,                &
         aveanet,avecrd,fcs,fca,cacn,cs,ci,gb,gs,                  &
         xwcl,xwel,xwsl,xwcs,xwes,xwss,                            &
         aveft1,aveft2, avehs )
!----------------------------------------------------------------------
!
       type(ssibdataios):: dataios
       INTEGER:: MONTH,JDAY, ITYPE, NPFT1,IL,ISTRT

       REAL DDTT
       real psurf,sig, zlat,zwind
       LOGICAL  L_VEG

       real                                                        &
        SUNANG, RLWDOWN,PPL1,PPC1,TM,UMM1,VMM1,QM, EM,             &
        RHOAIR,radn11,radn12,radn21,radn22,radnsum

      real,dimension(2)  ::CAPAC,RADT
      real,dimension(3)  ::WWW
      real,dimension(2,2)::SALB,RADFRAC
      real,dimension(3,2)::RADN
      real,dimension(2,4)::RSTFAC
      real,dimension(2)  :: zlt, green,vcover
!
      real xhlflx,xhsflx
      real fcs,fca,cs,ci,gb,gs,co2ca
      real &
         fsdown,fldown,fsup,flup,                    &
         umom,vmom,hlflx,hsflx,ulwsf1,evap,q2m,ei,rh,&  !Huilin Huang
         plantr,canopy,sheleg,cm,ch,fm,fh,rib,ZLT1,  &
         xwcl,xwel,xwsl,xwcs,xwes,xwss               
      real capac1,capac2
      real www1,www2,www3,pilphr                       !Huilin Huang output pilphr Feb 7 2018  
      real rstfac1,rstfac2,rstfac3,rstfac4
      real ra1,rb1,rd1,rsoil,shf,chf
      real TASAVE,TDSAVE,RASAVE
      real TG1,TC1,TD1,TGEFF
      real etmass,EVAPSOIL,EVAPWC
      real EVAPDC,EVAPSN,EVAPGX,GHTFLX,USTAR
      real DRAG,DRAGU,DRAGV,RADT1,RADT2,bedo,roff,soilm
      real salb11,salb12,salb21,salb22
      real aveft1,aveft2,avehs
      real aveanet,anetdgb,anetdgs,avecrd,cacn
!
      real,dimension(NTYPES):: xhlflxn,xhsflxn,gbn,gsn
      real,dimension(NTYPES)::                           &
           fsdownn,fldownn,fsupn,flupn,                  &
           umomn,vmomn,hlflxn,hsflxn,ulwsf1n,evapn,q2mn, &
           plantrn,cmn,chn,fmn,fhn,ribn,ZLT1n,           &
               xwcln,xweln,xwsln,xwcsn,xwesn,xwssn
!
      real,dimension(NTYPES)::                           &
           WWW1N,WWW2N,WWW3N,CAPAC1N,CAPAC2N,TGEFFN,     &
           TGSN,TCN,TDN,TASAVN,TDSAVN,RASAVN,            & 
           RSTFAC11N,RSTFAC12N,RSTFAC13N,RSTFAC14N,      &
           rsoiln,shfn,chfn,ETMASSN,EVAPSOILN,           &
           EVAPWCN,EVAPGXN,EVAPDCN,EVAPSNN,              &
           GHTFLXN,ANETN,USTARN,DRAGN,DRAGUN,DRAGVN,     &
           RADT1N,RADT2N,ROFFN,SOILMN,BEDON,             &
           salb11n,salb12N,salb21N,salb22N,              & 
           PILPHRn                                   !Huilin output HR Feb 7 2018
!
      real,dimension(NTYPES)::                           &
           ft1n,ft2n,hsn,TAIRN,ran,rbn,Crdn
      real,dimension(NTYPES)::  veg_fracn
!
      real,dimension(2,3,2):: TRAN,REF
      real,dimension(2)    :: CHIL,TOPT,TL,TU,DEFAC,PH1,PH2,ROOTD
      real,dimension(3)    :: SOREF,ZDEPTH
      real,dimension(2,3)  :: RSTPAR
      real:: poros,vm0,fc3,z0,xdd,z2,z1,  &
          rdc,rbc,bee,phsat,satco,slope
!
!yliu18Oct2016
      real,dimension(NTYPES)::FPARBDn,KPARBDn
!
      real,dimension(2)    :: CAPAC01,CAPACM1
      real,dimension(3)    :: W01,WM1
      real:: TD01,TG01,TC01,TDM1,TGM1,TCM1
      real:: RA,RB,zfw1,filtr,ea,TA,TAx
      real:: xneglect

      xneglect =1.0e-4

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
!-------------------------------------
!rr
! **  Read in vegetation parameters
!
      CALL VEGOUT1(POROS,ZDEPTH,ITYPE)

      DO 1001 NPFT1=1,NTYPES
 
       if(npft1.eq.itype) then
            dataios % xveg_fracn(npft1)=1
       else
            dataios % xveg_fracn(npft1)=0
       endif

       veg_fracn(npft1) = dataios % xveg_fracn(npft1)
       if(veg_fracn(npft1).lt.xneglect) goto 1001

       L_VEG = .TRUE.
!       IF(NPFT1.GT.6) L_VEG=.FALSE.

      call VEGOUT2(TRAN,REF,GREEN,VCOVER,CHIL,              &
           RSTPAR,VM0,TOPT,TL,TU,FC3,DEFAC,PH1,PH2,         &
           ZLT,Z0,XDD,Z2,Z1,RDC,RBC,ROOTD,                  &
           SOREF,BEE, PHSAT,SATCO,SLOPE,                    &
           MONTH,ITYPE)

      IF (npft1.EQ.9) CALL CROPS(ZLAT,JDAY,CHIL,           &
             ZLT,GREEN,VCOVER,RSTPAR,TOPT,TL,TU,DEFAC,PH2,PH1)
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
!
        zFW1 = dataios % xZFW(npft1)
        RA   = dataios % xCO2RA(npft1)
        RB   = dataios % xCO2RB(npft1)
        EA   = dataios % xea(npft1)
        filtr= dataios % xfiltr(npft1)
        TAx  = dataios % xTA(npft1)
!
       call SFCBIO(                                           &
!in ......................................................
         DDTT,SUNANG, RADFRAC,RADN,RLWDOWN,                   &
         PPL1,PPC1,TM,UMM1,VMM1,QM,EM,RHOAIR,                 &
         PSURF,SIG,ZWIND,co2ca,L_VEG,ZLAT,MONTH,              &  
! vegetation parameters
         TRAN,REF,GREEN,VCOVER,CHIL,                          &
         RSTPAR,VM0,TOPT,TL,TU,fc3,DEFAC,PH1,PH2,             &
         ZLT,Z0,XDD,Z2,Z1,RDC,RBC,ROOTD,SOREF,BEE,            &
         PHSAT,SATCO,SLOPE,POROS,ZDEPTH,                      &
!inout.....................................................
         TD01,TG01,TC01,W01,CAPAC01,                          &
         TDM1,TGM1,TCM1,WM1,CAPACM1,                          &
         ra,rb,zfw1,filtr,ea,TAx,                             &
!out.......................................................
         TASAVN(npft1), TDSAVN(npft1),  RASAVN(npft1),        &
         Tdn(npft1),    Tgsn(npft1),    Tcn(npft1),           &
         WWW,    CAPAC,                                       &
         ETMASSN(npft1),                                      &
         EVAPSOILn(npft1), EVAPWCn(npft1), EVAPDCn(npft1),    &
         EVAPSNn(npft1),   EVAPGXn(npft1), GHTFLXn(npft1),    &
         Xhlflxn(NPFT1),   Xhsflxn(NPFT1), USTARn(npft1),     &
         Dragn(npft1),     Dragun(npft1),  Dragvn(npft1),     &
         Tgeffn(npft1),    Bedon(npft1),                      &
         SALB,    RADT,                                       &
         Roffn(npft1),   soilmn(npft1), Anetn(npft1),         &
         Gbn(npft1),     Gsn(npft1),    Crdn(npft1),          &
         RSTFAC,                                              &
         rsoiln(NPFT1),  shfn(NPFT1),    chfn(NPFT1),         &
         Fsdownn(NPFT1), Fldownn(NPFT1), Fsupn(NPFT1),        &
         Flupn(NPFT1),   Umomn(NPFT1),   Vmomn(NPFT1),        &
         Hlflxn(NPFT1),  Hsflxn(NPFT1),  Ulwsf1n(NPFT1),      &
         Evapn(NPFT1),   Q2mn(NPFT1),    Plantrn(NPFT1),      &
         Cmn(NPFT1),     Chn(NPFT1),                          &
         Fmn(NPFT1),     Fhn(NPFT1),                          &
         RIbn(NPFT1),                                         &
         Zlt1n(NPFT1),   Xwcln(NPFT1),   Xweln(NPFT1),        &
         Xwsln(NPFT1),                                        &
         Xwcsn(NPFT1),   Xwesn(NPFT1),   Xwssn(NPFT1),        &
!yliu09May2017 PAR --         Ft1n(npft1),    Ft2n(npft1),    hsn(npft1) )    
         PILPHRn(npft1),                                      &   !Huilin for hr output  Feb 7th 2018   
         Ft1n(npft1),    Ft2n(npft1),    hsn(npft1),          &
         FPARBDn(npft1), KPARBDn(npft1) )    

         dataios % xZFW(npft1)  = zFW1
         dataios % xCO2RA(npft1)= RA
         dataios % xCO2RB(npft1)= RB
         dataios % xTA(npft1)  = TAx
         dataios % xea(npft1)   = EA
         dataios % xfiltr(npft1)= filtr
!
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
         ran(NPFT1)=ra
         rbn(NPFT1)=rb
         CAPAC1N(NPFT1)=CAPAC(1)
         CAPAC2N(NPFT1)=CAPAC(2)
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
         TAIRN(NPFT1)=TAx

!yliu CROP 30Sep2016
         dataios%xVEG_FRACN(NPFT1)  = VCOVER(1)
         if ( VCOVER(1).ne.0 ) then
           dataios%xLAI_TRIF(NPFT1)  = zlt(1)/VCOVER(1)
         else
           dataios%xLAI_TRIF(NPFT1)  = 0.0001
         endif
         dataios%xGREEN_TRIF(NPFT1) = GREEN(1)
!
 1001    continue
!
         RA1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RAN)
         RB1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RBN)
         RD1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,CRDN)
         TA     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TAIRN)
         TC1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TCN)
         TG1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TGSN)
         TD1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TDN)
         TASAVE = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TASAVN)
         TDSAVE = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TDSAVN)
         RASAVE = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RASAVN)
         CAPAC1 = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,CAPAC1N)
         CAPAC2 = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,CAPAC2N)
         WWW1   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,WWW1N)
         WWW2   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,WWW2N)
         WWW3   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,WWW3N)
         ETMASS    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ETMASSN)
         EVAPSOIL  = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,EVAPSOILN)
         EVAPWC    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,EVAPWCN)
         EVAPDC    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,EVAPDCN)
         EVAPSN    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,EVAPSNN)
         EVAPGX    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,EVAPGXN)
         GHTFLX    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,GHTFLXN)
         xhlflx    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xhlflxN)
         xhsflx    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xhsflxN)
         USTAR     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,USTARN)
         DRAG      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,DRAGN)
         DRAGU     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,dRAGUN)
         DRAGV     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,DRAGVN)
         RADT1     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RADT1N)
         RADT2     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RADT2N)
         TGEFF     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,TGEFFN)
         bedo      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,bedon)
         roff      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ROFFN)
         soilm     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,SOILMN)
         fsdown    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,fsdownN)
         fldown    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,fldownN)
         fsup      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,fsupN)
         flup      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,flupN)
         umom      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,umomN)
         vmom      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,vmomN)
         hlflx     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,hlflxN)
         hsflx     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,hsflxN)
         ulwsf1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ulwsf1N)
         evap      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,evapN)
         q2m       = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,q2mN)
         PILPHR    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,PILPHRn)        !Huilin Feb 7 2018
!
!  --- ... canopy water and snow water equivalent
          canopy = capac1 * 1000.0
          sheleg = capac2 * 1000.0
!
         cm        = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,cmN)
         ch        = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,chN)
         fm        = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,fmN)
         fh        = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,fhN)
         rib       = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ribN)
         ZLT1      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ZLT1N)
         xwcl      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xwclN)
         xwel      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xwelN)
         xwsl      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xwslN)
         xwcs      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xwcsN)
         xwes      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xwesN)
         xwss      = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,xwssN)
         salb11    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,salb11N)
         salb12    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,salb12N)
         salb21    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,salb21N)
         salb22    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,salb22N)
         gb        = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,gbN)
         gs        = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,gsN)
         aveanet   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,anetN)
         avecrd    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,crdN)
         aveft1    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ft1N)
         aveft2    = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,ft2N)
         avehs     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,hsN)
         RSTFAC1   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RSTFAC11N)
         RSTFAC2   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RSTFAC12N)
         RSTFAC3   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RSTFAC13N)
         RSTFAC4   = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,RSTFAC14N)
         rsoil     = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,rsoiln)
         shf       = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,shfn)
         chf       = summ(1,NTYPES,NTYPES,veg_fracn,xneglect,chfn)
!
        anetdgb = 0.0
        anetdgs = 0.0
        DO NPFT1 = 1, NTYPES
           if(veg_fracn(npft1).ge.xneglect) then
            anetdgb =anetdgb + veg_fracn(npft1)*anetn(npft1)/gbn(npft1)
            anetdgs =anetdgs + veg_fracn(npft1)*anetn(npft1)/gsn(npft1)
          endif
       ENDDO
!
      FCS = 6.46E-7*EXP((TDsave-273.2-24.3)/10.*LOG(2.3))
      FCA = FCS - aveANET
      CACN=CO2CA-(aveANET-FCS)*          &
         PSURF*1.4*RAsave/(44.6*273.2/TAsave*PSURF/1.013E5)
!
      cs  = cacn - 1.4*psurf*anetdgb
      ci  = cs   - 1.6*psurf*anetdgs
      ei  = 100.*exp(21.18123-5418./tm)/.622       !Huilin add ei
      rh  = em/ei*100                              !Huilin add rh

      if (ci.gt.60.) then
      print 777,cs,aveanet,psurf,gs,ci
      endif
!khs  cacn = amax1(cacn,1.e-3)
!khs  cs   = amax1(cs,  1.e-3)
!khs  ci   = amax1(ci,  1.e-3)
  777  format(1x,'sfcbio: cs an p gs ci = ',1x,5e12.5)
!
       END SUBROUTINE
END MODULE                                                           
