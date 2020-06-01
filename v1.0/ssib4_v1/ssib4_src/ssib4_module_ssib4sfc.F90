MODULE ssib4_module_ssib4sfc
       use ssib4_module_comsconst
       use ssib4_module_comssib4sfc
       use ssib4_module_sfcssib2
       use ssib4_module_sfctrif
       use ssib4_module_pressibtrif
!
!#define DEBUG_PRINT
!
       IMPLICIT NONE
contains
!=========================================================================
! GCM_ssib4 interface
!=========================================================================
      subroutine ssib4sfc(                                         &
! --- input:
        DDTT,SUNANG,RADN11,RADN12,RADN21,RADN22,RLWDOWN,PPC,PPL,   &
        SIG,TM,UMM,VMM,QM,PSURF,ZWIND,RHOAIR,TDD,KREC,             &
        LIGHT,DP,AGR,GDP,PEATF,                                    &  !Huilin add Light, DP, AGR, GDP, PEATF 
        zlon,zlat,IYEAR0,MONTH0,IDAY0,FHOUR0,THOUR,ISTRT,          &
! --- input/output:
        iosdat,                                                    & 
! --- output:
        xhsflx,xhlflx,GHTFLX,RIB,                                  &
        EVAPSOIL,EVAPWC,EVAPDC,EVAPSN,EVAPGX,                      &
        SALB11,SALB12,SALB21,SALB22,                               &
        fsdown,fldown,fsup,flup,                                   &
        TA,TGEFF,ROFF,q2m,USTAR,soilm,plantr,                      &!Huilin Huang
        cm,ch,fm,fh,                                               &
! --- other outputs:
        outdat                                                     &
        )
!
      integer,intent(in) :: IYEAR0,MONTH0,IDAY0
      real,   intent(in) :: FHOUR0, THOUR, DDTT
      INTEGER:: MONTH,MON_COR,JDAY,iyear,ITYPE,K,ISTRT,KREC
      real::  psurf,sig,tdd
      real:: umm1,vmm1,RCL,ZWIND,wind,LIGHT,DP,AGR,GDP,PEATF               !Huilin add light dp agr gdp
!
      real::                                                    &
        SUNANG,RLWDOWN,radn11,radn12,radn21,radn22,RADNSUM,     &
        TM,UMM,VMM,QM,EM,rhoair,co2ca,zlon,zlat
      real:: PPL,PPC,PREC,PREC60_AVG       !(mm) Huilin add PREC Jul 2019
!
      real:: xhlflx,xhsflx
      real:: cs,ci,gb,gs
      real::                                                    &
         fsdown,fldown,fsup,flup,                               &
         umom,vmom,hlflx,hsflx,ulwsf1,evap,q2m,                 &!Huilin Huang
         plantr,canopy,sheleg,cm,ch,fm,fh,rib,cacn,ZLT1,        &
         xwcl,xwel,xwsl,xwcs,xwes,xwss     
      real:: capac1,capac2
      real:: www1,www2,www3
      real:: rstfac1,rstfac2,rstfac3,rstfac4
      real:: ra1,rb1,rd1,rsoil,shf,chf
      real:: TA,TGS,TC,TD,TGEFF
      real:: etmass,EVAPSOIL,EVAPWC
      real:: EVAPDC,EVAPSN,EVAPGX,GHTFLX,USTAR
      real:: DRAG,DRAGU,DRAGV,RADT1,RADT2,bedo,roff,soilm
      real:: salb11,salb12,salb21,salb22
      real:: aveft1,aveft2, avehs
      real:: aveanet,avecrdn,fca,fcs
!
      real,dimension(2):: tvcover,tzlt,greenfrac
!yliu09May2017 addtype      real,dimension(6)::                                        &
      real,dimension(NTYPET)::                                        &
         GPP,NPP,RESP_P, LEAF1,ROOT1,WOOD1,DCVEG1, LIT_C
!yliu09May2017
      real,dimension(NTYPET):: RESP_P_M,RESP_P_G,ANETC,RDC_TRIF
      real:: RESP_S,LIT_C_T,NEP,RH,PILPHR
      real:: AGB_AVG,FEO_AVG,FES_AVG,FE1_AVG,FE2_AVG,FD_AVG
      real:: FCLI,FSAT                                                !Huilin add fire output
      real:: BAF_AVG,NFIRE_AVG,BURN_AVG                               !Huilin add fire output
      real,dimension(NTYPET):: BURN_FIRE                              !Huilin add fire output  Mar.2019
      real:: EM_CO2,EM_CO,EM_CH4,EM_NHMC,EM_H2,EM_NOX,EM_N2O          !Huilin add fire output  May.2019
      real:: EM_PM25,EM_TPM,EM_TC,EM_OC,EM_BC

      type(ssibdataios) :: iosdat                                                         
      type(ssibdataout) :: outdat
!.......
      real:: cam
      real,dimension(2,2)  :: RADFRAC
      real,dimension(3,2)  :: RADN
!rr
       itype=int(iosdat % xvtype)
       call qtime(IYEAR0,MONTH0,IDAY0,FHOUR0,THOUR,IYEAR,MONTH,JDAY)
       MON_COR = correct_month(MONTH,ZLAT)
!
      UMM1=UMM 
      VMM1=VMM

!zzq080321  UM=SQRT(UMM1**2+VMM1**2)

!rr ratko June 29,2007
!rr   UM = AMAX1(UM,2.)

!zzq080321   WIND=SQRT(RCL*UMM1*UMM1+VMM1*VMM1)
!zzq080321   WIND=AMAX1(WIND,1.0)
!zzq080321 the wind variable is not used by ssib4, only for GCM.
!zzq11/2013      WIND=RCL*SQRT(UMM1*UMM1+VMM1*VMM1)
!
!     END OF EMPIRICAL EQUATIONS
!     *********************************************************
!
!=========================================================================
! CONVERT TO VAPOR PRES. TO MB,(psurf, pa)
!=========================================================================
        EM=(psurf*qm)/62.2

!Huilin Jul 2019
       PREC = PPL + PPC
!zzq2012OCT08    CALL CO2ATM(IYEAR,CAM)
       CALL CO2ATM(IYEAR,CAM)   !Huilin real
!      CALL CO2ATM(1948,CAM)
! overwrite the input co2ca
       CO2CA   =CAM/10.
!
       if(istrt.eq.1) then
          call init_ssib(iosdat,MON_COR,TM,EM,tdd,PREC)
       endif
!
!!!![yliu 16Jun2016] set zero for ssib4 output data
       DO K=1,2
         outdat % qtvcover(K)   = 0.0
         outdat % qtzlt(K)      = 0.0
         outdat % qgreenfrac(K) = 0.0
       ENDDO
!
       outdat % QNEP      = 0.0
       outdat % QRESP_S   = 0.0
       outdat % QLIT_C_T  = 0.0
       outdat % qhs       = 0.0
       outdat % QRH       = 0.0 !Huilin Huang
       outdat % QHR       = 0.0 !Huilin Huang Feb 7 2018
       outdat % QPREC60_AVG = 0.0 !Huilin Huang
       outdat % QAGB_AVG  = 0.0 !Huilin Huang
       outdat % QFEO_AVG  = 0.0 !Huilin Huang
       outdat % QFCLI     = 0.0 !Huilin Huang
       outdat % QFSAT     = 0.0 !Huilin Huang
       outdat % QFES_AVG  = 0.0 !Huilin Huang
       outdat % QFE1_AVG  = 0.0 !Huilin Huang
       outdat % QFE2_AVG  = 0.0 !Huilin Huang
       outdat % QFD_AVG   = 0.0 !Huilin Huang
       outdat % QBAF_AVG  = 0.0 !Huilin Huang
       outdat % QBURN_AVG = 0.0 !Huilin Huang Mar 2019
       outdat % QNFIRE_AVG= 0.0 !Huilin Huang
       outdat % QPEATF    = 0.0 !Huilin Huang
       outdat % QEM_CO2   = 0.0 
       outdat % QEM_CO    = 0.0 
       outdat % QEM_CH4   = 0.0 
       outdat % QEM_NHMC  = 0.0 
       outdat % QEM_H2    = 0.0 
       outdat % QEM_NOX   = 0.0 
       outdat % QEM_N2O   = 0.0 
       outdat % QEM_PM25  = 0.0 
       outdat % QEM_TPM   = 0.0 
       outdat % QEM_TC    = 0.0 
       outdat % QEM_OC    = 0.0 
       outdat % QEM_BC    = 0.0 
!
       DO K=1,NTYPET
         outdat % QWOOD1(K)    = 0.0
         outdat % QROOT1(K)    = 0.0
         outdat % QLEAF1(K)    = 0.0
         outdat % QDCVEG1(K)   = 0.0
         outdat % QRESP_P(K)   = 0.0
         outdat % QLIT_C(K)    = 0.0
       ENDDO
!yliu CROP 30Sep2016
       DO K=1,NTYPES
         outdat % QNPP(K)      = 0.0
         outdat % QGPP(K)      = 0.0
!yliu18Oct2016
         outdat % QRESP_P_M(K) = 0.0
         outdat % QRESP_P_G(K) = 0.0
         outdat % QRDC         = 0.0
         outdat % QANET        = 0.0
         outdat % QBURN_FIRE(K)= 0.0
       ENDDO
       IF(iosdat%xschsel.eq.1.0) THEN

      CALL trif_prescribed(DDTT)
      IF(ISTRT.EQ.1) THEN
        CALL init_trif(iosdat,zlon,zlat,DDTT,MON_COR,JDAY,ITYPE)
      ENDIF
         call sfctrif(                                                 & 
         iosdat,                                                       & 
         ISTRT,MON_COR,JDAY,DDTT,                                      &
         SUNANG,radn11,radn12,radn21,radn22,RLWDOWN,                   &
         PPL,PPC,TM,UMM1,VMM1,QM,EM,rhoair,KREC,                       &
         PSURF,SIG,ITYPE,Co2Ca,ZLAT,ZWIND,                             &
         LIGHT,DP,AGR,GDP,PEATF,                                       &!Huilin add all input in this line 
!out.......................................................
         TD,TGS,TC,TA,                                                 &
         WWW1,WWW2,WWW3,CAPAC1,CAPAC2,                                 &
         SALB11,SALB12,SALB21,SALB22,                                  &
         radt1,radt2,                                                  &
         rstfac1,rstfac2,rstfac3,rstfac4,                              &
         ra1,rb1,rd1,rsoil,shf,chf,                                    &
         ETMASS,EVAPSOIL,                                              &
         EVAPWC,EVAPDC,EVAPSN,EVAPGX,                                  &
         GHTFLX,xhlflx,xhsflx,                                         &
         USTAR,DRAG,DRAGU,DRAGV,                                       &
         TGEFF,BEDO,roff,soilm,                                        &
         fsdown,fldown,fsup,flup,                                      &
         umom,vmom,                                                    &
         hlflx,hsflx,ulwsf1,evap,q2m,rh,pilphr,prec60_avg,             &  !Huilin Feb 7 2018
         plantr,canopy,sheleg,cm,ch,fm,fh,rib,ZLT1,                    &
         aveanet,avecrdn,fcs,fca,cacn,cs,ci,gb,gs,                     &
         xwcl,xwel,xwsl,xwcs,xwes,xwss,                                &
         aveft1,aveft2,avehs,                                          &
!TRIF output to GCM
         tvcover,tzlt,greenfrac,                                       &
         GPP,NEP,NPP,RESP_P,RESP_S,                                    &
         LEAF1,ROOT1,WOOD1,DCVEG1,                                     &
         LIT_C,LIT_C_T,                                                &
         AGB_AVG,FEO_AVG,FES_AVG,FE1_AVG,FE2_AVG,FD_AVG,FCLI,FSAT,     &
         BAF_AVG,NFIRE_AVG,BURN_AVG,BURN_FIRE,                         &  !Huilin add for fire output       
         EM_CO2,EM_CO,EM_CH4,EM_NHMC,EM_H2,EM_NOX,EM_N2O,              &
         EM_PM25,EM_TPM,EM_TC,EM_OC,EM_BC,                             &   
!yliu09May2017 
         RESP_P_M,RESP_P_G,ANETC,RDC_TRIF)
!
       DO K=1,2
         outdat % qtvcover(K)   = tvcover(K)
         outdat % qtzlt(K)      = tzlt(K)
         outdat % qgreenfrac(K) = greenfrac(K)
       ENDDO
!
       outdat % QNEP      = NEP
       outdat % QRESP_S   = RESP_S
       outdat % QLIT_C_T  = LIT_C_T
       outdat % qhs       = avehs
!Huilin Huang add fire variables
       outdat % qrh       = rh
       outdat % qhr       = pilphr
       outdat % qagb_avg  = agb_avg
       outdat % qfeo_avg  = feo_avg
       outdat % qfes_avg  = fes_avg
       outdat % qfe1_avg  = fe1_avg
       outdat % qfe2_avg  = fe2_avg
       outdat % qfd_avg   = fd_avg
       outdat % qfcli     = fcli 
       outdat % qfsat     = fsat 
       outdat % qpeatf    = peatf
       outdat % qbaf_avg  = baf_avg
       outdat % qburn_avg = burn_avg
       outdat % qnfire_avg= nfire_avg
       outdat % qem_co2   = em_co2
       outdat % qem_co    = em_co
       outdat % qem_ch4   = em_ch4
       outdat % qem_nhmc  = em_nhmc
       outdat % qem_h2    = em_h2
       outdat % qem_nox   = em_nox
       outdat % qem_n2o   = em_n2o
       outdat % qem_pm25  = em_pm25
       outdat % qem_tpm   = em_tpm
       outdat % qem_tc    = em_tc
       outdat % qem_oc    = em_oc
       outdat % qem_bc    = em_bc
       outdat % qprec60_avg = prec60_avg
!
       DO K=1,NTYPET
         outdat % QWOOD1(K)    = WOOD1(K)
         outdat % QROOT1(K)    = ROOT1(K)
         outdat % QLEAF1(K)    = LEAF1(K)
         outdat % QDCVEG1(K)   = DCVEG1(K)
         outdat % QNPP(K)      = NPP(K)
         outdat % QGPP(K)      = GPP(K)
         outdat % QRESP_P(K)   = RESP_P(K)
         outdat % QLIT_C(K)    = LIT_C(K)
!yliu09May2017 --
         outdat % QRESP_P_M(K) = RESP_P_M(K)
         outdat % QRESP_P_G(K) = RESP_P_G(K)
         outdat % QRDC(K)      = RDC_TRIF(K)
         outdat % QANET(K)     = ANETC(K)
         outdat % QBURN_FIRE(K)= BURN_FIRE(K) 
!!!--
       ENDDO
!
#ifdef DEBUG_PRINT
       print*,'trif',MONTH,JDAY,TGS,ETMASS,xHsFLX
#endif
      ELSE
!
      call sfcssib2(                                            &
! for ssib2
         iosdat,                                                &  
         ISTRT,MON_COR,JDAY,DDTT,                               &
         SUNANG,radn11,radn12,radn21,radn22,RLWDOWN,            &
         PPL,PPC,TM,UMM1,VMM1,QM,EM,rhoair,                     &
         PSURF,SIG,ITYPE,co2ca,ZLAT,ZWIND,                      &
!out.......................................................
         TD,TGS,TC,TA,                                          &
         WWW1,WWW2,WWW3,CAPAC1,CAPAC2,                          &
         SALB11,SALB12,SALB21,SALB22,                           &
         radt1,radt2,                                           &
         rstfac1,rstfac2,rstfac3,rstfac4,                       &
         ra1,rb1,rd1,rsoil,shf,chf,                             &
         ETMASS,EVAPSOIL,                                       &
         EVAPWC,EVAPDC,EVAPSN,EVAPGX,                           &
         GHTFLX,xhlflx,xhsflx,                                  &
         USTAR,DRAG,DRAGU,DRAGV,                                &
         TGEFF,BEDO,roff,soilm,                                 &
         fsdown,fldown,fsup,flup,                               &
         umom,vmom,                                             &
         hlflx,hsflx,ulwsf1,evap,q2m,rh,pilphr,                 &  !Huilin Feb 7 2018
         plantr,canopy,sheleg,cm,ch,fm,fh,rib,ZLT1,             &
         aveanet,avecrdn,fcs,fca,cacn,cs,ci,gb,gs,              &
         xwcl,xwel,xwsl,xwcs,xwes,xwss,                         &
         aveft1,aveft2,avehs)
!
#ifdef DEBUG_PRINT
       print*,'ssib2',MONTH,JDAY,TGS,ETMASS,xHsFLX
#endif
!yliu09May2017 CROP --
         outdat % qtvcover(1)   = iosdat%xVEG_FRACN(itype)
         outdat % qtzlt(1)      = iosdat%xLAI_TRIF(itype)
         outdat % qgreenfrac(1) = iosdat%xGREEN_TRIF(itype)
         outdat % QNPP(itype)   = aveanet*0.012
         outdat % QGPP(itype)   = (avecrdn+aveanet)*0.012
!!!0--        
       ENDIF
!
       outdat % tgs      = TGS
       outdat % tc       = tc
       outdat % td       = td
       outdat % ulwsf1   = ulwsf1
       outdat % hlflx    = hlflx
       outdat % hsflx    = hsflx
       outdat % evap     = evap
       outdat % umom     = umom
       outdat % vmom     = vmom
       outdat % zlt1     = zlt1
       outdat % canopy   = canopy
       outdat % sheleg   = sheleg
       outdat % etmass   = etmass
       outdat % radt1    = radt1
       outdat % radt2    = radt2
       outdat % ra       = ra1
       outdat % rb       = rb1
       outdat % rd       = rd1
       outdat % rsoil    = rsoil
       outdat % bedo     = bedo
       outdat % shf      = shf
       outdat % chf      = chf
       outdat % www1     = www1
       outdat % www2     = www2
       outdat % www3     = www3
       outdat % capac1   = capac1
       outdat % capac2   = capac2
       outdat % rstfac1  = rstfac1
       outdat % rstfac2  = rstfac2
       outdat % rstfac3  = rstfac3
       outdat % rstfac4  = rstfac4
       outdat % qhs      = avehs
       outdat % aveanet  = aveanet
       outdat % fcs      = fcs
       outdat % fca      = fca
       outdat % cacn     = cacn
       outdat % cs       = cs
       outdat % ci       = ci
       outdat % gb       = gb
       outdat % gs       = gs
       outdat % xwcl     = xwcl
       outdat % xwel     = xwel
       outdat % xwsl     = xwsl
       outdat % xwcs     = xwcs
       outdat % xwes     = xwes
       outdat % xwss     = xwss
       outdat % aveft1   = aveft1
       outdat % aveft2   = aveft2
       outdat % qrh      = rh           !Huilin huang add rh
       outdat % qhr      = pilphr
!
       END SUBROUTINE
!=========================================================================
!
      function correct_month(month,zlat) result(f_result)
!
!=========================================================================
!
      real:: zlat
      integer,intent(in)::month
      integer:: mon_new,f_result

       IF(ZLAT.LT.0.0) THEN
          MON_NEW=MONTH+6
          IF(MON_NEW.GT.12) MON_NEW=MON_NEW-12
       ELSE
          MON_NEW=MONTH
       ENDIF
!
       f_result=mon_new
!
       end function
!=========================================================================
!
      subroutine qtime(nbgnyr,nbgnmn,nbgndy,fbgnhr,thour,iyear,month,jday)
!
!=========================================================================
        IMPLICIT NONE
        integer,intent(in)::nbgnyr,nbgnmn,nbgndy
        real,intent(in)::thour,fbgnhr
        integer,intent(out)::iyear,month,jday
        integer::maxdays,ndays,totdays,im
        integer,dimension(12)::IDAYS
        DATA IDAYS/0,31,59,90,120,151,181,212,243,273,304,334/
!
        totdays = IDAYS(nbgnmn)+nbgndy+(int(fbgnhr+thour))/24   !the number of days from the year start
        iyear = nbgnyr
        if(mod(iyear,4).eq.0.and.nbgnmn.ge.3) then
            totdays =totdays+1
        endif
!
       do
          maxdays=365
          if(MOD(iyear,4).eq.0) maxdays=366
          if(totdays.le.maxdays) exit
          totdays =totdays-maxdays
          iyear = iyear +1
       enddo
!
        jday=totdays
        MONTH=1
!
        DO im=1,12
          ndays = idays(im)
          if(mod(iyear,4).eq.0.and.im.ge.3) then
               ndays = ndays+1
          ENDIF
          IF(JDAY.GT.ndays) MONTH=im
        ENDDO
!
      end subroutine
END MODULE
