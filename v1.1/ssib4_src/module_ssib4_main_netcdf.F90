!     THE MAIN PROGRAM TO TEST SSIB
      program main
!
      use ssib4_module_comsconst
      use ssib4_module_comsveg
      use ssib4_module_trifparms
      use ssib4_module_ssib4sfc
      use module_records
      use module_writenc
      use module_time
      use module_ssibio
      use module_radc2
      use module_read1deg_clim
!
      IMPLICIT NONE
!-----------------------------------------------------------------------
      real:: sunang,tdayhr,zlon,zlat,xlat
      real,dimension(lonf2s:lonf2e,latg2s:latg2e) ::  &
          SWDOWN,RLWDOWN,PPL,PPC,TM,UMM,VMM,QM,PSURF, &
          LIGHT,PREC60,PREC10,PRECD,DP,AGR,GDP,PEATF
      real,dimension(lonf2s:lonf2e,latg2s:latg2e) ::  &
          xvisb,xvisd,xnirb,xnird
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)   :: ta
      real,dimension(lonf2s:lonf2e,latg2s:latg2e,2) :: RADT,capac
      real,dimension(lonf2s:lonf2e,latg2s:latg2e,3) :: WWW
      real,dimension(lonf2s:lonf2e,latg2s:latg2e) ::  &
         ETMASS,EVAPSOIL,EVAPWC,EVAPDC,               &
         EVAPGX,EVAPSN,GHTFLX,USTAR,                  &
         SALB11,SALB12,SALB21,SALB22,                 &
         BEDO,TGEFF,TGS,tc,TD,                        &
         rstfac1,rstfac2,rstfac3,rstfac4,xft1,xft2
      real,dimension(lonf2s:lonf2e,latg2s:latg2e) ::  &
         ra,rb,rd,rsoil,shf,chf,ROFF,hlflx,hsflx,     &
         ulwsf1,evap,q2m,umom,vmom,plantr,sheleg,cm,  &
         ch,fm,fh,ZLT1,soilm,fsdown,fldown,           &
         fsup,flup,xhsflx,xhlflx,                     &
         sibnet,anet,fcs,fca,cacn,                    &!Huilin add rh to output
         cs,ci,gb,gs,rib
      real,dimension(lonf2s:lonf2e,latg2s:latg2e) ::  &
         xwcl,xwel,xwsl,xwcs,xwes,xwss
!
! output variables
              
      real, dimension(lonf2s:lonf2e,latg2s:latg2e,2)      ::   &
          qtvcover, qtzlt,qgreenfrac
      real,dimension(lonf2s:lonf2e,latg2s:latg2e,NTYPET)  ::   &
          QROOT1,QWOOD1,QLEAF1,QDCVEG1,                        &
          QNPP,QGPP,QRESP_P,QLIT_C,xphen,qburn_fire,           &    !Huilin add qburn_fire Mar.2019 
          xleaf_loss_dr,xwood_loss_dr,xroot_loss_dr
          
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)         ::   &
          QRESP_S, QNEP,QLIT_C_T,QHS                          
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)         ::   &
          QRH,QHR,QPREC60_AVG,QPEATF                                !Huilin add fire variable
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)         ::   &
          QAGB_AVG,QFEO_AVG,QFES_AVG,QFE1_AVG,                 &
          QFE2_AVG,QFD_AVG,QFCLI,QFSAT,                        &
          QLH,QIGN,QFM_AVG,QFT_AVG,QFB_AVG,QFM_FIRE,QAGR,      &
          QBAF_AVG,QNFIRE_AVG,QBURN_AVG,QBURN_PEAT,QBAF_PEAT,  &
          QBAF_AVG60                                                !Huilin add fire variable
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)         ::   &
          qem_co2,qem_co,qem_ch4,qem_nhmc,qem_h2,qem_nox,qem_n2o,&  !Huilin add fire output  May.2019
          qem_pm25,qem_tpm,qem_tc,qem_oc,qem_bc
      real, dimension(lonf2s:lonf2e,latg2s:latg2e,NTYPES) ::   &
          xht_trif,xgreen_trif,xlai_trif,xsai_trif,            &
          xgreen_ssib,xlaid_trif,xveg_fracn
      real, dimension(lonf2s:lonf2e,latg2s:latg2e,NTYPES) ::   &
          QRESP_P_M,QRESP_P_G
!-----------------------------------------------------------------------
      INTEGER:: irestart
      REAL:: RADFRAC11,RADFRAC12,RADFRAC21,RADFRAC22
      real:: CLOUD,DIFRAT,VNRAT,RHOAIR,SIG,TDD
!
!-----------------------------------------------------------------------
! Input and loop variables
!-----------------------------------------------------------------------
      real:: day,thour,dtt
      INTEGER:: i,j,k,it,ITIME,irec,krec,kyr
      INTEGER:: IYEAR,IMONTH,IDAY,jday,nthour,ndom
      integer:: ibgnyr,ibgnmn,ibgndy
      real:: fbgnhr
      integer:: krecst,kyrst,ierr,nrec,datalen

      real,dimension(lonf2s:lonf2e,latg2s:latg2e):: &
          zwind,pmr
      character(len=255):: outdir,expnum,strfile1,strfile2
      character(len=4):: skyr
      NAMELIST /nctl/expnum,outdir
!
      print*,'starting ...'
      print*,outdir
!
      rewind(1)
      read(1,nml=nctl)
      close(1)
!-----------------------------------------------------------------------
!
!initial output variables
      call initial
! ............. ssib constant start
      call vegin
! ............. ssib constant end
!
      call param_ini(ssibios,dtt,  &
          ibgnyr,ibgnmn,ibgndy,fbgnhr,xlongit,xlatit)

      thour=0.0
      iyear = ibgnyr
      iday  = ibgndy
      imonth= ibgnmn
      tdayhr= fbgnhr
!
      itime = 1 
!
!Judge restarting and loading previous data
!
! irestart = 0 for starting the model.
! irestart = 1 for starting the model and read ssib.inp.
! irestart = 2 for restart
      irestart  = 0        !Huilin revise to 2
!
      call ssibin(ssibios,irestart,outdir,dtt, &
            iyear,imonth,iday,tdayhr,itime,thour)   
!
      call datetojday(iyear,imonth,iday,jday,ierr)
      day = jday
      kyrst=iyear
      krecst=int(((day-1)*24+tdayhr)*3600/dtt)+1
!
      datalen =nstepofyear(iyear,dtt)
      if(krecst.gt.datalen) then
          krecst=krecst-datalen
          kyrst =kyrst +1
      endif
!
      call seton_reg_global(1,"monthly average for each timestep")
      call seton_reg_global(2,"monthly average")

!=========Reading forcing data ===================
!
! Ending year is 2008.
!
!yliu add test output
!     open(330, file=outdir(1:len_trim(outdir))//'diag.asc')
       print*, 'starting doloop',kyrst
!      do 2000 kyr = kyrst,2048
      do 2000 kyr = kyrst,2014
!yliu add to keep out ssib.inp
      if (mod(kyr,1) .eq. 0.) then
!     if (mod(kyr,10) .eq. 0.) then
        call sfromi(skyr,4,kyr)
        strfile1 = outdir(1:len_trim(outdir))//'/ssib.inp'
        strfile2 = outdir(1:len_trim(outdir))//'/ssib.inp_'//skyr
        call system('cp '//trim(strfile1)//' '//trim(strfile2))
      end if

!READ IN THE ATMOSPHERIC data
      datalen =nstepofyear(kyr,dtt)
      do 999 krec = krecst,datalen
!
      call readforcing(dtt,iyear,krec,krecst,datalen,   &
        swdown,rlwdown,psurf,pmr,ppl,ppc,tm,qm,umm,vmm, &
        zwind,light,prec60,prec10,precd,dp,agr,gdp,peatf)
!Huilin add light prec60 prec10 precd dp agr gdp peatf
!
      if(krec.eq.datalen) krecst = 1 
!
      JDAY =DAY
!
      do 1000 i = LONF2S, LONF2E
      do 1001 j = LATG2S, LATG2E
!
        zlon = xlongit(i)
        zlat = xlatit(j)
        CALL RADC2_cosz(1,SUNANG,IYEAR,DAY,tdayhr,DTT,zlon,zlat)
!
        CLOUD = (1160.*SUNANG-SWDOWN(i,j))/(963.*SUNANG)
        CLOUD = AMAX1(CLOUD,0.)
        CLOUD = AMIN1(CLOUD,1.)
        CLOUD = AMAX1(0.58,CLOUD)
!
!zzq sept26,2013   DIFRAT = 0.0604 / ( SUNANG-0.0223 )  +  0.0683
!zzq sept26,2013   IF ( DIFRAT .LT. 0. ) DIFRAT = 0.
!zzq sept26,2013   IF ( DIFRAT .GT. 1. ) DIFRAT = 1.
!since DIFRAT is negative for all SUNANG less than 0.0223
!
        IF(SUNANG.LE.0.0223) THEN
           DIFRAT =0.0
        ELSE
           DIFRAT = 0.0604 / ( SUNANG-0.0223 )  +  0.0683
           IF ( DIFRAT .GT. 1. ) DIFRAT = 1.
        ENDIF
!
        DIFRAT = DIFRAT  +  ( 1. - DIFRAT ) * CLOUD
        VNRAT = ( 580. - CLOUD*464. ) / ( ( 580. - CLOUD*499. )   &
                 +  ( 580. - CLOUD*464. ) )

        RADFRAC11 = (1.-DIFRAT)*VNRAT
        RADFRAC12 = DIFRAT*VNRAT
        RADFRAC21 = (1.-DIFRAT)*(1.-VNRAT)
        RADFRAC22 = DIFRAT*(1.-VNRAT)

        xvisb(i,j) = RADFRAC11*SWDOWN(i,j)
        xvisd(i,j) = RADFRAC12*SWDOWN(i,j)
        xnirb(i,j) = RADFRAC21*SWDOWN(i,j)
        xnird(i,j) = RADFRAC22*SWDOWN(i,j)
!
        landmsk(i,j) = 0
        if(ssibios(i,j) % xvtype.gt.0) then
          landmsk(i,j)=1
!
        RHOAIR = (PSURF(I,J)+PMR(I,J))/GASR/(TA(I,J)+TM(I,J))
        SIG = 1.0
        TDD = TM(I,J)
!
        CALL ssib4sfc(                                               & 
!input:
        DTT,SUNANG,                                                  &
        xvisb(i,j), xvisd(i,j), xnirb(i,j), xnird(i,j),              &
        RLWDOWN(i,j), PPC(i,j),PPL(i,j),                             &
        SIG, TM(i,j), UMM(i,j), VMM(i,j), QM(i,j),PSURF(i,j),        &
        ZWIND(i,j),  RHOAIR, TDD, KREC,                              &
        LIGHT(i,j), PREC60(i,j), PREC10(i,j), PRECD(i,j),            &
        DP(i,j), AGR(i,j), GDP(i,j), PEATF(i,j),                     &  !Huilin add Light, DP 
        zlon,zlat,ibgnyr,ibgnmn,ibgndy,fbgnhr,thour,itime,           &
!input/output:
        ssibios(i,j),                                                & 
!output:
        xhsflx(i,j),xhlflx(i,j),GHTFLX(i,j),rib(i,j),                &
        EVAPSOIL(i,j),EVAPWC(i,j),EVAPDC(i,j),                       &
        EVAPSN(i,j), EVAPGX(i,j),                                    &
        SALB11(i,j),SALB12(i,j),SALB21(i,j),SALB22(i,j),             &
        fsdown(i,j),fldown(i,j),fsup(i,j), flup(i,j),                &
        TA(i,j),TGEFF(i,j), ROFF(i,j),q2m(i,j),USTAR(i,j),           &!Huilin Hunang
        soilm(i,j),plantr(i,j),                                      &
        cm(i,j),ch(i,j),fm(i,j),fh(i,j),                             &
!other output:
        dataout(i,j)                                                 &
         )
!
        ETMASS(i,j) = dataout(i,j)%etmass 
        RADT(i,j,1) = dataout(i,j)%radt1
        RADT(i,j,2) = dataout(i,j)%radt2
        BEDO(i,j)   = dataout(i,j)%bedo
        TGS(i,j)    = dataout(i,j)%tgs
        TD(i,j)     = dataout(i,j)%td
        tc(i,j)     = dataout(i,j)%tc
        WWW(i,j,1)  = dataout(i,j)%www1
        www(i,j,2)  = dataout(i,j)%www2
        www(i,j,3)  = dataout(i,j)%www3
        capac(i,j,1)= dataout(i,j)%capac1
        capac(i,j,2)= dataout(i,j)%capac2
        rstfac1(i,j)= dataout(i,j)%rstfac1
        rstfac2(i,j)= dataout(i,j)%rstfac2
        rstfac3(i,j)= dataout(i,j)%rstfac3
        rstfac4(i,j)= dataout(i,j)%rstfac4
        ra(i,j)     = dataout(i,j)%ra
        rb(i,j)     = dataout(i,j)%rb
        rd(i,j)     = dataout(i,j)%rd
        rsoil(i,j)  = dataout(i,j)%rsoil
        shf(i,j)    = dataout(i,j)%shf
        chf(i,j)    = dataout(i,j)%chf
        hlflx(i,j)  = dataout(i,j)%hlflx
        hsflx(i,j)  = dataout(i,j)%hsflx
        ulwsf1(i,j) = dataout(i,j)%ulwsf1
        evap(i,j)   = dataout(i,j)%evap
        umom(i,j)   = dataout(i,j)%umom
        vmom(i,j)   = dataout(i,j)%vmom
        sheleg(i,j) = dataout(i,j)%sheleg
        ZLT1(i,j)   = dataout(i,j)%zlt1
        anet(i,j)   = dataout(i,j)%aveanet
        fcs(i,j)    = dataout(i,j)%fcs
        fca(i,j)    = dataout(i,j)%fca
        cacn(i,j)   = dataout(i,j)%cacn
        cs(i,j)     = dataout(i,j)%cs
        ci(i,j)     = dataout(i,j)%ci
        gb(i,j)     = dataout(i,j)%gb
        gs(i,j)     = dataout(i,j)%gs
        xwcl(i,j)   = dataout(i,j)%xwcl
        xwel(i,j)   = dataout(i,j)%xwel
        xwsl(i,j)   = dataout(i,j)%xwsl
        xwcs(i,j)   = dataout(i,j)%xwcs
        xwes(i,j)   = dataout(i,j)%xwes
        xwss(i,j)   = dataout(i,j)%xwss
        xft1(i,j)   = dataout(i,j)%aveft1
        xft2(i,j)   = dataout(i,j)%aveft2
        fca(i,j)    = fca(i,j)*1e6
        fcs(i,j)    = fcs(i,j)*1e6
        anet(i,j)   = anet(i,j)*1e6
        qem_co2(i,j) = dataout(i,j)%qem_co2
        qem_co(i,j)  = dataout(i,j)%qem_co
        qem_ch4(i,j) = dataout(i,j)%qem_ch4
        qem_nhmc(i,j)= dataout(i,j)%qem_nhmc
        qem_h2(i,j)  = dataout(i,j)%qem_h2
        qem_nox(i,j) = dataout(i,j)%qem_nox
        qem_n2o(i,j) = dataout(i,j)%qem_n2o
        qem_pm25(i,j)= dataout(i,j)%qem_pm25
        qem_tpm(i,j) = dataout(i,j)%qem_tpm
        qem_tc(i,j)  = dataout(i,j)%qem_tc
        qem_oc(i,j)  = dataout(i,j)%qem_oc
        qem_bc(i,j)  = dataout(i,j)%qem_bc
       SIBNET(I,J)  = RADT(I,J,1)+RADT(I,J,2)
!
! for output
!
         do k=1,2
           qtvcover(i,j,k)   = dataout(i,j) % qtvcover(k)
           qtzlt(i,j,k)      = dataout(i,j) % qtzlt(k)
           qgreenfrac(i,j,k) = dataout(i,j) % qgreenfrac(k)
         enddo
!
         do k=1,NTYPET
           QROOT1(i,j,k)  = dataout(i,j) % QROOT1(k)
           QWOOD1(i,j,k)  = dataout(i,j) % QWOOD1(k)
           QLEAF1(i,j,k)  = dataout(i,j) % QLEAF1(k)
           QDCVEG1(i,j,k) = dataout(i,j) % QDCVEG1(k)
           QNPP(i,j,k)    = dataout(i,j) % QNPP(k)
           QGPP(i,j,k)    = dataout(i,j) % QGPP(k)
           QRESP_P(i,j,k) = dataout(i,j) % QRESP_P(k) 
           QLIT_C(i,j,k)  = dataout(i,j) % QLIT_C(k)
           QBURN_FIRE(i,j,k) = dataout(i,j) % QBURN_FIRE(k)           !Huilin add fire variable Mar.2019
           xphen(i,j,k)   = ssibios(i,j) % xphen(k) 
           xleaf_loss_dr(i,j,k) = ssibios(i,j) % xleaf_loss_dr(k)     !Huilin add fire variable
           xwood_loss_dr(i,j,k) = ssibios(i,j) % xwood_loss_dr(k)     !Huilin add fire variable
           xroot_loss_dr(i,j,k) = ssibios(i,j) % xroot_loss_dr(k)     !Huilin add fire variable
         enddo
!
          QRESP_S(i,j) = dataout(i,j) % QRESP_S
          QNEP(i,j)    = dataout(i,j) % QNEP
          QLIT_C_T(i,j)= dataout(i,j) % QLIT_C_T
!         QHS(i,j)     = dataout(i,j) % QHS
          QHS(i,j)     = dataout(i,j) % QRH                           !Huilin add fire variable
          QHR(i,j)     = dataout(i,j) % QHR                           !Huilin add fire variable
          QAGB_AVG(i,j)= dataout(i,j) % QAGB_AVG                      !Huilin add fire variable
          QFEO_AVG(i,j)= dataout(i,j) % QFEO_AVG                      !Huilin add fire variable
          QFES_AVG(i,j)= dataout(i,j) % QFES_AVG                      !Huilin add fire variable
          QFE1_AVG(i,j)= dataout(i,j) % QFE1_AVG                      !Huilin add fire variable
          QFE2_AVG(i,j)= dataout(i,j) % QFE2_AVG                      !Huilin add fire variable
          QFD_AVG(i,j) = dataout(i,j) % QFD_AVG                       !Huilin add fire variable
          QLH(i,j)     = dataout(i,j) % QLH       
          QIGN(i,j)    = dataout(i,j) % QIGN      
          QFM_AVG(i,j) = dataout(i,j) % QFM_AVG   
          QFT_AVG(i,j) = dataout(i,j) % QFT_AVG   
          QFB_AVG(i,j) = dataout(i,j) % QFB_AVG   
          QFM_FIRE(i,j)= dataout(i,j) % QFM_FIRE       
          QAGR(i,j)    = dataout(i,j) % QAGR      
          QFCLI(i,j)   = dataout(i,j) % QFCLI                         !Huilin add fire variable
          QFSAT(i,j)   = dataout(i,j) % QFSAT                         !Huilin add fire variable
          QPEATF(i,j)  = dataout(i,j) % QPEATF                        !Huilin add fire variable
          QBAF_AVG(i,j)= dataout(i,j) % QBAF_AVG                      !Huilin add fire variable
          QBAF_AVG60(i,j)= dataout(i,j) % QBAF_AVG60                  !Huilin add fire variable
          QBAF_PEAT(i,j)= dataout(i,j) % QBAF_PEAT
          QBURN_AVG(i,j)= dataout(i,j) % QBURN_AVG                    !Huilin add fire variable Mar.2019
          QBURN_PEAT(i,j)= dataout(i,j) % QBURN_PEAT                  !Huilin add fire variable Mar.2019
          QNFIRE_AVG(i,j)= dataout(i,j) % QNFIRE_AVG                  !Huilin add fire variable
          QPREC60_AVG(i,j)= dataout(i,j) % QPREC60_AVG                !Huilin add fire variable
!
         do k=1,NTYPES
           xht_trif(i,j,k)    = ssibios(i,j) % xht_trif(k)
           xgreen_trif(i,j,k) = ssibios(i,j) % xgreen_trif(k)
           xlai_trif(i,j,k)   = ssibios(i,j) % xlai_trif(k)
           xsai_trif(i,j,k)   = ssibios(i,j) % xsai_trif(k)
           xgreen_ssib(i,j,k) = ssibios(i,j) % xgreen_ssib(k)
           xlaid_trif(i,j,k)  = ssibios(i,j) % xlaid_trif(k)
           xveg_fracn(i,j,k)  = ssibios(i,j) % xveg_fracn(k) 
         enddo
!
       if(krec.ne.0) then
!zzq           write(1000+krec,*) i,j, fca(i,j)
        endif
!
       endif      !========================
!
 1001  continue
 1000  continue
       ITIME =2
!
       nthour=int(thour)
       call seton_reg_records(1)
       call recsh(1,dtt,nthour,1,tgs,"tg","soil Surface Temp","degC")
       call recsh(1,dtt,nthour,1,td,"td","deep soil Temperature","degC")
       call recsh(1,dtt,nthour,1,tc,"tc","canopy temperature","degC")
       call recsh(1,dtt,nthour,1,tgeff,"tgeff","radiative Temp","degC")
       call recsh(1,dtt,nthour,1,etmass,"etmass", &
                  "eveporation","mm/per time step")
       call recsh(1,dtt,nthour,1,evapdc,"ect","transpriation","j/m2")
       call recsh(1,dtt,nthour,1,evapwc,"eci","evap intercep","j/m2")
       call recsh(1,dtt,nthour,1,evapsoil,"egs","evap soil","j/m2")
       call recsh(1,dtt,nthour,1,evapsn,"egi","evap snow","j/m2")
       call recsh(1,dtt,nthour,1,ra,"ra","aerodynamic res","s/m")
       call recsh(1,dtt,nthour,1,rstfac2,"rstfac2","rstfac2","-")
       call recsh(1,dtt,nthour,1,anet,"anet",  &
                   "net primary product","umol/(m*m*sec)")
       call recsh(1,dtt,nthour,1,fcs,"fcs",    &
                   "Soil surface co2 flux","umol/(m*m*sec)")
       call recsh(1,dtt,nthour,1,fca,"fca",    &
                   "simulated Co2 Flux","umol/(m*m*sec)")
       call recsh(1,dtt,nthour,1,sibnet,"sibnet","net radiation","W/m2")
       call recsh(1,dtt,nthour,1,xhsflx,"hflux", &
                  "sensible heat flux","W/m2")
       call recsh(1,dtt,nthour,1,xhlflx,"siblh", &
                  "latent heat flux","W/m2")
       call recsh(1,dtt,nthour,1,ghtflx,"ghtflx","CHF+SHF","W/m2")
       call recsh(1,dtt,nthour,1,bedo,"bedo","total surface albedo","-")
       call recsh(1,dtt,nthour,1,ustar,"ustar",  &
                 "friction velocity","m/s")
       call recsh(1,dtt,nthour,1,plantr,"rc","canopy res","s/m")
       call recsh(1,dtt,nthour,1,xft1,"ft1","temp. ctrl factor1","-")
       call recsh(1,dtt,nthour,1,xft2,"ft2","temp. ctrl factor2","-")
       call recsh(1,dtt,nthour,1,qhs,"hs","relative humidity","-")
       call recsh(1,dtt,nthour,1,swdown,"swdown","DSW radiation","W/m2")
       call recsh(1,dtt,nthour,1,rlwdown,"rlwdown",  &
                 "DLW radiation","W/m2")
       call recsh(1,dtt,nthour,1,psurf,"psur","surface pressure","bar")
       call recsh(1,dtt,nthour,1,tm,"tm","reference tempearature","k")
       call recsh(1,dtt,nthour,1,qm,"qm","reference moisture","kg/kg")
       call recsh(1,dtt,nthour,1,umm,"um", &
                 "wind at reference height","m/s")
       call recsh(1,dtt,nthour,1,ppl,"ppl", &
                 "large scale precipitation","mm/per time step")
       call setoff_reg_records(1)
       call netcdf_inital_monthlyh(expnum,outdir, &
                 1,dtt,ibgnyr,ibgnmn,ibgndy,fbgnhr,irestart)
       call setoff_reg_global(1)
!
       call seton_reg_records(2)
       call records(2,3,www,"www","Soil wetness in layer","-")
       call records(2,1,tgs,"tg","soil Surface Temp","degC")
       call records(2,1,td,"td","deep soil Temperature","degC")
       call records(2,1,tc,"tc","canopy temperature","degC")
       call records(2,1,tgeff,"tgeff","radiative Temp","degC")
       call records(2,1,etmass,"etmass","eveporation","mm/per tstep")
       call records(2,1,evapdc,"ect","transpriation","j/m2")
       call records(2,1,evapwc,"eci","evap intercep","j/m2")
       call records(2,1,evapsoil,"egs","evap soil","j/m2")
       call records(2,1,evapsn,"egi","evap snow","j/m2")
       call records(2,1,sibnet,"sibnet","net radiation","W/m2")
       call records(2,1,xhsflx,"hflux","sensible heat flux","W/m2")
       call records(2,1,xhlflx,"siblh","latent heat flux","W/m2")
       call records(2,1,ghtflx,"ghtflx","CHF+SHF","W/m2")
       call records(2,1,bedo,"bedo","total surface albedo","-")
       call records(2,1,qhs,"hs","relative humidity","-")
       call records(2,1,swdown,"swdown","DSW radiation","W/m2")
       call records(2,1,rlwdown,"rlwdown","DLW radiation","W/m2")
       call records(2,1,fsup,"swup","USW radiation","W/m2")
       call records(2,1,flup,"lwup","ULW radiation","W/m2")
       call records(2,1,psurf,"psur","surface pressure","bar")
       call records(2,1,tm,"tm","reference tempearature","k")
       call records(2,1,qm,"qm","reference moisture","kg/kg")
       call records(2,1,umm,"um","wind at reference height","m/s")
       call records(2,1,ppl,"ppl","large scale prec","mm/per tstep")
!!     call records(2,1,salb11,"salb11","surface albedo11","-")
!!     call records(2,1,salb12,"salb12","surface albedo12","-")
!!     call records(2,1,salb21,"salb21","surface albedo21","-")
!!     call records(2,1,salb22,"salb22","surface albedo22","-")
       call records(2,1,roff,"roff","Simulated runoff","m")
       call records(2,2,capac,"capac","water landcover","m")
       call records(2,1,rb,"rb","aerodynamic res rb","s/m")
       call records(2,1,rd,"rd","aerodynamic res rd","s/m")
       call records(2,1,ra,"ra","aerodynamic res ra","s/m")
       call records(2,1,plantr,"rc","aerodynamic res rc","s/m")
       call records(2,1,rsoil,"rsoil","soil res","s/m")
       call records(2,1,cs,"cs","carbon at soil","-")
!!     call records(2,1,xwcs,"wcs","saturation limited rate","-")
!!     call records(2,1,xwes,"wes","electron trans.limited rate","-")
!!     call records(2,1,xwss,"wss","The sink limited rate","-")
!!     call records(2,1,xwcl,"wcl","saturation limited rate","-")
!!     call records(2,1,xwel,"wel","electron trans.limited rate","-")
!!     call records(2,1,xwsl,"wsl","The sink limited rate","-")
!yliu09May2017 modified for add type
       call records(2,8,xveg_fracn,"frac",&
            "fract. for plant functype","-")
       call records(2,7,xlai_trif,"lai","Leaf Area Index","-")
!!     call records(2,7,xsai_trif,"sai","Stem Area Index","-")
       call records(2,7,xgreen_ssib,"green","green fraction","-")
       call records(2,7,xht_trif,"vht","veg height","m")
       call records(2,7,xphen,"phen","veg phenology","-")
       call records(2,7,qleaf1,"leaf","leaf1","-")
       call records(2,7,qroot1,"root","rootl","-")
       call records(2,7,qwood1,"wood","wood1","-")
       call records(2,7,QDCVEG1,"DCVEG","DCVEG1","-")
       call records(2,7,qnpp,"npp","npp","-")
       call records(2,7,qgpp,"gpp","gpp","-")
       call records(2,7,QLIT_C,"litc","LIT_C","-")
       call records(2,7,qresp_p,"resp","resp_p","-")
!!
       call records(2,1,qresp_s,"resps","resp_s","-")
!Huilin add output variables
       call records(2,1,rstfac2,"rstfac2","rstfac2","-")
!      call records(2,1,qhr,"hr","top layer soil moisture","-")
!      call records(2,1,qprec60_avg,"prec60","60-day prec","-")
       call records(2,1,qagb_avg,"agb_avg","fuel avalibility","-")
       call records(2,1,qfeo_avg,"feo_avg","human sup on fire occ","-")
       call records(2,1,qfes_avg,"fes_avg","human sup on fire sprd","-")
       call records(2,1,qfe1_avg,"fe1_avg","fe1_avg","-")
       call records(2,1,qfe2_avg,"fe2_avg","fe2_avg","-")
       call records(2,1,qfcli,"fcli","fcli","-")
       call records(2,1,qfsat,"fsat","fsat","-")
       call records(2,1,qpeatf,"peatf","peatf","-")
       call records(2,1,qfd_avg,"fd_avg","fd_avg","-")
       call records(2,1,qlh,"qlh","h ignition","-")
       call records(2,1,qign,"ignition","t ignition","-")
       call records(2,1,qfm_avg,"qfm_avg","sm mos factor","-")
       call records(2,1,qfm_fire,"qfm","mos fac","-")
       call records(2,1,qft_avg,"qft_avg","temp fac","-")
       call records(2,1,qfb_avg,"qfb_avg","fuel fac","-")
       call records(2,1,qagr,"qagr","agr frac","-")
       call records(2,1,qbaf_avg,"baf_avg","burn fraction","-")
       call records(2,1,qbaf_avg60,"baf_avg60","60-day bf","-")
       call records(2,1,qbaf_peat,"baf_peat","burn fraction","-")
       call records(2,1,qburn_avg,"burn_avg","c combustion","-")
       call records(2,1,qburn_peat,"burn_peat","cc from peat","-")
       call records(2,1,qnfire_avg,"nfire_avg","fire count","-")     
       call records(2,7,qburn_fire,"burn_fire","type-depend CC","-")
       call records(2,1,qem_co2,"em_co2","co2 emission","-")
       call records(2,1,qem_co,"em_co","co emission","-")
       call records(2,1,qem_ch4,"em_ch4","ch4 emission","-")
       call records(2,1,qem_nhmc,"em_nhmc","nhmc emission","-")
       call records(2,1,qem_h2,"em_h2","h2 emission","-")
       call records(2,1,qem_nox,"em_nox","nox emission","-")
       call records(2,1,qem_n2o,"em_n2o","n2o emission","-")
       call records(2,1,qem_pm25,"em_pm25","pm25 emission","-")
       call records(2,1,qem_tpm,"em_tpm","tpm emission","-")
       call records(2,1,qem_tc,"em_tc","tc emission","-")
       call records(2,1,qem_oc,"em_oc","oc emission","-")
       call records(2,1,qem_bc,"em_bc","bc emission","-")
!      call records(2,1,xleaf_loss_dr,"leaf_loss", &
!                   "leaf corbon loss accumulated","-")
!      call records(2,1,xwood_loss_dr,"wood_loss", & 
!                   "wood corbon loss accumulated","-")
!      call records(2,1,xroot_loss_dr,"root_loss", & 
!                   "root corbon loss accumulated","-")

       call setoff_reg_records(2)
       call netcdf_inital_monthly(expnum,outdir,2,&
              ibgnyr,ibgnmn,ibgndy,fbgnhr,irestart)
       call setoff_reg_global(2)
!
       print*,'step=',iyear,imonth,iday,krec,int(thour/24.)

       if(nthour.eq.21) then
          write(*,*) 'writing initial step'
           call netcdf_writeouth(1,dtt,1)
       endif
!
       call daynumofmonth(iyear,imonth,ndom)
!
       if((iday.eq.ndom).and.(mod(nthour,24).eq.21)) then
           call ssibout(ssibios,outdir,ibgnyr,ibgnmn,ibgndy,fbgnhr,thour)
           nrec = (iyear-ibgnyr)*12+imonth
!
           call netcdf_writeouth(1,dtt,nrec)
           call netcdf_writeoutm(2,nrec)
           call initial
       endif

       CALL RADC2_time(IYEAR, imonth,DAY,iday,tdayhr,DTT)
       THOUR =THOUR + DTT/3600.
!
  999  continue
 2000  continue
!
       call netcdf_writeouth_close(1,dtt)
       call netcdf_writeoutm_close(2)
!
        write(*,*) 'Have finished'
      END program
