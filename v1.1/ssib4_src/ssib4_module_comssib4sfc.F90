MODULE ssib4_module_comssib4sfc
!
    use ssib4_module_comsveg
!
    IMPLICIT none
!yliu09May2017 addtype    integer, parameter :: NTYPET=6
    integer, parameter :: NTYPET=7
!---------------------------------------------------------------------------
! define the inout variable type
!---------------------------------------------------------------------------
TYPE ssibdataios
    real,dimension(NTYPES):: tc0,tg0,td0,tcm,tgm,tdm,xta
    real,dimension(NTYPES):: w01,w02,w03,wm1,wm2,wm3      ! unit: without multiplied by porosity of soil 
    real,dimension(NTYPES):: cap01,capm1,cap02,capm2
    real,dimension(NTYPES)::                            &
       xveg_fracn,                                      &
       xco2ra,                                          &
       xco2rb,                                          & 
       xzfw,                                            &
       xfiltr,                                          &
       xea
!yliu add for ssib4 initial condition
    real :: xtc0,xtg0,xtd0,xsmc1,xsmc2,xsmc3,xcap1,xcap2
!-----------------------------------------------------------------------
    real,dimension(NTYPES) ::                           &
       xht_trif,                                        &
       xgreen_trif,                                     &
       xgreen_ssib,                                     &
       xlai_trif,                                       &
       xsai_trif,                                       &
       xlaid_trif
!
    real ::                                             &
       xitrifstep,                                      &
       xktrifday,                                       &
       xschsel
!
    real::                                              &
       xresp_s_dr,                                      &
       xveg_frac,                                       &
       xfrac_vs
!
    real:: xco2cs
    real:: xvtype
!Huilin add for fire
    real,dimension(480)    :: xbaf60
    real,dimension(NTYPES) :: xleaf_loss_dr
    real,dimension(NTYPES) :: xwood_loss_dr
    real,dimension(NTYPES) :: xroot_loss_dr
!
    real,dimension(NTYPET)  ::                          &
       xg_leaf_day,                                     &
       xg_leaf_dr,                                      &
       xgpp_dr,                                         &
       xnpp_dr,                                         &
       xresp_w_dr,                                      &
!Huilin add for fire
       xleaf1, xwood1, xroot1
!
    real,dimension(NTYPET) ::                           &
       xphen,                                           &
       xfsmc,                                           &
       xfsmc_ave
!
END TYPE
!
!---------------------------------------------------------------------------
! define the variable type for output
!---------------------------------------------------------------------------
!
TYPE ssibdataout
! ********
!for triffid output sept29,2013
!
    real,dimension(NTYPET)::                            &
       QROOT1,QWOOD1,QLEAF1,QDCVEG1,                    &
       QRESP_P,QLIT_C
!yliu
    real,dimension(NTYPES)::QNPP,QGPP
!
    real::                                              & 
       QRESP_S, QNEP,QLIT_C_T,QHS,QRH,QHR               !Huilin add QRH, QHR for fire
!Huilin add for fire
    real:: QAGB_AVG, QFEO_AVG, QFES_AVG                 
    real:: QFE1_AVG, QFE2_AVG, QFD_AVG
    real:: QPREC, QPEATF, QPREC60_AVG, QBAF_AVG60
    real:: QLH,QIGN,QFM_AVG,QFT_AVG,QFB_AVG,QFM_FIRE,QAGR
    real:: QFSAT, QFCLI, QBAF_PEAT                
    real:: QBAF_AVG, QNFIRE_AVG, QBURN_AVG, QBURN_PEAT  !Huilin add QBURN_AVG Mar. 2019
    real,dimension(NTYPES):: QBURN_FIRE                 !Huilin add QBURN_FIRE Mar. 2019
    real:: QEM_CO2,QEM_CO,QEM_CH4,QEM_NHMC,QEM_H2,QEM_NOX,QEM_N2O
    real:: QEM_PM25,QEM_TPM,QEM_TC,QEM_OC,QEM_BC 
!yliu09May2017 --
    real,dimension(NTYPES)::                            &
       QRESP_P_M,QRESP_P_G,QANET,QRDC
!!!--
!
    real,dimension(2)::                                 &
       qtvcover,                                        &
       qtzlt,                                           &
       qgreenfrac
!
    real::                                              & 
       aveanet,                                         &
       fcs,fca,                                         &
       cacn,cs,ci,                                      &
       gb,gs,                                           &
       xwcl,xwel,xwsl,xwcs,xwes,xwss,                   &
       aveft1,aveft2 
!
!output for model use
    real::                                              &
       evap,                                            &
       canopy,                                          & 
       sheleg
!
    real:: www1,www2,www3
    real:: capac1,capac2
    real::                                              &
       tgs,tc,td,                                       &
       hlflx,hsflx,                                     &
       ulwsf1,                                          & !upward longwave radiation
       radt1,radt2,etmass,                              &
       bedo,shf,chf,                                    &
       rstfac1,rstfac2,rstfac3,rstfac4,                 &
       ra,rb,rd,rsoil,                                  &
       umom,vmom,                                       &
       zlt1 
!
END TYPE

END MODULE
