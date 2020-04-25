MODULE ssib4_module_fireimpact 
      use ssib4_module_comsconst
      use ssib4_module_comsveg
      use ssib4_module_trifparms
      use ssib4_module_fireparms

      IMPLICIT NONE

CONTAINS 
!=========================================================================
! FIRE MODULE (CALLED IN SFCTRIF)
!=========================================================================
  SUBROUTINE FIREIMP(TIMESTEP,ZLAT,DAY_TRIF,L_TRIF,                     &
             FRAC_TRIF,BAF_AVG1,BAF_PEAT,                               &
             LEAF1,PHEN,WOOD1,ROOT1,CS,                                 & 
             BURN_AVG,BURN_PEAT,BURN_FIRE,                              &
             LEAF_LOSS_DR,WOOD_LOSS_DR,ROOT_LOSS_DR,                    &
             EM_CO2,EM_CO,EM_CH4,EM_NMHC,EM_H2,EM_NOX,EM_N2O,           &
             EM_PM25,EM_TPM,EM_TC,EM_OC,EM_BC)

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: FireImp changed by Huilin. Huang based on (Li et al., 2012 Biogeosciences)
!
! !INTERFACE:
!
! !DESCRIPTION:


! !called from sfctrif:   
    implicit none
    real :: timestep                      !in model timestep (s).  
    real :: zlat                          !in grid latitude (degree) 
    logical:: l_trif                      !in .t. if vegetation to be updated.  
    integer :: day_trif                   !in number of days between triffid. 
    real,dimension(ntypes)::  frac_trif   !in vegetation coverage
    real :: baf_avg1                      !in grib-averaged burned area fraction (/s)
    real :: baf_peat                      !in grib-averaged peat burned area fraction (/s)
    real,dimension(npft) :: leaf1         !in leaf carbon (kg C/m2) 
    real,dimension(npft) :: phen          !leaf*phen is the actual leaf carbon (kg C/m2)
    real,dimension(npft) :: wood1         !in wood carbon (kg C/m2)
    real,dimension(npft) :: root1         !in root carbon (kg C/m2)
    real :: cs                  !in soil carbon (kg C/m2)

! !inout to sfctrif:   
    real,dimension(npft) :: leaf_loss_dr  !inout accumulated leaf loss due to fire (kg C/m2) 
    real,dimension(npft) :: wood_loss_dr  !inout accumulated wood loss due to fire (kg C/m2)
    real,dimension(npft) :: root_loss_dr  !inout accumulated root loss due to fire (kg C/m2)

! !output to sfctrif:   
    real :: burn_avg                      !out carbon emission due to fire (kg C/m2)
    real :: burn_peat                     !out carbon emission due to peat fire (kg C/m2)
    real,dimension(npft) :: burn_fire     !out carbon emission by different types due to fire
    real,dimension(npft) :: abgl_burn     !out carbon emission by different types due to fire
    real :: em_co2, em_co, em_ch4, em_nmhc, em_h2, em_nox
    real :: em_n2o,em_pm25,em_tpm,em_tc,em_oc,em_bc

! !local variables     
    real,dimension(npft) :: baf_avg       !grib-averaged burned area fraction per sec (%/s) 
    real,dimension(npft) :: leaf1_fire    !leaf update due to fire (kg C/m2) 
    real,dimension(npft) :: wood1_fire    !wood update due to fire (kg C/m2) 
    real,dimension(npft) :: root1_fire    !root update due to fire (kg C/m2)
    real,dimension(npft) :: leaf_loss     !leaf loss due to fire (kg C/m2/s)
    real,dimension(npft) :: wood_loss     !wood loss due to fire (kg C/m2/s)
    real,dimension(npft) :: root_loss     !root loss due to fire (kg C/m2/s)
    real,dimension(npft) :: leaf_burn     !leaf burn due to fire (kg C/m2/s)
    real,dimension(npft) :: wood_burn     !wood burn due to fire (kg C/m2/s)
    real,dimension(npft) :: root_burn     !root burn due to fire (kg C/m2/s)
    real,dimension(npft) :: leaf_mort     !leaf mortality due to fire (kg C/m2/s)
    real,dimension(npft) :: wood_mort     !wood mortality due to fire (kg C/m2/s)
    real,dimension(npft) :: root_mort     !root mortality due to fire (kg C/m2/s) 
    real,dimension(npft) :: leaf_min      !leaf minimum from lai_min (kg C/m2) 
    real,dimension(npft) :: wood_min      !wood minimum from lai_min (kg C/m2) 
    real,dimension(npft) :: root_min      !root minimum from lai_min (kg C/m2)
    real,dimension(npft) :: em_co2_trif   !emission of co2 
    real,dimension(npft) :: em_co_trif    !emission of co 
    real,dimension(npft) :: em_ch4_trif   !emission of ch4 
    real,dimension(npft) :: em_nmhc_trif  !emission of nmhc
    real,dimension(npft) :: em_h2_trif    !emission of h2 
    real,dimension(npft) :: em_nox_trif   !emission of nox 
    real,dimension(npft) :: em_n2o_trif   !emission of n2o 

    real,dimension(npft) :: em_pm25_trif  !emission of pm25
    real,dimension(npft) :: em_tpm_trif   !emission of tpm 
    real,dimension(npft) :: em_tc_trif    !emission of tc 
    real,dimension(npft) :: em_oc_trif    !emission of oc  
    real,dimension(npft) :: em_bc_trif    !emission of bc  

    real :: ftime         !Weighting factor for accumulations
    real :: sec_year      !Weighting factor for accumulations
    integer :: npft1

    real, parameter :: day_year = 365.0        ! number of days in a year (days).
    real, parameter :: sec_day  = 86400.       ! number of seconds in a day (s).
    real, parameter :: borealat=40.0           ! latitude threshold for boreal peat fire     
                                             
    ftime=timestep/real(sec_day*day_trif)  !factors for accummulation variables 
    sec_year = sec_day*day_year

!-------------------------------------------------------------------
!       This is the start of fire impact program
!-------------------------------------------------------------------

    burn_avg                    = 0.0 
    burn_peat                   = 0.0 
    if(zlat.lt.borealat) then
       burn_peat                = max(0.0,cs)*baf_peat*6.0/33.9
    else
       burn_peat                = baf_peat*2.2
    end if
    do npft1=1,npft  
        leaf_loss(npft1)        = 0.0
        wood_loss(npft1)        = 0.0
        root_loss(npft1)        = 0.0
        burn_fire(npft1)        = 0.0
        baf_avg(npft1)          = 0.0
        em_co2_trif(npft1)      = 0.0 
        em_co_trif(npft1)       = 0.0 
        em_ch4_trif(npft1)      = 0.0 
        em_nmhc_trif(npft1)     = 0.0 
        em_h2_trif(npft1)       = 0.0 
        em_nox_trif(npft1)      = 0.0 
        em_n2o_trif(npft1)      = 0.0 
        em_pm25_trif(npft1)     = 0.0 
        em_tpm_trif(npft1)      = 0.0 
        em_tc_trif(npft1)       = 0.0 
        em_oc_trif(npft1)       = 0.0 
        em_bc_trif(npft1)       = 0.0 
        em_co2                  = 0.0
        em_co                   = 0.0
        em_ch4                  = 0.0
        em_nmhc                 = 0.0
        em_h2                   = 0.0
        em_nox                  = 0.0
        em_n2o                  = 0.0
        em_pm25                 = 0.0
        em_tpm                  = 0.0
        em_tc                   = 0.0
        em_oc                   = 0.0
        em_bc                   = 0.0
        leaf1_fire(npft1)       = leaf1(npft1)-leaf_loss_dr(npft1)
        wood1_fire(npft1)       = wood1(npft1)-wood_loss_dr(npft1)
        root1_fire(npft1)       = root1(npft1)-root_loss_dr(npft1)

        wood_min(npft1)         = a_wl(npft1)*lai_min(npft1)**b_wl(npft1)
        leaf_min(npft1)         = sigl(npft1)*lai_min(npft1) 
        root_min(npft1)         = sigl(npft1)*lai_min(npft1) 
!-------------------------------------------------------------------- 
        baf_avg(npft1)          = min(baf_avg1,1.0)     
        leaf_burn(npft1)        = leaf1_fire(npft1)*cc_leaf(npft1)*baf_avg(npft1)
        wood_burn(npft1)        = wood1_fire(npft1)*cc_wood(npft1)*baf_avg(npft1)
        root_burn(npft1)        = root1_fire(npft1)*cc_root(npft1)*baf_avg(npft1)
        burn_fire(npft1)        = leaf_burn(npft1)*phen(npft1)+wood_burn(npft1)+root_burn(npft1)

! update the pft fration 
        leaf_mort(npft1)        = leaf1_fire(npft1)*(1-cc_leaf(npft1))*mf_leaf(npft1)*baf_avg(npft1)
        wood_mort(npft1)        = wood1_fire(npft1)*(1-cc_wood(npft1))*mf_wood(npft1)*baf_avg(npft1)
        root_mort(npft1)        = root1_fire(npft1)*(1-cc_root(npft1))*mf_root(npft1)*baf_avg(npft1)

        leaf_loss(npft1)        = leaf_burn(npft1) + leaf_mort(npft1)
        wood_loss(npft1)        = wood_burn(npft1) + wood_mort(npft1)
        root_loss(npft1)        = root_burn(npft1) + root_mort(npft1)
     
        em_co2_trif(npft1)      = (burn_fire(npft1)+abgl_burn(npft1))*efg_co2(npft1)/0.45
        em_co_trif(npft1)       = (burn_fire(npft1)+abgl_burn(npft1))*efg_co(npft1)/0.45
        em_ch4_trif(npft1)      = (burn_fire(npft1)+abgl_burn(npft1))*efg_ch4(npft1)/0.45
        em_nmhc_trif(npft1)     = (burn_fire(npft1)+abgl_burn(npft1))*efg_nmhc(npft1)/0.45
        em_h2_trif(npft1)       = (burn_fire(npft1)+abgl_burn(npft1))*efg_h2(npft1)/0.45
        em_nox_trif(npft1)      = (burn_fire(npft1)+abgl_burn(npft1))*efg_nox(npft1)/0.45
        em_n2o_trif(npft1)      = (burn_fire(npft1)+abgl_burn(npft1))*efg_n2o(npft1)/0.45
        em_pm25_trif(npft1)     = (burn_fire(npft1)+abgl_burn(npft1))*efa_pm25(npft1)/0.45
        em_tpm_trif(npft1)      = (burn_fire(npft1)+abgl_burn(npft1))*efa_tpm(npft1)/0.45
        em_tc_trif(npft1)       = (burn_fire(npft1)+abgl_burn(npft1))*efa_tc(npft1)/0.45
        em_oc_trif(npft1)       = (burn_fire(npft1)+abgl_burn(npft1))*efa_oc(npft1)/0.45
        em_bc_trif(npft1)       = (burn_fire(npft1)+abgl_burn(npft1))*efa_bc(npft1)/0.45    
    end do

    do npft1=1,npft  
        leaf_loss_dr(npft1)     = leaf_loss_dr(npft1)+leaf_loss(npft1)*timestep
        wood_loss_dr(npft1)     = wood_loss_dr(npft1)+wood_loss(npft1)*timestep
        root_loss_dr(npft1)     = root_loss_dr(npft1)+root_loss(npft1)*timestep
    end do

    do npft1=1,npft  
        burn_avg                = burn_avg+(burn_fire(npft1)+abgl_burn(npft1))*frac_trif(npft1)
        em_co2                  = em_co2+em_co2_trif(npft1)*frac_trif(npft1)
        em_co                   = em_co+em_co_trif(npft1)*frac_trif(npft1)
        em_ch4                  = em_ch4+em_ch4_trif(npft1)*frac_trif(npft1)
        em_nmhc                 = em_nmhc+em_nmhc_trif(npft1)*frac_trif(npft1)
        em_h2                   = em_h2+em_h2_trif(npft1)*frac_trif(npft1)
        em_nox                  = em_nox+em_nox_trif(npft1)*frac_trif(npft1)
        em_n2o                  = em_n2o+em_n2o_trif(npft1)*frac_trif(npft1)
        em_pm25                 = em_pm25+em_pm25_trif(npft1)*frac_trif(npft1)
        em_tpm                  = em_tpm+em_tpm_trif(npft1)*frac_trif(npft1)
        em_tc                   = em_tc+em_tc_trif(npft1)*frac_trif(npft1) 
        em_oc                   = em_oc+em_oc_trif(npft1)*frac_trif(npft1)
        em_bc                   = em_bc+em_bc_trif(npft1)*frac_trif(npft1)
    end do

  END SUBROUTINE 
END MODULE
