MODULE ssib4_module_fireseason 
      use ssib4_module_comsconst
      use ssib4_module_comsveg
      use ssib4_module_trifparms
      use ssib4_module_fireparms

      IMPLICIT NONE

CONTAINS 
!=========================================================================
! FIRE MODULE (CALLED IN SFCTRIF)
!=========================================================================
  SUBROUTINE FIRESEASON (TIMESTEP,ZLAT,VFRAC,DP,LIGHT,AGR,  &
             GDP,PEATF,ZWIND,TM,EM,TG1,WWW1N,WWW2N,PREC60,  &
             FEO_AVG,FES_AVG,FE1_AVG,FE2_AVG,FD_AVG,        &
             LH,IGN,FM1_AVG,FT_AVG,FM2_AVG,                 &
             FRAC_TRIF,LEAF1,PHEN,WOOD1,                    &
             LEAF_LOSS_DR,WOOD_LOSS_DR,                     & 
             RH,FB_AVG,AGB_AVG,FCLI,FSAT,                   &
             BAF_AVG,BAF_PEAT,NFIRE_AVG)

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: FireSeason changed by Huilin Huang based on (Li et al., 2012 Biogeosciences)
!
! !DESCRIPTION:
! fire code was called 3-hourly 
!
! !USES:
! !nfire_avg: grib-averaged the number of fires count/km^2/s
! !baf_avg: grib-averaged burned area fraction /s

! !CALLED FROM SFCTRIF:   
    implicit none
    real :: timestep                      !IN model timestep (s).  
    real :: zlat                          !IN grid latitude (degree)
    real :: dp                            !IN population density
    real :: light                         !IN lightning frequency (count/km2/hr *1e6) 
    real :: agr                           !IN agriculture frac
    real :: gdp                           !IN GDP data K 1995US$/capita
    real :: peatf                         !in peatland fraction (need to devide scale factor 1e6) 
    real :: vfrac                         !IN fraction of veg 
    real :: zwind                         !IN air wind
    real :: tm                            !IN air tmp 
    real :: em                            !IN air vapor pressure 
    real :: tg1                           !IN ground temperature
    real :: prec60                        !IN 60-days of prec (mm/s)
    real,dimension(nnpft):: www1n         !IN 1nd layer soil moisture
    real,dimension(nnpft):: www2n         !IN 2nd layer soil moisture
    real,dimension(npft) ::   leaf1       !IN leaf carbon
    real,dimension(npft) ::   phen        !IN leaf1*phen is current leaf carbon 
    real,dimension(npft) ::   wood1       !IN wood carbon
    real,dimension(ntypes)::  frac_trif   !IN vegetation coverage. 

! !inout to sfctrif:
    real,dimension(npft) :: leaf_loss_dr  !inout accumulated leaf loss due to fire
    real,dimension(npft) :: wood_loss_dr  !inout accumulated wood loss due to fire

! !OTHER LOCAL VARIABLES: 
    real :: ei            ! vapor pressure
    real :: lh            ! anthropogenic ignition sources (/person/month) 
    real :: fs            ! fire unsupressed by human 
    real :: rh            ! relative humidity
    real :: www1n_avg     ! grib-averaged 1nd layer soil moisture
    real :: www2n_avg     ! grib-averaged 2nd layer soil moisture
    real :: ign           ! presence of a source of  ignition in 1 hours
    real :: agb_avg       ! above-ground biomass (leaf+wood) in grib level(gC/m2)
    real :: fm1_avg       ! dependence of fuel combustibility on surface soil wetness 
    real :: fm2_avg       ! dependence of fuel combustibility on rh 
    real :: ft_avg        ! dependence of fuel combustibility on surface temperature 
    real :: fa_avg        ! grib-averaged umax (m/s), parameter table 

    real :: ign_f         ! presence of a source of ignition/s
    real :: lb_lf         ! length-to-breadth ratio   
    real :: fd_lf         ! duration time
    real :: fb_avg        ! available of fuel(output for check) 
    real :: feo_avg       ! socialeconomic influence on fire occurrence
    real :: fe1_avg       ! economic influence on fire occurrence
    real :: fd_avg        ! demographic conditions on the average spread area 
    real :: fe2_avg       ! economic conditions on the average spread area  
    real :: fes_avg       ! socialeconomic conditions on the average spread area  
    real :: fm            ! dependence of fuel combustibility on RH (%) and on surface soil wetness (fm1_avg)
    real :: fcli          ! dependence of peat combusitbility on long-term precipitation
    real :: fsat          ! dependence of peat combusitbility on water table 
    real :: baf_peat      ! peat fire fraction 
    integer :: npft1
    integer :: day  
    real,dimension(npft) :: leaf1_fire !leaf carbon used in fire model
    real,dimension(npft) :: wood1_fire !wood carbon used in fire model
    real,dimension(npft) :: litt1_fire !litter carbon used in fire model

    real, parameter :: lfuel=100.0  ! lower threshold of fuel mass (gC/m2) for ignition
    real, parameter :: ufuel=800.0 ! upper threshold of fuel mass(gC/m2) for ignition 
    real, parameter :: g0=0.05
    real, parameter :: borealat=40.0! latitude threshold for boreal peat fire
    real, parameter :: tf=273.15    ! 273K for frozen peat land
    real, parameter :: c1=0.17*10.0**(-3.0)/3600.0 ! parameter for tropical peat fire (s-1)
    real, parameter :: c2=0.9*10.0**(-5.0)/3600.0 ! parameter for boreal peat fire (s-1)
    real,dimension(npft) :: umax, fm_up, fm_low 

    data umax/0.10,0.13,0.16,0.16,0.12,0.12,0.10/  
    data fm_up/0.70,0.7,0.70,0.70,0.6,0.60,0.70/
    data fm_low/0.3,0.3,0.30,0.30,0.3,0.3,0.3/


! !OUTPUT TO SFCTRIF:   
    real :: nfire_avg      ! OUT grib-averaged the number of fires count/km^2/s
    real :: baf_avg        ! OUT grib-averaged burned area fraction /s
!-----------------------------------------------------------------------
    ! Get model step size

    ei                          = 100.*exp(21.18123-5418./tm)/.622    !Huilin add ei for vapor pressure at saturation
    lh                          = 0.01*6.8*dp**(0.4)/30.0/24.0
    light                       = light/1000000.0
    light                       = light/(5.16+2.16*cos(3*min(60.0,abs(zlat))/90.0*PIE))*0.22
    ign                         = ( lh + light)*vfrac    
    rh                          = em/ei*100
    ft_avg                      = tg1 - 273.15  
    prec60                      = prec60*3600*24.0
    peatf                       = peatf/1000000.0

    agb_avg                     = 0.0
    fm1_avg                     = 0.0
    fm2_avg                     = 0.0
    fa_avg                      = 0.0
    feo_avg                     = 0.0
    fes_avg                     = 0.0
    fd_avg                      = 0.0
    fe1_avg                     = 0.0
    fe2_avg                     = 0.0

    do npft1=1,npft
       leaf1_fire(npft1)        = leaf1(npft1)-leaf_loss_dr(npft1) 
       wood1_fire(npft1)        = wood1(npft1)-wood_loss_dr(npft1) 
       if ((npft1.eq.1).or.(npft1.eq.2).or.(npft1.eq.7)) then
          litt1_fire(npft1)     = (leaf1_fire(npft1)+wood1_fire(npft1))*0.10
       else
          litt1_fire(npft1)     = (leaf1_fire(npft1)+wood1_fire(npft1))*0.30
       end if
       agb_avg                  = agb_avg + (leaf1_fire(npft1)*phen(npft1) + wood1_fire(npft1)+litt1_fire(npft1)) * frac_trif(npft1)*1000.0   !kg->g c/m^2

       fm1_avg                  = fm1_avg + ((1.0-max(0.0,min(1.0,(www2n(npft1)/0.42 &
                                  -fm_low(npft1))/(fm_up(npft1)-fm_low(npft1)))))**1.3*frac_trif(npft1))        !soil www2 on fire
       fm2_avg                  = (1.0-max(0.0,min(1.0,(min(max(rh,0.0),1.0)-0.3)/(0.7-0.3))))**1.1

       if (dp .gt. 0.1) then
           fs                   = 1-(0.01+0.98*exp(-0.025*dp))
           if((npft1.eq.1).or.(npft1.eq.2).or.(npft1.eq.7)) then               ! Tree PFTs
               if (gdp.gt.20) then
                   fe1_avg      = fe1_avg + 0.39*frac_trif(npft1)/vfrac        !GDP on fire occurrence
                   fe2_avg      = fe2_avg + 0.62*frac_trif(npft1)/vfrac        !GDP on fire spread 
               else if ((gdp.le.20.0) .and. (gdp.gt.8.0)) then
                   fe1_avg      = fe1_avg + 0.79*frac_trif(npft1)/vfrac        !GDP on fire occurrence
                   fe2_avg      = fe2_avg + 0.83*frac_trif(npft1)/vfrac        !GDP on fire spread 
               else
                   fe1_avg      = fe1_avg + 1.0*frac_trif(npft1)/vfrac         !GDP on fire occurrence
                   fe2_avg      = fe2_avg + 1.0*frac_trif(npft1)/vfrac         !GDP on fire spread 
               end if
               fd_avg           = fd_avg + (0.4+0.6*exp(-1.0*PIE*(dp/125.0)))*frac_trif(npft1)/vfrac 
           else 
               fe1_avg          = fe1_avg + (0.1+0.9*exp(-1.0*PIE*(gdp/8.0)**0.5))*frac_trif(npft1)/vfrac   !GDP on fire occurrence
               fd_avg           = fd_avg + (0.2+0.8*exp(-1.0*PIE*(dp/450.0)**0.5))*frac_trif(npft1)/vfrac 
               fe2_avg          = fe2_avg + (0.2+0.8*exp(-1.0*PIE*(gdp/7.0)**0.5))*frac_trif(npft1)/vfrac
           end if   
       else
           fs                   = 0.0
           fe1_avg              = 1.0 
           fd_avg               = 1.0 
           fe2_avg              = 1.0 
       end if
       feo_avg                  = fe1_avg*(1-fs)
       fes_avg                  = fe2_avg*fd_avg
       fa_avg                   = fa_avg + umax(npft1)*frac_trif(npft1)
    end do


! ********fire occurrence**********
    ign_f                       = ign /3600.0
    fb_avg                      = max(0.0,min(1.0,(agb_avg-lfuel)/(ufuel-lfuel)))
    if (ft_avg.ge.-10.0) then
       fm                       = fm1_avg*fm2_avg 
    else
       fm                       = 0.0
    end if
    nfire_avg                   = ign_f*fb_avg*fm*feo_avg*(1.0-agr)*(1.0-agr)
     
! ********fire spread **********
    fd_lf                       = 24*3600
    lb_lf                       = 1.0+10.0*(1.0-exp(-0.06*zwind))
      
! ********burn area fraction****
    baf_avg                     = (g0*sqrt(fm)*fa_avg *fd_lf/1000)**2*fes_avg*nfire_avg*PIE*lb_lf

! ********Peat fire*************

    www1n_avg                   = 0.0
    www2n_avg                   = 0.0
    do npft1=1,npft
      www1n_avg                 = www1n_avg + www1n(npft1)*frac_trif(npft1)         
      www2n_avg                 = www2n_avg + www2n(npft1)*frac_trif(npft1)         
    end do
    if(zlat.lt.borealat) then
      fcli                      = max(0.0,min(1.0,(4.0-prec60)/4.0))**2 
      fsat                      = max(0.0,(1.0-www1n_avg/0.42)) 
      baf_peat                  = c1*fcli*peatf*fsat
    else
      fcli                      = exp(-PIE*www2n_avg/0.3)*max(0.0,min(1.0,(tg1-tf)/10.0)) 
      fsat                      = max(0.0,(1.0-www1n_avg/0.42)) 
      baf_peat                  = c2*fcli*peatf*fsat
    end if

    baf_avg                     = baf_avg+baf_peat 

  END SUBROUTINE 
END MODULE
