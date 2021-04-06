MODULE ssib4_module_fireparms
      use ssib4_module_trifparms
!
! INTEGER,parameter:: NPFT=6,NNPFT=7                                                   
!----------------------------------------------------------------------- 
! Functional Type dependent parameters                                  
!-----------------------------------------------------------------------
      REAL,DIMENSION(NPFT)::     &                                                           
       UMAX,                     &! Average miximum fire spread rate in natural vegetation region
       FLAM,                     &! Plant flammability threshold  
       CC_leaf,                  &! Combustion Completeness for leave
       CC_wood,                  &! Combustion Completeness for wood 
       CC_root,                  &! Combustion Completeness for root 
       CC_abgl,                  &! Combustion Completeness for abg litter 
       MF_leaf,                  &! Mortality Factor for leave
       MF_wood,                  &! Mortality Factor for wood 
       MF_root,                  &! Mortality Factor for root 
       EFG_CO2,                  &! Emission Factor for CO2
       EFG_CO,                   &! Emission Factor for CO 
       EFG_CH4,                  &! Emission Factor for CH4 
       EFG_NMHC,                 &! Emission Factor for NMHC
       EFG_H2,                   &! Emission Factor for H2 
       EFG_NOx,                  &! Emission Factor for NOx 
       EFG_N2O,                  &! Emission Factor for N2O 
       EFA_PM25,                 &! Emission Factor for PM25
       EFA_TPM,                  &! Emission Factor for TPM
       EFA_TC,                   &! Emission Factor for TC 
       EFA_OC,                   &! Emission Factor for OC 
       EFA_BC                     ! Emission Factor for BC 
!                                                              
!----------------------------------------------------------------------  
!              rainforest boreal C3   C4 shrub tundra deciduous_broadleaf  
!---------------------------------------------------------------------- 
      DATA FLAM/    0.69, 0.69, 0.69,0.69,0.69,0.69,0.69 /
      DATA CC_leaf/ 0.60, 0.60, 0.70,0.70,0.70,0.70,0.60 /    
      DATA CC_wood/ 0.08, 0.13, 0.50,0.50,0.50,0.25,0.08 /
      DATA CC_root/ 0.00, 0.00, 0.00,0.00,0.00,0.00,0.00 / 
      DATA CC_abgl/ 0.70, 0.70, 0.80,0.80,0.60,0.80,0.70 / 
      DATA MF_leaf/ 0.60, 0.60, 0.60,0.60,0.60,0.60,0.60 /
      DATA MF_wood/ 0.12, 0.20, 0.50,0.50,0.40,0.20,0.12 /
      DATA MF_root/ 0.10, 0.15, 0.30,0.30,0.15,0.15,0.10 /
      DATA EFG_CO2/ 1613, 1549, 1647,1647,1647,1647,1566 /
      DATA EFG_CO/   108,  124,   70,  70,  70,  70, 108 / 
      DATA EFG_CH4/  6.3,  5.1,  2.5, 2.5, 2.5, 2.5, 5.8 /
      DATA EFG_NMHC/ 7.1,  5.3,  5.7, 5.7, 5.7, 5.7,14.6 /
      DATA EFG_H2/  3.11, 1.66, 0.97,0.97,0.97,0.97,2.09 /
      DATA EFG_NOx/ 2.55, 1.69, 2.58,2.58,2.58,2.58,2.09 /
      DATA EFG_N2O/ 0.20, 0.25, 0.18,0.18,0.18,0.18,0.25 /
      DATA EFA_PM25/ 8.3, 15.3,  7.5, 7.5, 7.5, 7.5,18.1 /
      DATA EFA_TPM/ 10.9, 15.3,  8.5, 8.5, 8.5, 8.5,18.1 /
      DATA EFA_TC/   6.0, 10.6,  3.4, 3.4, 3.4, 3.4, 8.4 /
      DATA EFA_OC/   4.5, 10.1,  3.1, 3.1, 3.1, 3.1, 8.9 /
      DATA EFA_BC/  0.49, 0.50, 0.51,0.51,0.51,0.51,0.66 /
!
END MODULE
