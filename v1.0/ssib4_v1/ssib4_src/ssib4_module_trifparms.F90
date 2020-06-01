MODULE ssib4_module_trifparms
!
!yliu DBL      INTEGER,parameter:: NPFT=6,NNPFT=7                                                   
      INTEGER,parameter:: NPFT=7,NNPFT=8                                                   
!----------------------------------------------------------------------- 
! Functional Type dependent parameters                                  
!-----------------------------------------------------------------------
      INTEGER,DIMENSION(NPFT) :: &                                                           
            C3,                  &! 1 for C3 Plants, 0 for C4 Plants.  
            CROP                  ! 1 for crop type, 0 for non-crop.  
                                                                      
      REAL,DIMENSION(NPFT)::     &                                                           
       ALPHA,                    &! Quantum efficiency (mol CO2/mol PAR photons).       
       A_WL,                     &! Allometric coefficient relating 
                                  ! the target woody biomass to the leaf area index (kg C/m2).     
       A_WS,                     &! Woody biomass as a multiple of live stem biomass.               
       B_WL,                     &! Allometric exponent relating    
                                  ! the target woody biomass to the leaf area index.          
       DGL_DM,                   &! Rate of change of leaf turnover   
                                  ! rate with moisture availability. 
       DGL_DT,                   &! Rate of change of leaf turnover rate with temperature (/K)    
       DQCRIT,                   &! Critical humidity deficit(kg H2O/kg air).                     
       DSAI_DLAI,                &! Rate of change of stem area index with balanced LAI (for SSiB).             
       ETA_SL,                   &! Live stemwood coefficient (kg C/m/LAI).                     
       FSMC_OF,                  &! Moisture availability below which leaves are dropped.       
       F0,                       &! CI/CA for DQ = 0.               
       GLMIN,                    &! Minimum leaf conductance for H2O 
       G_AREA,                   &! Disturbance rate (/360days).   
       G_GROW,                   &! Rate of leaf growth (/360days).  
       G_LEAF_0,                 &! Minimum turnover rate for leaves (/360days).                     
       G_ROOT,                   &! Turnover rate for root biomass (/360days).                   
       G_WOOD,                   &! Turnover rate for woody biomass (/360days).                    
       KPAR,                     &! PAR Extinction coefficient(m2 leaf/m2 ground).        
       LAI_MAX,                  &! Maximum projected LAI.    
       LAI_MIN,                  &! Minimum projected LAI.            
       NL0,                      &! Top leaf nitrogen concentration(kg N/kg C).                    
       NR_NL,                    &! Ratio of root nitrogen concentration to leaf nitrogen concentration.  
       NS_NL,                    &! Ratio of stem nitrogen concentration to leaf nitrogen concentration.       
       OMEGA,                    &! Leaf scattering coefficient for PAR.                   
       R_GROW,                   &! Growth respiration fraction.   
       SIGL,                     &! Specific density of leaf carbon (kg C/m2 leaf).                
       TLEAF_OF,                 &! Temperature below which leaves are dropped.                          
       TLOW,                     &! Lower temperature for photosynthesis (deg C)          
       TUPP                       ! Upper temperature for photosynthesis (deg C)      
!                                                              
!----------------------------------------------------------------------  
!                        BT     NT    C3G    C4G     S                 
!---------------------------------------------------------------------- 
      DATA C3      /   1, 1, 1, 0, 1 ,1 ,1/                               
      DATA CROP    /   0, 0, 1, 1, 0 ,0 ,0/                                  
      DATA ALPHA   /   0.08,  0.08,  0.08, 0.040,  0.08,0.04 ,0.08/           
      DATA A_WL    /   0.65,  0.65, 0.005, 0.005,  0.10,0.10 ,0.65/          
      DATA A_WS    /  10.00, 10.00,  1.00,  1.00, 10.00 ,10.0,10.0/         
      DATA B_WL    /  1.667, 1.667, 1.667, 1.667, 1.667,1.667,1.667 /       
!sav      DATA DGL_DM  /  0.0,  0.0,  0.0, 90.0, 0.0,90.0 /          
!sav      DATA DGL_DT  /    9.0,   9.0,   0.0, 0.0, 9.0,9.0 /       
      DATA DQCRIT  /  0.090, 0.060, 0.100, 0.075, 0.100,0.100,0.090 /     
      DATA DSAI_DLAI /  0.50,  0.50, 0.000, 0.000,  0.25,0.25,0.50 /
      DATA ETA_SL  /   0.01,  0.01,  0.01,  0.01,  0.01,0.01,0.01 /     
      DATA F0      /  0.875, 0.875, 0.900, 0.800, 0.900,0.80,0.875 /    
!sav      DATA FSMC_OF /   0.0,  0.0,  0.0,  0.0,  0.5,0.4 /      
      DATA GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6,1.0E-6,10E-6,1.0E-6 /  
!Huilin      DATA G_AREA  /  0.004, 0.004,  0.10,  0.10,  0.05,0.05,0.004 /   
!     DATA G_AREA  /  0.004, 0.004,  0.0,  0.0,  0.04,0.0,0.004 /   
      DATA G_AREA  /  0.004, 0.004,  0.02,  0.02,  0.04,0.01,0.004 /   
      DATA G_GROW  /  20.00, 20.00, 20.00, 20.00, 20.00,20.00,20.00 / 
      DATA G_LEAF_0/   0.25,  0.25,  0.25,  0.25,  0.25,0.25, 0.25 / 
      DATA G_ROOT  /   0.25,  0.25,  0.25,  0.25,  0.25,0.25, 0.25 /   
      DATA G_WOOD  /   0.01,  0.01,  0.20,  0.20,  0.05,0.05, 0.01 /  
      DATA KPAR    /   0.50,  0.50,  0.50,  0.50,  0.50,0.50, 0.50 /  
!zzq20090313      DATA LAI_MAX /   9.00,  9.00,  4.0,  4.0,  4.0,4.0 / 
!zzq20090313      DATA LAI_MIN /   3.00,  3.00,  1.0,  1.0,  1.0,1.0 / 
!zzq2012nov30      DATA LAI_MAX /   9.00,  9.00,  4.0,  4.0,  4.0,4.0 / 
!zzq2012dec08      DATA LAI_MAX /   9.00,  9.00,  4.0,  4.0,  2.0,2.0 /
!zzq2013-1-14      DATA LAI_MAX /   9.00,  9.00,  4.0,  4.0,  4.0,2.0 / 
!yliu23Sep2014      DATA LAI_MAX /   8.00,  8.00,  4.0,  4.0,  3.0,3.0 /         
      DATA LAI_MAX /   8.00,  8.00,  4.0,  4.0,  3.0,3.0, 8.00 /         
!zzq20093-23   DATA LAI_MIN /   3.00,  3.00,  0.7,  1.0,  0.25, 0.35/    
!zzq2012nov30      DATA LAI_MIN /   0.3,  0.3,  0.1,  0.1,  0.1, 0.1/ 
      DATA LAI_MIN /   2.00,  2.00,  0.5,  0.5,  0.5, 0.1,1.00/          
      DATA NL0     /  0.040, 0.030, 0.060, 0.030, 0.030,0.060,0.040 /   
      DATA NR_NL   /   2.00,  2.00,  2.00,  1.00,  1.00,2.00,2.00 /  
      DATA NS_NL   /   0.10,  0.10,  1.00,  1.00,  0.10,0.10,0.10 /
      DATA OMEGA   /   0.15,  0.15,  0.15,  0.17,  0.15 ,0.17,0.15/
      DATA R_GROW  /   0.25,  0.25,  0.25,  0.25,  0.25,0.25,0.25 / 
      DATA SIGL    / 0.0375,0.1000,0.0250,0.0500,0.0500,0.075,0.0375 / 
!sav      DATA TLEAF_OF/ 273.15,243.15,258.15,258.15,243.15,273.15 /
      DATA TLOW    /    0.0,  -5.0,   0.0,  13.0,   0.0,13.0,-0.5 /      
      DATA TUPP    /   36.0,  31.0,  36.0,  45.0,  36.0,45.0,36.0 /    
!
      DATA DGL_DM  /    90.0,  90.0,  90.0, 90.0,  90.0,90.0,90.0 /      
      DATA DGL_DT  /    9.0,   9.0,   9.0, 9.0, 9.0,9.0,9.0 /          
!zzq2012nov30      DATA FSMC_OF / 0.6,  0.2,  0.2,  0.2,  0.2, 0.2 /    
      DATA FSMC_OF / 0.6,  0.2,  0.2,  0.2,  0.1, 0.2, 0.4 /            
!zzq2012nov30      DATA TLEAF_OF/ 283.15,265.15,268.15,273.15,273.15,263.15 /  
!zzq2012dec04      DATA TLEAF_OF/ 280.15,265.15,270.15,270.15,273.15,258.15 /
!zzq2012dec05      DATA TLEAF_OF/ 280.15,265.15,268.15,268.15,273.15,268.15 / 
!zzq2012dec09      DATA TLEAF_OF/ 280.15,265.15,268.15,268.15,267.15,267.15 /
      DATA TLEAF_OF/ 284.15,265.15,268.15,268.15,273.15,268.15,278.15 /   
!
END MODULE
