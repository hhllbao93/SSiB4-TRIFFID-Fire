!================================================================================    
! Note:This version was changed by Jack Zhang from original one. 
!   The main change are: 
!        1) Dimension change; 
!        2) Inclusion of the trifparms.h file.  01/03/2008
!        3) TS1,STHU,V_SAT and V_WILT are input for each type. 05/12/2008
!        4) re-constructed from F77 to F90 by zzq 09/24/2013
!================================================================================    
!#define DEBUG_PRINT
!*DECK SITRIFFID
MODULE  ssib4_module_sitriffid
      use ssib4_module_trifparms
      use ssib4_module_comsveg
!
      IMPLICIT NONE                                                     
!
!--------------------------------------------------------------------------------
contains
      SUBROUTINE SITRIFFID(   &
             TIMESTEP,     & !IN (1) Model timestep 
             L_PHEN,       & !IN (1) Logical
             DAY_PHEN,     & !IN (1) Time
             L_TRIF,       & !IN (1) Logical
             DAY_TRIF,     & !IN (1) Time
             VEG_EQUIL,    & !IN (1) Logical
             ANETC,        & !IN (NPFT) Anet fraction
             RDC,          & !IN (NPFT) RDC fraction
             FSMC,         & !IN (NPFT) RSTFAC(1,2)
             STHU,         & !IN (npft) www(2) 
             TS1,          & !IN (npft) TD 
             TSTAR,        & !IN (NPFT) TC average
             V_SAT,        & !IN (npft) SPOROS prescribed
             V_WILT,       & !IN (npft) wilting point prescribed
             FPARn,        & !IN FPAR
             KPARn,        & !IN KPAR
             FRACA,        & !IN (1) Agriculture factor prescribed
             G_ANTH,       & !IN (1) Anthropogenic factor prescribed
             AGR,          & !IN Huilin read agriculture fraction 
             CS,           & !INOUT (1)     diagnose
             FRAC,         & !INOUT (NTYPES)diagnose
             HT,           & !INOUT (NPFT)  diagnose
             LAI,          & !INOUT (NPFT)  diagnose
             G_LEAF_DAY,   & !INOUT (NPFT)  diagnose
             G_LEAF_DR,    & !INOUT (NPFT)  diagnose
             GPP_DR,       & !INOUT (NPFT)  accumulate and initialize
             NPP_DR,       & !INOUT (NPFT)  accumulate and initialize
             RESP_S_DR,    & !INOUT (1)     accumulate and initialize
             RESP_W_DR,    & !INOUT (NPFT)  accumulate and initialize
             VEG_FRAC,     & !INOUT (1)     accumulate and initialize
             FRAC_VS,      & !INOUT (1)     accumulate and initialize
             LEAF_LOSS_DR, & !INOUT(NPFT) leaf loss due to fire Huilin
             WOOD_LOSS_DR, & !INOUT(NPFT) wood loss due to fire Huilin
             ROOT_LOSS_DR, & !INOUT(NPFT) rood loss due to fire Huilin
             GPP,          & !OUT (NPFT)
             NEP,          & !OUT (1)
             NPP,          & !OUT (NPFT)
             RESP_P,       & !OUT (NPFT)
             RESP_S,       & !OUT (1)
             LEAF1,        & !OUT (NPFT) 
             ROOT1,        & !OUT (NPFT) 
             WOOD1,        & !OUT (NPFT) 
             DCVEG1,       & !OUT (NPFT)
             LIT_C,        & !OUT (NPFT)
             LIT_C_T,      & !OUT (1)
             PHEN,         & !OUT (NPFT) 
             LAID,         & !OUT (NPFT)
             green,        & !OUT (NPFT)
             RESP_P_M,     & !OUT (NPFT)
             RESP_P_G      & !OUT (NPFT)
             )

!-----------------------------------------------------------------------
!
!      Subroutine for coupling TRIFFID to SSiB for Yongkang (Dec 2001).
!
!-----------------------------------------------------------------------
      INTEGER,PARAMETER::        &                                      
       SOIL=9                      ! Index of the surface type 'Soil'       
!
      INTEGER,PARAMETER::        &                                      
       ITER_EQ =10                ! Number of TRIFFID iterations for         
                                  ! gradient descent to equilibrium.         
      REAL,PARAMETER::           &                                      
      GAMMA_EQ =1.0E-1            ! Inverse timestep for gradient     
                                  ! descent to equilibrium (/360days)

      INTEGER::                  &                                  
       N,KITER,II                 ! Loop counters

!-----------------------------------------------------------------------
! Timesteps and logical flags
!-----------------------------------------------------------------------
      REAL::                    &
       TIMESTEP                   ! IN Model timestep (s).
      
      INTEGER::                 &
       DAY_PHEN,                & ! IN Number of days between phenology.
       DAY_TRIF                   ! IN Number of days between TRIFFID.

      LOGICAL::                 &
       L_PHEN,                  & ! IN .T. if phenology to be updated.
       L_TRIF,                  & ! IN .T. if vegetation to be updated.
       VEG_EQUIL                  ! IN .T. if the vegetation equilibrium 
!                                 !    is required.

!-----------------------------------------------------------------------
! Inputs from SSiB                             
!-----------------------------------------------------------------------
      REAL,DIMENSION(NPFT)::    &
       ANETC,                   & ! IN Net canopy photosynthesis(mol CO2/m2/s).                        
       RDC,                     & ! IN Canopy dark respiration  (mol CO2/m2/s).                       
       FSMC,                    & ! IN Moisture availability factor.
       STHU,                    & ! IN Unfrozen soil moisture as a fraction
                                  !    of saturation.
       TS1,                     & ! IN Sub-surface temperature (K).
       TSTAR,                   & ! IN Surface temperature on veg tiles (K).                  
       V_SAT,                   & ! IN Volumetric soil moisture                  
                                  !    concentration at saturation (m3 H2O/m3 soil).                       
       V_WILT,                  & ! IN Volumetric soil moisture   
                                  !    concentration below which stomata             
                                  !    close (m3 H2O/m3 soil).                    
       FPARn,                   &
       KPARn
!-----------------------------------------------------------------------
! Disturbance prescription.
!-----------------------------------------------------------------------
      REAL::                    &
       FRACA,                   & ! IN Areal fraction of agriculture.
       G_ANTH                     ! IN Anthropogenic disturbance rate (/yr).

!-----------------------------------------------------------------------
! Agriculture 
!-----------------------------------------------------------------------
      REAL::                    &
       AGR                        ! IN Agriculture fraction

!-----------------------------------------------------------------------
! Prognostic variables
!-----------------------------------------------------------------------
      REAL::                    &
       CS                         ! INOUT Soil carbon (kg C/m2).

      REAL,DIMENSION(NTYPES)::  &
       FRAC                       ! INOUT Areal coverage.

      REAL,dimension(NPFT)::    &
       HT,                      & ! INOUT Canopy height (m)
       LAI,                     & ! INOUT Leaf area index.
       PHEN,                    & ! INOUT Leaf area index.
       LEAF1,                   & ! INOUT Canopy height (m)
       ROOT1,                   & ! INOUT Leaf area index.
       WOOD1,                   & ! INOUT Canopy height (m)
       DCVEG1,                  & ! INOUT Leaf area index.
       LAID,                    & ! OUT Dead leaf area index.
       green                      ! OUT green fraction.

!-----------------------------------------------------------------------
! Other variables which need to be saved between calls
!-----------------------------------------------------------------------
      REAL,DIMENSION(NPFT)::    &
       G_LEAF_DAY,              & ! INOUT Daily mean leaf turnover rate (/yr).
       G_LEAF_DR,               & ! INOUT Accumulated leaf turnover rate(/yr).
       GPP_DR,                  & ! INOUT Accumulated Gross Primary Productivity(kg C/m2/yr). 
       NPP_DR,                  & ! INOUT Accumulated Net Primary Productivity (kg C/m2/yr). 
       RESP_W_DR,               & ! INOUT Accumulated wood respiration rate (kg C/m2/yr).
       LEAF_LOSS_DR,            & ! INOUT ACCUMULATED leaf loss due to fire Huilin 
       WOOD_LOSS_DR,            & ! INOUT ACCUMULATED wood loss due to fire Huilin 
       ROOT_LOSS_DR               ! INOUT ACCUMULATED root loss due to fire Huilin

      REAL::                    &
       RESP_S_DR,               & ! INOUT Accumulated soil respiration rate(kg C/m2/yr).
       VEG_FRAC,                & ! INOUT Vegetated fraction.
       FRAC_VS                    ! INOUT Total fraction of gridbox covered by
!                                 !       vegetation and soil.
!-----------------------------------------------------------------------
! Output Carbon Fluxes.
!-----------------------------------------------------------------------
      REAL,DIMENSION(NPFT)::    &
       GPP,                     & ! OUT Gross Primary Productivity (kg C/m2/s). 
       NPP,                     & ! OUT Net Primary Productivity (kg C/m2/s). 
       RESP_P                     ! OUT Plant respiration rate (kg C/m2/s).

      REAL::                    &
       NEP,                     & ! OUT Net Ecosystem Productivity (kg C/m2/s). 
       RESP_S,                  & ! OUT Soil respiration rate (kg C/m2/s).
       RESP_S1                    ! WORK Soil respiration rate (kg C/m2/s).

!-----------------------------------------------------------------------
! Work Carbon variables.
!-----------------------------------------------------------------------
      REAL,DIMENSION(NPFT)::    &
       C_VEG,                   & ! Vegetation carbon (kg C/m2).
       G_LEAF,                  & ! Leaf turnover rate (/yr).
       G_LEAF_PHEN,             & ! Daily leaf turnover rate including phenology (/yr).
       LIT_C,                   & ! Carbon litter (kg C/m2/yr).
       RESP_W                     ! Wood maintenance respiration rate (kg C/m2/s).

      REAL::                    &
       CV,                      & ! Gridbox mean vegetation carbon (kg C/m2).
       FTIME,                   & ! Weighting factor for accumulations.
       FTIME_PHEN,              & ! Weighting factor for accumulations.
       LIT_C_T                    ! Gridbox mean carbon litter (kg C/m2/yr).
          
!-----------------------------------------------------------------------
! Additional variables required to calculate plant respiration.
!-----------------------------------------------------------------------
      REAL::                    &                                      
       FPAR,                    & ! PAR absorption factor. 
       LAI_BAL,                 & ! Leaf area index in balanced growth state.                         
       NL,                      & ! Mean leaf nitrogen concentration (kg N/kg C).            
       NL_BAL,                  & ! Mean leaf nitrogen concentration in balanced             
                                  ! growth state (kg N/kg C).             
       N_LEAF,                  & ! Nitrogen contents of the leaf,        
       N_ROOT,                  & ! root,                                 
       N_STEM,                  & ! and stem (kg N/m2).                   
       ROOT                       ! Root carbon (kg C/m2).               
       REAL,DIMENSION(NPFT)::                    &
        RESP_P_G,                & ! Plant growth respiration rate(kg C/m2/sec).
        RESP_P_M                   ! Plant maintenance respiration rate (kg C/m2/sec).
!!!--
 
!-----------------------------------------------------------------------
! Time parameters
!-----------------------------------------------------------------------
      REAL::                    &
       DTIME_PHEN,              & ! Timestep for phenology (/yr).
       FORW,                    & ! Forward timestep weighting.
       GAMMA                      ! Inverse timestep (/yr).
      
      REAL, PARAMETER::         & 
       DAY_YEAR = 365.0,        & ! Number of days in a year (days).
       SEC_DAY  = 86400.,       & ! Number of seconds in a day (s).
       SEC_YEAR =               & 
              DAY_YEAR*SEC_DAY    ! Number of seconds in a year (s).

!-----------------------------------------------------------------------
! Factors for accummulation variables
!-----------------------------------------------------------------------
      FTIME=TIMESTEP/REAL(SEC_DAY*DAY_TRIF)
      FTIME_PHEN=TIMESTEP/REAL(SEC_DAY*DAY_PHEN)

!-----------------------------------------------------------------------
! Calculate fluxes required by TRIFFID from those given by SSiB
!-----------------------------------------------------------------------
      DO N=1,NPFT
!-----------------------------------------------------------------------     
! Assume that root biomass is equal to balanced growth leaf biomass         
!-----------------------------------------------------------------------      
          LAI_BAL = (A_WS(N)*ETA_SL(N)*HT(N)/A_WL(N))**(1.0/(B_WL(N)-1))                                
          ROOT = SIGL(N) * LAI_BAL                                     
                                                                       
!-----------------------------------------------------------------------      
! Calculate the actual and balanced mean leaf nitrogen concentration          
! assuming perfect light acclimation                                          
!-----------------------------------------------------------------------    
          FPAR = FPARn(N)
          NL = (FPAR / LAI(N)) * NL0(N)                                 
          NL = (FPAR / LAI(N)) * NL0(N)                                 
          NL_BAL = FPARn(N)/LAI_BAL * NL0(N)
                                                                        
!-----------------------------------------------------------------------      
! Calculate the total nitrogen content of the leaf, root and stem             
!----------------------------------------------------------------------      
          N_LEAF = NL * SIGL(N) *LAI_BAL
          N_ROOT = NR_NL(N) * NL * ROOT
          N_STEM = NS_NL(N) * NL * ETA_SL(N) * HT(N) * LAI_BAL
!!!--
!-----------------------------------------------------------------------      
! Calculate the Gross Primary Productivity, the plant maintenance             
! respiration rate, and the wood maintenance respiration rate                 
! in kg C/m2/se
! NB. This assumes that all respiration fluxes are reduced as FSMC drops
!     (this is different from implementation in MOSES)                                            
!-----------------------------------------------------------------------     
          GPP(N) = 12.0E-3 * (ANETC(N) + RDC(N))                       
          RESP_P_M(N) = 12.0E-3 * RDC(N) * (N_LEAF + N_STEM + N_ROOT) / N_LEAF          
          RESP_W(N) = 12.0E-3 * RDC(N) * N_STEM / N_LEAF               
                                                                        
!-----------------------------------------------------------------------      
! Calculate the total plant respiration and the Net Primary Productivity    
!-----------------------------------------------------------------------      
          RESP_P_G(N) = R_GROW(N) * (GPP(N) - RESP_P_M(N))                   
          RESP_P(N) = RESP_P_M(N) + RESP_P_G(N)                             
          NPP(N) = GPP(N) - RESP_P(N)                                   
                                                                        
      ENDDO
                                                                       
!-----------------------------------------------------------------------
! Calculate the leaf turnover rate for each PFT.
!-----------------------------------------------------------------------
      DO N=1,NPFT

        CALL LEAF_LIT (N, FSMC,TSTAR,G_LEAF(N))

      ENDDO ! FT Loop

!-----------------------------------------------------------------------
! Calculate the soil respiration
!-----------------------------------------------------------------------

      resp_s=0 
      veg_frac=0

      DO N=1,NPFT
        CALL MICROBE (CS,STHU(N),V_SAT(N),V_WILT(N),TS1(N),RESP_S1)
        resp_s=resp_s+FRAC(n)*resp_s1
        veg_frac =veg_frac+frac(n)
      enddo

      if(veg_frac.gt.0.0) then
        resp_s = resp_s/veg_frac
      else
        resp_s =0.0
      endif
!      
!-----------------------------------------------------------------------
! Update phenology accumulation variable.
!-----------------------------------------------------------------------
      DO N=1,NPFT
          G_LEAF_DAY(N) = G_LEAF_DAY(N) + G_LEAF(N) * FTIME_PHEN
      ENDDO
      DO N=1,NPFT
          G_LEAF_PHEN(N)=0.0
      ENDDO
!-----------------------------------------------------------------------
! Update leaf phenological state
!-----------------------------------------------------------------------
      IF (L_PHEN) THEN 
        DTIME_PHEN = REAL(DAY_PHEN)/DAY_YEAR
        DO N=1,NPFT
          CALL PHENOL (N,G_LEAF_DAY(N),HT(N),DTIME_PHEN,   &
                 G_LEAF_PHEN(N),LAI(N))           
          G_LEAF_DAY(N)=0.0
          PHEN(N)   = LAI(N)/(LEAF1(N)/SIGL(N))
        ENDDO

      ENDIF ! End of PHENOL call

!-----------------------------------------------------------------------
! Accumulate TRIFFID driving variables
!-----------------------------------------------------------------------
        RESP_S_DR = RESP_S_DR+RESP_S*FTIME*SEC_YEAR 

        DO N=1,NPFT
          GPP_DR(N) = GPP_DR(N) + GPP(N)*FTIME*SEC_YEAR
          NPP_DR(N) = NPP_DR(N) + NPP(N)*FTIME*SEC_YEAR
          G_LEAF_DR(N) = G_LEAF_DR(N) + G_LEAF_PHEN(N)*FTIME
          RESP_W_DR(N) = RESP_W_DR(N) + RESP_W(N)*FTIME*SEC_YEAR
        ENDDO

!----------------------------------------------------------------------
! Update the vegetation areal coverages, structural parameters,
! and soil carbon.
!----------------------------------------------------------------------
      IF (L_TRIF) THEN 

        IF (VEG_EQUIL) THEN
          FORW = 1.0
          GAMMA = GAMMA_EQ
          KITER=ITER_EQ
        ELSE
           FORW = 0.0
           GAMMA = DAY_YEAR/DAY_TRIF
           KITER=1
        ENDIF

        DO II=1,KITER

          CALL TRIFFID (FORW,GAMMA,                              &
                        FRAC_VS,FRACA,G_LEAF_DR,                 &
                        NPP_DR,RESP_S_DR,RESP_W_DR,              &
                        CS,FRAC,HT,LAI,C_VEG,CV,LIT_C,LIT_C_T,   &
                    AGR,LEAF_LOSS_DR,WOOD_LOSS_DR,ROOT_LOSS_DR,  &  !Huilin add for fire
                 LEAF1,ROOT1,WOOD1,DCVEG1,PHEN,LAID,green,       & 
                 FPARn,KPARn)
        ENDDO

!---------------------------------------------------------------------
! Zero the accumulated driving variables
!----------------------------------------------------------------------
        RESP_S_DR = 0.0 
        DO N=1,NPFT
          GPP_DR(N) = 0.0
          NPP_DR(N) = 0.0 
          G_LEAF_DR(N) = 0.0 
          RESP_W_DR(N) = 0.0 
          LEAF_LOSS_DR(N) = 0.0         !Huilin clear LEAF_LOSS
          WOOD_LOSS_DR(N) = 0.0         !Huilin clear WOOD_LOSS
          ROOT_LOSS_DR(N) = 0.0         !Huilin clear ROOT_LOSS  
        ENDDO

!----------------------------------------------------------------------
! Define gridbox mean diagnostics
!----------------------------------------------------------------------
        VEG_FRAC = 0.0
        NEP=-RESP_S

        DO N=1,NPFT
            VEG_FRAC = VEG_FRAC + FRAC(N)
            NEP = NEP + FRAC(N)*NPP(N)
        ENDDO
     
#ifdef DEBUG_PRINT
        print*,FRAC
#endif
        FRAC_VS = VEG_FRAC + FRAC(SOIL)

      ENDIF  ! End of TRIFFID call

      END SUBROUTINE
! *****************************COPYRIGHT******************************    
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.   
!                                                                         
! Use, duplication or disclosure of this code is subject to the           
! restrictions as set forth in the contract.                              
!                                                                         
!                Meteorological Office                                    
!                London Road                                              
!                BRACKNELL                                                
!                Berkshire UK                                             
!                RG12 2SZ                                                 
!                                                                         
! If no contract has been raised with this copy of the code, the use,     
! duplication or disclosure of it is strictly prohibited.  Permission     
! to do so must first be obtained in writing from the Head of Numerical   
! Modelling at the above address.                                         
! ******************************COPYRIGHT******************************   
!!! Subroutine COMPETE ------------------------------------------------   
!!!                                                                       
!!! Purpose : Updates fractional coverage of each functional type.        
!!!           Requires a dominance hierachy as input.                     
!!!                                                                       
!!!                                                                       
!!!  Model            Modification history:                               
!!! version  Date                                                         
!!!  4.4     10/97     New Deck. Peter Cox                                
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.    
!!!                    Richard Betts                                      
!!!                                                                       
!!!END ----------------------------------------------------------------   
      SUBROUTINE COMPETE (DOM,B,DB_DFRAC,FORW,GAMMA,NOSOIL, FRAC,DFRAC)
!
      INTEGER,PARAMETER::         &
        SOIL = 9                    ! Index of the surface type 'Soil'      

      REAL, PARAMETER::           &
       FRAC_MIN = 1.0E-6,         & ! Minimum areal fraction for PFTs.      
       FRAC_SEED =  0.01            ! "Seed" fraction for PFTs.             
!
      REAL,PARAMETER::            &
       DENOM_MIN = 1.0E-6,        & ! Minimum value for the denominator     
       GAMMA_EQ  = 1.0E-1,        & ! of the update equation. Ensures
       ITER_EQ   = 10               ! that gradient descent does not
                                    ! lead to an unstable solution.      
                                    ! descent to equilibrium (/360days)
      INTEGER::                   &
!                                   ! TRIFFID may operate                
        K,M,N                       ! WORK Loop counters.                   
                                                                       
      INTEGER,DIMENSION(NPFT)::   &
       DOM                          ! IN Dominance hierachy.                
                                                                       
      REAL,DIMENSION(NPFT)::      &
       B                            ! IN Mean rate of change of vegetation fraction over           
                                    !    the timestep (kg C/m2/360days).    
      REAL,DIMENSION(NPFT,NPFT):: & 
       DB_DFRAC                     ! IN Rate of change of B                
                                    !    with vegetation fraction.          
      
      REAL::                      &
       FORW,                      & ! IN Forward weighting factor.          
       GAMMA,                     & ! IN Inverse timestep (/360days).       
       NOSOIL                       ! IN Fractional area not available      
                                    !    to vegetation.                     
      REAL,DIMENSION(NTYPES)::    &
       FRAC                         ! INOUT Updated areal fraction.         
 
      REAL,DIMENSION(NPFT)::      &
       DFRAC                        ! OUT Increment to areal fraction.      

      REAL::                      &
       DENOM,                     & ! WORK Denominator of update equation.                        
       FRACN,FRACM,               & ! WORK Fractions used in the spreading  
                                    !      calculation.                     
       NUMER,                     & ! WORK Numerator of the update equation.                        
       SPACE,                     & ! WORK Available space.                 
       P1,P2,Q1,Q2,R1,R2            ! WORK Coefficients in simultaneous equations.                       
!----------------------------------------------------------------------   
! Local parameters                                                        
!----------------------------------------------------------------------   
                                                                       
!----------------------------------------------------------------------   
! Initialisations. Set increments to zero and define the space            
! available to the dominant type leaving space for the seeds of others.   
!----------------------------------------------------------------------   
        DO N=1,NPFT                                                    
          DFRAC(N) = 0.0                                                
        ENDDO                                                          
        SPACE = 1-NOSOIL-FRAC_MIN*(NPFT-1)                          
                                                                       
!----------------------------------------------------------------------    
! Calculate the increments to the tree fractions                           
!----------------------------------------------------------------------    
        N = DOM(1)                                                     
        M = DOM(2)                                                     
                                                                        
        FRACN=FRAC(N)                                                  
        FRACN=MAX(FRACN,FRAC_SEED)                                      
                                                                        
        FRACM=FRAC(M)                                                  
        FRACM=MAX(FRACM,FRAC_SEED)                                     
                                                                        
        P1 = GAMMA/FRACN-FORW*DB_DFRAC(N,N)                            
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(M,M)                            
        Q1 = -FORW*DB_DFRAC(N,M)                                       
        Q2 = -FORW*DB_DFRAC(M,N)                                       
        R1 = B(N)                                                      
        R2 = B(M)                                                      
        DO K=1,NPFT                                                     
          R1 = R1+FORW*(DB_DFRAC(N,K)*DFRAC(K))                        
          R2 = R2+FORW*(DB_DFRAC(M,K)*DFRAC(K))                        
        ENDDO                                                           
                                                                        
        NUMER = R1-(Q1/P2)*R2                                          
        DENOM = P1-(Q1/P2)*Q2                                          
        DENOM = MAX(DENOM,DENOM_MIN)                                    
        DFRAC(N) = NUMER/DENOM                                         
        FRAC(N) = FRAC(N)+DFRAC(N)                                   
                                                                        
        IF (FRAC(N).LT.FRAC_MIN) THEN                                  
          DFRAC(N) = DFRAC(N)+(FRAC_MIN-FRAC(N))                     
          FRAC(N) = FRAC_MIN                                           
        ELSEIF (FRAC(N).GT.SPACE) THEN                                
          DFRAC(N) = DFRAC(N)+(SPACE-FRAC(N))                     
          FRAC(N) = SPACE                                             
        ENDIF                                                           
                                                                        
        SPACE = SPACE-FRAC(N)+FRAC_MIN                             
                                                                        
        NUMER = R2-Q2*DFRAC(N)                                         
        DENOM = P2                                                      
        DENOM = MAX(DENOM,DENOM_MIN)                                    
        DFRAC(M) = NUMER/DENOM                                         
        FRAC(M) = FRAC(M)+DFRAC(M)                                   
                                                                        
        IF (FRAC(M).LT.FRAC_MIN) THEN                                  
          DFRAC(M) = DFRAC(M)+(FRAC_MIN-FRAC(M))                     
          FRAC(M) = FRAC_MIN                                           
        ELSEIF (FRAC(M).GT.SPACE) THEN                                
          DFRAC(M) = DFRAC(M)+(SPACE-FRAC(M))                     
          FRAC(M) = SPACE                                             
        ENDIF                                                           
                                                                        
        SPACE = SPACE-FRAC(M)+FRAC_MIN                             
                                                                        
!----------------------------------------------------------------------    
!!! Calculate the increment to the shrub fraction                            
! Calculate the increment to the third kind of tree
!----------------------------------------------------------------------    
        N = DOM(3)                                                       
                                                                        
        FRACN=FRAC(N)                                                   
        FRACN=MAX(FRACN,FRAC_SEED)                                        
                                                                       
        DENOM = GAMMA/FRACN-FORW*DB_DFRAC(N,N)                          
        DENOM = MAX(DENOM,DENOM_MIN)                                       
                                                                        
        NUMER = B(N)                                                     
        DO K=1,NPFT                                                        
          NUMER = NUMER+FORW*(DB_DFRAC(N,K)*DFRAC(K))                  
        ENDDO                                                              
                                                                       
        DFRAC(N) = NUMER/DENOM                                           
        FRAC(N) = FRAC(N)+DFRAC(N)                                   
                                                                      
        IF (FRAC(N).LT.FRAC_MIN) THEN                                    
          DFRAC(N) = DFRAC(N)+(FRAC_MIN-FRAC(N))                     
          FRAC(N) = FRAC_MIN                                             
        ELSEIF (FRAC(N).GT.SPACE) THEN                                
          DFRAC(N) = DFRAC(N)+(SPACE-FRAC(N))                     
          FRAC(N) = SPACE                                             
        ENDIF                                                              
                                                                       
        SPACE = SPACE-FRAC(N)+FRAC_MIN                             
                                                                        
!----------------------------------------------------------------------    
! Calculate the increment to the grass fraction                            
!----------------------------------------------------------------------    
                                                                        
        N = DOM(6)
        M = DOM(7)
                                                                        
        FRACN=FRAC(N)                                                  
        FRACN=MAX(FRACN,FRAC_SEED)                                      
                                                                        
        FRACM=FRAC(M)                                                  
        FRACM=MAX(FRACM,FRAC_SEED)                                      
                                                                        
        P1 = GAMMA/FRACN-FORW*DB_DFRAC(N,N)                            
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(M,M)                            
        Q1 = -FORW*DB_DFRAC(N,M)                                       
        Q2 = -FORW*DB_DFRAC(M,N)                                       
        R1 = B(N)                                                      
        R2 = B(M)                                                      
        DO K=1,NPFT                                                     
          R1 = R1+FORW*(DB_DFRAC(N,K)*DFRAC(K))                        
          R2 = R2+FORW*(DB_DFRAC(M,K)*DFRAC(K))                        
        ENDDO                                                           
                                                                        
        NUMER = R1-(Q1/P2)*R2                                           
        DENOM = P1-(Q1/P2)*Q2                                           
        DENOM = MAX(DENOM,DENOM_MIN)                                    
        DFRAC(N) = NUMER/DENOM                                         
        FRAC(N) = FRAC(N)+DFRAC(N)                                   
                                                                        
        IF (FRAC(N).LT.FRAC_MIN) THEN                                  
          DFRAC(N) = DFRAC(N)+(FRAC_MIN-FRAC(N))                     
          FRAC(N) = FRAC_MIN                                            
        ELSEIF (FRAC(N).GT.SPACE) THEN                                
          DFRAC(N) = DFRAC(N)+(SPACE-FRAC(N))                     
          FRAC(N) = SPACE                                             
        ENDIF                                                           
                                                                        
        SPACE = SPACE-FRAC(N)+FRAC_MIN                             
                                                                        
        NUMER = R2-Q2*DFRAC(N)                                         
        DENOM = P2                                                      
        DENOM = MAX(DENOM,DENOM_MIN)                                    
        DFRAC(M) = NUMER/DENOM                                         
        FRAC(M) = FRAC(M)+DFRAC(M)                                   
                                                                        
        IF (FRAC(M).LT.FRAC_MIN) THEN                                  
          DFRAC(M) = DFRAC(M)+(FRAC_MIN-FRAC(M))                     
          FRAC(M) = FRAC_MIN                                           
        ELSEIF (FRAC(M).GT.SPACE) THEN                                
          DFRAC(M) = DFRAC(M)+(SPACE-FRAC(M))                     
          FRAC(M) = SPACE                                             
        ENDIF                                                           
                                                                        
        SPACE = SPACE-FRAC(M)+FRAC_MIN                             
                                                                        
!----------------------------------------------------------------------    
! Calculate the increments to the shrub fractions                          
!----------------------------------------------------------------------    
        N = DOM(4)
        M = DOM(5)
                                                                        
        FRACN=FRAC(N)                                                  
        FRACN=MAX(FRACN,FRAC_SEED)                                      
                                                                        
        FRACM=FRAC(M)                                                  
        FRACM=MAX(FRACM,FRAC_SEED)                                      
                                                                        
        P1 = GAMMA/FRACN-FORW*DB_DFRAC(N,N)                            
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(M,M)                            
        Q1 = -FORW*DB_DFRAC(N,M)                                       
        Q2 = -FORW*DB_DFRAC(M,N)                                       
        R1 = B(N)                                                      
        R2 = B(M)                                                      
        DO K=1,NPFT                                                     
          R1 = R1+FORW*(DB_DFRAC(N,K)*DFRAC(K))                        
          R2 = R2+FORW*(DB_DFRAC(M,K)*DFRAC(K))                        
        ENDDO                                                           
                                                                        
        NUMER = R1-(Q1/P2)*R2                                           
        DENOM = P1-(Q1/P2)*Q2                                           
        DENOM = MAX(DENOM,DENOM_MIN)                                    
        DFRAC(N) = NUMER/DENOM                                         
        FRAC(N) = FRAC(N)+DFRAC(N)                                   
                                                                       
        IF (FRAC(N).LT.FRAC_MIN) THEN                                   
          DFRAC(N) = DFRAC(N)+(FRAC_MIN-FRAC(N))                     
          FRAC(N) = FRAC_MIN                                           
        ELSEIF (FRAC(N).GT.SPACE) THEN                                
          DFRAC(N) = DFRAC(N)+(SPACE-FRAC(N))                     
          FRAC(N) = SPACE                                             
        ENDIF                                                           
                                                                        
        SPACE = SPACE-FRAC(N)+FRAC_MIN                             
                                                                        
        NUMER = R2-Q2*DFRAC(N)                                          
        DENOM = P2                                                      
        DENOM = MAX(DENOM,DENOM_MIN)                                    
        DFRAC(M) = NUMER/DENOM                                          
        FRAC(M) = FRAC(M)+DFRAC(M)                                   
                                                                        
        IF (FRAC(M).LT.FRAC_MIN) THEN                                  
          DFRAC(M) = DFRAC(M)+(FRAC_MIN-FRAC(M))                     
          FRAC(M) = FRAC_MIN                                           
        ELSEIF (FRAC(M).GT.SPACE) THEN                                
          DFRAC(M) = DFRAC(M)+(SPACE-FRAC(M))                     
          FRAC(M) = SPACE                                             
        ENDIF                                                           
                                                                       
        SPACE = SPACE-FRAC(M)+FRAC_MIN                             
                                                                       
!----------------------------------------------------------------------     
! Diagnose the new bare soil fraction                                       
!----------------------------------------------------------------------     
        FRAC(SOIL) = 1.0-NOSOIL                                        
        DO N=1,NPFT                                                     
          FRAC(SOIL) = FRAC(SOIL)-FRAC(N)                             
        ENDDO                                                          
                                                                       
      END SUBROUTINE                                                             
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine DECAY --------------------------------------------------    
!!!                                                                        
!!! Purpose : Updates carbon contents of the soil.                         
!!!                                                                        
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.      
!!!                    Richard Betts                                        
!!!                                                                        
!!!END ----------------------------------------------------------------    
      SUBROUTINE DECAY (DPC_DCS,FORW,GAMMA,PC,CS)                      
                                                                        
      REAL::                    &
       DPC_DCS,                 & !IN Rate of change of PC with soil carbon (yr).                   
       FORW,                    & ! IN Forward timestep weighting.         
       GAMMA,                   & ! IN Inverse timestep (/360days).         
       PC,                      & ! IN Net carbon flux into the soil (kg C/m2/360days).              
       CS,                      & ! INOUT Soil carbon (kg C/m2).           
       DENOM,                   & ! WORK Denominator of update equation.                         
       NUMER                      ! WORK Numerator of the update equation.                         
!----------------------------------------------------------------------    
! Local parameters                                                         
!----------------------------------------------------------------------    
      REAL,PARAMETER::          &
       CS_MIN = 1.0E-6            ! Minimum soil carbon (kg C/m2).         

      REAL,PARAMETER::          &
       DENOM_MIN = 1.0E-6         ! Minimum value for the denominator     
!                                 ! of the update equation. Ensures      
!                                 ! that gradient descent does not      
!                                 ! lead to an unstable solution.      

      NUMER = PC                                                      
      DENOM = GAMMA+FORW*DPC_DCS                                      
      DENOM = MAX(DENOM,DENOM_MIN)                                     
                                                                        
      CS = CS+NUMER/DENOM                                          
      CS = MAX(CS_MIN,CS)                                          
                                                                        
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine GROWTH -------------------------------------------------    
!!!                                                                        
!!! Purpose : Increments leaf, root and wood carbon.                       
!!!                                                                        
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.      
!!!                    Richard Betts                                        
!!!                                                                        
!!!END ----------------------------------------------------------------    
      SUBROUTINE GROWTH(N,DPCG_DLAI,FORW,GAMMA,PC_G,LEAF,ROOT,WOOD,   &
              DLEAF,DROOT,DWOOD)

      INTEGER::                 &
       N                          ! IN Vegetation type.                    
                                                                        
      REAL::                    &
       DPCG_DLAI,               & ! IN Rate of change of PC_G with         
!                                 !    leaf area index (kg C/m2/360days/LAI).               
       FORW,                    & ! IN Forward timestep weighting.         
       GAMMA,                   & ! IN Inverse timestep (/360days).         
       PC_G,                    & ! IN Net carbon flux available           
!                                 !    for growth (kg C/m2/360days).        
       LEAF,                    & ! INOUT Leaf biomass (kg C/m2).          
       ROOT,                    & ! INOUT Root biomass (kg C/m2).          
       WOOD                       ! INOUT Woody biomass (kg C/m2).         
!                                                                        
      REAL::                    &
       DENOM,                   & ! WORK Denominator of update equation.                         
       DLEAF,DROOT,DWOOD,       & ! WORK Increments to leaf, root          
!                                 !      and woody biomass (kg C/m2).      
       DL_DW,                   & ! WORK Rate of change of leaf            
!                                 !      carbon with wood carbon.          
       DLAI_DW,                 & ! WORK Rate of change of leaf area       
!                                 !      index with wood carbon (LAI m2/kg C).                    
       DR_DW,                   & ! WORK Rate of change of root            
!                                 !      carbon with wood carbon.          
       NUMER,                   & ! WORK Numerator of the update equation.                         
       WOOD_MAX,                & ! WORK Maximum wood carbon (kg C/m2).    
       WOOD_MIN                   ! WORK Minimum wood carbon (kg C/m2).    
!
       real,parameter:: DENOM_MIN = 1.0e-6

!----------------------------------------------------------------------    
! Calculate the increment to the wood carbon                               
!----------------------------------------------------------------------    
        DL_DW = LEAF/(B_WL(N)*WOOD)                                  
        DR_DW = DL_DW                                                   
        DLAI_DW = DL_DW/SIGL(N)                                         
                                                                        
        NUMER = PC_G                                                    
        DENOM = (1+DL_DW+DR_DW)*GAMMA-FORW*DLAI_DW*DPCG_DLAI            
        DENOM = MAX(DENOM,DENOM_MIN)                                    
                                                                        
        DWOOD = NUMER/DENOM                                             
                                                                        
!----------------------------------------------------------------------    
! Ensure that the local leaf area index does not drop below its            
! minimum value or exceed its maximum value.                               
!----------------------------------------------------------------------    
        WOOD_MIN =A_WL(N)*LAI_MIN(N)**B_WL(N)
        WOOD_MAX = A_WL(N)*LAI_MAX(N)**B_WL(N)                          
        DWOOD = MAX((WOOD_MIN-WOOD),DWOOD)                              
        DWOOD = MIN((WOOD_MAX-WOOD),DWOOD)                              
                                                                        
!----------------------------------------------------------------------    
! Diagnose the increments to leaf and root carbon                          
!----------------------------------------------------------------------    
        DLEAF = SIGL(N)*((WOOD+DWOOD)/A_WL(N))**(1.0/B_WL(N))  &        
               -LEAF                                                  
        DROOT = DLEAF                                                   
                                                                        
!----------------------------------------------------------------------    
! Update carbon contents                                                   
!----------------------------------------------------------------------    
        LEAF = LEAF+DLEAF                                            
        ROOT = ROOT+DROOT                                            
        WOOD = WOOD+DWOOD                                            

      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!***********************************************************************   
! Calculates the leaf turnover rate as a function of temperature and       
! soil water availability                                                  
!***********************************************************************   
      SUBROUTINE LEAF_LIT(N, FSMC,TSTAR, G_LEAF)                       
!
      INTEGER::                 &
       N                          ! IN Plant functional type.              
                                                                        
      REAL, DIMENSION(NPFT)::   &
       FSMC,                    & ! IN Soil moisture availability factor.                                        
       TSTAR                      ! IN Surface temperature (K).                       

      REAL::                    &
       G_LEAF,                  & ! OUT Rate of leaf turnover(/360days).                                    
       FM,FT                      ! WORK Soil moisture and leaf                       
                                  !  temperature amplifiers of leaf turnover.                               
                                                                        
!-----------------------------------------------------------------------   
! Calculate the leaf turnover rate                                         
!-----------------------------------------------------------------------   
                                                                        
      FT = 1.0                                                         
      FM = 1.0                                                         

      IF (TSTAR(N) .LT. TLEAF_OF(N)) THEN                              
          FT = 1.0 + DGL_DT(N)*(TLEAF_OF(N)-TSTAR(N))                   
      ENDIF
!
      IF (FSMC(N) .LT. FSMC_OF(N)) THEN                              
          FM = 1.0 + DGL_DM(N)*(FSMC_OF(N)-FSMC(N))                     
      ENDIF                                                            
                                                                        
      G_LEAF = G_LEAF_0(N)*FT*FM                                      
                                                                        
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine LOTKA --------------------------------------------------    
!!!                                                                        
!!! Purpose : Updates fractional coverage of each functional type.         
!!!           Based on the Lotka-Volterra equations of interspecies        
!!!           competition.                                                 
!!!                                                                        
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!  4.5  12/05/98     Operate only on points indexed with TRIF_INDEX      
!!!                    and correct calculation of NOSOIL.  Richard Betts   
!!!                                                                        
!!!END ----------------------------------------------------------------    
      SUBROUTINE LOTKA (C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI,PC_S, &
                      FRAC,DFRAC)                                     
                                                                        
      INTEGER,PARAMETER::          &
       NNVG=4,SOIL=9                 ! Number of plant functional types.      
!
      REAL,PARAMETER::             &
       FRAC_MIN = 1.0E-6,          & ! Minimum areal fraction for PFTs.        
       FRAC_SEED = 0.01
!
      INTEGER::                    &
       K,M,N                         !    which TRIFFID may operate           
                                     ! WORK Loop counters.                    
!
      INTEGER,DIMENSION(NPFT)::    &
       DOM                           ! WORK Dominance hierachy.               
!                                                                        
      REAL,DIMENSION(NPFT)::       &
       C_VEG                         ! IN Carbon content of vegetation        
                                     !    (kg C/m2).                          
!
      REAL::                       &
       FORW,                       & ! IN Forward timestep weighting.         
       FRAC_VS,                    & ! IN Total fractional cover of           
!                                    !    vegetation and soil.                
       FRAC_AGRIC,                 & ! IN Fraction of agriculture.            
       GAMMA                         ! IN Inverse timestep (/360days).        
!
      REAL,DIMENSION(NPFT)::       &
       LAI,                        & ! IN Leaf area index.                    
       PC_S                          ! IN Net carbon flux available for       
!                                    !    spreading (kg C/m2/360days).        
!
      REAL,DIMENSION(NTYPES)::     &
       FRAC                          ! INOUT Fractional cover of each         
                                     !       Functional Type.                 
!
      REAL,DIMENSION(NPFT)::       &
       DFRAC,                      & ! OUT Increment to the areal fraction    
!                                    !     during the timestep (/timestep).   
       B                             ! WORK Mean rate of change of            
                                     !      vegetation fraction over          
                                     !      the timestep (kg C/m2/360days).   
!
      REAL,DIMENSION(NPFT,NPFT)::  & 
       DB_DFRAC,                   & ! WORK Rate of change of B               
!                                    !      with vegetation fraction.         
       COM                           ! WORK Coefficients representing         
                                     !      the influence of one type         
                                     !      (second argument) on another      
                                     !      (first argument).                 
!
      REAL::                       &
       DIFF_SUM,                   & ! WORK Difference divided by sum         
!                                    !      for competing canopy heights.     
       HC1,HC2,HC3,HC4,HC5,HC6,HC7,& ! WORK Competing canopy heights (m).     
       NOSOIL                        ! WORK Fractional area not available     
                                     !      to vegetation.                    
!
      REAL,DIMENSION(NPFT)::       &
       SPACE                         ! WORK Space available for invasion.     
!
      REAL,PARAMETER::             &
       POW  =20.0                    ! Power in sigmoidal function.           

!----------------------------------------------------------------------    
! Define competition coefficients and the dominance hierachy               
!----------------------------------------------------------------------    
                                                                        
      DO N=1,NPFT                                                       
        DO M=1,NPFT                                                     
            COM(N,M) = 1.0                                             
        ENDDO                                                           
      ENDDO                                                             
                                                                        
                                                                        
      HC1 = A_WL(1)/(A_WS(1)*ETA_SL(1))*(LAI(1)**(B_WL(1)-1))          
      HC2 = A_WL(2)/(A_WS(2)*ETA_SL(2))*(LAI(2)**(B_WL(2)-1))          
      DIFF_SUM = (HC1-HC2)/(HC1+HC2)                                    
                                                                        
      COM(1,2) = 1.0/(1+EXP(POW*DIFF_SUM))    ! BT vs NT               
      COM(1,3) = 0.0                          ! BT vs C3G              
      COM(1,4) = 0.0                          ! BT vs C4G              
      COM(1,5) = 0.0                          ! BT vs C3S              
      COM(1,6) = 0.0                          ! BT vs C4S              
                                                                        
      COM(2,1) = 1.0-COM(1,2)                 ! NT vs BT               
      COM(2,3) = 0.0                          ! NT vs C3G              
      COM(2,4) = 0.0                          ! NT vs C4G              
      COM(2,5) = 0.0                          ! NT vs C3S              
      COM(2,6) = 0.0                          ! NT vs C4S              
                                                                        
      HC3 = A_WL(3)/(A_WS(3)*ETA_SL(3))*(LAI(3)**(B_WL(3)-1))          
      HC4 = A_WL(4)/(A_WS(4)*ETA_SL(4))*(LAI(4)**(B_WL(4)-1))          
      DIFF_SUM = (HC3-HC4)/(HC3+HC4)                                   
                                                                        
      COM(3,4) = 1.0/(1+EXP(POW*DIFF_SUM))    ! C3G vs C4G             
      COM(4,3) = 1.0-COM(3,4)                 ! C4G vs C3G             
                                                                        
      HC5 = A_WL(5)/(A_WS(5)*ETA_SL(5))*(LAI(5)**(B_WL(5)-1))          
      HC6 = A_WL(6)/(A_WS(6)*ETA_SL(6))*(LAI(6)**(B_WL(6)-1))          

      DIFF_SUM = (HC5-HC6)/(HC5+HC6)                                   

      COM(5,6)=1/(1+EXP(POW*DIFF_SUM))        ! C3S VS C4S 
      COM(6,5)=1.0-COM(5,6)                   ! C4S VS C3S
      COM(5,3) = 0.0                          ! C3S vs C3G             
      COM(5,4) = 0.0                          ! C3S vs C4G             
      COM(6,3) = 0.0                          ! C4S vs C3G             
      COM(6,4) = 0.0                          ! C4S vs C4G             
!
!
      DIFF_SUM = (HC3-HC5)/(HC3+HC5)                                   
      COM(3,5) = 1/(1+EXP(POW*DIFF_SUM))
      COM(5,3) =1-COM(3,5)
!
      DIFF_SUM = (HC3-HC6)/(HC3+HC6)                                    
      COM(3,6) = 1/(1+EXP(POW*DIFF_SUM))
      COM(6,3) =1-COM(3,6)
!
      HC7 = A_WL(7)/(A_WS(7)*ETA_SL(7))*(LAI(7)**(B_WL(7)-1))
!
      DIFF_SUM = (HC1-HC7)/(HC1+HC7)
      COM(1,7) = 1.0/(1+EXP(POW*DIFF_SUM))
      COM(7,1) = 1-COM(1,7)
!
      DIFF_SUM = (HC2-HC7)/(HC2+HC7)
      COM(2,7) = 1.0/(1+EXP(POW*DIFF_SUM))
      COM(7,2) = 1-COM(2,7)
!
      COM(7,3) = 0.0
      COM(7,4) = 0.0
      COM(7,5) = 0.0
      COM(7,6) = 0.0

      IF (HC1 .GE. HC2 .and. HC1.GE.HC7) THEN
         DOM(1) = 1
         IF (HC2 .GE. HC7) THEN
             DOM(2) = 2
             DOM(3) = 7
         ELSEIF ( HC2 .LT. HC7) THEN
             DOM(2) = 7
             DOM(3) = 2
         ENDIF
      ELSEIF (HC2 .GE. HC1 .and. HC2.GE.HC7) THEN
         DOM(1) = 2
         IF (HC1 .GE. HC7) THEN
             DOM(2) = 1
             DOM(3) = 7
         ELSEIF ( HC1 .LT. HC7) THEN
             DOM(2) = 7
             DOM(3) = 1
         ENDIF
      ELSEIF (HC7 .GE. HC1 .and. HC7.GE.HC2) THEN
         DOM(1) = 7
         IF (HC1 .GE. HC2) THEN
             DOM(2) = 1
             DOM(3) = 2
         ELSEIF ( HC1 .LT. HC2) THEN
             DOM(2) = 2
             DOM(3) = 1
         ENDIF
      ENDIF
!
      IF (HC5 .GE. HC6) THEN
          DOM(4) = 5
          DOM(5) = 6
      ELSE
          DOM(4) = 6
          DOM(5) = 5
      ENDIF
!      
      IF (HC3 .GE. HC4) THEN
          DOM(6) = 3
          DOM(7) = 4
      ELSE
          DOM(6) = 4
          DOM(7) = 3
      ENDIF
!
!----------------------------------------------------------------------    
! Calculate the space available for the expansion of each FT               
!----------------------------------------------------------------------    
      NOSOIL = 1 - FRAC_VS                                         
                                                                        
!----------------------------------------------------------------------    
! Exclude non-crop types from agricultural regions                         
!----------------------------------------------------------------------    
      DO K=1,NPFT                                                       
          N=DOM(K)                                                     
          SPACE(N)=1.0-NOSOIL-FRAC_AGRIC*(1-CROP(N)) &             
                                  -FRAC_MIN*(NPFT-K)                   
      ENDDO                                                             
                                                                        
      DO N=1,NPFT                                                       
        DO M=1,NPFT                                                     
            SPACE(N)=SPACE(N)-COM(N,M)*FRAC(M)                     
        ENDDO                                                           
      ENDDO                                                             
                                                                        
!----------------------------------------------------------------------    
! Calculate the variables required for the implicit calculation.           
! Divide the update equation by FRAC to eliminate the (unstable)           
! bare soil solution.                                                      
!----------------------------------------------------------------------    
      DO N=1,NPFT                                                       
          B(N) = PC_S(N)*SPACE(N)/C_VEG(N)-G_AREA(N)

      if (N.eq.1) then
         if ((FRAC(3).GT.0.3).OR.(FRAC(4).GT.0.3).OR.((FRAC(3)+FRAC(4)).GT.0.5)) then
           if (FRAC(N).GT.0.10) then
             B(N) = PC_S(N)*SPACE(N)/C_VEG(N)-0.03
           endif
         endif
      endif
      if (N.eq.7) then
         if ((FRAC(3).GT.0.3).OR.(FRAC(4).GT.0.3).OR.((FRAC(3)+FRAC(4)).GT.0.5)) then
           if (FRAC(N).GT.0.10) then
             B(N) = PC_S(N)*SPACE(N)/C_VEG(N)-0.04
           endif
         endif
      endif
!
      if (N.eq.2) then
         if ((FRAC(3).GT.0.3).OR.(FRAC(4).GT.0.3).OR.(FRAC(6).GT.0.3).OR.((FRAC(3)+FRAC(4)+FRAC(6)).GT.0.5)) then
           if (FRAC(N).GT.0.05) then
           B(N) = PC_S(N)*SPACE(N)/C_VEG(N)-0.02
           endif
         endif
      endif
!
          DO M=1,NPFT                                                   
            DB_DFRAC(N,M) = -COM(N,M)*PC_S(N)/C_VEG(N)             
          ENDDO                                                        
      ENDDO                                                            
                                                                       
!----------------------------------------------------------------------    
! Update the areal fractions                                               
!----------------------------------------------------------------------    
      CALL COMPETE (DOM,B,DB_DFRAC,FORW,GAMMA,NOSOIL,FRAC,DFRAC)        
                                                                        
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!***********************************************************************   
! Calculates the soil respiration based on a simplified version of the     
! model of Raich et al. (1991).                                            
!***********************************************************************   
      SUBROUTINE MICROBE (CS,STH_SOIL,V_SAT,V_WILT,TSOIL,RESP_S)        
                                                                        
      REAL::                    &
       CS,                      & ! IN Soil carbon (kg C/m2).              
       STH_SOIL,                & ! IN Top layer soil moisture as a        
                                  !    fraction of saturation (m3/m3).     
       V_SAT,                   & ! IN Volumetric soil moisture            
                                  !    concentration at saturation         
                                  !    (m3 H2O/m3 soil).                   
       V_WILT,                  & ! IN Volumetric soil moisture            
                                  !    concentration below which           
                                  !    stomata close (m3 H2O/m3 soil).     
                                  !    as a fraction of saturation.        
       TSOIL,                   & ! IN Soil temperature (K).               
       RESP_S,                  & ! OUT Soil respiration (kg C/m2/s).      
       FSTH,FTEMP,              & ! WORK Factors describing the            
                                  !      influence of soil moisture and    
                                  !      soil temperature respectively     
                                  !      on the soil respiration.          
       STH_OPT                    ! WORK Fractional soil moisture at       
                                  !      which respiration is maximum.     
                                                                       
!-----------------------------------------------------------------------   
! Local parameters                                                         
!-----------------------------------------------------------------------   
      REAL,PARAMETER::          &
       KAPS = 0.5E-8,           & ! Specific soil respiration rate         
                                  ! at 25 deg ! and optimum soil moisture (/s).                         
       Q10 = 2.0                  ! Q10 factor for soil respiration.       

        IF (V_SAT .GT. 0.0) THEN  
                                                                        
!          STH_WILT = V_WILT / V_SAT                                 
          STH_OPT = 0.5 * (1 +V_WILT)    
                                                                        
          IF (STH_SOIL .LE. V_WILT) THEN                             
            FSTH = 0.2                                                  
          ELSEIF (STH_SOIL .GT. V_WILT .AND.        &                         
                  STH_SOIL .LE. STH_OPT) THEN                           
            FSTH = 0.2 + 0.8 * ((STH_SOIL - V_WILT) &                 
                              / (STH_OPT - V_WILT))                     
          ELSEIF (STH_SOIL .GT. STH_OPT) THEN                           
            FSTH = 1 - 0.8 * (STH_SOIL - STH_OPT)                       
          ENDIF                                                         
                                                                        
          FTEMP = Q10 ** (0.1 * (TSOIL - 298.15))                       
                                                                        
          RESP_S = KAPS * CS * FSTH * FTEMP                          
                                                                        
        ELSE                                                           
                                                                       
          RESP_S = 0.0                                                  
                                                                       
        ENDIF                                                          
                                                                       
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine PHENOL -------------------------------------------------    
!!!                                                                        
!!!                                                                        
!!! Purpose :  Parametrizes leaf phenological changes and updates the      
!!!            leaf area index and the leaf turnover rate.                 
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!                                                                        
!!!END -----------------------------------------------------------------   
      SUBROUTINE PHENOL (N,G_LEAF,HT,DTIME_PHEN,G_LEAF_PHEN,LAI)    

      INTEGER::                 &
       N                          ! IN Plant functional type.              
                                                                        
      REAL::                    &
       G_LEAF,                  & ! IN Rate of leaf turnover (/360days).   
       HT,                      & ! IN Canopy height (m).                  
       DTIME_PHEN,              & ! IN Timestep (years).                   
       G_LEAF_PHEN,             & ! OUT Rate of leaf turnover              
!                                 !     including leaf phenology           
!                                 !     (/360days).                        
       LAI,                     & ! INOUT Leaf area index.                 
       DPHEN,                   & ! WORK Increment to phenological         
!                                 !      state.                            
       LAI_BAL,                 & ! WORK Balanced growth LAI.              
       PHEN                       ! WORK Phenological state.               
!
       real,dimension(NPFT):: xphen_max
       data xphen_max/0.05,0.05,0.15,0.15,0.05,0.15,0.05/

!-----------------------------------------------------------------------   
! Diagnose the phenological state                                          
!-----------------------------------------------------------------------   
        LAI_BAL = (A_WS(N)*ETA_SL(N)*HT          &                   
                     /A_WL(N))**(1.0/(B_WL(N)-1))                       
        PHEN = LAI/LAI_BAL                                        
                                                                        
!-----------------------------------------------------------------------   
! Update the phenological state and output the leaf turnover rate in       
! terms of the balanced growth LAI                                         
!-----------------------------------------------------------------------   
                                                                        
        IF (G_LEAF.GT.2*G_LEAF_0(N)) THEN                               
          DPHEN = -DTIME_PHEN*G_GROW(N)                                 
          DPHEN = MAX(DPHEN,(xphen_max(n)-PHEN))                        
          DPHEN = MIN(DPHEN,0.)
          G_LEAF_PHEN = -DPHEN/DTIME_PHEN                               
        ELSE                                                           
          DPHEN = DTIME_PHEN*G_GROW(N)*(1.0-PHEN)                       
          DPHEN = MIN(DPHEN,(1.0-PHEN))                                 
          G_LEAF_PHEN = PHEN*G_LEAF                               
        ENDIF                                                           
                                                                        
!-----------------------------------------------------------------------   
! Update the leaf area index                                               
!-----------------------------------------------------------------------   
        PHEN = PHEN + DPHEN                                          
        LAI = PHEN*LAI_BAL                                        
                                                                       
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine SOILCARB -----------------------------------------------    
!!!                                                                        
!!! Purpose : Updates carbon contents of the soil.                         
!!!                                                                        
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.     
!!!                    Richard Betts                                       
!!!                                                                        
!!!END ----------------------------------------------------------------    
      SUBROUTINE SOILCARB (FORW,GAMMA,LIT_C_T,RESP_S,CS)                
                                                                        
      REAL::                    &
       FORW,                    & ! IN Forward timestep weighting.         
       GAMMA,                   & ! IN Inverse timestep (/360days).        
       LIT_C_T,                 & ! IN Total carbon litter                 
!                                 !    (kg C/m2/360days).                  
       RESP_S,                  & ! INOUT Soil respiration                 
!                                 !    (kg C/m2/360days).                  
       CS,                      & ! INOUT Soil carbon (kg C/m2).           
       DCS,                     & ! WORK Increment to the soil carbon      
!                                 !      (kg C/m2).                        
       DPC_DCS,                 & ! WORK Rate of change of PC with         
!                                 !      soil carbon (/360days).           
       PC                         ! WORK Net carbon accumulation in        
!                                 !      the soil (kg C/m2/360days).       
                                                                        
!----------------------------------------------------------------------    
! Diagnose the net local carbon flux into the soil                         
!----------------------------------------------------------------------    
        PC = LIT_C_T-RESP_S                                       
                                                                        
!----------------------------------------------------------------------    
! Variables required for the implicit and equilibrium calculations         
!----------------------------------------------------------------------    
        DPC_DCS = RESP_S/CS                                       
                                                                        
!----------------------------------------------------------------------    
! Save current value of soil carbon                                        
!----------------------------------------------------------------------    
        DCS = CS                                                     
                                                                        
!----------------------------------------------------------------------    
! Update soil carbon                                                       
!----------------------------------------------------------------------    
      CALL DECAY (DPC_DCS,FORW,GAMMA,PC,CS)                             
                                                                        
!----------------------------------------------------------------------    
! Apply implicit correction to the soil respiration rate.                  
!----------------------------------------------------------------------    
                                                                        
        DCS = CS - DCS                                            
        RESP_S = RESP_S + FORW*DPC_DCS*DCS                     
                                                                        
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine TRIFFID ------------------------------------------------    
!!!                                                                        
!!!                     Top-down                                           
!!!                     Representation of                                  
!!!                     Interactive                                        
!!!                     Foliage and                                        
!!!                     Flora                                              
!!!                     Including                                          
!!!                     Dynamics                                           
!!!                                                                        
!!! Purpose : Simulates changes in vegetation structure, areal             
!!!           coverage and the carbon contents of vegetation and soil.     
!!!           can be used to advance these variables dynamically           
!!!           (GAMMA=1/TIMESTEP) or to iterate towards  equilibrium        
!!!           (GAMMA --> 0.0, FORW=1.0).                                   
!!!                                                                        
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.     
!!!                    Richard Betts                                       
!!!                                                                        
!!!END ----------------------------------------------------------------    
      SUBROUTINE TRIFFID (FORW,GAMMA,                   & 
       FRAC_VS,FRAC_AGRIC,G_LEAF,NPP,RESP_S,RESP_W,     &              
       CS,FRAC,HT,LAI,C_VEG,CV,LIT_C,LIT_C_T,AGR,       &               
       LEAF_LOSS,WOOD_LOSS,ROOT_LOSS,                   &  !Huilin add for fired
       LEAF1,ROOT1,WOOD1,DCVEG1,PHEN,LAID,green,        &
       FPARn,KPARn )
!
      REAL,PARAMETER:: FRAC_MIN = 1.0E-6

      INTEGER::                 &
       N                          ! WORK Loop counters                     
!                                                                           
      REAL::                    &
       FORW,                    &  ! IN Forward timestep weighting.         
       FRAC_VS,                 & ! IN Total fraction of gridbox           
!                                 !    covered by veg or soil.             
       GAMMA,                   & ! IN Inverse timestep (/360days).        
       FRAC_AGRIC,              & ! IN Fraction of agriculture.            
       AGR                        ! IN Agriculture fraction
      
       REAL,DIMENSION(NPFT)::   &
       G_LEAF                     ! IN Turnover rate for leaf and          
!                                 !    fine root biomass (/360days).       
       REAL,DIMENSION(NPFT)::   &
       NPP,                     & ! INOUT Net primary productivity         
!                                 !       (kg C/m2/360days).               
       RESP_W,                  & ! INOUT Wood maintenance respiration     
!                                 !       (kg C/m2/360days).               
       HT,                      & ! INOUT Vegetation height (m).           
       LAI,                     & ! INOUT Leaf area index.                 
       LEAF_LOSS,               & ! INOUT ACCUMULATED leaf loss due to fire Huilin (kg C/m2/360 days) 
       WOOD_LOSS,               & ! INOUT ACCUMULATED wood loss due to fire Huilin (kg C/m2/360 days) 
       ROOT_LOSS                  ! INOUT ACCUMULATED root loss due to fire Huilin (kg C/m2/360 days)

      REAL::                    &
       RESP_S,                  & ! INOUT Soil respiration (kg C/m2/360days).               
       CS                         ! INOUT Soil carbon (kg C/m2). 

      REAL,DIMENSION(NTYPES)::  & 
       FRAC                       ! INOUT Fractional cover of each Functional Type.                 

      REAL,DIMENSION(NPFT)::    &
       LAID,                    & ! OUT Dead leaf area index.
       GREEN,                   & ! OUT green fraction.
       C_VEG,                   & ! OUT Total carbon content of            
!                                 !     the vegetation (kg C/m2).          
       LIT_C                      ! OUT Carbon Litter (kg C/m2/360days).   

      REAL::                    &
       CV,                      & ! OUT Gridbox mean vegetation carbon (kg C/m2).                  
       LIT_C_T                    ! OUT Gridbox mean carbon litter         
!                                 !     (kg C/m2/360days).                 
      REAL::                    &
       FRAC_FLUX                  ! WORK PFT fraction to be used           
!                                 !      in the calculation of             
!                                 !      the gridbox mean fluxes.          
      REAL,DIMENSION(NPFT)::    &
       DCVEG,                   & ! WORK Change in vegetation carbon       
!                                 !      during the timestep               
!                                 !      (kg C/m2/timestep).               
       DFRAC,                   & ! WORK Change in areal fraction          
!                                 !      during the timestep               
!                                 !      (/timestep).                      
       LAI_BAL,                 & ! WORK Leaf area index in balanced       
!                                 !      growth state.                     
       LEAF,                    & ! WORK Leaf biomass (kg C/m2).           
       LEAF1,                   & ! WORK Leaf biomass (kg C/m2).          
       PC_S,                    & ! WORK Net carbon flux available         
!                                 !      for spreading                     
!                                 !      (kg C/m2/yr).                     
       PHEN,                    & ! WORK Phenological state.               
       ROOT,                    & ! WORK Root biomass (kg C/m2).           
       WOOD,                    & ! WORK Woody biomass (kg C/m2).          
       ROOT1,                   & ! WORK Root biomass (kg C/m2).           
       WOOD1,                   & ! WORK Woody biomass (kg C/m2).          
       DCVEG1,                  & ! WORK Woody biomass (kg C/m2).          
       DLEAF,                   & ! WORK DLEAF due to NPP Nov 2018
       DWOOD,                   & ! WORK DLEAF due to NPP Nov 2018
       DROOT                      ! WORK DLEAF due to NPP Nov 2018

      REAL,DIMENSION(NPFT)::    &
       FPARn,                   & ! 
       KPARn                      ! 
      REAL,DIMENSION(NPFT)::    &
       DFRAC_FIRE                 ! WORK, Fire disturbance on Veg frac.
                                                                        
!----------------------------------------------------------------------    
! Loop through Functional Types                                            
!----------------------------------------------------------------------    
      DO N=1,NPFT                                                       
                                                                        
!----------------------------------------------------------------------    
! Diagnose the balanced-growth leaf area index and the associated leaf,    
! wood, root and total vegetation carbon                                   
!----------------------------------------------------------------------    
          LAI_BAL(N) = (A_WS(N)*ETA_SL(N)*HT(N)   &                     
                    /A_WL(N))**(1.0/(B_WL(N)-1))                       
          LEAF(N) = SIGL(N)*LAI_BAL(N)                                 
          ROOT(N) = LEAF(N)                                            
          WOOD(N) = A_WL(N)*(LAI_BAL(N)**B_WL(N))                      
          C_VEG(N) = LEAF(N) + ROOT(N) + WOOD(N)                   
!----------------------------------------------------------------------    
! Diagnose the phenological state                                          
!----------------------------------------------------------------------    
!        PHEN(N) = LAI(N)/LAI_BAL(N)                                
! LAI_BAL is calculated based on HT, contradict with carbon pool
! We use LAI_BAL calculated cased on LEAF1, different from LAI_BAL here 
! PHEN is input based on updates from PHENOL
!----------------------------------------------------------------------    
! Update vegetation carbon contents                                        
!----------------------------------------------------------------------    
        CALL VEGCARB (N,FORW,GAMMA,G_LEAF(N),NPP(N),RESP_W(N),  &
                      FPARn(N),KPARn(N),                        &
                      LEAF(N),ROOT(N),WOOD(N),                  &
                      DCVEG(N),PC_S(N))                                
 
        LEAF1(N)= LEAF1(N)+DLEAF(N)
        ROOT1(N)= ROOT1(N)+DROOT(N)
        WOOD1(N)= WOOD1(N)+DWOOD(N)
        DCVEG1(N)= DCVEG(N) 
      
      ENDDO                                                             
                                                                        
!-----------------------------------------------------------------------   
! Diagnose the new value of Canopy Height, Leaf Area Index and Total       
! Vegetation Carbon                                                        
!-----------------------------------------------------------------------   
      DO N=1,NPFT                                                       
                                                                        
          HT(N) = WOOD1(N) / (A_WS(N) * ETA_SL(N))      &               
                  * (A_WL(N)/WOOD1(N))**(1.0/B_WL(N))                   
          LAI_BAL(N) = LEAF1(N) / SIGL(N)                               
          LAI(N) = PHEN(N) * LAI_BAL(N)                              
          C_VEG(N) = LEAF1(N) + ROOT1(N) + WOOD1(N)                   
!
!
          LAID(N) = G_LEAF(N)*LAI_BAL(N)
          IF(LAID(N).LT.0.001) LAID(N) = 0.001
          green(N)=LAI(N)/(LAI(N)+2*LAID(N))
      ENDDO                                                             
!----------------------------------------------------------------------    
! Fire disturbence Ye && Huilin Mar.28 
!----------------------------------------------------------------------    
       DO N = 1, NPFT
           DFRAC_FIRE(N) = (LEAF_LOSS(N) + ROOT_LOSS(N) + WOOD_LOSS(N)) /(1.0-AGR) / C_VEG(N)*FRAC(N)
           DFRAC_FIRE(N) = min(DFRAC_FIRE(N), FRAC(N) - FRAC_MIN)
           FRAC(N)       = FRAC(N) - DFRAC_FIRE(N)
       ENDDO
!----------------------------------------------------------------------    
! Update the areal coverage of each functional type                        
!----------------------------------------------------------------------    
      CALL LOTKA (C_VEG,FORW,FRAC_VS,FRAC_AGRIC,GAMMA,LAI_BAL,PC_S, &
                  FRAC,DFRAC)                                           
                                                                        
!----------------------------------------------------------------------    
! Diagnose the litterfall from the carbon balance of each vegetation       
! type                                                                     
!----------------------------------------------------------------------    
                                                                        
        LIT_C_T = 0.0                                                   
                                                                        
        DO N=1,NPFT                                                     
          FRAC_FLUX=FRAC(N)-(1.0-FORW)*DFRAC(N)                        
          LIT_C(N) = NPP(N)-GAMMA/FRAC_FLUX*(C_VEG(N)*FRAC(N)  &    
                     -(C_VEG(N)-DCVEG(N))*(FRAC(N)-DFRAC(N)))      
          LIT_C_T = LIT_C_T+FRAC_FLUX*LIT_C(N)                     
        ENDDO                                                           

!----------------------------------------------------------------------    
! Call SOIL_C to update the soil carbon content                            
!----------------------------------------------------------------------    
      CALL SOILCARB (FORW,GAMMA,LIT_C_T,RESP_S,CS)                      
                                                                        
!----------------------------------------------------------------------    
! Diagnose the gridbox mean vegetation carbon                              
!----------------------------------------------------------------------    
        CV = 0.0                                                        
        DO N=1,NPFT                                                     
          CV = CV+FRAC(N)*C_VEG(N)                               
        ENDDO                                                          
                                                                       
      END SUBROUTINE
! *****************************COPYRIGHT******************************     
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.    
!                                                                          
! Use, duplication or disclosure of this code is subject to the            
! restrictions as set forth in the contract.                               
!                                                                          
!                Meteorological Office                                     
!                London Road                                               
!                BRACKNELL                                                 
!                Berkshire UK                                              
!                RG12 2SZ                                                  
!                                                                          
! If no contract has been raised with this copy of the code, the use,      
! duplication or disclosure of it is strictly prohibited.  Permission      
! to do so must first be obtained in writing from the Head of Numerical    
! Modelling at the above address.                                          
! ******************************COPYRIGHT******************************    
!!! Subroutine VEGCARB ------------------------------------------------    
!!!                                                                        
!!! Purpose : Updates carbon contents of the vegetation.                   
!!!                                                                        
!!!                                                                        
!!!  Model            Modification history:                                
!!! version  Date                                                          
!!!  4.4     10/97     New Deck. Peter Cox                                 
!!!  4.5   12/05/98    Operate only on points indexed with TRIF_INDEX.     
!!!                    Richard Betts                                       
!!!                                                                        
!!!END ----------------------------------------------------------------    
       SUBROUTINE VEGCARB (N,FORW,GAMMA,G_LEAF,NPP,RESP_W,  &
                           FPARn,KPARn,                     &
                           LEAF,ROOT,WOOD,                  &
                           DCVEG,PC_S)                   
                                                                        
!                                                                           
      INTEGER::                 &
       N                          ! IN Plant functional type.              
!                                                                           
      REAL::                    &
       FORW,                    & ! IN Forward timestep weighting.         
       GAMMA,                   & ! IN Inverse timestep (/360days).        
       G_LEAF,                  & ! IN Turnover rate for leaf and          
!                                 !    fine root biomass (/360days).       
       NPP,                     & ! INOUT Net primary productivity         
!                                 !       (kg C/m2/360days).               
       RESP_W,                  & ! INOUT Wood maintenance respiration     
!                                 !       (kg C/m2/360days).               
       LEAF,                    & ! INOUT Leaf biomass (kg C/m2).          
       ROOT,                    & ! INOUT Root biomass (kg C/m2).          
       WOOD,                    & ! INOUT Woody biomass (kg C/m2).         
       DCVEG,                   & ! OUT Change in vegetation carbon        
!                                 !     during the timestep                
!                                 !     (kg C/m2/timestep).                
       PC_S,                    & ! OUT Net carbon flux available          
!                                 !     for spreading (kg C/m2/360days).   
       DFPAR_DLAI,              & ! WORK Rate of change of FPAR            
!                                 !      with leaf area index.             
       DLAI,                    & ! WORK Increment to the leaf area        
!                                 !      index.                            
       DLAMG_DLAI,DLIT_DLAI,    & ! WORK Required for calculation          
!                                 !      of the equilibrium increments.    
       DNPP_DLAI,               & ! WORK Rate of change of NPP             
!                                 !      with leaf area index              
!                                 !      (kg C/m2/360days/LAI).            
       DPC_DLAI,                & ! WORK Rate of change of PC              
!                                 !      with leaf area index              
!                                 !      (kg C/m2/360days/LAI).            
       DPCG_DLAI,               & ! WORK Rate of change of PC_G            
!                                 !      with leaf area index              
!                                 !      (kg C/m2/360days/LAI).            
       DRESPW_DLAI,             & ! WORK Rate of change of RESP_W          
!                                 !      with leaf area index              
       FPAR,                    & ! WORK PAR interception factor.          
       LAI,                     & ! WORK Leaf area index.                  
       LAMBDA_G,                & ! WORK Fraction of NPP available         
!                                 !      for spreading.                    
       LIT_C_L,                 & ! WORK Local rate of Carbon Litter       
!                                 !      production (kg C/m2/360days).     
       PC,                      & ! WORK Net carbon flux available         
!                                 !      to vegetation (kg C/m2/360days)   
       PC_G                       ! WORK Net carbon flux available         
!                                 !      for growth (kg C/m2/360days).     
      REAL::                    &
       FPARn,                   &
       KPARn
!
        LAI = LEAF/SIGL(N)                                           

!----------------------------------------------------------------------    
! Calculate the local production rate for carbon litter                    
!----------------------------------------------------------------------    
        LIT_C_L = G_LEAF*LEAF+G_ROOT(N)*ROOT + G_WOOD(N)*WOOD

!----------------------------------------------------------------------    
! Diagnose the net local carbon flux into the vegetation                   
!----------------------------------------------------------------------    
        PC = NPP - LIT_C_L                                        
                                                                        
!----------------------------------------------------------------------    
! Variables required for the implicit and equilibrium calculations         
!----------------------------------------------------------------------    
        DLIT_DLAI = (G_LEAF*LEAF+G_ROOT(N)*ROOT)/LAI  &          
                  + B_WL(N)*G_WOOD(N)*WOOD/LAI                       
                                                                        
        FPAR = FPARn
        DFPAR_DLAI = EXP(-KPARn*LAI)                                  
                                                                        
        DNPP_DLAI = NPP*DFPAR_DLAI/FPAR + (1-R_GROW(N))*RESP_W   &
                     *(DFPAR_DLAI/FPAR-B_WL(N)/LAI)                     
                                                                        
        LAMBDA_G = 1 - (LAI - LAI_MIN(N)) /(LAI_MAX(N) - LAI_MIN(N))
        DLAMG_DLAI = -1.0/(LAI_MAX(N) - LAI_MIN(N))                    
                                                                       
        PC_G = LAMBDA_G * NPP - LIT_C_L                           
        DPCG_DLAI = LAMBDA_G*DNPP_DLAI   &                            
                     + DLAMG_DLAI*NPP - DLIT_DLAI
        DPC_DLAI = DNPP_DLAI - DLIT_DLAI                             
                                                                       
!----------------------------------------------------------------------    
! Update vegetation carbon contents                                        
!----------------------------------------------------------------------    
        DCVEG = LEAF+ROOT+WOOD                                 
                                                                        
        CALL GROWTH (N,DPCG_DLAI,FORW,GAMMA,PC_G,LEAF,ROOT,WOOD,   &
                    DLEAF,DROOT,DWOOD)        
                                                                        
        DCVEG = LEAF+ROOT+WOOD-DCVEG                        
                                                                       
!----------------------------------------------------------------------    
! Diagnose the carbon available for spreading and apply implicit           
! corrections to the driving fluxes.                                       
!----------------------------------------------------------------------    
        DLAI = LEAF/SIGL(N) - LAI                                    
        PC_S = PC + FORW*DPC_DLAI*DLAI - DCVEG*GAMMA           
                                                                        
        FPAR = FPARn
        DFPAR_DLAI = EXP(-KPARn*LAI)                                  
        DRESPW_DLAI = RESP_W*B_WL(N)/LAI                             
                                                                        
        NPP = NPP + FORW*DNPP_DLAI*DLAI                           
        RESP_W = RESP_W + FORW*DRESPW_DLAI*DLAI                      
                                                                        
      END SUBROUTINE
END MODULE
