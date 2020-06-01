MODULE  ssib4_module_comsveg
!========================================================================       
!    Zhengqiu Zhang        12 August 2013       
!========================================================================       
!
     use ssib4_module_comsconst
     use ssib4_module_veg9table
!
     implicit none
!    vegetation and soil parameters                                         
!------------------------------------------------------------------------
      integer, parameter :: NTYPES=NTYY
!yliu9May2017 for addtype      integer, parameter :: JBTYP=7        !bare soil index
      integer, parameter :: JBTYP=8        !bare soil index
      real, dimension(NTYPES,2,3,2):: tran0,ref0
      real, dimension(NTYPES,12,2) :: green0,vcover0
      real, dimension(NTYPES,2)    :: chil0,topt0,tl0,tu0,defac0,ph10,ph20
      real, dimension(NTYPES,2,3)  :: rstpar0
      real, dimension(NTYPES)      :: vm00, fc30
      real, dimension(NTYPES,12,2) :: zlt0
      real, dimension(NTYPES,12)   :: zz0,d0,z20,z10,rdc0,rbc0 
      real, dimension(NTYPES,2)    :: rootd0
      real, dimension(NTYPES,3)    :: soref0,depth0
      real, dimension(NTYPES)      :: bee0,phsat0,poros0,satco0,slope0 
!
!-----for bare soil--------------------------------------------------------
!
      real, dimension(2,3,2)       :: soil_tran0,soil_ref0
      real, dimension(12,2)        :: soil_green0,soil_vcover0,soil_zlt0
      real, dimension(2,3)         :: soil_rstpar0
      real, dimension(2)           :: soil_chil0
      real, dimension(2)           :: soil_topt0,soil_tl0,soil_tu0,soil_defac0,soil_ph10,soil_ph20 
      real                         :: soil_vm00,soil_fc30
      real, dimension(12)          :: soil_zz0,soil_d0,soil_z20,soil_z10,soil_rdc0,soil_rbc0
      real, dimension(2)           :: soil_rootd0
      real, dimension(3)           :: soil_soref0,soil_depth0
      real                         :: soil_bee0,soil_phsat0,soil_poros0,soil_satco0,soil_slope0 
!
!------------------------------------------------------------------------
!zzq...
      logical::                 &
        phen_on,                & ! .t. for interactive phenology.
        trif_on,                & ! .t. for dynamic vegetation. 
        veg_equil                 ! .t. if the vegetation equilibrium is required.
      integer::                 &
         day_phen,              & ! number of days between phenology.
         day_trif,              & ! number of days between triffid.
         step_day                 ! number of timesteps. in a day.
!zzq05-09-2008 
      real,dimension(NTYPES) :: &
         v_sat0,                & ! volumetric soil moisture concentration 
                                  !       at saturation (m3 h2o/m3 soil).
         v_wilt0                  ! volumetric soil moisture concentration below
                                  !       which stomata close (m3 h2o/m3 soil).
      real::                    &
         fraca,                 & ! areal fraction of agriculture.
         g_anth,                & ! anthropogenic disturbance rate (/yr).
         anetc_def_kgc            ! default anetc in kgc.  
!zzq
!------------------------------------------------------------------------       
contains
!=====================================================================
!
      SUBROUTINE VEGIN
!                                                         August 2000
!=====================================================================
!     READ VEGETATION PHYSIOLOGY
!
!    SURFACE PARAMETERS ARE READ IN SAME ORDER AS IN GCM
!    SUBROUTINE SIBINP. ONLY EXCEPTION IS THAT 1-D VERSION READS IN
!    SITE SPECIFIC PARAMETERS CORB1 ... ZMET .
!
!     VARIABLES THAT ENTER THROUGH COMSIB:
!        SUBSCRIPTS (IV, IW, IL) :
!              IV = VEGETATION STORY; 1 = TOP AND 2 = BOTTOM
!              IW = RADIATION WAVELENGTH; 1 = VISIBLE, 2 = NEAR
!                   INFRARED AND 3 = THERMAL INFRARED
!              IL = VEGETATION STATE; 1 = LIVE (GREEN) AND
!                   2 = DEAD (STEMS AND TRUNK)
!
!   TRAN(IV,IW,IL): LEAF TRANSMITTANCE
!   REF (IV,IW,IL): LEAF REFLECTANCE
!   RSTPAR(IV,IW) : PAR-DEPENDENT LEAF STOMATAL RESISTANCE COEFFICIENTS
!                          A =(J/M**3) B = 2(W/M**2) C = 3(S/M)
!   SOREF(IW)     : SOIL REFLECTANCE
!   CHIL(IV)      : LEAF ANGLE DISTRIBUTION FACTOR
!   TOPT(IV)      : OPTIMUM TEMPERATURE FOR STOMATAL FUNCTIONING
!   TL(IV)        : LOWER TEMPERATURE LIMIT FOR STOMATAL FUNCTIONING
!   TU(IV)        : UPPER TEMPERATURE LIMIT FOR STOMATAL FUNCTIONING
!   DEFAC(IV)     : VAPOR PRESSURE DEFICIT PARAMETER
!   PH1(IV)       :
!   PH2(IV)       :
!   ROOTD(IV)     : ROOTING DEPTH
!   BEE           : SOIL WETNESS EXPONENT
!   PHSAT         : SOIL TENSION AT SATURATION
!   SATCO         : HYDRAULIC CONDUCTIVITY AT SATURATION
!   POROS         : SOIL POROSITY
!   DEPTH         : DEPTH OF 3 SOIL MOISTURE LAYERS
!   Z0            : ROUGHNESS LENGTH
!   D             : ZERO PLANE DISPLACEMENT
!   ZLT(IV)       : LEAF AREA INDEX
!   GREEN(IV)     : GREEN LEAF FRACTION
!   VCOVER(IV)    : VEGETATION COVER FRACTION
!
!
!***  READ VEGETATION MORPHOLOGICAL AND PHYSIOLOGICAL DATA
!
!
      INTEGER:: JTYP,IV,IW,IL,K,IM,JMON,JCG
      real:: f0001

      DO 100 JTYP=1,NTYPES
         VM00(JTYP)=VM0x(JTYP)
         fc30(JTYP)=fc3x(JTYP)
         BEE0(JTYP)=BEEx(JTYP)
         PHSAT0(JTYP)= PHSATx(JTYP)
         SATCO0(JTYP)=SATCOx(JTYP)
         POROS0(JTYP)=POROSx(JTYP)
         SLOPE0(JTYP)= SLOPEx(JTYP)

         DO IV=1,2
           CHIL0(JTYP,IV)=CHILx(JTYP,IV)
           TOPT0(JTYP,IV)=TOPTx(JTYP,IV)
           TL0(JTYP,IV)=TLx(JTYP,IV)
           TU0(JTYP,IV)=TUx(JTYP,IV)
           DEFAC0(JTYP,IV)=DEFACx(JTYP,IV)
           PH10(JTYP,IV)=PH1x(JTYP,IV)
           PH20(JTYP,IV)= PH2x(JTYP,IV)
           ROOTD0(JTYP,IV)=ROOTDx(JTYP,IV)

           DO IW=1,3
             RSTPAR0(JTYP,IV,IW)=RSTPARx(JTYP,IV,IW)
             DO IL=1,2
               TRAN0(JTYP,IV,IW,IL)=TRANx(JTYP,IV,IW,IL)
               REF0(JTYP,IV,IW,IL)=REFx(JTYP,IV,IW,IL)
             ENDDO
           ENDDO
         ENDDO

         DO IW=1,3
           SOREF0(JTYP,IW)=SOREFx(JTYP,IW)
         ENDDO
!
         DO K=1,3
           DEPTH0(JTYP,K)=DEPTHx(JTYP,K)
         ENDDO
!
         DO IM=1,12
           DO IV=1,2
              ZLT0(JTYP,IM,IV)=ZLTx(JTYP,IM,IV)
              GREEN0(JTYP,IM,IV)=GREENx(JTYP,IM,IV)
              VCOVER0(JTYP,IM,IV)=VCOVERx(JTYP,IM,IV)
           ENDDO

           Z20(JTYP,IM)=Z2x(JTYP,IM)
           Z10(JTYP,IM)=Z1x(JTYP,IM)
           ZZ0(JTYP,IM)=Z0x(JTYP,IM)
           D0(JTYP,IM)=Dx(JTYP,IM)
           RBC0(JTYP,IM)=RBCx(JTYP,IM)
           RDC0(JTYP,IM)=RDCx(JTYP,IM)
         ENDDO
!        write(*,120) JTYP,(ZLT0(JTYP,IM,1)*GREEN0(JTYP,IM,1),im=1,12)
!  120     format('data zlt',I1,'/',                  &
!               f7.5,',',f7.5,',',f7.5,',',f7.5,',',  &
!               f7.5,',',f7.5,',',f7.5,',',f7.5,',',  &
!               f7.5,',',f7.5,',',f7.5,',',f7.5,'/')

  100    CONTINUE
!
      F0001=0.0001
      DO 110 JCG =1, 2
      DO 110 JMON=1,12
      DO 110 JTYP=1,NTYPES
        GREEN0(JTYP,JMON,JCG)=AMAX1(F0001,GREEN0(JTYP,JMON,JCG))
  110   CONTINUE
!
!
! save bare ground parameters
!
         soil_vm00=VM00(JBTYP)
         soil_fc30=fc30(JBTYP)
         soil_bee0=BEE0(JBTYP)
         soil_PHSAT0=PHSAT0(JBTYP)
         soil_SATCO0=SATCO0(JBTYP)
         soil_POROS0=POROS0(JBTYP)
         soil_SLOPE0=SLOPE0(JBTYP)

         DO IV=1,2
           soil_CHIL0(IV)=CHIL0(JBTYP,IV)
           soil_TOPT0(IV)=TOPT0(JBTYP,IV)
           soil_TL0(IV)=TL0(JBTYP,IV)
           soil_TU0(IV)=TU0(JBTYP,IV)
           soil_DEFAC0(IV)=DEFAC0(JBTYP,IV)
           soil_PH10(IV)=PH10(JBTYP,IV)
           soil_PH20(IV)=PH20(JBTYP,IV)
           soil_ROOTD0(IV)=ROOTD0(JBTYP,IV)

           DO IW=1,3
             soil_RSTPAR0(IV,IW)=RSTPAR0(JBTYP,IV,IW)
             DO IL=1,2
               soil_TRAN0(IV,IW,IL)=TRAN0(JBTYP,IV,IW,IL)
               soil_REF0(IV,IW,IL)=REF0(JBTYP,IV,IW,IL)
             ENDDO
           ENDDO
         ENDDO

         DO IW=1,3
           soil_SOREF0(IW)=SOREF0(JBTYP,IW)
         ENDDO
!
         DO K=1,3
           soil_DEPTH0(K)=DEPTH0(JBTYP,K)
         ENDDO
!
         DO IM=1,12
           DO IV=1,2
              soil_ZLT0(IM,IV)=ZLT0(JBTYP,IM,IV)
              soil_GREEN0(IM,IV)=GREEN0(JBTYP,IM,IV)
              soil_VCOVER0(IM,IV)=VCOVER0(JBTYP,IM,IV)
           ENDDO

           soil_Z20(IM)=Z20(JBTYP,IM)
           soil_Z10(IM)=Z10(JBTYP,IM)
           soil_ZZ0(IM)=ZZ0(JBTYP,IM)
           soil_D0(IM)=D0(JBTYP,IM)
           soil_RBC0(IM)=RBC0(JBTYP,IM)
           soil_RDC0(IM)=RDC0(JBTYP,IM)
         ENDDO
!
      END SUBROUTINE
!=====================================================================
!
      SUBROUTINE BARELAND
!                                                         March 2010
!=====================================================================
      integer:: IW,IV,IL,K,IM

        VM00(JBTYP)  = soil_vm00
        fc30(JBTYP)  = soil_fc30
        BEE0(JBTYP)  = soil_bee0
        PHSAT0(JBTYP)= soil_PHSAT0
        SATCO0(JBTYP)= soil_SATCO0
        POROS0(JBTYP)= soil_POROS0
        SLOPE0(JBTYP)= soil_SLOPE0

        DO IV=1,2
          CHIL0(JBTYP,IV) = soil_CHIL0(IV)
          TOPT0(JBTYP,IV) = soil_TOPT0(IV)
          TL0(JBTYP,IV)   = soil_TL0(IV)
          TU0(JBTYP,IV)   = soil_TU0(IV)
          DEFAC0(JBTYP,IV)= soil_DEFAC0(IV)
          PH10(JBTYP,IV)  = soil_PH10(IV)
          PH20(JBTYP,IV)  = soil_PH20(IV)
          ROOTD0(JBTYP,IV)= soil_ROOTD0(IV)

          DO IW=1,3
             RSTPAR0(JBTYP,IV,IW) = soil_RSTPAR0(IV,IW)
             DO IL=1,2
               TRAN0(JBTYP,IV,IW,IL)= soil_TRAN0(IV,IW,IL)
               REF0(JBTYP,IV,IW,IL) = soil_REF0(IV,IW,IL)
             ENDDO
           ENDDO
        ENDDO

        DO IW=1,3
           SOREF0(JBTYP,IW) = soil_SOREF0(IW)
        ENDDO
!
        DO K=1,3
           DEPTH0(JBTYP,K) = soil_DEPTH0(K)
        ENDDO
!
        DO IM=1,12
          DO IV=1,2
             ZLT0(JBTYP,IM,IV)   = soil_ZLT0(IM,IV)
             GREEN0(JBTYP,IM,IV) = soil_GREEN0(IM,IV)
             VCOVER0(JBTYP,IM,IV)= soil_VCOVER0(IM,IV)
          ENDDO

          Z20(JBTYP,IM) = soil_Z20(IM)
          Z10(JBTYP,IM) = soil_Z10(IM)
          ZZ0(JBTYP,IM) = soil_ZZ0(IM)
          D0(JBTYP,IM)  = soil_D0(IM)
          RBC0(JBTYP,IM)= soil_RBC0(IM)
          RDC0(JBTYP,IM)= soil_RDC0(IM)
        ENDDO
!
      END SUBROUTINE
!
!==================================================================================
!
      SUBROUTINE CO2ATM( IYR,CAM)
!                                                         NOV 24, 2010
!==============================================================================
      integer IYR,IDY
      real cam
! Atmospheric CO2 Concentration from Xue,
! reference: (http://www.esrl.noaa.gov/gmd/ccgg/trends)

       real,dimension(107):: ycam
       data ycam /295.800,                                   &! 1900
         296.125,  296.475,  296.825,  297.200,  297.625,    &! 1901-1905
         298.075,  298.500,  298.900,  299.300,  299.700,    &! 1906-1910
         300.075,  300.425,  300.775,  301.100,  301.400,    &! 1911-1915
         301.725,  302.075,  302.400,  302.700,  303.025,    &! 1916-1920
         303.400,  303.775,  304.125,  304.525,  304.975,    &! 1921-1925
         305.400,  305.825,  306.300,  306.775,  307.225,    &! 1926-1930
         307.700,  308.175,  308.600,  309.000,  309.400,    &! 1931-1935
         309.750,  310.000,  310.175,  310.300,  310.375,    &! 1936-1940
         310.375,  310.300,  310.200,  310.125,  310.100,    &! 1941-1945
         310.125,  310.200,  310.325,  310.500,  310.750,    &! 1946-1950
         311.100,  311.500,  311.925,  312.425,  313.000,    &! 1951-1955
         313.600,  314.225,  314.8475, 315.500,  316.2725,   &! 1956-1960
         317.075,  317.795,  318.3975, 318.925,  319.6475,   &! 1961-1965
         320.6475, 321.605,  322.635,  323.9025, 324.985,    &! 1966-1970
         325.855,  327.140,  328.6775, 329.7425, 330.585,    &! 1971-1975
         331.7475, 333.2725, 334.8475, 336.525,  338.360,    &! 1976-1980
         339.7275, 340.7925, 342.1975, 343.7825, 345.2825,   &! 1981-1985
         346.7975, 348.645,  350.7375, 352.4875, 353.855,    &! 1986-1990
         355.0175, 355.885,  356.7775, 358.1275, 359.8375,   &! 1991-1995
         361.4625, 363.155,  365.3225, 367.3475, 368.865,    &! 1996-2000
         370.4675, 372.5225, 374.760,  376.8125, 378.8125,   &! 2001-2005
         380.8275/                                            ! 2006
!
       IF(IYR.GE.1900.AND.IYR.LE.2006) THEN
          IDY = IYR-1900+1
          cam = YCAM(IDY)
       ELSE IF(IYR.GT.2006) THEN
          cam=381.
       ELSE
          cam=296.
       ENDIF
!
       END SUBROUTINE
!------------------------------------------------------------------------------
END MODULE
