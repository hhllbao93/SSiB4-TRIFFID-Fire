MODULE ssib4_module_ssibsub
       use ssib4_module_comsconst
       use ssib4_module_comsveg
!
contains
!=======================================================================
!                                                                       
      SUBROUTINE CROPS(XLAT,JDAY,CHIL,ZLT,GREEN,XCOVER,   &
              RSTPAR,TOPT,TL,TU,DEFAC,PH2,PH1)
!
!=======================================================================
!
!     A NEW CROP VERSION BY XUE.                            AUG., 1998
!
!     XLAT IS FROM -90 TO 90 DEGREES  FROM S. TO N.
!
!----------------------------------------------------------------------
      DIMENSION GREEN (2),XCOVER(2),                      &
                CHIL  (2),ZLT   (2),                      &
                RSTPAR(2,3), TOPT(2),                     &
                TL(2), TU(2), DEFAC(2),                   &
                PH1(2), PH2(2)
!
      DIMENSION PHENST(9),WLAI(9),WGRN(9)
!
!C-----------------------------------------------------------------
!**               E   J    H    SD   R   HRV  CUT  PRE-E   E
!     SAVE WLAI,WGRN,IHEAD,IEND,DEND,IWHEAT,SYR
      DATA WLAI/1.0, 2.0, 6.0, 4.0, 3.0, 1.0, 0.01, 0.01, 1.0/
      DATA WGRN/0.6, 0.9, 0.8, 0.5, 0.2, 0.1, 0.01, 0.01, 0.6/
      DATA IHEAD,IEND,DEND,IWHEAT/3,9,244.,12/,SYR/365.25E0/
      IF (XLAT.LT.0.) THEN
      RDAY= JDAY+184
      IF (RDAY.GT.365) RDAY=RDAY-365
      ELSE
      RDAY= JDAY
      END IF
      JULDAY=INT(RDAY+0.2)
!rr   PHI=XLAT
      PHI=90.0 - 180.0E0/3.14159 * XLAT
      APHI = ABS(PHI)
      IF (APHI.GT.55.) PHI=SIGN(55.,PHI)
      IF (APHI.LT.20.) PHI=SIGN(20.,PHI)
!
      FLIP =   0.0
!     IF(PHI.LT.0.0 E0 )FLIP = 182.5
!
! ** DETERMINE WHEAT PHENOLOGY FOR LATITUDE AND JULIAN DAY
       PHENST(2) = 4.50    *ABS(PHI) - 64.0     + FLIP
       PHENST(3) = 4.74    *ABS(PHI) - 46.2     + FLIP
       PHENST(4) = 4.86    *ABS(PHI) - 30.8     + FLIP
       PHENST(5) = 4.55    *ABS(PHI) -  3.0     + FLIP
       PHENST(6) = 4.35    *ABS(PHI) + 11.5     + FLIP
       PHENST(7) = PHENST(6) + 3.0
       DEMG      = ABS( 5.21    *ABS(PHI) - 0.3    )
       PHENST(1) = PHENST(2) - DEMG
       PHENST(9) = PHENST(1)
       PHENST(8) = PHENST(9) - 5.0
!
       DO 10 NS = 1,9
       IF(PHENST(NS) .LT. 0.0E0)PHENST(NS) = PHENST(NS) + 365.
       IF(PHENST(NS) .GT. 365. )PHENST(NS) = PHENST(NS) - 365.
   10  CONTINUE
!
       ROOTGC = 1.0
       CHILW  =-0.02
       TLAI   = 0.5
       GRLF   = 0.6
!
! ** FIND GROWTH STAGE GIVEN LATITUDE AND DAY
       DO 50 NS = 1,8
       TOP = PHENST(NS+1)
       BOT = PHENST(NS)
       DIFF1 = TOP-BOT
       DIFF2 = RDAY-BOT
       IF(RDAY.GE. BOT .AND. RDAY .LE. TOP ) GO TO 40
       IF(BOT .LT. TOP ) GO TO 50
!
! ** PHENOLOGY STAGES OVERLAP THE END OF YEAR?
       ICOND = 0
       IF(RDAY .GE. BOT   .AND. RDAY .LE. 365.) ICOND = 1
       IF(RDAY .GE. 0.0   .AND. RDAY .LE. TOP ) ICOND = 2
!
       IF(ICOND .EQ. 0)GO TO 50
       IF(ICOND .EQ. 2)GO TO 35
           DIFF1 = 365.    - BOT + TOP
           DIFF2 = RDAY     - BOT
           GO TO 40
!
   35  CONTINUE
           DIFF1 = 365.   - BOT + TOP
           DIFF2 = 365.   - BOT + RDAY
!
! ** DATE FOUND IN PHENOLOGY STAGE
   40  CONTINUE
       IF ((RDAY.GT.PHENST(IHEAD)).AND.(RDAY.LE.DEND)) THEN      
           TLAI=WLAI(IHEAD)
           GRLF=WGRN(IHEAD)
           GO TO 77
       END IF
       IF ((RDAY.GT.DEND).AND.(RDAY.LE.PHENST(IEND))) THEN
          DIFF1=PHENST(IEND)-DEND
          DIFF2=RDAY-DEND
          PERC =  DIFF2/DIFF1
          TLAI =  PERC*(WLAI(IEND)-WLAI(IHEAD)) + WLAI(IHEAD)
          GRLF =  PERC*(WGRN(IEND)-WGRN(IHEAD)) + WGRN(IHEAD)
          GO TO 77
       END IF
       PERC =  DIFF2/DIFF1
       TLAI =  PERC*(WLAI(NS+1)-WLAI(NS)) + WLAI(NS)
       GRLF =  PERC*(WGRN(NS+1)-WGRN(NS)) + WGRN(NS)
   77     CONTINUE
       GO TO  95
   50  CONTINUE
   95  CONTINUE
       XCOVER(1)=0.90*(1.0 - EXP(-TLAI))
       ZLTGMX = WLAI(IHEAD)
       ROOTGC = 2910.0    * (0.5    +0.5    *TLAI/ZLTGMX * GRLF)
       IF (NS.NE.1.AND.NS.NE.2) CHILW=-0.2
!
       ZLT   (1) = TLAI
       GREEN (1) = GRLF
       CHIL  (1) = CHILW
!
       END SUBROUTINE
!=======================================================================
!                                                                       
      SUBROUTINE INTERC(                                              &
          DTT,VCOVER,ZLT,TM,TC,TGS,CAPAC,WWW,PPC,PPL,ROFF,            &
          ZDEPTH,POROS,CCX,CG,SATCO,SATCAP,SPWET,EXTK,RNOFFS,FILTR,   &
          SMELT)
!                                                         12 AUGUST 2000
!=======================================================================
!                                                                       
!     CALCULATION OF (1) INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW 
!                    (2) SPECIFIC HEAT TERMS FIXED FOR TIME STEP        
!                                                                       
!     MODIFICATION 30 DEC 1985 : NON-UNIFORM PRECIPITATION             
!     ------------      CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY  
!                       RELATIONSHIP :-                                 
!                                                                       
!                                        F(X) = A*EXP(-B*X)+C           
!                                                                       
!                       THROUGHFALL, INTERCEPTION AND INFILTRATION      
!                       EXCESS ARE FUNCTIONAL ON THIS RELATIONSHIP      
!                       AND PROPORTION OF LARGE-SCALE PPN.              
!---------------------------------------------------------------------- 
!                                                                       
      DIMENSION VCOVER(2),ZLT(2),WWW(3),CAPAC(2),SATCAP(2),EXTK(2,3,2)
      DIMENSION ZDEPTH(3),SNOWW(2)
      DIMENSION CAPACP(2),SNOWP(2),PCOEFS(2,2)
      DATA PCOEFS(1,1)/ 20. /, PCOEFS(1,2)/ .206E-8 /,                & 
           PCOEFS(2,1)/ 0.0001 /, PCOEFS(2,2)/ 0.9999 /, BP /20. /      
!                                                                       
      AP = PCOEFS(2,1)                                                  
      CP = PCOEFS(2,2)                                                  
      TOTALP = PPC + PPL                                                
      IF(TOTALP.LT.1.E-8)GO TO 6000                                     
      AP = PPC/TOTALP * PCOEFS(1,1) + PPL/TOTALP * PCOEFS(2,1)          
      CP = PPC/TOTALP * PCOEFS(1,2) + PPL/TOTALP * PCOEFS(2,2)          
 6000 CONTINUE                                                          
!                                                                       
      ROFF = 0.                                                         
      THRU = 0.                                                         
      FPI  = 0.                                                         
!                                                                       
!---------------------------------------------------------------------- 
!     THERMAL CONDUCTIVITY OF THE SOIL, TAKING INTO ACCOUNT POROSITY    
!---------------------------------------------------------------------- 
!                                                                       
      THETA=WWW(1)*POROS                                                
      CHISL=( 9.8E-4+1.2E-3*THETA )/( 1.1-0.4*THETA )                   
      CHISL=CHISL*4.186E2                                               
!                                                                       
!---------------------------------------------------------------------- 
!     THERMAL DIFFUSIVITY AND HEAT CAPACITYOF THE SOIL                  
!---------------------------------------------------------------------- 
!                                                                       
      DIFSL=5.E-7                                                       
!                                                                       
      ROCS =CHISL/DIFSL                                                 
      D1   =SQRT(DIFSL*86400.0)                                         
      CSOIL=ROCS*D1/SQRT(PIE)/2.0 
      THALAS=0.
      OCEANS=0.
      POLAR=0.                                 
      CSOIL=CSOIL*(1.0-THALAS)+10.E10*OCEANS+POLAR*3.6*4.2E4            
!                                                                       
      P0 = TOTALP * 0.001                                               
!                                                                       
!---------------------------------------------------------------------- 
!     INPUT PRECIPITATION IS GIVEN IN MM, CONVERTED TO M TO GIVE P0.    
!---------------------------------------------------------------------- 
!                                                                       
      DO 1000 IVEG = 1, 2                                               
!                                                                       
      SPWET1 = AMIN1 ( 0.05, CAPAC(IVEG))*CW                            
!                                                                       
      TS = TC                                                           
      SPECHT = ZLT(1) * CLAI                                            
      IF ( IVEG .EQ. 1 ) GO TO 1100                                     
      TS = TGS                                                          
      SPECHT = CSOIL                                                    
 1100  CONTINUE                                                         
!                                                                       
      XSC = AMAX1(0., CAPAC(IVEG) - SATCAP(IVEG) )                      
      IF(IVEG.EQ.2 .AND. TS.LE.TF )GO TO 1170                           
      CAPAC(IVEG) = CAPAC(IVEG) - XSC                                   
      ROFF = ROFF + XSC                                                 
      RNOFFS = XSC*1000. + RNOFFS
 1170  CONTINUE                                                         
      CAPACP(IVEG) = 0.                                                 
      SNOWP(IVEG) = 0.                                                  
!                                                                       
      IF( TS .GT. TF ) CAPACP(IVEG) = CAPAC(IVEG)                       
      IF( TS .LE. TF ) SNOWP(IVEG) = CAPAC(IVEG)                        
      CAPAC(IVEG) = CAPACP(IVEG)                                        
      SNOWW(IVEG) = SNOWP(IVEG)                                         
      ZLOAD = CAPAC(IVEG) + SNOWW(IVEG)                                 
!                                                                       
      FPI = ( 1.-EXP( - EXTK(IVEG,3,1) * ZLT(IVEG)/VCOVER(IVEG) ) )  & 
            * VCOVER(IVEG)                                              
      TTI = P0 * ( 1.-FPI )                                            
!                                                                       
!---------------------------------------------------------------------- 
!    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)            
!---------------------------------------------------------------------- 
!                                                                       
      XS = 1.                                                           
      IF ( P0 .LT. 1.E-9 ) GO TO 1150                                   
      ARG =  ( SATCAP(IVEG)-ZLOAD )/( P0*FPI*AP ) -CP/AP                
      IF ( ARG .LT. 1.E-9 ) GO TO 1150                                  
      XS = -1./BP * ALOG( ARG )                                         
      XS = AMIN1( XS, 1. )                                              
      XS = AMAX1( XS, 0. )                                              
 1150  TEX = P0*FPI * ( AP/BP*( 1.- EXP( -BP*XS )) + CP*XS ) -       &   
            ( SATCAP(IVEG) - ZLOAD ) * XS                               
      TEX = AMAX1( TEX, 0. )                                            
!                                                                       
!---------------------------------------------------------------------- 
!    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION                    
!---------------------------------------------------------------------- 
!                                                                       
      THRU = TTI + TEX                                                  
      IF(IVEG.EQ.2.AND.TGS.LE.TF)THRU = 0.                              
!                                                                       
      PINF = P0 - THRU                                                  
      IF( TM .GT. TF ) CAPAC(IVEG) = CAPAC(IVEG) + PINF                 
      IF( TM .LE. TF ) SNOWW(IVEG) = SNOWW(IVEG) + PINF                 
!                                                                       
      IF( IVEG .EQ. 1 ) GO TO 1300                                      
      IF( TM .GT. TF ) GO TO 1200                                      
      SNOWW(IVEG) = SNOWP(IVEG) + P0                                    
      THRU = 0.                                                         
      GO TO 1300                                                        
!                                                                       
!---------------------------------------------------------------------- 
!    INSTANTANEOUS OVERLAND FLOW CONTRIBUTION ( ROFF )                  
!---------------------------------------------------------------------- 
!                                                                       
 1200  EQUDEP = SATCO * DTT                                            
!                                                                       
      XS = 1.                                                           
      IF ( THRU .LT. 1.E-9 ) GO TO 1250                                 
      ARG = EQUDEP / ( THRU * AP ) -CP/AP                               
      IF ( ARG .LT. 1.E-9 ) GO TO 1250                                  
      XS = -1./BP * ALOG( ARG )                                         
      XS = AMIN1( XS, 1. )                                              
      XS = AMAX1( XS, 0. )                                              
 1250  ROFFO = THRU * ( AP/BP * ( 1.-EXP( -BP*XS )) + CP*XS )      &        
             -EQUDEP*XS                                                 
      ROFFO = AMAX1 ( ROFFO, 0. )                                       
      ROFF = ROFF + ROFFO                                               
      RNOFFS = RNOFFS + ROFFO*1000.
      FILTR =  FILTR + (THRU - ROFFO)             
      WWW(1) = WWW(1) + (THRU - ROFFO) / ( POROS*ZDEPTH(1) )            
 1300  CONTINUE                                                        
!                                                                       
!---------------------------------------------------------------------- 
!    TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION                
!---------------------------------------------------------------------- 
!                                                                       
      DIFF = ( CAPAC(IVEG)+SNOWW(IVEG) - CAPACP(IVEG)-SNOWP(IVEG) )*CW 
      CCP = SPECHT + SPWET1                                             
      CCT = SPECHT + SPWET1 + DIFF                                      
!                                                                       
      TSD = ( TS * CCP + TM * DIFF ) / CCT                              
!                                                                       
      FREEZE = 0.                                                       
      IF ( TS .GT. TF .AND. TM .GT. TF ) GO TO 2000                     
      IF ( TS .LE. TF .AND. TM .LE. TF ) GO TO 2000                     
!                                                                       
      TTA = TS                                                          
      TTB = TM                                                          
      CCA = CCP                                                         
      CCB = DIFF                                                        
      IF ( TSD .GT. TF ) GO TO 2100                                     
!                                                                       
!---------------------------------------------------------------------- 
!    FREEZING OF WATER ON CANOPY OR GROUND                              
!---------------------------------------------------------------------- 
!                                                                       
      CCC = CAPACP(IVEG) * SNOMEL                                       
      IF ( TS .LT. TM ) CCC = DIFF * SNOMEL / CW                        
      TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT                       
!                                                                       
      FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )                 
      FREEZE = (AMIN1 ( CCC, FREEZE )) / SNOMEL                         
      IF(TSD .GT. TF)TSD = TF - 0.1                                     
!                                                                       
      GO TO 2000                                                        
!                                                                       
 2100  CONTINUE                                                      
!                                                                       
!---------------------------------------------------------------------- 
!    MELTING OF SNOW ON CANOPY OR GROUND                                
!---------------------------------------------------------------------- 
!                                                                       
      CCC = - SNOWW(IVEG) * SNOMEL                                      
      IF ( TS .GT. TM ) CCC = - DIFF * SNOMEL / CW                      
!                                                                       
      TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT                       
!                                                                       
      FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )                 
      FREEZE = (AMAX1( CCC, FREEZE )) / SNOMEL                          
      IF(TSD .LE. TF)TSD = TF - 0.1                                     
!                                                                       
 2000  CONTINUE
      SMELT = FREEZE
      SNOWW(IVEG) = SNOWW(IVEG) + FREEZE                                
      CAPAC(IVEG) = CAPAC(IVEG) - FREEZE                                
!                                                                       
      IF( IVEG .EQ. 1 ) TC = TSD                                        
      IF( IVEG .EQ. 2 ) TGS = TSD                                       
      IF( SNOWW(IVEG) .LT. 0.0000001 ) GO TO 3000                       
      ZMELT = 0.                                                        
!     modified to force water into soil. Xue Feb. 1994
      ZMELT = CAPAC(IVEG)                             
!     IF ( TD .GT. TF ) ZMELT = CAPAC(IVEG)                             
!     IF ( TD .LE. TF ) ROFF = ROFF + CAPAC(IVEG)                       
      CAPAC(IVEG) = 0.                                                  
      WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) )                   
      FILTR = FILTR + ZMELT                    
!                                                                       
 3000  CONTINUE                                                       
!                                                                       
      CAPAC(IVEG) = CAPAC(IVEG) + SNOWW(IVEG)                           
      SNOWW(IVEG) = 0.                                                  
!
!     **** LOAD PILPS PARAMETER
!
!     if (freeze.lt.0) snm(istat)=snm(istat)-freeze
      freeze=0.0
!                                                                       
      P0 = THRU                                                         
!                                                                       
 1000  CONTINUE                                                        
!                                                                       
!---------------------------------------------------------------------- 
!    CALCULATION OF CANOPY AND GROUND HEAT CAPACITIES.                  
!    N.B. THIS SPECIFICATION DOES NOT NECESSARILY CONSERVE ENERGY WHEN  
!    DEALING WITH VERY LATGE SNOWPACKS.                                 
!---------------------------------------------------------------------- 
!                                                                       
      CCX = ZLT(1) * CLAI + CAPAC(1) * CW                               
      SPWET = AMIN1 ( 0.05, CAPAC(2))*CW                                
      CG = (CSOIL + SPWET)                                              
!                                                                       
      END  SUBROUTINE
!====================================================================   
!zzq..dec 28,2007
!zzq   SUBROUTINE NEWTON(A1,Y,FINC,NOX,NONPOS,IWOLK,L)
       SUBROUTINE NEWTON(A1,Y,FINC,NOX,NONPOS,IWOLK,L,ZINC,A2,Y1,ITER)
!                                                                       
!====================================================================== 
! ** VERSION ACQUIRED FROM EROS 2/19/86.                                
!                                                                       
!=======================================================================
!                                                                       
! ** THE NEWTON RAPHSON ITERATIVE ROUTINE WILL BE USED TO GENERATE NEW  
! ** VALUES OF A1 IF DABSOLUTE VALUE OF Y IS GREATER THAN ERTOL;        
! ** A1 IS ESTIMATE, Y IS RESULTANT ERROR                               
! ** NEX IS EXIT CONDITION  (0=NO EXIT) OR (1 WHEN DABS(Y) LT ERTOL)    
! ** ERTOL IS THE DABSOLUTE VALUE OF Y NECESSARY TO OBTAIN AN EXIT      
! ** FINC IS INITIAL INCREMENT SIZE FOR SECOND ESTIMATE OF A1           
! ** NONPOS=0 IF QUANTITY TO BE MINIMIZED CAN BE LESS THAN ZERO;        
! ** NONPOS=1 IF QUANTITY CAN ONLY BE POSITIVE                          
! ** L IDENTIFIES WHICH QUANTITY IS BEING CALCULATED.                   
!                                                                       
! ** CONTROL VALUES: FINC,ERTOL,NOX,NONPOS,L:MUST BE SET BY USER        
!-----------------------------------------------------------------------
!                                                                       
       DIMENSION   IWALK(3), NEX(3)                                     
!zzq       COMMON/TONNEW/  ZINC(3), A2(3), Y1(3)                            
!zzq       COMMON/NEWT/  ITER(3)                                            
       real  ZINC(3), A2(3), Y1(3)                            
       integer ITER(3)                                            
       DATA CONS/1.0/                                                   
!                                                                       
       ERTOL = 0.05 * FINC                                              
       IWALK(L) = IWOLK                                                 
       NEX(L)=NOX                                                       
!                                                                       
       IF ( ITER(L) .GE. 490 ) GO TO 160                                
       IF (ERTOL .LT. 0.00000001) ERTOL=0.000001                        
       IF (ABS(Y) .LE. ERTOL) GO TO 150                                 
       IF((ABS(Y-Y1(L))).LE.0.01*ERTOL .AND. IWALK(L).EQ.0 ) GO TO 8    
!                                                                     
       IF(ABS(Y1(L)).GT.ERTOL) GO TO 1                                  
       A2(L)=A1                                                         
       A1=A1-Y                                                          
       NEX(L)=0                                                         
       Y1(L)=Y                                                          
       ITER(L)=1                                                        
       IF (IWALK(L) .EQ. 3) GO TO 101                                   
       IWALK(L)=0                                                       
       GO TO 101                                                        
    1   ITER(L)=ITER(L)+1                                            
       IF(ITER(L) .EQ. 10) IWALK(L)=1                                   
       IF(IWALK(L) .NE. 0) GO TO 2                                      
       IF(ABS(Y) .GT. ERTOL) GO TO 3                                    
       NEX(L)=1                                                         
       GO TO 150                                                        
    3   A=A1-Y*(A1-A2(L))/(Y-Y1(L))                                     
       IF(ABS(A-A1).GT.(10.0*FINC))                    &                           
                  A=A1+10.0*FINC*SIGN(CONS,(A-A1))                      
       A2(L)=A1                                                         
       A1=A                                                             
       Y1(L)=Y                                                          
       GO TO 101                                                        
    2   IF(IWALK(L).EQ.2)GO TO 4                                        
       IF(IWALK(L).EQ.3) GO TO 6                                        
       IF(SIGN(CONS,Y).EQ.SIGN(CONS,Y1(L))) GO TO  3                    
       ZINC(L)=(A1-A2(L))/4.0                                           
       A1=A2(L)+ZINC(L)                                                 
       IWALK(L)=2                                                       
       NEX(L)=0                                                         
       GO TO 101                                                        
    4   IF(SIGN(CONS,Y) .EQ.SIGN(CONS,Y1(L))) GO TO 5                   
       ZINC(L)=-ZINC(L)/4.0                                             
       A2(L)=A1                                                         
       A1=A1+ZINC(L)                                                    
       NEX(L)=0                                                         
       Y1(L)=Y                                                          
       GO TO 101                                                        
    5   A2(L)=A1                                                       
       A1=A1+ZINC(L)                                                    
       Y1(L)=Y                                                          
       NEX(L)=0                                                         
       GO TO 101                                                        
    6   IF(SIGN(CONS,Y).EQ.SIGN(CONS,Y1(L))) GO TO 7                   
       IWALK(L)=1                                                       
       GO TO 2                                                          
    7   A2(L) = A1                                                     
       A1 = A1+FINC                                                     
       Y1(L)=Y                                                          
       NEX(L) = 0                                                       
       GO TO 101                                                        
    8   A1 = A1 + FINC*2.0                                             
       NEX(L)=0                                                         
       GO TO 101                                                        
  160    CONTINUE                                                      
       WRITE(7,900) Y, A1                                               
!      STOP                                                             
  900   FORMAT ( 3X,' FAILURE TO CONVERGE AFTER 490 ITERATIONS',  &       
       /, 3X,' Y = ',2G12.5,2X,I14)                                     
!                                                                       
  150   NEX(L) = 1                                                     
       ZINC(L)=0.0                                                      
       ITER(L) = 0                                                      
       IWALK(L)=0                                                       
       Y1(L)=0.0                                                        
       Y=0.0                                                            
       A2(L)=0.0                                                        
  101   CONTINUE                                                       
       IF(NONPOS.EQ.1.AND.A1.LT.0.0) A1=A2(L)/2.0                       
       NOX = NEX(L)                                                     
       IWOLK = IWALK(L)                                                 
!                                                                       
       END SUBROUTINE
!=======================================================================
!                                                                       
!song 9.2.2011 SUBROUTINE RADAB(TRAN,REF,GREEN,VCOVER,CHIL,ZLT,Z2,Z1,SOREF,
      SUBROUTINE RADAB(TRAN,REF,GREEN,VCOVER,CHIL,ZLT,Z2,Z1,SOREF,  &
       WWW1,POROS,                                                  &
           TC,TGS,SATCAP,EXTK,CLOSS,GLOSS,THERMK,P1F,P2F,           &
           RADT,PAR,PD,SALB,ALBEDO,TGEFF,SUNANG,XADJ,CAPAC,         &
           RADN,bedo,ZLWUP,RADFRAC,SWDOWN,SCOV2,                    &
           fsdown,fldown,fsup,flup,FPARBD,KPARBD)
!rr  2     RADN,bedo,ZLWUP,RADFRAC,SWDOWN,SCOV2)
!                                                         11 AUGUST 2000
!=======================================================================
!                                                                       
!     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT       
!     AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY                     
!                                                                       
!-----------------------------------------------------------------------
      DIMENSION sr(2)
!
      DIMENSION TRANC1(2), TRANC2(2), TRANC3(2),RADFRAC(2,2)
      DIMENSION CAPAC(2), SATCAP(2)
      DIMENSION TRAN(2,3,2), REF(2,3,2),SOREF(3)
      DIMENSION GREEN(2), VCOVER(2), ZLT(2), CHIL(2)
      DIMENSION RADN(3,2), RADT(2), PAR(2), PD(2)
      DIMENSION RADFAC(2,2,2), RADSAV(12)
      DIMENSION SALB(2,2), ALBEDO(2,3,2), EXTK(2,3,2)
!9.2.2011 song for soil reflectence
      dimension bdcnsta(2),bdcnstb(2),bdcnstc(2)
      real FPARBD,KPARBD
      data sr/0.85,0.65/
      bdcnsta(1)=0.01
      bdcnsta(2)=0.01
      bdcnstb(1)=15
      bdcnstb(2)=20
! bdcnstc(1)=40
! bdcnstc(2)=40
      bdcnstc(1)=80
      bdcnstc(2)=80
!
      f=max(sunang,0.01746)
!
!----------------------------------------------------------------------
!     CALCULATION OF MAXIMUM WATER STORAGE VALUES.
!----------------------------------------------------------------------
!
      FMELT = 1.
      IF ( ABS(TF-TGS) .LT. 0.5 ) FMELT = 0.6
      SATCAP(1) =  ZLT(1) * 0.0001                                      
      SATCAP(2) =  ZLT(2) * 0.0001                                      
      DEPCOV = AMAX1( 0., (CAPAC(2)*5.-Z1) )                            
      DEPCOV = AMIN1( DEPCOV, (Z2-Z1)*0.95 )                            
      SATCAP(1) = SATCAP(1) * ( 1. - DEPCOV / ( Z2 - Z1 ) )             
!rr - thermal part is in use in temrs1
      albedo(1,3,1)=0.
      albedo(1,3,2)=0.
      albedo(2,3,1)=0.
      albedo(2,3,2)=0.
!rr
!                                                                       
!---------------------------------------------------------------------- 
      DO 1000 IWAVE = 1,2 
!
!song 9.2.2011
!zzq09/25/2013      soref0=soref(iwave)
      sorefxx=soref(iwave)
!
      if (soref(iwave) .lt. 0.25) then
            delsoref=bdcnsta(iwave)*(bdcnstb(iwave)-bdcnstc(iwave)*www1*poros)
           if (delsoref .le. -0.05) delsoref=-0.05
           soref(iwave)=sorefxx+ delsoref
      endif
!
      DO 2000 IVDUM = 1,2                                              
      IF ( IVDUM .EQ. 1 ) IVEG = 2                                      
      IF ( IVDUM .EQ. 2 ) IVEG = 1                                      
!---------------------------------------------------------------------- 
!     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO            
!         SNOW REFLECTANCE   = 0.80, 0.40 . MULTIPLY BY 0.6 IF MELTING  
!         SNOW TRANSMITTANCE = 0.20, 0.54                               
!         SNOW REFLECTANCE   = 0.85, 0.65 . MULTIPLY BY 0.6 IF MELTING  
!                                                                       
!---------------------------------------------------------------------- 
      SCOV = 0.                                                         
      IF( IVEG .EQ. 2 ) GO TO 100                                       
      IF( TC .LE. TF ) SCOV =  AMIN1( 0.5, CAPAC(1) / SATCAP(1) )       
  100   CONTINUE                                                       
      REFF1 = ( 1. - SCOV ) * REF(IVEG,IWAVE,1) + SCOV * ( 1.2 -     & 
              IWAVE * 0.4 ) * FMELT                                     
      REFF2 = ( 1. - SCOV ) * REF(IVEG,IWAVE,2) + SCOV * ( 1.2 -     & 
              IWAVE * 0.4 ) * FMELT                                     
      TRAN1 = TRAN(IVEG,IWAVE,1) * ( 1. - SCOV )                     & 
              + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT )         & 
              * TRAN(IVEG,IWAVE,1)                                      
      TRAN2 = TRAN(IVEG,IWAVE,2) * ( 1. - SCOV )                     & 
              + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT ) * 0.9   & 
              * TRAN(IVEG,IWAVE,2)                                      
                                                                        
!---------------------------------------------------------------------- 
!                                                                       
      SCAT = GREEN(IVEG)*( TRAN1 + REFF1 ) +( 1. - GREEN(IVEG) ) *   & 
             ( TRAN2 + REFF2)                                           
      CHIV = CHIL(IVEG)                                                 
!                                                                       
      IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01                            
      AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV                      
      BB = 0.877 * ( 1. - 2. * AA )                                     
!                                                                       
      PROJ = AA + BB * F                                                
      EXTKB = ( AA + BB * F ) / F                                       
      ZMEW = 1. / BB * ( 1. - AA / BB * ALOG ( ( AA + BB ) / AA ) )     
      ACSS = SCAT / 2. * PROJ / ( PROJ + F * BB )                       
      ACSS = ACSS * ( 1. - F * AA / ( PROJ + F * BB) * ALOG ( ( PROJ  & 
             +   F * BB + F * AA ) / ( F * AA ) ) )                     
!                                                                       
      EXTK( IVEG, IWAVE, 1 ) = PROJ / F * SQRT( 1.-SCAT )               
      EXTK( IVEG, IWAVE, 2 ) = 1. / ZMEW * SQRT( 1.-SCAT )              
      EXTK( IVEG, 3, 1 ) = AA + BB                                      
      EXTK( IVEG, 3, 2 ) = 1./ZMEW                                      
!                                                                       
      UPSCAT = GREEN(IVEG) * TRAN1 + ( 1. - GREEN(IVEG) ) * TRAN2       
      UPSCAT = 0.5 * ( SCAT + ( SCAT - 2. * UPSCAT ) *                &  
               (( 1. - CHIV ) / 2. ) ** 2 )                             
!                                                                       
      BETAO = ( 1. + ZMEW * EXTKB ) / ( SCAT * ZMEW * EXTKB ) * ACSS    
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     DICKINSON'S VALUES 
!                                                                       
      BE = 1. - SCAT + UPSCAT                                           
      CE = UPSCAT                                                       
      BOT = ( ZMEW * EXTKB ) ** 2 + ( CE**2 - BE**2 )                   
      IF ( ABS(BOT) .GT. 1.E-10) GO TO 200                              
      SCAT = SCAT* 0.98                                                 
      BE = 1. - SCAT + UPSCAT                                           
      BOT = ( ZMEW * EXTKB ) ** 2 + ( CE**2 - BE**2 )                   
  200   CONTINUE                                                       
      DE = SCAT * ZMEW * EXTKB * BETAO                                  
      FE = SCAT * ZMEW * EXTKB * ( 1. - BETAO )                         
!---------------------------------------------------------------------- 
!                                                                       
      CCE = DE * BE - ZMEW * DE * EXTKB + CE * FE                       
      FFE = BE * FE + ZMEW * FE * EXTKB + CE * DE                       
!                                                                       
      TORE = -CCE / BOT                                                 
      SIGE = -FFE / BOT                                                 
!                                                                       
      PSI = SQRT(BE**2 - CE**2)/ZMEW                                    
!                                                                       
!---------------------------------------------------------------------- 
!     REDUCTION IN EXPOSED HEIGHT OF UPPER STOREY AS SNOW ACCUMULATES   
!                                                                       
      SDEP = CAPAC(2) * 5.                                              
      FAC = ( SDEP - Z1 ) / ( Z2 - Z1 )                                 
      FAC = AMAX1( 0., FAC )                                            
      FAC = AMIN1( 0.99, FAC )                                          
!                                                                       
      ZAT = ZLT(IVEG) / VCOVER(IVEG)                                    
      IF ( IVEG .EQ. 1 ) ZAT = ZAT * (1.-FAC)                           
!                                                                       
      POWER1 = AMIN1( PSI*ZAT, 50. )                                    
      POWER2 = AMIN1( EXTKB*ZAT, 50. )                                  
      EPSI = EXP( - POWER1 )                                            
      EK = EXP ( - POWER2 )                                             
!                                                                       
      ROSB = SOREF(IWAVE)                                               
      ROSD = SOREF(IWAVE)                                               
      IF ( IVEG .EQ. 2 ) GO TO 300                                      
      ROSB = ALBEDO(2,IWAVE,1)                                          
      ROSD = ALBEDO(2,IWAVE,2)                                          
  300   CONTINUE                                                       
!                                                                       
      GE = ROSB / ROSD                                                  
!                                                                       
!-----------------------------------------------------------------------
!     CALCULATION OF DIFFUSE ALBEDOS                                    
!-----------------------------------------------------------------------
!                                                                       
      F1 = BE - CE / ROSD                                               
      ZP = ZMEW * PSI                                                   
!                                                                       
      DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -       &                    
            ( BE - ZP ) * ( F1 + ZP ) * EPSI                            
      ALPHA = CE * ( F1 - ZP ) / EPSI / DEN                             
      BETA = -CE * ( F1 + ZP ) * EPSI / DEN                             
      F1 = BE - CE * ROSD                                               
      DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI                     
!                                                                       
      GAMMA = ( F1 + ZP ) / EPSI / DEN                                  
      DELTA = - ( F1 - ZP ) * EPSI / DEN                                
!                                                                       
      ALBEDO(IVEG,IWAVE,2) =  ALPHA + BETA                              
!     XQQ(IVEG,IWAVE,2) = ALBEDO(IVEG, IWAVE, 2)                        
!                                                                       
      IF ( IVEG .EQ. 1 ) GO TO 400                                      
      SCOV2 = 0.                                                        
      IF ( TGS .LE. TF ) SCOV2 = AMIN1( 1., CAPAC(2) / 0.004 )          
      ALBEDO(2,IWAVE,2) =                                           &             
       ROSD * ( 1. - VCOVER(2) ) + ALBEDO(2,IWAVE,2) * VCOVER(2)        
      ALBEDO(2,IWAVE,2) =                                           &                
       ( 1. - SCOV2 ) * ALBEDO(2,IWAVE,2) + SCOV2 *                 &                
       ( 1.2-IWAVE*0.4 ) * FMELT
  400   CONTINUE                                                      
!                                                                       
      TRANC2(IWAVE) = GAMMA * EPSI + DELTA / EPSI                       
!                                                                       
!-----------------------------------------------------------------------
!     CALCULATION OF DIRECT ALBEDOS                                     
!-----------------------------------------------------------------------
!                                                                       
      F1 = BE - CE / ROSD                                               
      ZMK = ZMEW * EXTKB                                                
!                                                                       
      DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -                         &  
            ( BE - ZP ) * ( F1 + ZP ) * EPSI                            
      ALPHA = ( DE - TORE * ( BE + ZMK ) ) * ( F1 - ZP ) / EPSI -      &  
              ( BE - ZP ) * ( DE - CE*GE - TORE * ( F1 + ZMK ) ) * EK   
      ALPHA = ALPHA / DEN                                               
      BETA = ( BE + ZP ) * (DE - CE*GE - TORE * ( F1 + ZMK ))* EK -    &   
             ( DE - TORE * ( BE + ZMK ) ) * ( F1 + ZP ) * EPSI          
      BETA = BETA / DEN                                                 
      F1 = BE - CE * ROSD                                               
      DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI                     
      GAMMA = - SIGE * ( F1 + ZP ) / EPSI -                         &   
              ( FE + CE * GE * ROSD + SIGE * ( ZMK - F1 ) ) * EK        
      GAMMA = GAMMA / DEN                                               
      DELTA = ( CE * GE * ROSD + FE + SIGE * ( ZMK - F1 ) ) * EK    &   
              + SIGE * ( F1 - ZP ) * EPSI                               
      DELTA = DELTA / DEN                                               
!                                                                       
      ALBEDO(IVEG,IWAVE,1) = TORE + ALPHA + BETA                        
!     XQQ(IVEG,IWAVE,1) = ALBEDO(IVEG, IWAVE, 1)                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      IF( IVEG .EQ. 1 ) GO TO 500                                       
      ALBEDO(2,IWAVE,1) = ROSB * ( 1. - VCOVER(2) )               &   
                          + ALBEDO(2,IWAVE,1) * VCOVER(2)              
      ALBEDO(2,IWAVE,1) = ( 1. - SCOV2 ) * ALBEDO(2,IWAVE,1) +    &   
                          SCOV2 * ( 1.2-IWAVE*0.4 ) * FMELT             
!                                                                       
  500   CONTINUE                                                       
!                                                                       
      TRANC1(IWAVE) = EK                                                
      TRANC3(IWAVE) = SIGE * EK + GAMMA * EPSI + DELTA / EPSI           
!                                                                       
 2000  CONTINUE                                                        
!                                                                       
!---------------------------------------------------------------------- 
!     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES    
!     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND              
!---------------------------------------------------------------------- 
!                                                                       
      RADFAC(2,IWAVE,1) = ( 1.-VCOVER(1) ) * ( 1.-ALBEDO(2,IWAVE,1) )  & 
             + VCOVER(1) * ( TRANC1(IWAVE) * ( 1.-ALBEDO(2,IWAVE,1) )  & 
             + TRANC3(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )               
!                                                                       
      RADFAC(2,IWAVE,2) = ( 1.-VCOVER(1) ) * ( 1.-ALBEDO(2,IWAVE,2) )  & 
             + VCOVER(1) *  TRANC2(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) )    
!                                                                       
      RADFAC(1,IWAVE,1) = VCOVER(1) * ( ( 1.-ALBEDO(1,IWAVE,1) )   & 
             - TRANC1(IWAVE) * ( 1.-ALBEDO(2,IWAVE,1) )            & 
             - TRANC3(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )               
!                                                                       
      RADFAC(1,IWAVE,2) = VCOVER(1) * ( ( 1.-ALBEDO(1,IWAVE,2) )   & 
             - TRANC2(IWAVE) * ( 1.-ALBEDO(2,IWAVE,2) ) )               
!                                                                       
!     XQQ(1,IWAVE,1) = RADFAC(1,IWAVE,1)                                
!     XQQ(1,IWAVE,2) = RADFAC(1,IWAVE,2)                                
!     XQQ(2,IWAVE,1) = RADFAC(2,IWAVE,1)                                
!     XQQ(2,IWAVE,2) = RADFAC(2,IWAVE,2)                                
!                                                                       
!---------------------------------------------------------------------- 
!     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB )                     
!---------------------------------------------------------------------- 
!                                                                       
      DO 3000 IRAD = 1, 2                                               
      SALB(IWAVE,IRAD) = ( 1.-VCOVER(1) ) * ALBEDO(2,IWAVE,IRAD) +  & 
                         VCOVER(1) * ALBEDO(1,IWAVE,IRAD)               
 3000  CONTINUE                                                      
!                                                                       
!---------------------------------------------------------------------- 
!     SAVING OF EXTINCTION COEFFICIENTS ( PAR ) FOR STOMAT CALCULATION  
!---------------------------------------------------------------------- 
      IF ( IWAVE .EQ. 2 ) GO TO 600                                     
      RADSAV(1) = 1. - VCOVER(1)      & 
                + VCOVER(1) * ( TRANC1(IWAVE) + TRANC3(IWAVE) )         
      RADSAV(2) = 1. - VCOVER(1) + VCOVER(1) * TRANC2(IWAVE)            
!     XQQ(1,1,1) = RADSAV(1)                                            
!     XQQ(1,2,1) = RADSAV(2)                                            
  600   CONTINUE                                                       
!                                                                       
 1000  CONTINUE                                                        

!yliu18Oct2016
      IF ((RADN(1,1)+RADN(1,2)).lt.0.001)THEN
        FPARBD = 0.000001
      ELSE
        FPARBD = RADFAC(1,1,1)*(RADN(1,1)/(RADN(1,1)+RADN(1,2)))+  &
                 RADFAC(1,1,2)*(RADN(1,2)/(RADN(1,1)+RADN(1,2)))
      ENDIF

      IF ((RADN(1,1)+RADN(1,2)).lt.0.001)THEN
        KPARBD = 0.000001
      ELSE
        KPARBD = RADSAV(1)*(RADN(1,1)/(RADN(1,1)+RADN(1,2)))+      &
                 RADSAV(2)*(RADN(1,2)/(RADN(1,1)+RADN(1,2)))
      ENDIF
!
!     albedo adjustment ==============================================
!
      if (xadj.eq.0.) go to 730
      xx = radfac(1,1,2) + radsav(2)
      xy = radfac(1,1,1) + radsav(1)
      ssum = salb(1,1)*radfrac(1,1) + salb(1,2)*radfrac(1,2)+      &
             salb(2,1)*radfrac(2,1) + salb(2,2)*radfrac(2,2)
!     for diffuse albedo
      do 650 iwave = 1, 2
      salb(iwave,2) = salb(iwave,2) + xadj * salb(iwave,2) / ssum
      x0 = 1. - salb(iwave,2)
      x1 = radfac(1,iwave,2) + radfac(2,iwave,2)
      x2 = radfac(1,iwave,2) / x1
      x3 = radfac(2,iwave,2) / x1
      radfac(1,iwave,2) = x0 * x2
      radfac(2,iwave,2) = x0 * x3
      if (salb(iwave,2).gt.1..or.radfac(1,iwave,2).gt.1..or.         &
          radfac(2,iwave,2).gt.1..or.salb(iwave,2).lt.0..or.         &
          radfac(1,iwave,2).lt.0..or.radfac(2,iwave,2).lt.0.) then
          stop 999
      end if
  650  continue
  640  format(1x,'unrealistic value, dif',2i12,4e11.4)
!     for direct albedo
      do 750 iwave = 1, 2
      salb(iwave,1) = salb(iwave,1) + xadj * salb(iwave,1) / ssum
      x0 = 1. - salb(iwave,1)
      x1 = radfac(1,iwave,1) + radfac(2,iwave,1)
      x2 = radfac(1,iwave,1) / x1
      x3 = radfac(2,iwave,1) / x1
      radfac(1,iwave,1) = x0 * x2
      radfac(2,iwave,1) = x0 * x3
      radsav(1) =  xy - radfac(1,1,1)
      radsav(2) =  xx - radfac(1,1,2)
      if (salb(iwave,1).gt.1..or.radfac(1,iwave,1).gt.1..or.         &
          radfac(2,iwave,1).gt.1..or.salb(iwave,1).lt.0..or.         &
          radfac(1,iwave,1).lt.0..or.radfac(2,iwave,1).lt.0.) then
!rr       write(7,740) nymdh,iwave,salb(iwave,1),radfac(1,iwave,1),
!rr  &                  radfac(2,iwave,1)
          stop 999
      end if
  750  continue
  740  format(1x,'unrealistic value',2i12,4e11.4)
  730  continue
!!***************** end adjustment *******************************
      sibsu = radn(1,1)*salb(1,1) + radn(1,2)*salb(1,2)   &
                           + radn(2,1)*salb(2,1) + radn(2,2)*salb(2,2)
      if ((swdown.gt.0.01).and.(sibsu.gt.0.01)) then
!rr   if ((swdown.gt.0.1).and.(sibsu.gt.0.1)) then
         bedo = sibsu / swdown
         if (bedo.gt.1.) then
            sibsu =  0.
            bedo = 999.
!rr      write (7, *) 'albebo incorrect',nymdh,bedo
         endif
      else
         sibsu = 0.0
         bedo = 999.
      endif
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      TC4 = TC * TC * TC * TC                                           
      TG4 = TGS * TGS * TGS * TGS                                       
!                                                                       
      ZKAT = EXTK(1,3,2) * ZLT(1) / VCOVER(1)                           
      ZKAT = AMIN1( 50. , ZKAT )                                        
      ZKAT = AMAX1( 1.E-5, ZKAT )                                       
      THERMK = EXP(-ZKAT)                                               
!                                                                       
      FAC1 =  VCOVER(1) * ( 1.-THERMK )                                 
      FAC2 =  1.                                                        
      CLOSS =  2. * FAC1 * STEFAN * TC4                                 
      CLOSS =  CLOSS - FAC2 * FAC1 * STEFAN * TG4                       
      GLOSS =  FAC2 * STEFAN * TG4                                      
      GLOSS =  GLOSS - FAC1 * FAC2 * STEFAN * TC4                       
!                                                                       
      ZLWUP =  FAC1 * STEFAN * TC4 + (1. - FAC1 ) * FAC2 * STEFAN * TG4 
      TGEFF = SQRT( SQRT ( ( ZLWUP / STEFAN ) ) )                       
!                                                                       
      RADSAV(3) = EXTK(1,1,1)                                           
      RADSAV(4) = EXTK(1,1,2)                                           
      RADSAV(5) = EXTK(2,1,1)                                           
      RADSAV(6) = EXTK(2,1,2)                                           
      RADSAV(7) = THERMK                                                
      RADSAV(8) = EXTK(1,3,1)                                           
      RADSAV(9) = EXTK(2,3,1)                                           
      RADSAV(10)= CLOSS                                                 
      RADSAV(11)= GLOSS                                                 
      RADSAV(12)= TGEFF                                                 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!     CALL LONGRN( TRANC1, TRANC2, TRANC3)                              
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!     CALL RADUSE                                                       
!                                                                       
!=======================================================================
!                                                                       
!     CALCULATION OF ABSORPTION OF RADIATION BY SURFACE                 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      P1F         = RADSAV(1)                                           
      P2F         = RADSAV(2)                                           
      EXTK(1,1,1) = RADSAV(3)                                           
      EXTK(1,1,2) = RADSAV(4)                                           
      EXTK(2,1,1) = RADSAV(5)                                           
      EXTK(2,1,2) = RADSAV(6)                                           
      THERMK      = RADSAV(7)                                           
      EXTK(1,3,1) = RADSAV(8)                                           
      EXTK(2,3,1) = RADSAV(9)                                           
      CLOSS       = RADSAV(10)                                          
      GLOSS       = RADSAV(11)                                          
      TGEFF       = RADSAV(12)                                          
!                                                                       
!---------------------------------------------------------------------- 
!     SUMMATION OF SHORT-WAVE RADIATION ABSORBED BY CANOPY AND GROUND   
!---------------------------------------------------------------------- 
!                                                                       
      RADT(1) = 0.                                                      
      RADT(2) = 0.                                                      
!                                                                       
      DO 7000 IVEG  = 1, 2                                              
      DO 7000 IWAVE = 1, 2                                              
      DO 7000 IRAD  = 1, 2                                              
!                                                                       
      RADT(IVEG) = RADT(IVEG)+RADFAC(IVEG,IWAVE,IRAD)*RADN(IWAVE,IRAD)  
!                                                                       
 7000  CONTINUE                                                        
!KHS====================================================================
      fsdown = radn(1,1)+radn(1,2)+radn(2,1)+radn(2,2)
      fsup   = fsdown-radt(1)-radt(2)
!KHS====================================================================
!                                                                       
      SWCAN=RADT(1)
      SWGND=RADT(2)
!                                                                       
      RADT(1) = RADT(1) + RADN(3,2)*VCOVER(1)*(1.- THERMK)       &      
              - CLOSS                                                  
      RADT(2) = RADT(2) + RADN(3,2)*( 1.-VCOVER(1)*(1-THERMK) )  &       
              - GLOSS                                                   
!KHS====================================================================
      fldown = radn(3,2)
      flup   = closs+gloss
!KHS====================================================================
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      PAR(1) = RADN(1,1) + RADN(1,2) + 0.001                            
      PD(1) = ( RADN(1,1) + 0.001 ) / PAR(1)                            
      P1 = P1F * RADN(1,1) + 0.001                                      
      P2 = P2F * RADN(1,2)                                              
      PAR(2) = P1 + P2                                                  
      PD(2) = P1 / PAR(2)                                               
!                                                                       
      END SUBROUTINE
!=====================================================================  
!                                                                       
      SUBROUTINE RASIT5(TRIB,CTNI,CUNI,RA,Z2,Z0,XDD,ZWIND,UMM,   &
                 RHOAIR,TM,U2,USTAR,DRAG,TA,XAKMS,bps,rib,CU,CT)
!                                                                       
!=======================================================================
!     CUU AND CTT ARE LINEAR  (A SIMPLIFIED VERSION, XUE ET AL. 1991)  
!                                                                      
!----------------------------------------------------------------------
      FS(X) = 66.85 * X                                              
      FT(X) = 0.904 * X                                          
      FV(X) = 0.315 * X                                              
!                                                                      
!     CU AND CT ARE THE FRICTION AND HEAT TRANSFER COEFFICIENTS.    
!     CUN AND CTN ARE THE NEUTRAL FRICTION AND HEAT TRANSFER         
!     COEFFICIENTS.                                                     
!                                                                       
      G2= 0.75                                                          
      G3= 0.75                                                          
      Z22 = Z2                                                          
      ZL = Z2 + 11.785 * Z0                                             
      if(zwind.le.xdd.or.zl.le.xdd) xdd=min(zwind,zl)-0.1
!rr Ratko June 29,2007
!rr   Z2 = XDD + Z0                                                       
!zzq2012Sep09
      Z2 = XDD + Z0
      CUNI = ALOG((ZWIND-XDD)/Z0)/VKC                                 
      IF (ZL.LT.ZWIND) THEN                                             
         XCT1 = ALOG((ZWIND-XDD)/(ZL-XDD))                             
         XCT2 = ALOG((ZL-XDD)/(Z2-XDD))                             
         XCTU2 = ALOG((ZL-XDD)/(Z22-XDD))                             
         CTNI = (XCT1 + G3 * XCT2) / VKC                                
      ELSE                                                              
         XCT2 =  ALOG((ZWIND-XDD)/(Z2-XDD))                           
         XCTU2 =  ALOG((ZWIND-XDD)/(Z22-XDD))                         
         CTNI = G3 * XCT2 /VKC                                          
      END IF                                                            
!        NEUTRAL VALUES OF USTAR AND VENTMF                             
!                                                                       
!rr Ratko June 29,2007 - limit UM
         UM=AMAX1(UMM,2.0)
         USTARN=UM/CUNI                                                 
         VENTN =RHOAIR   /CTNI*USTARN                                   
      IF (ZL.LT.ZWIND) THEN                                             
         U2 = UM - 1. / VKC * USTARN * (XCT1 + G2 *  XCTU2)                                                 
      ELSE                                                              
         U2 = UM - 1. / VKC * USTARN * G2 * XCTU2                      
      END IF                                                            
      if(u2.lt.0.01) u2=0.01
!                                                                       
!     STABILITY BRANCH BASED ON BULK RICHARDSON NUMBER.                 
!                                                                       
      THM=TM*bps
      THVGM= TRIB-THM
      IF (TA.EQ.0.) THVGM = 0.                                          
      RIB  = -THVGM*GRAV*(ZWIND-XDD) / (THM*(UM-U2)**2)
      RIB  = MAX(-10.E0 ,RIB)                                          
      RIB  = MIN( .1643E0 ,RIB)                                        
!                                                                       
!     NON-NEUTRON CORRECTION  (SEE XUE ET AL(1991))                     
      IF(RIB.LT.0.0)THEN                                                
         GRIB = +RIB                                                    
         GRZL = +RIB*(ZL-XDD)/(ZWIND-XDD)                           
         GRZ2 = +RIB*(Z2-XDD)/(ZWIND-XDD)                            
         FVV =  FV(GRIB)                                                
         IF (ZL.LT.ZWIND) THEN                                          
             FTT = FT(GRIB) + (G3-1.) * FT(GRZL) - G3 * FT(GRZ2)        
         ELSE                                                           
             FTT = G3*(FT(GRIB) - FT(GRZ2))                             
         END IF                                                         
         CUI = CUNI + FVV                                               
         CTI = CTNI + FTT                                               
      ELSE                                                              
         RZL = RIB/(ZWIND-XDD)*(ZL-XDD)                               
         RZ2 = RIB/(ZWIND-XDD)*(Z2-XDD)                               
         FVV = FS(RIB)                                                  
         IF (ZL.LT.ZWIND) THEN                                          
             FTT = FS(RIB) + (G3-1) * FS(RZL) - G3 * FS(RZ2)            
         ELSE                                                           
             FTT = G3 * (FS(RIB) - FS(RZ2))                             
         END IF                                                         
  312     CUI = CUNI + FVV                                            
         CTI = CTNI + FTT                                               
      ENDIF                                                             
  310  CONTINUE                                                       
!                                                                       
      CU=1./CUI
      USTAR =UM*CU
!
      XAKMS=USTAR*CU
!
      RAF = CTI / USTAR                                                 
!yliu09May2017 from Qianli --    IF (RAF.LT.0.80) RAF = 0.80
      IF (RAF.LT.0.80) THEN
          RAF = 0.80
          CTI = RAF*USTAR
      ENDIF
      CT = 1./CTI
!!!--
!
      RA  = RAF                                                         
!                                                                       
      UEST  = USTAR                                                     
      DRAG = RHOAIR * UEST*UEST                                         
      Z2 = Z22                                                          
 1010 FORMAT(1X,'RIB,CTI,CUI,CTN,CUN',I10,7E10.3)                       
 1011 FORMAT(1X,'RIB,RAF,USTAR,UM,U2',7E10.3)                           
      END SUBROUTINE
!=======================================================================
!                                                                       
      SUBROUTINE ROOT1(PHSAT,BEE,WWW,PHSOIL)
!                                                         12 AUG 2000   
!=======================================================================
!                                                                       
!    CALCULATION OF SOIL MOISTURE POTENTIALS IN ROOT ZONE OF EACH       
!    VEGETATION LAYER AND SUMMED SOIL+ROOT RESISTANCE                   
!                                                                       
!-----------------------------------------------------------------------
      DIMENSION WWW(3),PHSOIL(3)
!                                                                       
      DO 1000 IL = 1, 3                                                 
      PHSOIL(IL) = PHSAT * AMAX1( 0.05, WWW(IL) ) ** ( - BEE )          
 1000 CONTINUE                                                          
!                                                                       
!-----------------------------------------------------------------------
!     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE USED FOR SOURCE      
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!     PHROOT(1) = PHSOIL(1)-0.01                                        
!                                                                       
!     DO 1200 I = 2 ,3                                                  
!1200 PHROOT(1) = AMAX1( PHROOT(1), PHSOIL(I) )                         
!     PHROOT(2) = PHROOT(1)                                             
!                                                                       
!                                                                       
      END SUBROUTINE
!                                                             
!=======================================================================
!
      SUBROUTINE SFILT(TD  ,TGS  ,TC  ,WWW  ,CAPAC  ,   &
          TD0 ,TGS0 ,TC0 ,WWW0 ,CAPAC0 ,                &
          TDM ,TGSM ,TCM ,WWWM ,CAPACM ,INTG)
!                                                         22 AUGUST 2000
!=======================================================================
!
!     LEAP-FROG IMPLICIT SCHEME
!
!=======================================================================
!
      DIMENSION WWW(3), CAPAC(2)
      DIMENSION WWW0(3),WWWM(3),CAPAC0(2),CAPACM(2)
      real TC0,TGS0,TD0,TCM,TGSM,TDM
      DATA FILTA/.92/
      EPSFLT=0.5 *(1. - FILTA )
!
      IF(INTG.EQ.2) THEN
      TD0 =TD0 +EPSFLT*(TD +TDM -2.0  *TD0)
      TGS0=TGS0+EPSFLT*(TGS+TGSM-2.0  *TGS0)
      TC0 =TC0 +EPSFLT*(TC +TCM -2.0  *TC0)
       IF(WWW0(1).GT.0.0 ) THEN
       WWW0(1)=WWW0(1)+EPSFLT*(WWW(1)+WWWM(1)-2.0*WWW0(1))
          END IF
       IF(WWW0(2).GT.0.0 ) THEN
       WWW0(2)=WWW0(2)+EPSFLT*(WWW(2)+WWWM(2)-2.0*WWW0(2))
          END IF
       IF(WWW0(3).GT.0.0 ) THEN
       WWW0(3)=WWW0(3)+EPSFLT*(WWW(3)+WWWM(3)-2.0*WWW0(3))
          END IF
       IF(CAPAC0(1).GT.0.0 ) THEN
       CAPAC0(1)=CAPAC0(1)+EPSFLT*(CAPAC(1)+CAPACM(1)   &
       -2.0*CAPAC0(1))
          END IF
       IF(CAPAC0(2).GT.0.0 ) THEN
       CAPAC0(2)=CAPAC0(2)+EPSFLT*(CAPAC(2)+CAPACM(2)   &
       -2.0*CAPAC0(2))
         END IF
!       
        TDM=TD0   
        TGSM=TGS0   
        TCM =TC0   
        WWWM(1) =WWW0(1)    
        WWWM(2) =WWW0(2)    
        WWWM(3) =WWW0(3)    
        CAPACM(1)=CAPAC0(1)
        CAPACM(2)=CAPAC0(2)
!
        TD0=TD    
        TGS0=TGS    
        TC0 =TC    
        WWW0(1) =WWW(1)    
        WWW0(2) =WWW(2)    
        WWW0(3) =WWW(3)    
        CAPAC0(1)=CAPAC(1)
        CAPAC0(2)=CAPAC(2)
!
      ELSE
!
        TDM=TD   
        TGSM=TGS   
        TCM =TC   
        WWWM(1) =WWW(1)    
        WWWM(2) =WWW(2)    
        WWWM(3) =WWW(3)    
        CAPACM(1)=CAPAC(1)
        CAPACM(2)=CAPAC(2)
!
        TD0 =TD   
        TGS0=TGS   
        TC0 =TC   
        WWW0(1) =WWW(1)    
        WWW0(2) =WWW(2)    
        WWW0(3) =WWW(3)    
        CAPAC0(1)=CAPAC(1)
        CAPAC0(2)=CAPAC(2)
      END IF
!
      END SUBROUTINE
!=======================================================================
!                                                                       
      SUBROUTINE STOMA1(GREEN,VCOVER,CHIL,ZLT,PAR,PD,EXTK,SUNANG,RST, &
                        RSTPAR,CTLPA)
!                                                         12 AUG 2000
!=======================================================================
!                                                                       
!     CALCULATION OF PAR-LIMITED STOMATAL RESISTANCE                    
!                                                                       
!-----------------------------------------------------------------------
      DIMENSION GREEN(2), VCOVER(2), ZLT(2), CHIL(2)
      DIMENSION PAR(2), PD(2), EXTK(2,3,2), RST(2), RSTPAR(2,3)
!                                                                       
      DO 1000 IVEG = 1, 2                                               
!                                                                       
      AT = ZLT(IVEG) / VCOVER(IVEG)                                     
!                                                                       
      IF (SUNANG .LE. 0.02) THEN                                        
         XABC = RSTPAR(IVEG,1) / RSTPAR(IVEG,2) + RSTPAR(IVEG,3)        
         RST(IVEG) = 0.5 / XABC * AT                                  
         IF (RST(IVEG) .LT. 0.) RST(IVEG) = 0.00001                     
         GO TO 1010                                                     
      END IF                                                            
!                                                                       
      GAMMA = ( RSTPAR(IVEG,1) + RSTPAR(IVEG,2) * RSTPAR(IVEG,3) ) /  & 
                RSTPAR(IVEG,3)                                          
!                                                                       
      POWER1 = AMIN1( 50., AT * EXTK(IVEG,1,1) )                        
      POWER2 = AMIN1( 50., AT * EXTK(IVEG,1,2) )                        
!                                                                       
!-----------------------------------------------------------------------
!     ROSS INCLINATION FUNCTION                                         
!-----------------------------------------------------------------------
!                                                                       
      AA = 0.5 - 0.633 * CHIL(IVEG)- 0.33 * CHIL(IVEG)* CHIL(IVEG)      
      BB = 0.877 * ( 1. - 2. * AA )                                     
!                                                                       
!-----------------------------------------------------------------------
!     COMBINED ESTIMATE OF K-PAR USING WEIGHTS FOR DIFFERENT COMPONENTS 
!-----------------------------------------------------------------------
!                                                                       
      ZAT = ALOG( ( EXP(-POWER1) + 1. )/2. ) * PD(IVEG) &              
            / ( POWER1/AT )                                             
      ZAT = ZAT + ALOG( ( EXP(-POWER2) + 1. )/2. )      &               
       * ( 1. - PD(IVEG) ) / ( POWER2/AT )                              
!                                                                       
      POW1 = AMIN1( 50., (POWER1*ZAT/AT) )                              
      POW2 = AMIN1( 50., (POWER2*ZAT/AT) )                              
!                                                                       
      ZK = 1. / ZAT * ALOG( PD(IVEG) * EXP ( POW1 )     &               
            + ( 1. - PD(IVEG) ) * EXP ( POW2 ) )                        
!                                                                       
!                                                                       
      POW = AMIN1( 50., ZK*AT )                                         
      EKAT = EXP ( POW )                                                
!                                                                       
      AVFLUX = PAR(IVEG) * ( PD(IVEG) / SUNANG * ( AA + BB * SUNANG ) & 
            + ( 1. - PD(IVEG) )*( BB / 3. + AA * 1.5                  & 
            + BB / 4. * PIE ))                                          
!                                                                       
      RHO4 = GAMMA / AVFLUX                                             
!                                                                       
      RST(IVEG) = RSTPAR(IVEG,2)/GAMMA * ALOG(( RHO4 * EKAT + 1. ) /  & 
                    ( RHO4 + 1. ) )                                     
      RST(IVEG) = RST(IVEG) - ALOG (( RHO4 + 1. / EKAT ) /            & 
                    ( RHO4 + 1. ) )                                     
      RST(IVEG) = RST(IVEG) / ( ZK * RSTPAR(IVEG,3) )                   
!                                                                       
!---------------------------------------------------------------------- 
!     MODIFICATIONS FOR GREEN FRACTION : RST UPRIGHT                    
!---------------------------------------------------------------------- 
!                                                                       
 1010  RST(IVEG) = 1. / ( RST(IVEG) * GREEN(IVEG) + 0.0000001)         
 1000  CONTINUE                                                        
!                                                                       
      RST(1) = RST(1) * CTLPA
      END SUBROUTINE
!======================================================================
!
!zzq20110523      SUBROUTINE STRES1 (IFIRST, RSTM,ROOTP,
      SUBROUTINE STRES1 (IFIRST, RSTM,ROOTP,www,phsat,bee,          &
           RSTFAC,RST,TC,ETC,RB,TGS,ETGS,RD,TU,TL,TOPT,EA,          &
           DEFAC,PH1,PH2,NROOT,ZDEPTH,PHSOIL,ROOTD,VCOVER,DROP,zfw, &
           rst1,rst2,rst3,rst4)
!khs &     DEFAC,PH1,PH2,NROOT,ZDEPTH,PHSOIL,ROOTD,VCOVER,DROP)
!
!======================================================================
!
!     CALCULATION OF ADJUSTMENT TO LIGHT DEPENDENT STOMATAL RESISTANCE
!     BY TEMPERATURE, HUMIDITY AND STRESS FACTORS
!     SIMPLIFIED SEE XUE ET AL(1991)
!
!         RSTFAC(IVEG,1) = FD
!         RSTFAC(IVEG,2) = FP
!         RSTFAC(IVEG,3) = FT
!         RSTFAC(IVEG,4) = FTPD
!
!----------------------------------------------------------------------
      DIMENSION RSTM(2), DEP(3), VCOVER(2)
      DIMENSION TOPT(2), TL(2), TU(2), DEFAC(2)
      DIMENSION PH1(2), PH2(2), RST(2), RSTFAC(2,4)
      DIMENSION ROOTD(2), ROOTP(3), ZDEPTH(3),PHSOIL(3)
!zzq Aug19,2009
      DIMENSION XDRR(3),WWW(3)
      REAL, DIMENSION (3)   :: PHSOILA
!
!----------------------------------------------------------------------
!     HUMIDITY, TEMPERATURE AND TRANSPIRATION FACTORS
!----------------------------------------------------------------------
!
      DO 1000 IVEG = 1, 2
!
      TV = TC
      ETV = ETC
      RAIR = RB * 2.
      IF ( IVEG .EQ. 1 ) GO TO 100
      TV = TGS
      ETV = ETGS
      RAIR = RD
  100   CONTINUE
!
      TV = AMIN1 ( ( TU(IVEG) - 0.1 ), TV )
      TV = AMAX1 ( ( TL(IVEG) + 0.1 ), TV )
!
      IF( IFIRST .EQ. 0 ) GO TO 200
      RSTM(IVEG) = RST(IVEG)
      D2 = ( TU(IVEG) - TOPT(IVEG) ) / ( TOPT(IVEG) - TL(IVEG) )
      D1 = 1. /(( TOPT(IVEG) - TL(IVEG) )*   &
              EXP( ALOG( TU(IVEG) - TOPT(IVEG))*D2))
      RSTFAC(IVEG,3) = D1*( TV-TL(IVEG)) * EXP(ALOG(TU(IVEG)-TV)*D2)
!
      IF (RSTFAC(IVEG,3).LT.0.) RSTFAC(IVEG,3) = 0.
      IF (RSTFAC(IVEG,3).GT.1.) RSTFAC(IVEG,3) = 1.
!
!----------------------------------------------------------------------
!      SIMPLIFIED CALCULATION OF LEAF WATER POTENTIAL FACTOR , FP
!----------------------------------------------------------------------
!
      IF (NROOT.EQ.1) THEN
      XROT = ROOTD(1)
      DO 7400 I = 1, 3
 7400 DEP(I) = 0.
      DO 7500 I = 1, 3
      DEP(I) = MIN(ZDEPTH(I), XROT)
      XROT = XROT - ZDEPTH(I)
      IF (XROT.LE.0.) GO TO 7410
 7500 CONTINUE
 7410 CONTINUE
!---------------------------------------------------------------------
!yliu modified based on Jiwoo's wilting point adjustment 25Jun2016
      DO 8000 IL=1,3
        PHSOILA(IL)=-AMIN1(abs(PHSOIL(IL)),abs(EXP(PH2(IVEG)*1.0)))
 8000 CONTINUE
      XDR = (PHSOILA(1) * DEP(1) + PHSOILA(2) * DEP(2) &
     &      +PHSOILA(3) * DEP(3)) /ROOTD(1)

      ELSE
      DO 8010 IL=1,3
      PHSOILA(IL)=-AMIN1(abs(PHSOIL(IL)),abs(EXP(PH2(IVEG)*1.0)))
 8010 CONTINUE
      XDR = (PHSOILA(1) * DEP(1) + PHSOILA(2) * DEP(2) &
     &      +PHSOILA(3) * DEP(3)) /ROOTD(1)
      END IF
!
      XDR = - XDR
      IF (XDR .LE. 0.001) XDR = 0.001
      XDR = ALOG (XDR)
!---------------------------------------------------------------------
!
      RSTFAC(IVEG,2) = 1. - EXP(- PH1(IVEG) * (PH2(IVEG) - XDR))
      IF (RSTFAC(IVEG,2).GT.1.) RSTFAC(IVEG,2) = 1.
      IF (RSTFAC(IVEG,2).LT.0.) RSTFAC(IVEG,2) = 0.
!
!khs..for CO2 routines
!zzq20110923      zfw = amax1(rstfac(iveg,2),1.e-5)
!khs
  200   RST(IVEG) = RSTM(IVEG)
!     EPOT = AMAX1(0.0001,(ETV-EA))
!
!               ***** PJS mod 10/9/92 *****
! ***** based on Verma FIFE-87 function for C4 grasses *****
!
!     RSTFAC(IVEG,1) = 1. - DROP * DEFAC(IVEG)
!     rstfac(iveg,1) = 1./ ( 1 + defac(iveg)*drop )
!     IF (RSTFAC(IVEG,1).GT.1.) RSTFAC(IVEG,1) = 1.
!
!----------------------------------------------------------------------
!     VALUE OF FP FOUND
!----------------------------------------------------------------------
!
!khs300   FTPD = RSTFAC(IVEG,1) * RSTFAC(IVEG,2) * RSTFAC(IVEG,3)
  300   FTPD = RSTFAC(IVEG,2)
      RSTFAC(IVEG,4) = AMAX1( FTPD, 0.00001 )
!----------------------------------------------------------------------
!
      RST(IVEG) = RST(IVEG) / RSTFAC(IVEG,4) / VCOVER(IVEG)
!
!khs  RST(IVEG) = AMIN1( RST(IVEG), 100000.)
      RST(IVEG) = AMIN1( RST(IVEG), 10000. )
!
 1000  CONTINUE
!zzq20110814
      zfw = amax1(rstfac(1,2),1.e-4)
      rst1 = rstfac(1,1)
      rst2 = rstfac(1,2)
      rst3 = rstfac(1,3)
      rst4 = rstfac(1,4)
!
      END SUBROUTINE
!======================================================================
!
      SUBROUTINE TEMRS1(                                               &
         DTT,TC,TGS,TD,TA,TM,QM,EM_IN,PSURF,ZLAT,MONTH,WWW,CAPAC,SATCAP,    &
         DTC,DTG,RA,RST,ZDEPTH,BEE,PHSAT,POROS,XDD,Z0,RDC,RBC,VCOVER,  &
         Z2,ZLT,DEFAC,TU,TL,TOPT,RSTFAC,NROOT,ROOTD,PHSOIL,ROOTP,      &
         PH1,PH2,ECT,ECI,EGT,EGI,EGS,HC,HG,EC,EG,EA,RADT,CHF,SHF,      &
         ALBEDO,ZLWUP,THERMK,RHOAIR,ZWIND,UM,USTAR,DRAG,CCX,CG,        &
!song 11.22.2011
         rd,                                                           &
         pilphr,                                                       & !Huilin for hr output  Feb 7th 2018
         BPS,XXT,XXQ,XAKMS,rcc,rib,CU,CT,flup,zfw,rb,                  &
         rst1,rst2,rst3,rst4,rsoil)
!khs &   BPS,XXT,XXQ,XAKMS,rcc,rib,CU,flup)
!rr  &   BPS,XXT,XXQ,XAKMS,rcc,rib,CU)
!                                                          11 AUG 2000 
!=======================================================================
!     A MODIFIED SIMPLIFIED VERSION (XUE ET AL. 1991)
!     FLUX COUPLING
!     CORE ROUTINE: CALCULATION OF CANOPY AND GROUND TEMPERATURE        
!     INCREMENTS OVER TIME STEP, FLUXES DERIVED.                        
!-----------------------------------------------------------------------
!                                                                       
!     SUBROUTINES IN THIS SUBROUTINE:
!                                 STRES1                                
!                                 STOMA1
!                                 NEWTON
!     SUBPARTS IN THIS SUBROUTINE:                                
!     -------------------------   DELRN                                 
!                                 DELHF                                 
!                                 DELEF                                 
!----------------------------------------------------------------------
      DIMENSION ZINC(3), A2(3), Y1(3)  , ITEX(3)
      DIMENSION RSTM(2)                                                 
!                                                                       
!---------------------------------------------------------------------- 
!     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNCTION OF TEMPERATURE     
!     GE(X) IS D E(X) / D ( TEMP )                                      
!---------------------------------------------------------------------- 
!                                                                       
      DIMENSION WWW(3), CAPAC(2), SATCAP(2), ZDEPTH(3)
      DIMENSION VCOVER(2), ZLT(2), RADT(2),ALBEDO(2,3,2)
      DIMENSION TOPT(2), TL(2), TU(2), DEFAC(2)
      DIMENSION PH1(2), PH2(2), RST(2), RSTFAC(2,4)
      DIMENSION ROOTD(2), ROOTP(3), PHSOIL(3)
      REAL ZLAT
      INTEGER MONTH
      REAL :: EM

      E(X) = EXP( 21.18123 - 5418. / X ) / .622                         
      GE(X) = EXP( 21.18123 - 5418. / X ) * 5418. / (X*X) / .622 

      ETC   = E(TC)                                                     
      ETGS  = E(TGS)                                                    
      GETC  = GE(TC)                                                    
      GETGS = GE(TGS)                                                   
      PSY      = CPAIR / HLAT * PSURF/100. / .622                    
!
      RCP = RHOAIR * CPAIR                                              
!     RADD = 44.                                                        
!                                                                       
      WC = AMIN1( 1., CAPAC(1)/SATCAP(1) ) 
      WG = AMIN1( 1., CAPAC(2)/SATCAP(2) )                              

      EM  = EM_IN
!                                                                       
!---------------------------------------------------------------------- 
!      RSOIL FUNCTION FROM FIT TO CAMILLO AND GURNEY (1984) DATA.       
!      WETNESS OF UPPER 0.5 CM OF SOIL CALCULATED FROM APPROXIMATION    
!      TO MILLY FLOW EQUATION WITH REDUCED (1/50 ) CONDUCTIVITY IN      
!      TOP LAYER.                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     WT = WWW(1) + 0.75 * ZDEPTH(1) / ( ZDEPTH(1) + ZDEPTH(2) )        
!    &     * (WWW(1) - (WWW(2)**2)/WWW(1) ) / 2. * 50.                  
!     FAC = AMIN1( WT, 0.99 )                                           
!     FAC = AMAX1( FAC, WWW(1) * 0.1 )                            
!
!------------------------------------------------------------
!     Y.K. Xue changed Jan.18,1994
!------------------------------------------------------------
!                                 
      FAC = AMIN1( www(1), 0.99 )                                      
      FAC = AMAX1( FAC, 0.02 )                            
      RSOIL =  101840. * (1. - fac ** 0.0027) 
! 
!------------------------------------------------------------  
!                                                                       
      FAC = AMIN1( WWW(1), 1. )
      FAC = AMAX1( FAC, 0.02  )  
      PSIT = PHSAT * FAC ** (- BEE )                                    
      ARGG = AMAX1(-10.,(PSIT*GRAV/461.5/TGS))                         
      HR = EXP(ARGG)                                                    
      pilphr=hr
!                                                                       
!---------------------------------------------------------------------- 
!     ALTERATION OF AERODYNAMIC TRANSFER PROPERTIES IN CASE OF SNOW     
!     ACCUMULATION.                                                     
!---------------------------------------------------------------------- 
!                                                                       
      RESD = XDD                                                       
      RESZ0 = Z0                                                        
      RESRDC = RDC                                                      
      RESRBC = RBC                                                      
      RESV2 = VCOVER(2)                                                 
!                                                                       
      IF ( TGS .GT. TF ) GO TO 100                                      
!                                                                       
      SDEP = CAPAC(2) * 5.                                              
      SDEP = AMIN1( SDEP, (Z2*0.95) )                                   
      XDD = Z2 - ( Z2-XDD ) / Z2 * ( Z2 - SDEP )                        
      Z0 = Z0 / ( Z2-RESD ) * ( Z2-XDD )                               
      RDC = RDC * ( Z2-SDEP ) / Z2                                      
      RBC = RBC * Z2 / ( Z2-SDEP )                                      
      VCOVER(2) = 1.                                                    
      WG = AMIN1( 1., CAPAC(2) / 0.004 )                                
      RST(2) = RSOIL                                                    
!                                                                       
  100   CONTINUE                                                     
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!      CALCULATION OF EA, TA, RA, RB, RD AND SOIL MOISTURE STRESS       
!      FOR THE BEGINNING OF THE TIME STEP                               
!                                                                       
!---------------------------------------------------------------------- 
      IFIRST = 1                                                        
      ICOUNT = 0                                                        
      IONCE = 1                                                         
!                                                                       
!     TA = TGS                                                          
      TGEN = TGS
      TCEN = TC
      FC = 1.
      FG = 1.
      TRIB = TA
      EA = EM                                                           
      HT = 0.                                                           
      IONCE = 0                                                         
!                                                                       
 1000  CONTINUE                                                        
      ICOUNT = ICOUNT + 1                                               
      
      CALL RASIT5(TRIB,CTNI,CUNI,RA,Z2,Z0,XDD,ZWIND,UM,   &
            RHOAIR,TM,U2,USTAR,DRAG,TA,XAKMS,bps,rib,CU,CT)
!     *******   IF ( IFIRST .EQ. 1 ) CALL RBRD1 ******
      IF ( IFIRST .EQ. 1 ) THEN
          TCTA = TC - TA
          RB  = 1.0/(SQRT(U2)/RBC+ZLT(1)*.004)

          X1 = TEMDIF
!
          TGTA = TGS- TA
          TEMDIF = ( TGTA + SQRT(TGTA*TGTA) ) / 2. + 0.1
          FIH = SQRT( 1. + 9.*GRAV*TEMDIF*Z2/TGS/( U2*U2) )
          RD  = RDC / U2 / FIH
      END IF
!     ******    END RBRD    ********** 
      D1 = 1./RA + 1./RB + 1./RD                                        
      TA = ( TGS/RD + TC/RB + TM/RA *bps) / D1                        
      HT = ( TA - TM ) * RCP / RA                                       
      RCC = RST(1)*FC + 2. * RB                                         
      COC = (1.-WC)/RCC + WC/(2.*RB)                                    
      RG = RST(2)*FG                                                    
      RSURF = RSOIL*FG                                                  
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HR  &  
             + VCOVER(2)/(RSURF+RD+44.)*HR
      COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)     &  
             + VCOVER(2)/(RSURF+RD+44.)
      COG1 = COG1 + WG/RD * VCOVER(2) 
      COG2 = COG2 + WG/RD * VCOVER(2)
      D2 = 1./RA + COC + COG2 
      TOP = COC * ETC + COG1 * ETGS + EM / RA                           
      EA = TOP / D2
      DROP = AMAX1( 0., (E(TA)-EA) )
!---------------------------------------------------------------------- 
!
!zzq20110523      CALL STRES1 (IFIRST, RSTM,ROOTP,
      CALL STRES1 (IFIRST, RSTM,ROOTP,www,phsat,bee,                &
           RSTFAC,RST,TC,ETC,RB,TGS,ETGS,RD,TU,TL,TOPT,EA,          &
           DEFAC,PH1,PH2,NROOT,ZDEPTH,PHSOIL,ROOTD,VCOVER,DROP,zfw  &
          ,rst1,rst2,rst3,rst4)
!khs &     DEFAC,PH1,PH2,NROOT,ZDEPTH,PHSOIL,ROOTD,VCOVER,DROP)
!---------------------------------------------------------------------- 
!                                                                       
      IFIRST = 0                                                        
      ERIB = EA                                                         
      TRIB = TA                                                         
!CC                                                                     
      IF ( ICOUNT .LE. 4 ) GO TO 1000                                   
!                                                                       
!---------------------------------------------------------------------- 
!
!     CALL DELRN ( RNCDTC, RNCDTG, RNGDTG, RNGDTC )                     
!
      TC3 = TC * TC * TC
      TG3 = TGS * TGS * TGS
      FAC1 = ( 1. - ALBEDO(1,3,2) ) * ( 1.-THERMK ) * VCOVER(1)
      FAC2 =   1. - ALBEDO(2,3,2)
!                
      RNCDTC = - 2. * 4. * FAC1 * STEFAN * TC3
      RNCDTG = 4. * FAC1 * FAC2 * STEFAN * TG3
!                
      RNGDTG = - 4. * FAC2 * STEFAN * TG3
      RNGDTC = 4. * FAC1 * FAC2 * STEFAN * TC3
!                
!---------------------------------------------------------------------- 
!                                                                       
!     DEW CALCULATION : DEW CONDITION IS SET AT BEGINNING OF TIME STEP. 
!     IF SURFACE CHANGES STATE DURING TIME STEP, LATENT HEAT FLUX IS    
!     SET TO ZERO.                                                      
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      IF ( EA .GT. ETC ) FC = 0.                                        
      IF ( EA .GT. ETGS) FG = 0.                                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     WET FRACTION EXHAUSTION TEST : IF CAPAC(X) IS EXHAUSTED IN        
!     A TIME STEP, INTERCEPTION LOSS IS LIMITED TO CAPAC(X).            
!                                                                       
!---------------------------------------------------------------------- 
!     START OF NON-NEUTRAL RESISTANCE CALCULATION LOOP                  
!---------------------------------------------------------------------- 
!                                                                       
      I = 0                                                             
!                                                                      
!    ----- INITIALIZE NEWTON-RAPHSON ITERATIVE ROUTINE FOR RASIT 3,5,8  
                    NOX = 0                                             
                 NONPOS = 1                                             
                  IWALK = 0                                             
                     LX = 2                                             
                   FINC = 1.                                            
                   ITEX(LX) = 0.                                        
                   ZINC(LX) = 0.                                        
                   A2(LX)   = 0.                                       
                   Y1(LX)   = 0.                                        
 2000  CONTINUE                                                       
!                                                                       
      CALL RASIT5(TRIB,CTNI,CUNI,RA,Z2,Z0,XDD,ZWIND,UM,   &
            RHOAIR,TM,U2,USTAR,DRAG,TA,XAKMS,bps,rib,CU,CT)
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     CALL DELHF ( HCDTC, HCDTG, HGDTG, HGDTC )                         
!
      RCP = RHOAIR * CPAIR
      D1 = 1./RA + 1./RB + 1./RD
      TA = ( TGS/RD + TC/RB + TM/RA *bps) / D1
!                
      HC = RCP * ( TC - TA ) / RB * DTT
      HG = RCP * ( TGS - TA ) / RD * DTT
!----------------------------------------------------------------------
!     N.B. FLUXES EXPRESSED IN JOULES M-2
!----------------------------------------------------------------------
!                
      HCDTC = RCP / RB * ( 1./RA + 1./RD ) / D1
      HCDTG = - RCP / ( RB * RD ) / D1
! FOR TM
      HCDTM = - RCP / ( RB * RA ) / D1 * BPS
!                
      HGDTG = RCP / RD * ( 1./RA + 1./RB ) / D1
      HGDTC = - RCP / ( RD * RB ) / D1
! FOR TM
      HGDTM = - RCP / ( RD * RA ) / D1 *BPS
!                
!     CALL DELEF ( ECDTC,ECDTG,EGDTG,EGDTC,DEADTC,DEADTG,EC,EG, 
!    &             WC, WG, FC, FG, HR )                                 
!
!     RCP = RHOAIR * CPAIR
!----------------------------------------------------------------------
!     MODIFICATION FOR SOIL DRYNESS : HR = REL. HUMIDITY IN TOP LAYER
!----------------------------------------------------------------------
!                
      HRR = HR   
      IF ( FG .LT. .5 ) HRR = 1.
!                
!khs...
      if (rst(1).ge.10000.) then
      rcc = 10000.
      else
      RCC = RST(1)*FC + 2. * RB
      endif
!khs...
      COC = (1.-WC)/RCC + WC/(2.*RB)
      RG = RST(2)*FG
      RSURF = RSOIL*FG
      COG1 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)*HRR  &
           + VCOVER(2)/(RSURF+RD+44.)*HRR
      COG2 = VCOVER(2)*(1.-WG)/(RG+RD)+(1.-VCOVER(2))/(RSURF+RD)      &
           + VCOVER(2)/(RSURF+RD+44.)
      COG1 = COG1 + WG/RD * VCOVER(2)
      COG2 = COG2 + WG/RD * VCOVER(2)
!                
      D2 = 1./RA + COC + COG2
      TOP = COC * ETC + COG1 * ETGS + EM/RA
      EA = TOP / D2
!                
      EC = ( ETC - EA ) * COC * RCP/PSY * DTT
!                
      EG = ( ETGS*COG1 - EA*COG2 ) * RCP/PSY * DTT
!                
      DEADTC = GETC * COC / D2
      DEADTG = GETGS * COG1 / D2
!                
      ECDTC = ( GETC - DEADTC ) * COC * RCP / PSY
      ECDTG = - DEADTG * COC * RCP / PSY
!                
      EGDTG = ( GETGS*COG1 - DEADTG*COG2 ) * RCP / PSY
      EGDTC = - DEADTC * COG2 * RCP / PSY
!   FOR QM
!      DEADQM = 0.622 * PSURF /( (0.622+QM)**2 * RA * D2 )    
!ljw - 6/3/2015, bugfix
      DEADQM = 0.622 * PSURF/100. /( (0.622+QM)**2 * RA * D2 )    
      ECDQM =        -DEADQM * COC * RCP / PSY
      EGDQM =        -DEADQM * COG2 * RCP / PSY
!   FOR YPDATING TM AND QM
      AK = 1/ RCP / BPS
      AH = 1/ (HLAT*RHOAIR)
      XXT = -DTT * AK * (HGDTM+HCDTM)
      XXQ = -DTT * AH * (EGDQM+ECDQM)
!             
!---------------------------------------------------------------------- 
!                                                                       
!     CALCULATION OF COEFFICIENTS OF TEMPERATURE TENDENCY EQUATIONS     
!        C - CANOPY                                                     
!        G - GROUND                                                     
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      CCODTC = CCX / DTT - RNCDTC + HCDTC + ECDTC                       
      CCODTG = - RNCDTG + HCDTG + ECDTG                                 
      CCORHS = RADT(1) - ( HC + EC ) / DTT                              
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      GCODTG = CG / DTT + TIMCON*CG*2. - RNGDTG + HGDTG + EGDTG         
      GCODTC = - RNGDTC + HGDTC + EGDTC                                 
      GCORHS = RADT(2) - TIMCON*CG*2. * ( TGS -TD ) - ( HG + EG ) / DTT 
!                                                                       
      DENOM = CCODTC * GCODTG - CCODTG * GCODTC                         
!                                                                       
      DTC = ( CCORHS * GCODTG - CCODTG * GCORHS ) / DENOM               
      DTG = ( CCODTC * GCORHS - CCORHS * GCODTC ) / DENOM               
!                                                                       
!---------------------------------------------------------------------- 
!     CHECK IF INTERCEPTION LOSS TERM HAS EXCEEDED CANOPY STORAGE       
!---------------------------------------------------------------------- 
!                                                                       
      ECPOT = ( (ETC - EA) + (GETC - DEADTC)*DTC - DEADTG*DTG )         
      ECI = ECPOT * WC /(2.*RB) * RCP/PSY * DTT                         
      ECIDIF=AMAX1(0.0,(ECI-CAPAC(1)*1.E3*HLAT))                        
      ECI   =AMIN1(ECI,(    CAPAC(1)*1.E3*HLAT))                        
!                                                                       
      EGPOT = ( (ETGS - EA) + (GETGS - DEADTG)*DTG - DEADTC*DTC )       
      EGI = EGPOT * VCOVER(2) * WG/RD * RCP/PSY * DTT                   
      EGIDIF=AMAX1(0.0,(EGI-CAPAC(2)*1.E3*HLAT))                        
      EGI   =AMIN1(EGI,(    CAPAC(2)*1.E3*HLAT))                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      TGEN = TGS + DTG                                                  
      TCEN = TC + DTC                                                   
      D1 = 1./RA + 1./RB + 1./RD                                        
      TAEN = ( TGEN / RD + TCEN / RB + TM / RA *bps) / D1            
!                                                                       
      HEND = ( TAEN - TM ) * RCP / RA + (ECIDIF + EGIDIF)/DTT           
      Y= TRIB - TAEN                                                    
      I = I + 1                                                         
      HT   = HEND                                                       
      IF ( I .GT. ITRUNK ) GO TO 200                                    
!                                                                       
!zzq      CALL NEWTON(TRIB,Y,FINC,NOX,NONPOS,IWALK,LX)                      
!zzq080321      CALL NEWTON(TRIB,Y,FINC,NOX,NONPOS,IWALK,LX,ZINC,A2,Y1,ITER)
      CALL NEWTON(TRIB,Y,FINC,NOX,NONPOS,IWALK,LX,ZINC,A2,Y1,ITEX)

      IF(NOX.NE.1)GO TO 2000                                            
!                                                                       
  200   CONTINUE                                                     
!     IQIN = IQIN + I                                                   
!     IF (I.GT.10) IQIN1 = IQIN1 + 1                                    
 1010  FORMAT(1X,I3,1X,'TR1B,Y,RA,RIB,EGDF',7E11.4)                    
 1011  FORMAT(1X,'HEND,HT,Y,TA,TC,TG,ECDF',8E11.4)                    
 1012  FORMAT(5X,I10,I5)                                              
!      WRITE(6,1014)  RIB                                                
 1014  FORMAT(5X,F12.5)                                               
!                                                                       
!---------------------------------------------------------------------- 
!     EXIT FROM NON-NEUTRAL CALCULATION                                 
!                                                                       
!     EVAPOTRANSPIRATION FLUXES CALCULATED FIRST ( J M-2 )              
!---------------------------------------------------------------------- 
!                                                                       
      HRR = HR                                                          
      IF ( FG .LT. .5 ) HRR = 1.                                        
      RSURF = RSOIL*FG                                                  
!                                                                       
      COCT = (1.-WC)/RCC                                                
      COGT = VCOVER(2) * (1.-WG)/( RG + RD )                            
      COGS1 = (1.-VCOVER(2)) / ( RD + RSURF ) * HRR       &             
              + VCOVER(2) / ( RD + RSURF + 44.) * HRR                   
      COGS2 = COGS1 / HRR                                               
!                                                                       
      ECT = ECPOT * COCT * RCP/PSY * DTT                                
!                                                                       
      EGT = EGPOT * COGT * RCP/PSY * DTT                                
      EGS = (ETGS + GETGS*DTG ) * COGS1                   &             
            - ( EA + DEADTG*DTG + DEADTC*DTC ) * COGS2                 
      EGS = EGS * RCP/PSY * DTT                                         
      EGSMAX = WWW(1) / 2. * ZDEPTH(1) * POROS * HLAT * 1000.           
      EGIADD = AMAX1( 0., EGS - EGSMAX )                                
      EGS = AMIN1 ( EGS, EGSMAX )                                       
      EGIDIF = EGIDIF + EGIADD                                          
!                                                                       
!---------------------------------------------------------------------- 
!     SENSIBLE HEAT FLUX CALCULATED WITH LATENT HEAT FLUX CORRECTION    
!---------------------------------------------------------------------- 
      HC = HC + (HCDTC*DTC + HCDTG*DTG)*DTT + ECIDIF                    
      HG = HG + (HGDTC*DTC + HGDTG*DTG)*DTT + EGIDIF                    
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     TEST OF DEW CONDITION. LATENT HEAT FLUXES SET TO ZERO IF SIGN     
!     OF FLUX CHANGES OVER TIME STEP.EXCESS ENERGY DONATED TO SENSIBLE  
!     HEAT FLUX.                                                        
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      ECF = SIGN( 1., ECPOT )                                           
      EGF = SIGN( 1., EGPOT )                                           
      DEWC = FC * 2. - 1.                                               
      DEWG = FG * 2. - 1.                                               
!                                                                       
      IF(DEWC*ECF.GT.0.0) GO TO 300                                     
      HC = HC + ECI + ECT                                               
      ECI = 0.                                                          
      ECT = 0.                                                          
  300   IF(DEWG*EGF.GT.0.0) GO TO 400                                 
      HG = HG + EGS + EGI + EGT                                         
      EGS = 0.                                                          
      EGI = 0.                                                          
      EGT = 0.                                                          
  400   CONTINUE                                                       
!                                                                       
      EC = ECI + ECT                                                    
      EG = EGT + EGS + EGI                                              
!                                                                       
!---------------------------------------------------------------------- 
!     ADJUSTMENT OF TEMPERATURES AND VAPOR PRESSURE , CALCULATION OF    
!     SENSIBLE HEAT FLUXES.                                             
!---------------------------------------------------------------------- 
!                                                                       
      TC  = TCEN                                                        
      TGS = TGEN                                                        
      TA  = TAEN                                                        
      EA = EA + DEADTC*DTC + DEADTG*DTG                                 
!                                                                       
      RADT(1) = RADT(1) + RNCDTC*DTC + RNCDTG*DTG                       
      RADT(2) = RADT(2) + RNGDTC*DTC + RNGDTG*DTG 
!KHS====================================================================
      FLUP = FLUP - (RNCDTC+RNGDTC)*DTC   &
                          - (RNCDTG+RNGDTG)*DTG
!KHS====================================================================
!
! ** simulated net all-wave radiation **
!
!     sibnet(nmm,ndd,nhh) = RADT(1) + RADT(2)                   
!                                                                       
      CHF = CCX / DTT * DTC                                             
      SHF = CG / DTT * DTG + TIMCON*CG*2. * ( TGS - TD )                
!                                                                       
      ZLWUP = ZLWUP - RNCDTC * DTC / 2.   &                             
                    - RNGDTG * DTG * (1.-VCOVER(1)*(1.-THERMK) )        
!                                                                       
      IF ( TGS .GT. TF ) GO TO 500                                      
      EGS = EG - EGI                                                    
      EGT = 0.                                                          
  500   CONTINUE                                                       
!                                                                       
      VCOVER(2) = RESV2                                                 
      XDD = RESD                                                       
      Z0 = RESZ0                                                        
      RDC = RESRDC                                                      
      RBC = RESRBC                                                      
!                                                                       
      END SUBROUTINE
!=======================================================================
!                                                                     
      SUBROUTINE UPDAT1(DTT,TC,TGS,TD,CAPAC,DTC,DTG,ECT,ECI,EGT,EGI, &
         EGS,EG,HC,HG,HFLUX,ETMASS,FILTR,SOILDIF,SOILDRA,ROFF,       &
         RNOFFB,RNOFFS,NROOT,ROOTD,ROOTP,POROS,BEE,SATCO,SLOPE,      &
         PHSAT,ZDEPTH,WWW,CCX,CG,CHF,SHF,SMELT)
!                                                         12 AUGUST 2000
!=======================================================================
!                                                                       
!     UPDATING OF SOIL MOISTURE STORES AND INTERCEPTION CAPACITY        
!                                                                       
!----------------------------------------------------------------------
      DIMENSION WWW(3), CAPAC(2),EF(3),SNOWW(2)
      DIMENSION ROOTD(2), ZDEPTH(3), ROOTP(3)
      dimension temw(3), temwp(3), temwpp(3),           &             
                aaa(2) , bbb(2)  , ccc(2)   , qqq(2)                
!                                                                       
!---------------------------------------------------------------------- 
!     EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY       
!     ( HLAT*1000.) LOSS IS IN M M-2                                    
!     MASS TERMS ARE IN KG M-2 DT-1                                     
!---------------------------------------------------------------------- 
!                                                                       
      SNOFAC = HLAT / ( HLAT + SNOMEL /1000. )                          
      FACKS = 1.                                                        
      IF ( (TC-DTC) .LE. TF ) FACKS = SNOFAC                            
      IF ( (ECT+ECI) .GT. 0.) GO TO 100                                 
      ECI = ECT + ECI                                                   
      ECT = 0.                                                          
      FACKS = 1. / FACKS                                                
  100   CAPAC(1)=CAPAC(1) - ECI*FACKS/HLAT/1000.                       
!                                                                       
      ECMASS = ( ECT + ECI * FACKS ) / HLAT                             
!                                                                       
      FACKS = 1.                                                        
      IF ( (TGS-DTG) .LE. TF ) FACKS = SNOFAC                           
      IF ( (EGT+EGI) .GT. 0. ) GO TO 200                                
      EGI = EGT + EGI                                                   
      EGT = 0.                                                          
      FACKS = 1. / FACKS                                                
  200   CAPAC(2)=CAPAC(2) - EGI*FACKS/HLAT/1000.                      
!                                                                       
      EGMASS = ( EGT + EGS + EGI * FACKS ) / HLAT                       
!                                                                       
      ETMASS = ECMASS + EGMASS                                          
!                                                                       
!rr   HFLUX = ( HC + HG ) / DTT ( output for Eta)
      HFLUX = ( HC + HG )       ! output for GCM
!rr
!                                                                       
!---------------------------------------------------------------------- 
!      DUMPING OF SMALL CAPAC VALUES ONTO SOIL SURFACE STORE            
!---------------------------------------------------------------------- 
!                                                                       
      DO 1000 IVEG = 1, 2                                               
      IF ( CAPAC(IVEG) .GT. 0.000001 ) GO TO 300                        
      FILTR = FILTR + CAPAC(IVEG)                
      WWW(1) = WWW(1) + CAPAC(IVEG) / ( POROS*ZDEPTH(1) )               
      CAPAC(IVEG) = 0.                                                  
  300   CONTINUE                                                      
 1000  CONTINUE                                                       
!---------------------------------------------------------------------- 
!     SNOWMELT / REFREEZE CALCULATION                                   
!---------------------------------------------------------------------- 
!                                                                       
!     CALL SNOWM                                                        
!
!=======================================================================
!                                                                       
!     CALCULATION OF SNOWMELT AND MODIFICATION OF TEMPERATURES          
!     N.B. THIS VERSION DEALS WITH REFREEZING OF WATER                  
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      DO 7000 IVEG = 1, 2                                               
!                                                                       
      CCT = CCX                                                         
      TS = TC                                                           
      DTS = DTC                                                         
      FLUX = CHF                                                        
      IF ( IVEG .EQ. 1 ) GO TO 7100                                    
      CCT = CG                                                          
      TS = TGS                                                          
      DTS = DTG                                                         
      FLUX = CCT * DTG / DTT                                            
                                                                        
 7100  CONTINUE                                                        
!                                                                       
      TTA = TS - DTS                                                    
      TTB = TS                                                          
      SNOWW(IVEG) = 0.                                                  
      IF ( TTA .LE. TF ) SNOWW(IVEG) = CAPAC(IVEG)                      
      CAPAC(IVEG) = CAPAC(IVEG) - SNOWW(IVEG)                           
      IF ( TTA .GT. TF .AND. TTB .GT. TF ) GO TO 7200                  
      IF ( TTA .LE. TF .AND. TTB .LE. TF ) GO TO 7200                 
!                                                                       
      DTF = TF - TTA                                                    
      DTIME1 = CCT * DTF /  FLUX                                        
      HF = FLUX*(DTT-DTIME1)                                            
      FCAP = - CAPAC(IVEG)  * SNOMEL                                    
      SPWET = AMIN1( 5. , SNOWW(IVEG) )                                 
      IF ( DTS .GT. 0. ) FCAP =  SPWET * SNOMEL                         
      DTIME2 = FCAP / FLUX                                              
      DTF2 =   FLUX * (DTT-DTIME1-DTIME2)/CCT                           
      TN = TF + DTF2                                                    
      TS = TF - 0.1                                                     
      IF (ABS(HF) .GE.ABS(FCAP) ) TS = TN                               
      CHANGE = HF                                                       
      IF (ABS(CHANGE) .GE.ABS(FCAP) ) CHANGE = FCAP                     
!                                                                       
      CHANGE = CHANGE / SNOMEL                                          
!
      IF (CHANGE.GT.0.0) SMELT=CHANGE+SMELT
!
      SNOWW(IVEG) = SNOWW(IVEG) - CHANGE                                
      CAPAC(IVEG) = CAPAC(IVEG) + CHANGE                                
!                                                                       
      IF ( IVEG .EQ. 1 ) TC = TS                                        
      IF ( IVEG .EQ. 2 ) TGS = TS                                       
      IF ( SNOWW(IVEG) .LT. 0.00001 ) GO TO 7200                       
      ZMELT = 0.                                                        
!     modified to force water into soil. Xue Feb. 1994
      ZMELT = CAPAC(IVEG)                             
!     IF ( TD .GT. TF ) ZMELT = CAPAC(IVEG)                             
      FILTR =  FILTR+ ZMELT                    
      WWW(1) = WWW(1) + ZMELT / ( POROS * ZDEPTH(1) )                   
!     IF ( TD .LE. TF )  ROFF = ROFF + CAPAC(IVEG)   
      CAPAC(IVEG) = 0.                                                  
 7200   CONTINUE                                                      
!                                                                       
      CAPAC(IVEG) = CAPAC(IVEG) + SNOWW(IVEG)                           
!                                                                       
 7000  CONTINUE                                                       
!                                                                       
      FLUXEF = SHF - CCT*DTG/DTT                                        
      TD = TD + FLUXEF / ( CG * 2. * SQRT ( PIE*365. ) ) * DTT          
!
!       *** LOAD PILPS DATA
!
!     if (change .gt. 0) snm(istat)=snm(istat)+(change*1000.)
      change=0.0
!                                                                       
!---------------------------------------------------------------------- 
!     BARE SOIL EVAPORATION LOSS                                        
!---------------------------------------------------------------------- 
!                                                                       
      FILTR = FILTR - EGS / HLAT / 1000.       
      WWW(1) = WWW(1) - EGS / HLAT / 1000. / ( POROS * ZDEPTH(1) )      
!                                                                       
!---------------------------------------------------------------------- 
!   EXTRACTION OF TRANSPIRATION LOSS FROM ROOT ZONE                     
!---------------------------------------------------------------------- 
!                                                                       
      DO 2000 IVEG = 1, 2                                               
!                                                                       
      IF ( IVEG .EQ. 1 ) ABSOIL = ECT / HLAT / 1000.                    
      IF ( IVEG .EQ. 2 ) ABSOIL = EGT / HLAT / 1000.                    
!                                                                       
      if (NROOT.EQ.1) then
      EF(2) = 0.                                                        
      EF(3) = 0.                                                        
      TOTDEP = ZDEPTH(1)                                                
!                                                                       
      DO 3000 IL = 2, 3                                                 
      TOTDEP = TOTDEP + ZDEPTH(IL)                                      
!                                                                       
!     DIV = AMAX1 ( 1., ( PHSOIL(IL) - PHL(IVEG) ) )                    
!                                                                       
      IF ( ROOTD(IVEG) .LT. TOTDEP ) GO TO 400                          
!                                                                       
      EF(IL) = ZDEPTH(IL) / ROOTD(IVEG)                                 
      GO TO 500                                                         
!                                                                       
  400   CONTINUE                                                          
      EF(IL) = ROOTD(IVEG) - TOTDEP + ZDEPTH(IL)                        
      EF(IL) = EF(IL) / ROOTD(IVEG)                                     
      GO TO 600                                                         
!                                                                       
  500   CONTINUE                                                          
 3000  CONTINUE                                                          
!                                                                       
  600   EFT = EF(2) + EF(3)                                               
!
      EFT = MAX(EFT,0.1E-5)
!
      EF(2) = EF(2) / EFT                                               
      EF(3) = EF(3) / EFT                                               
!
      DO 4000 IL = 2, 3                                                 
      WWW(IL) = WWW(IL) - ABSOIL * EF(IL) / ( POROS * ZDEPTH(IL) )      
 4000  CONTINUE                                                          
      else
      ef(1) = rootp(1)
      ef(2) = rootp(2)
      ef(3) = rootp(3)
      DO 4004 IL = 1, 3                                                 
      WWW(IL) = WWW(IL) - ABSOIL * EF(IL) / ( POROS * ZDEPTH(IL) )      
 4004  CONTINUE                                                          
      end if
!                                                                       
 2000  CONTINUE                                                          
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
!     CALCULATION OF INTERFLOW, INFILTRATION EXCESS AND LOSS TO         
!     GROUNDWATER .  ALL LOSSES ARE ASSIGNED TO VARIABLE 'ROFF' .       
!                                                                       
!---------------------------------------------------------------------- 
!                                                                       
      DO 5000 IL = 1, 2                                                 
      IF ( WWW(IL) .GT. 0. ) GO TO 700                                  
      WWW(IL+1) = WWW(IL+1) + WWW(IL) * ZDEPTH(IL)/ZDEPTH(IL+1)         
      WWW(IL) = 0.                                                      
  700   CONTINUE                                                          
 5000  CONTINUE                                                          
!                                                                       
!     IF ( TD .LT. TF ) GO TO 800                                       
!  
!     call run2
!                                                                               
!=======================================================================        
!    calculation of interflow, infiltration excess and loss to                  
!    groundwater .  all losses are assigned to variable 'roff' .                
!----------------------------------------------------------------------         
!                                                                               
      do 8000 i = 1, 3                                                          
!                                                                               
      temw(i)   = amax1( 0.03, www(i) )                                         
      temwp(i)  = temw(i) ** ( -bee )                                           
      temwpp(i) = amin1( 1., temw(i)) ** ( 2.*bee+ 3. )                         
 8000  continue                                                                  
!                                                                               
!-----------------------------------------------------------------------        
!                                                                               
!    calculation of gravitationally driven drainage from w(3) : taken           
!    as an integral of time varying conductivity.addition of liston             
!    baseflow term to original q3g to insure flow in                            
!    dry season. modified liston baseflow constant scaled                       
!    by available water.                                                        
!                                                                               
!     q3g (q3) : equation (62) , SE-86                                          
!                                                                               
!-----------------------------------------------------------------------        
!                                                                               
      pows = 2.*bee+2.                                                          
      q3g = temw(3)**(-pows) + satco/zdepth(3)/poros*slope*pows*dtt             
      q3g = q3g ** ( 1. / pows )                                                
      q3g = - ( 1. / q3g - www(3) ) * poros * zdepth(3) / dtt                   
      q3g = amax1( 0., q3g )                                                    
      q3g = amin1( q3g, www(3)*poros*zdepth(3)/dtt )                            
!                                                                               
      q3g = q3g + 0.002*poros*zdepth(3)*0.5 / 86400. * www(3)                   
!                                                                               
!----------------------------------------------------------------------         
!                                                                               
!    calculation of inter-layer exchanges of water due to gravitation           
!    and hydraulic gradient. the values of w(x) + dw(x) are used to             
!    calculate the potential gradients between layers.                          
!    modified calculation of mean conductivities follows ME-82 ), 
!    reduces recharge flux to top layer.                      
!                                                                               
!      dpdw           : estimated derivative of soil moisture potential         
!                       with respect to soil wetness. assumption of             
!                       gravitational drainage used to estimate likely          
!                       minimum wetness over the time step.                     
!                                                                               
!      qqq  (q     )  : equation (61) , SE-86                                   
!             i,i+1                                                             
!            -                                                                  
!      avk  (k     )  : equation (4.14) , ME-82                                 
!             i,i+1                                                             
!                                                                               
!----------------------------------------------------------------------         
!                                                                               
      wmax = amax1( www(1), www(2), www(3), 0.05 )                              
      wmax = amin1( wmax, 1. )                                                  
      pmax = wmax**(-bee)                                                       
      wmin = (pmax-2./( phsat*(zdepth(1)+2.*zdepth(2)+zdepth(3))))   &
              **(-1./bee)                                                       
      wmin = amin1( www(1), www(2), www(3), wmin )                              
      wmin = amax1( wmin, 0.02 )                                                
      pmin = wmin**(-bee)                                                       
      dpdw = phsat*( pmax-pmin )/( wmax-wmin )                                  
!                                                                               
      do 8200 i = 1, 2                                                          
!                                                                               
      rsame = 0.                                                                
      avk  = temwp(i)*temwpp(i) - temwp(i+1)*temwpp(i+1)                        
      div  = temwp(i+1) - temwp(i)                                              
      if ( abs(div) .lt. 1.e-6 ) rsame = 1.                                     
      avk = satco*avk / ( ( 1. + 3./bee ) * div + rsame )                       
      avkmin = satco * amin1( temwpp(i), temwpp(i+1) )                          
      avkmax = satco * amax1( temwpp(i), temwpp(i+1) )*1.01                     
      avk = amax1( avk, avkmin )                                                
      avk = amin1( avk, avkmax )                                                
!                                                                               
!-----------------------------------------------------------------------        
!     conductivities and base flow reduced when temperature drops below         
!     freezing.                                                                 
!-----------------------------------------------------------------------        
!                                                                               
      tsnow = amin1 ( tf-0.01, tgs ) 
      areas = amin1 (0.999,13.2*snoww(2))
      tgg = tsnow*areas + tgs*(1.-areas)
      ts    = tgg*(2-i) + td*(i-1)                                              
      props = ( ts-(tf-10.) ) / 10.                                             
!     props = 1.+5*(ts-tf)                                             
      props = amax1( 0.05, amin1( 1.0, props ) )                                
      avk  = avk * props                                                        
      q3g  = q3g * props                                                        
!                                                                               
!-----------------------------------------------------------------------        
!     backward implicit calculation of flows between soil layers.               
!-----------------------------------------------------------------------        
!                                                                               
      dpdwdz = dpdw * 2./( zdepth(i) + zdepth(i+1) )                            
      aaa(i) = 1. + avk*dpdwdz*( 1./zdepth(i)+1./zdepth(i+1) )  &       
                  *dtt/poros                                                    
      bbb(i) =-avk *   dpdwdz * 1./zdepth(2)*dtt/poros                          
      ccc(i) = avk * ( dpdwdz * ( www(i)-www(i+1) ) + 1. +      &      
                 (i-1)*dpdwdz*q3g*1./zdepth(3)*dtt/poros )                      
 8200  continue                                                                  
!                                                                               
      denom  = ( aaa(1)*aaa(2) - bbb(1)*bbb(2) )                                
      rdenom = 0.                                                               
      if ( abs(denom) .lt. 1.e-6 ) rdenom = 1.                                  
      rdenom = ( 1.-rdenom)/( denom + rdenom )                                  
      qqq(1)   = ( aaa(2)*ccc(1) - bbb(1)*ccc(2) ) * rdenom                     
      qqq(2)   = ( aaa(1)*ccc(2) - bbb(2)*ccc(1) ) * rdenom                     
!                                                                               
!-----------------------------------------------------------------------        
!     update wetness of each soil moisture layer due to layer interflow         
!        and base flow.                                                         
!-----------------------------------------------------------------------        
!                                                                               
      www(3) = www(3) - q3g*dtt/(poros*zdepth(3))                               
      roff = roff + q3g * dtt                                                   
!                                                                               
      do 8300 i = 1, 2                                                          
!                                                                               
      qmax   =  www(i)   * (poros*zdepth(i)  /dtt)                              
      qmin   = -www(i+1) * (poros*zdepth(i+1)/dtt)                              
      qqq(i) = amin1( qqq(i),qmax)                                              
      qqq(i) = amax1( qqq(i),qmin)                                              
      www(i)   =   www(i)   - qqq(i)/(poros*zdepth(i)  /dtt)                    
      www(i+1) =   www(i+1) + qqq(i)/(poros*zdepth(i+1)/dtt)                    
 8300  continue                                                                  
!
!       *** LOAD water flow & root-zone drainage PILPS DATA
      soildif=soildif+qqq(1)*dtt*1000.
      soildra=soildra+q3g*dtt*1000.
!
      do 8400 i = 1, 3                                                          
      excess = amax1(0.,(www(i) - 1.))                                          
      www(i) = www(i) - excess                                                  
      roff   = roff   + excess * poros*zdepth(i)                                
!
!       *** LOAD IN as root-drainage for PILPS
      if (i.lt.2) then
        RNOFFS= RNOFFS+ 1000.*excess*POROS*ZDEPTH(I)
      else
        RNOFFB= RNOFFB+ 1000.*excess*POROS*ZDEPTH(I)
      endif
 8400  continue                                                                  
!                                                                               
!-----------------------------------------------------------------------        
!     prevent negative values of www(i)                                         
!-----------------------------------------------------------------------        
!                                                                               
      do 8402 i = 1,2                                                           
      deficit   = amax1 (0.,(1.e-12 - www(i)))                                  
      if (i.eq.1) soildif=soildif-deficit*                 &
                  zdepth(1)*poros
      www (i)   = www(i) + deficit                                              
      www (i+1) = www(i+1) - deficit * zdepth(i) / zdepth (i+1)                 
 8402 continue                                                                  
      www(3)    = amax1 (www(3),1.e-12)                                         
! 
  800   CONTINUE                                                          
!                                                                       
      IF (WWW(1) .GT.1.) then 
          WWW(2) = WWW(2) + (WWW(1)-1.) * ZDEPTH(1)/ZDEPTH(2)
          soildif=soildif+(www(1)-1.)*ZDEPTH(1) *poros*1000.
          WWW(1) = 1.                                    
      end if
      If (WWW(2) .GT.1.) WWW(3) = WWW(3) + (WWW(2)-1.) *   &
                                  ZDEPTH(2) / ZDEPTH(3) 
!
!       *** LOAD IN AS PILP ROOT DRAINAGE
      IF (WWW(2) .GT.1.) WWW(2) = 1.                                    
      if (WWW(3) .GT.1.) then 
          ROFF   = ROFF + (WWW(3)-1.)*POROS*ZDEPTH(3)    
          RNOFFB=RNOFFB+((WWW(3)-1.)*ZDEPTH(3)*POROS*1000.)
          WWW(3) = 1.                                    
      end if
!                                                                       
      END SUBROUTINE
!=======================================================================
!                                                                       
      SUBROUTINE VEGOUT1(XPOROS,XDEPTH,ITYPE)
!zzq 2009Feb23   SUBROUTINE VEGOUT1(XPOROS,XDEPTH,MONTH,ITYPE)
!     &     XSOREF,
!     &     XBEE, XPHSAT, XPOROS, XSATCO,XSLOPE,
!     &     XDEPTH,MONTH,ITYPE)
!                                                     12 AUGUSTY 2000 
!=======================================================================
!                                                                       
!     ASSIGN VEGETATION PHYSIOLOGY                                        
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
!   ZDEPTH        : DEPTH OF 3 SOIL MOISTURE LAYERS                     
!   Z0            : ROUGHNESS LENGTH                                    
!   XDD           : ZERO PLANE DISPLACEMENT                             
!   ZLT(IV)       : LEAF AREA INDEX                                     
!   GREEN(IV)     : GREEN LEAF FRACTION                                 
!   VCOVER(IV)    : VEGETATION COVER FRACTION                           
!                                                                       
!     VARIABLES ( SPECIFIC TO SIB 1-D VERSION ONLY ) FROM COMSIB        
!                                                                       
!      ZWIND  : REFERENCE HEIGHT FOR WIND MEASUREMENT                   
!      ZMET   : REFERENCE HEIGHT FOR TEMPERATURE, HUMIDITY MEASUREMENT  
!        THE ABOVE ARE GENERATED FROM SIBX + MOMOPT OUTPUT              
!                                                                       
!----------------------------------------------------------------------
!the following not used. commented by zzq
!zzq      DIMENSION XTRAN(2,3,2), XREF(2,3,2),XSOREF(3),
!zzq     &          XGREEN(2), XVCOVER(2), XZLT(2), XCHIL(2),
!zzq     &          XRSTPAR(2,3), XTOPT(2), XTL(2), XTU(2), XDEFAC(2),
!zzq     &          XPH1(2), XPH2(2), XROOTD(2), XDEPTH(3)
         dimension  XDEPTH(3)
!-----------------------------------------------------------------------
!                                                                       
!       DO IW=1,3
!       XTRAN(1,IW,1)=TRAN0(ITYPE,1,IW,1)
!       XTRAN(1,IW,2)=TRAN0(ITYPE,1,IW,2)
!       XTRAN(2,IW,1)=TRAN0(ITYPE,2,IW,1)
!       XTRAN(2,IW,2)=TRAN0(ITYPE,2,IW,2)
!       XREF (1,IW,1)= REF0(ITYPE,1,IW,1)
!       XREF (1,IW,2)= REF0(ITYPE,1,IW,2)
!       XREF (2,IW,1)= REF0(ITYPE,2,IW,1)
!       XREF (2,IW,2)= REF0(ITYPE,2,IW,2)
!       XRSTPAR(1,IW)=RSTPAR0(ITYPE,1,IW)
!       XRSTPAR(2,IW)=RSTPAR0(ITYPE,2,IW)
!       XSOREF  (IW) =SOREF0(ITYPE,IW)
!       END DO
!       DO IV=1,2
!       XCHIL(IV)=CHIL0(ITYPE,IV) 
!       XTOPT(IV)=TOPT0(ITYPE,IV)
!       XTL(IV)=TL0(ITYPE,IV)
!       XTU(IV)=TU0(ITYPE,IV)
!       XDEFAC(IV)=DEFAC0(ITYPE,IV)
!       XPH1(IV)=PH10(ITYPE,IV)
!       XPH2(IV)=PH20(ITYPE,IV)
!       XROOTD(IV)=ROOTD0(ITYPE,IV)
!       XZLT(IV)=ZLT0(ITYPE,MONTH,IV)
!       XGREEN(IV)=GREEN0(ITYPE,MONTH,IV)
!       XVCOVER(IV)=VCOVER0(ITYPE,MONTH,IV)
!       END DO
       DO IDEP=1,3
       XDEPTH(IDEP)=DEPTH0(ITYPE,IDEP)
       END DO
!                   
!        xvm0=VM00(ITYPE)
!        Xfc3=fc30(ITYPE)
!       write(*,100)xvm0,xfc3
   36     format(2f8.2)
!       XBEE=BEE0(ITYPE)
!       XPHSAT=PHSAT0(ITYPE)
!       XSATCO=SATCO0(ITYPE)
        XPOROS=POROS0(ITYPE)
!       XSLOPE=SLOPE0(ITYPE)
!       XZ2=Z20(ITYPE,MONTH)
!       XZ1=Z10(ITYPE,MONTH)
!       XZ0= ZZ0(ITYPE,MONTH)
!       XDD= D0(ITYPE,MONTH)
!       XRBC=RBC0 (ITYPE,MONTH)
!       XRDC=RDC0 (ITYPE,MONTH)
!                        
!      rootp(idep),idep=1,3
!
      END SUBROUTINE
!
      SUBROUTINE VEGOUT2(XTRAN,XREF,XGREEN,XVCOVER,XCHIL,      &
           XRSTPAR,XVM0,XTOPT,XTL,XTU,Xfc3,XDEFAC,XPH1,XPH2,   &
           XZLT,XZ0,XDD,XZ2,XZ1,XRDC,XRBC,XROOTD,              &
           XSOREF,XBEE, XPHSAT,  XSATCO,XSLOPE,                &
           MONTH,NXTYPE)
!                                                     12 AUGUSTY 2000 
!=======================================================================
!                                                                       
!     ASSIGN VEGETATION PHYSIOLOGY                                        
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
!   ZDEPTH        : DEPTH OF 3 SOIL MOISTURE LAYERS                     
!   Z0            : ROUGHNESS LENGTH                                    
!   XDD           : ZERO PLANE DISPLACEMENT                             
!   ZLT(IV)       : LEAF AREA INDEX                                     
!   GREEN(IV)     : GREEN LEAF FRACTION                                 
!   VCOVER(IV)    : VEGETATION COVER FRACTION                           
!                                                                       
!     VARIABLES ( SPECIFIC TO SIB 1-D VERSION ONLY ) FROM COMSIB        
!                                                                       
!      ZWIND  : REFERENCE HEIGHT FOR WIND MEASUREMENT                   
!      ZMET   : REFERENCE HEIGHT FOR TEMPERATURE, HUMIDITY MEASUREMENT  
!        THE ABOVE ARE GENERATED FROM SIBX + MOMOPT OUTPUT              
!                                                                       
!----------------------------------------------------------------------
!
      DIMENSION XTRAN(2,3,2), XREF(2,3,2),XSOREF(3),                &
                XGREEN(2), XVCOVER(2), XZLT(2), XCHIL(2),           &
                XRSTPAR(2,3), XTOPT(2), XTL(2), XTU(2), XDEFAC(2),  &
                XPH1(2), XPH2(2), XROOTD(2)
!-----------------------------------------------------------------------
!                                                                       
       DO IW=1,3
       XTRAN(1,IW,1)=TRAN0(NXTYPE,1,IW,1)
       XTRAN(1,IW,2)=TRAN0(NXTYPE,1,IW,2)
       XTRAN(2,IW,1)=TRAN0(NXTYPE,2,IW,1)
       XTRAN(2,IW,2)=TRAN0(NXTYPE,2,IW,2)
       XREF (1,IW,1)= REF0(NXTYPE,1,IW,1)
       XREF (1,IW,2)= REF0(NXTYPE,1,IW,2)
       XREF (2,IW,1)= REF0(NXTYPE,2,IW,1)
       XREF (2,IW,2)= REF0(NXTYPE,2,IW,2)
       XRSTPAR(1,IW)=RSTPAR0(NXTYPE,1,IW)
       XRSTPAR(2,IW)=RSTPAR0(NXTYPE,2,IW)
       XSOREF  (IW) =SOREF0(NXTYPE,IW)
       END DO
       DO IV=1,2
       XCHIL(IV)=CHIL0(NXTYPE,IV) 
       XTOPT(IV)=TOPT0(NXTYPE,IV)
       XTL(IV)=TL0(NXTYPE,IV)
       XTU(IV)=TU0(NXTYPE,IV)
       XDEFAC(IV)=DEFAC0(NXTYPE,IV)
       XPH1(IV)=PH10(NXTYPE,IV)
       XPH2(IV)=PH20(NXTYPE,IV)
       XROOTD(IV)=ROOTD0(NXTYPE,IV)
       XZLT(IV)=ZLT0(NXTYPE,MONTH,IV)
       XGREEN(IV)=GREEN0(NXTYPE,MONTH,IV)
       XVCOVER(IV)=VCOVER0(NXTYPE,MONTH,IV)
       END DO
!       DO IDEP=1,3
!       XDEPTH(IDEP)=DEPTH0(ITYPE,IDEP)
!       END DO
!                   
        xvm0=VM00(NXTYPE)
        Xfc3=fc30(NXTYPE)
!       write(*,100)xvm0,xfc3
!36     format(2f8.2)
       XBEE=BEE0(NXTYPE)
       XPHSAT=PHSAT0(NXTYPE)
       XSATCO=SATCO0(NXTYPE)
!       XPOROS=POROS0(ITYPE)
       XSLOPE=SLOPE0(NXTYPE)
       XZ2=Z20(NXTYPE,MONTH)
       XZ1=Z10(NXTYPE,MONTH)
       XZ0= ZZ0(NXTYPE,MONTH)
       XDD= D0(NXTYPE,MONTH)
       XRBC=RBC0 (NXTYPE,MONTH)
       XRDC=RDC0 (NXTYPE,MONTH)
!                        
!      rootp(idep),idep=1,3
!
      END SUBROUTINE
END MODULE
