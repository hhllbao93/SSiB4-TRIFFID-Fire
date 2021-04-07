MODULE ssib4_module_ssibco2
    use ssib4_module_comsconst
contains
!=======================================================================
       SUBROUTINE GMUPI(SUNANG,FC3,CHIL,ZLT,VCOVER,EXTK,PD,GREEN,   &
           GMU1,GMU2,BPIDR3,BPIDR4,BPIDF3,BPIDF4)
!                                                         AUGUST 2003
!=======================================================================
!                                                                       
!     CALCULATION OF GMU and Scaling for the sunlid and sunshed potion
!                                                                       
!-----------------------------------------------------------------------
       DIMENSION CHIL(2),ZLT(2),VCOVER(2),EXTK(2,3,2),PD(2),GREEN(2)
!
       IVEG = 1
!
!-----------------------------------------------------------------------
!     ROSS INCLINATION FUNCTION
!-----------------------------------------------------------------------
!
      AA = 0.5 - 0.633 * CHIL(IVEG)- 0.33 * CHIL(IVEG)* CHIL(IVEG)
      BB = 0.877 * ( 1. - 2. * AA )
!-----------------------------------------------------------------------
!     ESTIMATE OF GMU FOR DIFFERENT COMPONENTS
!-----------------------------------------------------------------------
!
      GMU1 = ( AA + BB * SUNANG ) / SUNANG
      GMU2 = BB / 3. + AA * 1.5 + BB / 4. * PIE
!
      AT = ZLT(IVEG) / VCOVER(IVEG)
!     Fsl = 2*(1-exp(-0.5*AT/SUNANG))*SUNANG ! sunlit leaf area index
      POW = AMIN1 (10., AT * EXTK(IVEG,1,1) )
!      Fsl = 4*(1-exp(-POW/SUNANG))*SUNANG ! sunlit leaf area index
       Fsl = (1-exp(-POW/SUNANG))*SUNANG/EXTK(IVEG,1,1)!sunlit leaf area index
        if (Fsl.lt.0.0) Fsl=0.0
!zzq,jan5,2009 if (fsl.gt.at) fsl=0.95*at
        if (fsl.gt.0.95*at) fsl=0.95*at
        Fslf= AMAX1 (Fsl/AT, 0.001)
!
      DO 1000 ICO2 = 1, 2
      IF (ICO2.EQ.1) THEN
      ATO = ZLT(IVEG) * fc3 / VCOVER(IVEG)
      ATO=AMAX1 (ATO,0.0001)
      ELSE
      ATO = ZLT(IVEG) * (1-fc3) / VCOVER(IVEG)
      ATO=AMAX1 (ATO,0.0001)
      END IF
!
      AT = ATO * Fslf
!
      POWER1 = AMIN1( 50., AT * EXTK(IVEG,1,1) )
      POWER2 = AMIN1( 50., AT * EXTK(IVEG,1,2) )
!
!-----------------------------------------------------------------------
!     COMBINED ESTIMATE OF K-PAR USING WEIGHTS FOR DIFFERENT COMPONENTS
!-----------------------------------------------------------------------
!
      XX = PD (IVEG) + (1. - PD(IVEG)) * Fslf
      ZAT = ALOG( ( EXP(-POWER1) + 1. )/2. ) * (PD(IVEG) / XX)   &
            / ( POWER1/AT )
      ZAT = ZAT + ALOG( ( EXP(-POWER2) + 1. )/2. )               &
       * ((( 1. - PD(IVEG) ) * Fslf ) / XX) / ( POWER2/AT )

!
      POW1 = AMIN1( 50., (POWER1*ZAT/AT) )
      POW2 = AMIN1( 50., (POWER2*ZAT/AT) )
!
      ZK = 1. / ZAT * ALOG( ( PD(IVEG) / XX) * EXP ( POW1 )      &
            + (( 1. - PD(IVEG) ) *Fslf / XX) * EXP ( POW2 ) )
!
!kk
      ZK   = AMAX1(1.e-5,ZK)
!kk
      POW  = AMIN1( 50., ZK*AT )
      IF (ICO2.EQ.1) THEN
      BPIDR3 = VCOVER(IVEG)*GREEN(IVEG)*(1.- EXP ( -POW  ))/ZK
      ELSE
      BPIDR4 = VCOVER(IVEG)*GREEN(IVEG)*(1.- EXP ( -POW  ))/ZK
      END IF
!
      AT = ATO * (1-Fslf)
!
      POW22 = AMIN1( 50., AT * EXTK(IVEG,1,2) )
!
      ZK2 =  AMIN1( 50. / AT, EXTK(IVEG,1,2) )
      IF (ICO2.EQ.1) THEN
      BPIDF3 = VCOVER(IVEG)*GREEN(IVEG)*(1.- EXP ( -POW22 ))/ZK2
      ELSE
      BPIDF4 = VCOVER(IVEG)*GREEN(IVEG)*(1.- EXP ( -POW22 ))/ZK2
      END IF
!
 1000 CONTINUE
!
      END SUBROUTINE
!==============================================================================
!
      SUBROUTINE ANGS2(PAR,TAK,TCK,EA,RB,CA,ZLT,VCOVER,P,                 &
                 VM0,FC3,FW,FD,SE,GMU1,GMU2,BPIDR3,BPIDR4,                &
                 BPIDF3,BPIDF4,EXTK,TOPTM,THIGH,TLOW,AN,RST,ES,GB,GS,     &
                 WCL,WEL,WSL,WCS,WES,WSS,Crd,ft1,ft2,hs)
!==============================================================================
!------------------------------------------------------------------------------
!
!     This is to run the An-rst coupled model for a given environment.
!
!     INPUT variables:
!
!     PAR         ! PAR                                         [W m-2]
!     TaK         ! canopy air temperature                      [K]
!     TcK         ! canopy leaf temperature                     [K]
!     ea          ! vapor pressure at measurement height        [Pa]
!     rb          ! canopy to canopy air space resistance       [s m-1]
!     Ca          ! CO2 concentration                           [Pa]
!     LAI         ! leaf area index of the canopy               [m2 m-2]
!     p           ! air pressure                                [Pa]
!     fC3         ! fraction of C3 plant coverage
!                   (fC4 = 1 - fC3)
!
!     OUTPUT variables:
!
!     An          ! net photosynthetic rate                     [mol m-2 s-1]
!     rst         ! stomatal resistance to water vapor transfer [s m-1]
!     es          ! leaf surface vapor pressure                 [Pa]
!
!
!     Xiwu Zhan, University f Maryland Baltimore County, 10-2-2002
!     Yongkang Xue   8-15-03
!------------------------------------------------------------------------------
!kk   common /AnParam/fw
!kk   common /RSLT/wc,we,ws,al,gs,wc1,we1,ws1,al1
!kk   common /Tfac/fT1,fT2,Respd
!     common /se/se
      dimension Ac(4),gc(4),esnc(4),F34(4),a(4),zlt(2),vcover(2)
!khs
      dimension xwcl(4),xwsl(4),xwel(4)
      dimension xwcs(4),xwss(4),xwes(4)
!song
      dimension crdo(4)
      real crdsl, crdsh
!khs
      integer PT
      real LAI

      if (zlt(1).lt.0.001) then
          lai=0.001
      else
          lai=zlt(1)/vcover(1)
      endif
      if (LAI.le.0.001) then
!khs  An=0
!khs  rst=2000.
!khs  es=ea
!khs  print *,'LAI is too small'
!khs  stop 999
      lai = 0.001
      end if

      fC4=1.-fC3                ! C4 fraction
      F34(3) = LAI * fC3
      F34(4) = LAI * fC4
!khs
      gb = 0.5

!     fd = 0.23                 ! proportion of diffusive radiation
!     PARdf = PAR*fd            ! diffusive PAR
!     PARdr = PAR*(1-fd)        ! direct PAR beam
      PARdr = PAR*fd            ! direct PAR beam
      PARdf = PAR*(1-fd)        ! diffusive PAR
!     PARsh = PARdf*exp(-0.5*LAI**0.7)+0.07*PARdr*(1.1-
!    +     0.1*LAI)*exp(-se)    ! PAR on shaded leaves
!     PARsl = PARdr*0.5/se + PARsh ! PAR on sunlit leaves
!     Fsl = 2*(1-exp(-0.5*LAI/se))*se ! sunlit leaf area index
      POW = AMIN1 (10., LAI * EXTK)
!      Fsl = 4*(1-exp(-POW/se))*se ! sunlit leaf area index
       Fsl = (1-exp(-POW/se))*se/EXTK ! sunlit leaf area index
!khs    if (Fsl.lt.0.) Fs1=0.0
!khs    if (Fsl.gt.LAI) Fs1=0.95*LAI
        if (Fsl.lt.0.) Fsl=0.0
        if (Fsl.gt.LAI) Fsl=0.95*LAI
         Fslf= AMAX1(Fsl/LAI, 0.001)
       if (Fsl.lt.0..or.Fsl.gt.LAI) then
          write(6,*) 'Fsl is worng'
          stop 999
       end if

        if (se.lt.0.00001) then
         PARshf = 0.001
         PARsld = 0.001
         PARslf = 0.001
         Fsl   = 0
         Fslf  = 0
        end if

        do PT=3,4                 ! Plant Type (C3-3, C4-4)
         if (fC3.eq.0.and.PT.eq.3) then
            Ac(3)=0
            goto 90
         end if

         if (fC3.eq.1.and.PT.eq.4) then
            Ac(4)=0
            goto 90
         end if

         if (PT.eq.3) then
            bpi=bpidr3
            PARsld=fc3*PARdr
            PARslf=fc3*PARdf*Fslf
         else
            bpi=bpidr4
            PARsld=fc4*PARdr
            PARslf=fc4*PARdf*Fslf
         end if

         do LT=1,3              ! Limitation type (Wc-1, We-2, Ws-3)
            es=ea               ! es is initiallized to ea
         call angsa(vm0,toptm,thigh,tlow,PARsld,PARslf,Tak,ea,rb,Ca,   &
           Tck,es,p,F34(PT)*Fslf,PT,LT,An,gs,fw,gmu1,gmu2,bpi,         &
           Crd,ft1,ft2)
            a(LT) = An*F34(PT)*Fslf/LAI
         end do
!song 2.26.11
         crdsl= crd*F34(PT)*Fslf/LAI
         if (pt.eq.3) then
         atheta=0.98
         btheta=0.95
         else
         atheta=0.85
         btheta=0.85
         endif
         sqrtin=(a(1)+a(2))**2-4.*atheta*a(1)*a(2)
         a(4)=((a(1)+a(2))-sqrt(sqrtin))/(2.*atheta )

         if (sqrtin.lt.0.) then
         write(7,*) 'Problem in sunlid smoothing, a1,a2,atheta'
         write(7,*) a(1),a(2),atheta
         stop 777
         end if
         sqrtin=(a(3)+a(4))**2-4.*btheta*a(3)*a(4)
         if (sqrtin.lt.0.) then
         write(7,*) 'Problem in sunlid smoothing, a3,a4,atheta'
         write(7,*) a(3),a(4),btheta
         stop 777
         end if
         Ansl=((a(3)+a(4))-sqrt(sqrtin))/(2.*btheta )

         xwcl(pt) = a(1)
         xwel(pt) = a(2)
         xwsl(pt) = a(3)
!khs     al = Ansl
         if (PT.eq.3) then
            bpi=bpidf3
            PARshf=fc3*PARdf*(1.-Fslf)
         else
            bpi=bpidf4
            PARshf=fc4*PARdf*(1.-Fslf)
         end if

         do LT=1,3
            es=ea               ! es is initiallized to ea
       call angsa(vm0,toptm,thigh,tlow,0.001,PARshf,Tak,ea,rb,Ca,   &
         Tck,es,p,F34(PT)*(1-Fslf),PT,LT,An,gs,fw,gmu1,gmu2,bpi,    &
         Crd,ft1,ft2)
            a(LT) = An*F34(PT)*(1-Fslf)/LAI
         end do
!song 2.26.11
         crdsh= crd*F34(PT)*(1-Fslf)/LAI
         sqrtin=(a(1)+a(2))**2-4.*atheta*a(1)*a(2)
         if (sqrtin.lt.0.) then
         write(7,*) 'Problem in sunshed smoothing, a1,a2,atheta'
         write(7,*) a(1),a(2),atheta
         stop 777
         end if
         a(4)=((a(1)+a(2))-sqrt(sqrtin))/(2.*atheta )
         sqrtin=(a(3)+a(4))**2-4.*btheta*a(3)*a(4)
         if (sqrtin.lt.0.) then
         write(7,*) 'Problem in sunshed smoothing, a3,a4,atheta'
         write(7,*) a(3),a(4),btheta
         stop 777
         end if
         Ansh=((a(3)+a(4))-sqrt(sqrtin))/(2.*btheta )
         xwcs(pt) = a(1)
         xwes(pt) = a(2)
         xwss(pt) = a(3)
!khs     al1 = Ansh

         Ac(PT) = Ansl + Ansh
!song
         crdo(PT)=crdsh+crdsl
         ei=e(Tck)              ! satuation vapor pressure [Pa]
         hs=es/ei               ! relative humidity of leaf bndry lyr
         gb=(.5/rb)*44.6*273.16/Tak*p/1.013e5
          
         if (abs(gb).gt.100.) then
         print *,'angs2: rb tak p gb = ',rb,tak,p,gb
         endif
                                ! rb[s m-1] to gb [mol m-2 s-1]
         if (PT.eq.3) then
            m=9
            b=0.01
         end if

         if (PT.eq.4) then
            m=4
            b=0.04
         end if

   90  Csc=Ca-1.4*p*Ansl/gb
         gcsl=m*hs*p*Ansl/Csc+b*F34(PT)
         Csc=Ca-1.4*p*Ansh/gb
         gcsh=m*hs*p*Ansh/Csc+b*F34(PT)

         if (gcsl.le.0) gcsl=0.000001
         if (gcsh.le.0) gcsh=0.000001

!song 12.24.2010         gc(PT) = 1./(1./gcsl*Fslf + 1./gcsh*(1-Fslf))
         gc(PT) = gcsl*Fslf + gcsh*(1-Fslf)

         if (gc(PT).lt.0) then
            gc(PT) = 0.00001
         end if
         esnc(PT)=(gb*ea+gc(PT)*ei)/(gb+gc(PT))
      end do
                                ! composite for the whole canopy

!
!zzq080728  Notice that, 
!zzq080728  if fc3 or fc4 is not assigned, a NaN will take place. 
!
!zzq080728      An=fc3*Ac(3)+fc4*Ac(4)
!zzq080728      gs=fc3*gc(3)+fc4*gc(4)
!zzq080728      es=fc3*esnc(3)+fc4*esnc(4)
!khs..
!zzq080728      gs  = amax1(gs,0.01)
!zzq080728      wcl = fc3*xwcl(3)+fc4*xwcl(4)
!zzq080728      wsl = fc3*xwsl(3)+fc4*xwsl(4)
!zzq080728      wel = fc3*xwel(3)+fc4*xwel(4)
!zzq080728      wcs = fc3*xwcs(3)+fc4*xwcs(4)
!zzq080728      wss = fc3*xwss(3)+fc4*xwss(4)
!zzq080728      wes = fc3*xwes(3)+fc4*xwes(4)

        if(fc3.eq.1) then
          An=Ac(3)
          gs=gc(3)
          es=esnc(3)
          wcl = xwcl(3)
          wsl = xwsl(3)
          wel = xwel(3)
          wcs = xwcs(3)
          wss = xwss(3)
          wes = xwes(3)         
!song
          CRD=crdo(3)
        endif
!
        if(fc4.eq.1) then
          An=Ac(4)
          gs=gc(4)
          es=esnc(4)
          wcl = xwcl(4)
          wsl = xwsl(4)
          wel = xwel(4)
          wcs = xwcs(4)
          wss = xwss(4)
          wes = xwes(4)       
!song
          CRD=crdo(4)      
        endif
!
       gs  = amax1(gs,0.01)
!
!khs  print *,'angs2: co2cs p an gs = ',co2cs,p,an,gs,co2ci

      gsms=gs*0.0224*Tak/273.16*1.013e5/p
!fhs
      gsmst=gsms
!fhs
      if (gsms.lt.0) gsms=1.e-4
      rst=1/gsms

      if (rst.gt.10000.) rst=10000.
      if (rst.eq.10000. .and. par.ge.50) then
!zzq       write(67,55) PAR,gsmst,gs,Tak,se,an
      endif
!zzq   55   format(f7.2,2x,e9.2,2x,e8.2,2x,f6.2,2x,f6.4,2x,e10.3)
!fhs
  999  return
      end subroutine
!=======================================================================

      subroutine angsa(vm0,toptm,thigh,tlow,PARd,PARf,Ta,ea,rb,Ca,Tc,   &
          es,p,LAI,PT,LT,An,gs,fW,Gmu1,Gmu2,BPI,Rd,ft1,ft2)

!=======================================================================
!
!     This An-gs coupled model combines Farquhar's biochemical photosynthesis
!     model with a stomatal conductance model in order to estimate the
!     exchanges of CO2, water vapor and sensible heat between the leaf and
!     and its environment. 
!
!     The atmospheric conditions are 
!
!     PAR          [W m-2]            photosynthetically active radiation
!     Ta           [degree K]         canopy air space temperature
!     ea           [Pa]               c.a.s. vapor pressure
!     es           [Pa]               vapor pressure of leaf bndry lyr
!     rb           [s m-1]            aerodyn. resist. to LE of leaf bndry lyr
!     Ca           [Pa]               c.a.s. CO2 concentraton
!     p            [Pa]               air pressure
!     Tc           [degree K]         canopy leaf temperature
!     LAI          [m2 m-2]           leaf area index of the canopy
!     
!     Farquhar's photosynthesis model (Sellers et al., 1996, J. Climate):
!
! (1) An = (A - Rd)*BPI                                 [mol m-2 s-1]
!          A = min{Wc, We, Ws}
!              Wc=Vm0*(Ci-Gama)/(Ci+Kc*(1+O2/Ko))      for C3 
!                =Vm0                                  for C4 
!              We=PAR*epsn3*Gmu*(1-omg)*(Ci-Gama)/(Ci+2*Gama) for C3 
!                =PAR*epsn4*Gmu*(1-omg)                       for C4 
!              Ws=Vm0*0.5*BPI   for C3 
!                =Vm0*2e4*Ci    for C4 
!          Rd = Vm0*fd
!
!     CO2 transfer equations:
!
! (2) An = (Cs-Ci)*(gs/1.6)
! (3) An = (Ca-Cs)*(gb/1.4)
!
!     stomatal conductance model:
! 
! (4) gs = m*hs*p*An/Cs + b*LAI
!
!     Combining these equations, the following will be solved analytically:
!
!     An         [mol m-2 s-1]        leaf photosynthetic rate
!     gs         [mol m-2 s-1]        stomatal conductance for LE
!     Cs         [Pa]                 leaf surface CO2 concentration
!     Ci         [Pa]                 leaf inside CO2 concentration
!
!         Xiwu Zhan, University of Maryland, 10-20-1998
!         Yongkang Xue, 8-13-2003
!------------------------------------------------------------------------------

      integer PT
!zzq      real m,Kc,Ko,kaba,LAI,A(4,3,5),x(3)
      real m,Kc,Ko,LAI,x(3)

      PARdmol=4.6e-6*PARd         ! [W m-2] to [mol m-2 s-]
      PARfmol=4.6e-6*PARf         ! [W m-2] to [mol m-2 s-]

      gb = (.5/rb)*44.6*273.16/Ta*p/1.013e5 ! rb[s m-1] to gb [mol m-2 s-1]
!khs  print *,'rb ta p gb = ',rb,ta,p,gb
      ei = e(Tc)                ! satuation vapor pressure [Pa]
      hs = es/ei                ! relative humidity of leaf bndry lyr

      O2   = 20900.             ! [Pa]
!     kaba = 0.45               ! canopy atenuation extinction coefficient
!     BPI  = (1-exp(-kaba*LAI))/kaba ! leaf to canopy scaling factor
!     Gmu  = 1.                 ! [G(u)/u]

      Qt = 0.1*(Tc-Toptm)  
      
      Kc = 30.*2.1**Qt   
      Ko = 30000.*1.2**Qt
      S  = 2600.*0.57**Qt
      Gamma = 0.5*O2/S

      if (PT.eq.3) then
         m=9
         b=0.01*fW

!zzq24Feb2014         fT1=2.**Qt/(1+exp(.3*(Tc-Thigh)))
!zzq24Feb2014         fT2=2.**Qt/(1.+exp(1.3*(Tc-Thigh)))
!
         exph1 = exp(.3*(Tc-Thigh)/2.0)
         exph2 = exp(1.3*(Tc-Thigh)/2.0)
         fT1=2.**Qt/(1+exph1*exph1)
         fT2=2.**Qt/(1.+exph2*exph2)

         Vm=Vm0*fT1*fW
         Rd=0.015*Vm0*fT2*fW*BPI

         if (LT.eq.1) then
            A1=Vm*BPI
            A2=Gamma
            A3=1
            A4=Kc*(1+O2/Ko)
            A5=-Rd
         end if

         if (LT.eq.2) then
            A1=(PARdmol*Gmu1+PARfmol*Gmu2)*(1-0.15)*0.08*BPI
            A2=Gamma
            A3=1
            A4=2*Gamma
            A5=-Rd
         end if

         if (LT.eq.3) then
            A1=0.5*Vm*(PARdmol*Gmu1+PARfmol*Gmu2)/    &
              (PARdmol+PARfmol) * BPI
            A2=0
            A3=1
            A4=0
            A5=-Rd
         end if
      end if

      if (PT.eq.4) then
         m=4
         b=0.04*fW

!rr - by Huiping, June 15, 2007
!rr      fT1=2.**Qt*(1+exp(.2*(Tlow-Tc)))/(1+exp(.3*(Tc-Thigh)))
!zzq24Feb2014         fT1=2.**Qt/(1+exp(.2*(Tlow-Tc)))/(1+exp(.3*(Tc-Thigh)))
!zzq24Feb2014         fT2=2.**Qt/(1.+exp(1.3*(Tc-Thigh)))

         exph1=exp(.3*(Tc-Thigh)/2.0)
         exph2=exp(1.3*(Tc-Thigh)/2.0)

         fT1=2.**Qt/(1+exp(.2*(Tlow-Tc)))/(1+exph1*exph1)
         fT2=2.**Qt/(1.+exph2*exph2)

         Vm=Vm0*fT1*fW
         Rd=0.025*Vm0*fT2*fW*BPI

         if (LT.eq.1) then
            A1=Vm*BPI
            A2=0
            A3=1
            A4=0
            A5=-Rd
         end if

         if (LT.eq.2) then
            A1=(PARdmol*Gmu1+PARfmol*Gmu2)*(1-0.15)*0.05
            A2=0
            A3=1
            A4=0
            A5=-Rd
         end if

         if (LT.eq.3) then
            A1=2e4*Vm/p*BPI
            A2=0
            A3=0
            A4=1
            A5=-Rd
         end if
      end if
!     to solve the cubic equation for Ci
!
!zzq20081002
        jmpflag=0
!
      if (A1.ne.0.and.A3.ne.0) then
         A6=A1+A3*A5
         A7=A4*A5-A1*A2

         B1=Ca*gb*A3-1.4*p*A6
         B2=Ca*gb*A4-1.4*p*A7

!zzq20081002   aa=A3*A6*m*hs*gb*gb*p+A3*B1*gb*b*LAI
         aa1=A3*A6*m*hs*gb*gb*p
         aa2=A3*B1*gb*b*LAI
         aa = aa1+aa2
       
         bb=gb*p*(A6*((1.6-m*hs)*B1+A4*m*hs*gb)+A3*A7*m*hs*gb)+    &
              b*LAI*(A4*B1*gb+A3*B2*gb-B1*B1)

!zzq20081002 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!zzq20091203        if(aa1.eq.-aa2.or.abs(bb/aa).gt.1.0e9) then
        if(abs(aa).lt.1.0e-20.or.abs(bb/aa).gt.1.0e9) then
           jmpflag=1
!zzq27OCT2009
           print*,"warning: jmp due to big values:",aa 
           goto 59
        endif
!zzq20081002 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

         cc=gb*p*(A6*(1.6-m*hs)*B2+A7*((1.6-m*hs)*B1+A4*m*hs*gb))+ &
              b*LAI*(A4*B2*gb-2*B1*B2)
         dd=A7*B2*(1.6-m*hs)*gb*p-B2*B2*b*LAI
         
         pp=cc/aa-bb*bb/3./aa/aa
         qq=dd/aa+2./27.*bb*bb*bb/aa/aa/aa-bb*cc/3./aa/aa

!zzq aug 29,2008 delta=(0.5*qq)**2+(pp/3.)**3
!zzq aug 29,2008 delta is too big, get grid of it.
!zzq aug 29,2008 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         

          qqx=0.5*qq
          ppx=pp/3.0
          if(abs(qqx).gt.1.e6) then
              qqx = qqx*1.e-6
              ppx = ppx*1.e-4
              deltax=qqx**2+ppx**3
              if(deltax.ge.0) then
                 sqdelta=1.0e6*sqrt(deltax)
              endif                                    
          else
              deltax=qqx**2+ppx**3
              if(deltax.ge.0) then
                 sqdelta=sqrt(deltax)
              endif
          endif
            
!zzq aug 29,2008    if (delta.gt.0) then
!zzq aug 29,2008          x11=-.5*qq+sqrt(delta)
!zzq aug 29,2008          x12=-.5*qq-sqrt(delta)

         if (deltax.gt.0) then
            x11=-.5*qq+sqdelta
            x12=-.5*qq-sqdelta
!rr vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!rr - by Huiping, June 15, 2007
!rr         x1=x11/abs(x11)*abs(x11)**(1/3.)+
!rr  &           x12/abs(x12)*abs(x12)**(1/3.)
            if(x11.ne.0.0) then
            x1a=x11/abs(x11)*abs(x11)**(1.0/3.0)
            else
            x1a=0.0
            endif

            if(x12.ne.0.0) then
            x1b=x12/abs(x12)*abs(x12)**(1.0/3.0)
            else
            x1b=0.0
            endif
            x1=x1a+x1b
!rr ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            Ci=x1-bb/3./aa
         end if

!zzq aug 29,2008   if (delta.eq.0) then
         if (deltax.eq.0) then
            x11=-.5*qq
!zzq 20071223       x1=2*x11/abs(x11)*abs(x11)**(1/3.)
            if(x11.ne.0) then
               x1=2*x11/abs(x11)*abs(x11)**(1/3.)
            else
               x1=0.0
            endif
            Ci=x1-bb/3./aa
         end if

!zzq aug 29,2008    if (delta.lt.0.and.pp.lt.0) then
         if (deltax.lt.0.and.pp.lt.0) then
!zzq dec 8,2008   rr=sqrt(-(pp/3.)**3)
            rr=(-pp/3.)*sqrt(-pp/3.)
            qqrr= -qq/2./rr
            qqrr=amin1(qqrr,1.0)
            qqrr=amax1(qqrr,-1.0)

!zzq oct 09,2008
!zzq oct 09,2008  theta=acosd(-qq/2./rr)/3.
!zzq Sept19,2013            theta=acosd(qqrr)/3.
            theta=acosdeg(qqrr)/3.
            x(1)=2*rr**(1/3.)*cosdeg(theta)-bb/3./aa
            x(2)=2*rr**(1/3.)*cosdeg(theta+120.)-bb/3./aa
            x(3)=2*rr**(1/3.)*cosdeg(theta+240.)-bb/3./aa

            xmin=1e10           ! to choose a valid solution
            do ii=1,3
               if (x(ii).gt.0) then
                  if (abs(x(ii)-Ca).lt.abs(xmin-Ca)) then
                     xmin=x(ii)
                  end if
               end if
            end do
            Ci=xmin
         end if

!zzq...by huiping,done by zzq,Dec 24 2007
!zzq         An=(A6*Ci+A7)/(A3*Ci+A4)
            a3ci=A3*ci+A4
          if (a3ci.eq.0.0) a3ci=0.0000001
!         An=(A6*Ci+A7)/(A3*Ci+A4)
          An=(A6*Ci+A7)/a3ci
!zzq
         Cs=Ca-1.4*p*An/gb
          if (cs.eq.0.0) cs=0.001
         gs=m*hs*p*An/Cs+b*LAI
         end if

!     to solve the quadratic equation for Ci
!
!zzq20081002      if (A1.ne.0.and.A3.eq.0) then
!
   59   continue
      if (A1.ne.0.and.A3.eq.0.or.jmpflag.eq.1) then
         tmp1=1.6-m*hs
         tmp2=1.4*p*A1
         tmp3=Ca*gb-1.4*p*A5

         a0=gb*p*A1*(tmp2*tmp1-m*hs*gb)+b*LAI*(tmp2*(tmp2+gb))
         b0=-gb*p*A1*tmp1*tmp3+A5*gb*p*(tmp2*tmp1-m*hs*gb)-    &
              b*LAI*(2*tmp2+gb)*tmp3
         c0=-gb*p*A5*tmp1*tmp3+b*LAI*tmp3**2

         delta=b0*b0-4.*a0*c0
         
         if (delta.gt.0) then
            x1=(-b0-sqrt(delta))/2./a0
            x2=(-b0+sqrt(delta))/2./a0

            if (x1.gt.0) then
               Ci=x1
               if (x2.gt.0.and.(Ca-x2).lt.(Ca-x1)) Ci=x2
            else
               if (x2.gt.0) then
                  Ci=x2
               else
                  Ci=Ca
               end if
            end if
         end if
         
         if (delta.eq.0) then
            x1=-b0/2./a0
            if (x1.gt.0) then
               Ci=x1
            else
               Ci=Ca
            end if
         end  if

                                ! check the validity of the Ci
         xx=PARdmol+PARfmol
         if ((delta.lt.0.or.Ci.le.0).and.xx.gt.0) then
            Ci=Ca*0.5
         end if

         An=A1*Ci+A5
         Cs=Ca-1.4*p/gb*An
!zzq by huiping,done by zzq, Dec 24, 2007
           if(cs.eq.0.0)cs=0.001
         gs=m*hs*p*An/Cs+b*LAI
      end if

!     to solve Ci when An is constant

!zzq by huiping,done by zzq, Dec 24, 2007
!zzq      if (A1.eq.0) then
!zzq         An=A5
!zzq         Cs=Ca-1.4*p*An/gb
!zzq         gs=m*hs*p*An/Cs+b*LAI
!zzq         Ci=Cs-1.6*p*An/gs
!zzq      end if
      if (A1.eq.0) then
         An=A5
         Cs=Ca-1.4*p*An/gb
!zzq
          if(Cs.eq.0.0) Cs=0.001
         gs=m*hs*p*An/Cs+b*LAI
          if(gs.lt.0.00001)gs=0.00001
         Ci=Cs-1.6*p*An/gs
      end if

!-----------------------------------------------------------------------
! zzq Sept28,2013
! added by zzq, in agreement with that in subroutine ANGS2,
! since I found gs could be negative which leads to gs+gb=0 and
! causes floating overflow in the following formula. 
!-----------------------------------------------------------------------
!
       gs  = amax1(gs,0.00001)
!
!zzq sept28,2013,ended 
!-----------------------------------------------------------------------
!
      es=(gs*ei+gb*ea)/(gs+gb)
!khs  print *,'angsa: gs ei gb ea es = ',gs,ei,gb,ea,es

      Respd=Rd

      end subroutine

!-----------------------------------------------------------------------
      function e(x) result(f_result)
      f_result = 100.*exp(21.18123-5418./x)/.622 ! [Pa] with x in [degree K]
      end function
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ***************
!zzq sept25,2013: The following two functions cosdeg and acosdeg are added
!   to replace the two original functions cosd and acosd, since the two
!   are not recoginzed by some compilers or in different systems, such as 
!   LINUX and UNIX.
!----------------------------------------------------------------------
!
      function cosdeg(x) result(f_result)
!
!----------------------------------------------------------------------
!  The input for this function is degree rather than radian.
!----------------------------------------------------------------------
        real,intent(in)::x            !  degree angle
        real:: f_result
        f_result=cos(x*PIE/180.0)
!
      end function
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      function acosdeg(x) result(f_result)
!
!----------------------------------------------------------------------
!  This function is to convert the radian to degree after calculation.
!----------------------------------------------------------------------
        real,intent(in)::x
        real:: f_result
        f_result=(180.0/PIE)*acos(x)
!
      end function
!
END MODULE
