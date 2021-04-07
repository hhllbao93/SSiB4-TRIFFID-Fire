!
! For weighting variables from vegetation for bareland
!  zzq, sept 28,2013
! Yeliu fixed bugs Oct.18 2018
!
MODULE ssib4_module_weighting
      use ssib4_module_comsveg
      use ssib4_module_trifparms, only:NPFT
!
      implicit none
contains
!=======================================================================
!
      subroutine update_soil(m,FRAC,month,itype,istrt)
!                                           zzq, March 22,2010
!=======================================================================
!
      INTEGER it,IW,IV,m,month,itype,istrt
      real,dimension(M):: FRAC
!temporary varaibles
      real,dimension(13)::                                              &
           xTRAN11,xTRAN12,xTRAN21,xTRAN22,xREF11,xREF12,xREF21,xREF22, &
           xRSTPAR1,xRSTPAR2,xSOREF,xCHIL,xTOPT,xTL,xTU,xDEFAC,         &
           xPH1,xPH2,xROOTD, xZLT,xGREEN,xVCOVER, xvm0,xfc3,            &
           xBEE,xPHSAT,xSATCO,xSLOPE,xZ2,xZ1,xzZ0,xDD,xRBC,xRDC
      real ::vfrac
!
!yliu09May2017 addtype     if(itype.gt.7) return
      if(itype.gt.JBTYP) return
      if(istrt.eq.1) then
         call dominant_depth(M,FRAC,itype)
      endif
!
      call bareland
!
!yliu09May2017 addtype      if(itype.eq.7) return
      if(itype.eq.JBTYP) return
!
!yliu09May2017 addtype      vfrac=frac(1)+frac(2)+frac(3)+frac(4)+frac(5)+frac(6)
      vfrac = 0.
      do it = 1,NPFT
         vfrac=vfrac+frac(it)
      enddo
      if(vfrac.le.0.05) return
! 
      DO IW=1,3
!yliu09May2017         do it=1,6
         do it=1,NPFT
!          xTRAN11(it)=TRAN0(it,1,IW,1)
!          xTRAN12(it)=TRAN0(it,1,IW,2)
!          xTRAN21(it)=TRAN0(it,2,IW,1)
!          xTRAN22(it)=TRAN0(it,2,IW,2)
!          xREF11(it)= REF0(it,1,IW,1)
!          xREF12(it)= REF0(it,1,IW,2)
!          xREF21(it)= REF0(it,2,IW,1)
!          xREF22(it)= REF0(it,2,IW,2)
!          xRSTPAR1(it)=RSTPAR0(it,1,IW)
!          xRSTPAR2(it)=RSTPAR0(it,2,IW)
           xSOREF(it) =SOREF0(it,IW)
         enddo
!
!        call  avelinearly(m,FRAC,vfrac,xTRAN11)
!        call  avelinearly(m,FRAC,vfrac,xTRAN12)
!        call  avelinearly(m,FRAC,vfrac,xTRAN21)
!        call  avelinearly(m,FRAC,vfrac,xTRAN22)
!        call  avelinearly(m,FRAC,vfrac,xREF11)
!        call  avelinearly(m,FRAC,vfrac,xREF12)
!        call  avelinearly(m,FRAC,vfrac,xREF21)
!        call  avelinearly(m,FRAC,vfrac,xREF22)
!        call  avelinearly(m,FRAC,vfrac,xRSTPAR1)
!        call  avelinearly(m,FRAC,vfrac,xRSTPAR2)
         call  avelinearly(m,FRAC,vfrac,xSOREF)
!
!        TRAN0(7,1,IW,1)=xTRAN11(7)
!        TRAN0(7,1,IW,2)=xTRAN12(7)
!        TRAN0(7,2,IW,1)=xTRAN21(7)
!        TRAN0(7,2,IW,2)=xTRAN22(7)
!        REF0(7,1,IW,1)=xREF11(7)
!        REF0(7,1,IW,2)=xREF12(7)
!        REF0(7,2,IW,1)=xREF21(7)
!        REF0(7,2,IW,2)=xREF22(7)
!        RSTPAR0(7,1,IW)=xRSTPAR1(7)
!        RSTPAR0(7,2,IW)=xRSTPAR2(7)
!yliu09May2017 addtype         SOREF0(7,IW)=xSOREF(7)
         SOREF0(JBTYP,IW)=xSOREF(JBTYP)
       END DO
!
       DO IV=1,2
!yliu09May2017         do it=1,6
         do it=1,NPFT
!          xCHIL(it) = CHIL0(it,IV)
!          xTOPT(it) = TOPT0(it,IV)
!          xTL(it)   = TL0(it,IV)
!          xTU(it)   = TU0(it,IV)
!          xDEFAC(it)=DEFAC0(it,IV)
!          xPH1(it)  = PH10(it,IV)
!          xPH2(it)  = PH20(it,IV)
!          xROOTD(it)=ROOTD0(it,IV)
!          xZLT(it)  =ZLT0(it,MONTH,IV)
!          xGREEN(it)=GREEN0(it,MONTH,IV)
!          xVCOVER(it)=VCOVER0(it,MONTH,IV)
         enddo

!        call  avelinearly(m,FRAC,vfrac,xCHIL)
!        call  avelinearly(m,FRAC,vfrac,xTOPT)
!        call  avelinearly(m,FRAC,vfrac,xTL)
!        call  avelinearly(m,FRAC,vfrac,xTU)
!        call  avelinearly(m,FRAC,vfrac,xDEFAC)
!        call  avelinearly(m,FRAC,vfrac,xPH1)
!        call  avelinearly(m,FRAC,vfrac,xPH2)
!        call  avelinearly(m,FRAC,vfrac,xROOTD)
! song 4.14.2011
!         call  avelinearly(m,FRAC,vfrac,xZLT)
!         call  avelinearly(m,FRAC,vfrac,xGREEN)
!         call  avelinearly(m,FRAC,vfrac,xVCOVER)
!yliu09May2017 addtype --
!         xZLT(7)=0.0001
!         xVCOVER(7)=0.00001
!         xGREEN(7)=0.0001
         xZLT(JBTYP)=0.0001
         xVCOVER(JBTYP)=0.00001
         xGREEN(JBTYP)=0.0001
!!!--
!
!        CHIL0(7,IV)=xCHIL(7)
!        TOPT0(7,IV)=xTOPT(7)
!        TL0(7,IV)= xTL(7)
!        TU0(7,IV)=xTU(7)
!        DEFAC0(7,IV)=xDEFAC(7)
!        PH10(7,IV)= xPH1(7)
!        PH20(7,IV)= xPH2(7)
!        ROOTD0(7,IV)= xROOTD(7)
!yliu09May2017 addtype --
!         ZLT0(7,MONTH,IV)= xZLT(7)
!         GREEN0(7,MONTH,IV)= xGREEN(7)
!         VCOVER0(7,MONTH,IV)=xVCOVER(7)
         ZLT0(JBTYP,MONTH,IV)= xZLT(JBTYP)
         GREEN0(JBTYP,MONTH,IV)= xGREEN(JBTYP)
         VCOVER0(JBTYP,MONTH,IV)=xVCOVER(JBTYP)
!!!--
       END DO
!
!yliu09May2017 addtype       do it=1,6
       do it=1,NPFT
!         xvm0(it)=VM00(it)
!         xfc3(it)=fc30(it)
!        the following needs to test.
          xBEE(it)=BEE0(it)
          xPHSAT(it)=PHSAT0(it)
          xSATCO(it)=SATCO0(it)
          xSLOPE(it)=SLOPE0(it)
!         xZ2(it)=Z20(it,MONTH)
!         xZ1(it)=Z10(it,MONTH)
!         xzZ0(it)= ZZ0(it,MONTH)
!         xDD(it)= D0(it,MONTH)
!         xRBC(it)=RBC0(it,MONTH)
!         xRDC(it)=RDC0(it,MONTH)
       enddo

!      call  avelinearly(m,FRAC,vfrac,xvm0)
!      call  avelinearly(m,FRAC,vfrac,xfc3)
!      the following needs testing (XUE oct. 7. 2012)
       call  avelinearly(m,FRAC,vfrac,xBEE)
       call  avelinearly(m,FRAC,vfrac,xPHSAT)
       call  avelogarthmic(m,FRAC,vfrac,xSATCO)
       call  avelinearly(m,FRAC,vfrac,xSLOPE)
!song 4.14.2011
!       call  avelinearly(m,FRAC,vfrac,xZ2)
!       call  avelinearly(m,FRAC,vfrac,xZ1)
!       call  avelinearly(m,FRAC,vfrac,xzZ0)
!       call  avelinearly(m,FRAC,vfrac,xDD)
!       call  avelinearly(m,FRAC,vfrac,xRBC)
!       call  avelinearly(m,FRAC,vfrac,xRDC)
!      call  avelogarthmic(m,FRAC,vfrac,xZ2)
!      call  avelogarthmic(m,FRAC,vfrac,xZ1)
!      call  avelogarthmic(m,FRAC,vfrac,xzZ0)
!      call  avelogarthmic(m,FRAC,vfrac,xDD)
!      call  avelogarthmic(m,FRAC,vfrac,xRBC)
!      call  avelogarthmic(m,FRAC,vfrac,xRDC)
!
!      VM00(7)=xvm0(7)
!       fc30(7)=xfc3(7)
!yliu09May2017 addtype --
!       BEE0(7)=xBEE(7)
!       PHSAT0(7)=xPHSAT(7)
!       SATCO0(7)=xSATCO(7)
!       SLOPE0(7)=xSLOPE(7)
       BEE0(JBTYP)=xBEE(JBTYP)
       PHSAT0(JBTYP)=xPHSAT(JBTYP)
       SATCO0(JBTYP)=xSATCO(JBTYP)
       SLOPE0(JBTYP)=xSLOPE(JBTYP)
!!!--
!      Z20(7,MONTH)=xZ2(7)
!      Z10(7,MONTH)=xZ1(7)
!      ZZ0(7,MONTH)=xzZ0(7)
!      D0(7,MONTH)=xDD(7)
!      RBC0(7,MONTH)=xRBC(7)
!      RDC0(7,MONTH)=xRDC(7)
!
       END SUBROUTINE
!=======================================================================
!
      subroutine avelinearly(M,FRAC,vfrac,x)
!                                             March 03,2010 by zzq
!=======================================================================
!
      INTEGER:: M,K
      real,dimension(M)::frac
      real,dimension(13)::x
      real::vfrac,sum
!
      sum=0.
      DO K=1,NPFT
         sum = sum+x(k)*FRAC(K)
      ENDDO
      x(JBTYP)=sum/vfrac
!
      end subroutine
!=======================================================================
!
      subroutine avelogarthmic(M,FRAC,vfrac,x)
!                                             March 03,2010 by zzq
!=======================================================================
!
      INTEGER:: M,K
      real,dimension(M)::frac
      real,dimension(13)::x
      real::vfrac,sum

      sum=0.
      DO K=1,NPFT
         if(x(k).le.0) then 
           write(*,*) 'error in avelogarthmic due to negative'
         endif
         sum = sum+alog(x(k))*FRAC(K)
      ENDDO
      x(JBTYP)=exp(sum/vfrac)
!
      end subroutine

!=======================================================================
!
      subroutine dominant_depth(M,FRAC,itype)
!                                           March 22,2010
!=======================================================================
!
! Note: the returned type  is only for setting the dominant soil depths, 
!       not really represent the real dominant type. 
!
      INTEGER:: M,K,itype,ity
      real,dimension(M)::frac
      real,dimension(5)::fracx
      real::fracxm 

!yliu09May2017 addtype      fracx(1) = frac(1)+frac(2)          ! depths 0.02,1.48,2.00 for T1 and T2
      fracx(1) = frac(1)+frac(2)+frac(7)  ! depths 0.02,1.48,2.00 for T1 and T2
      fracx(2) = frac(3)+frac(4)          ! depths 0.02,0.47,1.00 for T3 and T4
      fracx(3) = frac(5)                  ! depths 0.02,0.47,1.00 for T5
      fracx(4) = frac(6)                  ! depths 0.02,0.17,1.00 for T6
!yliu09May2017 addtype      fracx(5) = frac(7)                  ! depths 0.02,0.17,0.30 for T7
      fracx(5) = frac(8)                  ! depths 0.02,0.17,0.30 for T7
!
      fracxm=fracx(1)+fracx(2)+fracx(3)+fracx(4)
      if(fracxm.le.0.05) then
!yliu09May2017 addtype        itype=7
        itype=8
        return
      endif

      fracxm=fracx(1)
      ity=1
      DO K=2,4
         if(fracx(K).gt.fracxm) then
              fracxm=fracx(k)
              ity=k
          endif
      ENDDO
      if(ity.eq.1) itype =1
      if(ity.eq.2) itype =3
      if(ity.eq.3) itype =1  !taking the soil depth of Type1 for Shrub
      if(ity.eq.4) itype =6
!
      end  subroutine
!=======================================================================
!
      subroutine dominant_type(M,FRAC,itype)
!                                           March 22,2010
!=======================================================================
!
      INTEGER:: M,K,itype
      real,dimension(M):: frac
      real:: fracx

      itype=1 
      fracx=FRAC(1)
      DO K=2,M
         if(FRAC(K).gt.fracx) then
              fracx=FRAC(k)
              itype=k
          endif
      ENDDO
!
      end subroutine
!=======================================================================
!
      function summ(N1,N2,M,FRAC,xneglect,AA) result(f_result)
!
!=======================================================================
!
      INTEGER:: N1,N2,M
      real,dimension(M):: FRAC       ! IN Areal coverage.
      real,dimension(M):: AA         ! IN Variable

      integer:: K
      real:: vfrac,xneglect,sum,f_result

      VFRAC=0.
      sum=0.
      f_result=0.

      DO K=N1,N2
        if(FRAC(K).GE.xneglect) THEN
           sum =sum +FRAC(K)*AA(K)
           VFRAC= VFRAC +FRAC(K)
        ENDIF
      ENDDO

      if(VFRAC.NE.0.) then
       f_result=sum/VFRAC
      endif

      end function
!=======================================================================
!
      function sumu(N1,N2,M,FRAC,navedex,AA) result(f_result)
!
!=======================================================================
!
      INTEGER:: N1,N2,M,K
      real,dimension(M):: FRAC,AA
      integer,dimension(M):: navedex
      real:: vfrac,sum,f_result

      VFRAC=0.
      sum=0.
      f_result=0.

      DO K=N1,N2
        if(navedex(K).eq.1.and.AA(K).ne.-999) THEN
           sum =sum +FRAC(K)*AA(K)
           VFRAC= VFRAC +FRAC(K)
        ENDIF
      ENDDO

      if(VFRAC.NE.0.) then
         f_result=sum
       else
         f_result=-999.0
      endif
!
      end function
END MODULE
