!
! mainly used to convert vegetation types from 13 types to 9 types
! zzq, sept28,2013
!
MODULE ssib4_module_vegcnvt
    use ssib4_module_trifparms
    use ssib4_module_comsveg
    use ssib4_module_comssib4sfc
!
    implicit none
contains
     subroutine changeVT8(datafs,xlon,xlat,MONTH,ITYPE,FRAC,LAI) 
      type(ssibdataios)::datafs
      integer:: ITYPE,N,MONTH
      real:: xlon,xlat,rnum
      integer:: itypexx
      REAL,dimension(NTYPES),intent(inout):: frac    ! Areal coverage
      REAL,dimension(NPFT)  ,intent(inout):: LAI     ! Leaf area index
!
!yliu DBL        if(ITYPE.ne.8) return

        if(ITYPE.eq.9) then
        rnum=num_section(xlon,xlat,9)

        DO N=1,NPFT
!zzq 2009Feb23
!zzq 2009Feb23:taking the itypexx as the type with the deepest root.
          if(rnum.eq.1) then
              if(N.eq.5) then
                FRAC(N) = VCOVER0(N,MONTH,1)
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
              itypexx= 5
          endif
!
          if(rnum.eq.2) then
              if(N.eq.4) then
                FRAC(N) = VCOVER0(N,MONTH,1)
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
              itypexx= 4
          endif
!
          if(rnum.eq.6) then
!yliu DBL              if(N.eq.1) then
              if(N.eq.7) then
                FRAC(N) = 0.35
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else if(N.eq.2) then
                FRAC(N) = 0.35
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else if(N.eq.4) then
                FRAC(N) = 0.30
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
!yliu09May2017              itypexx= 1
              itypexx= 4
          endif
!yliu add 3 
          if(rnum.eq.3) then
!yliu DBL              if(N.eq.1) then
              if(N.eq.7) then
                FRAC(N) = 0.10
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else if(N.eq.2) then
                FRAC(N) = 0.10
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else if(N.eq.4) then
                FRAC(N) = 0.80
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
!yliu09May2017              itypexx= 1
              itypexx= 4
          endif
!
          if(rnum.eq.4.or.rnum.eq.8) then
              if(N.eq.4) then
                FRAC(N) = VCOVER0(N,MONTH,1)
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
              itypexx= 4
          endif
!
          if(rnum.eq.5) then
              if(N.eq.4) then
                FRAC(N) = 0.35
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else if(N.eq.5) then
                FRAC(N) = 0.35
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else if(N.eq.1) then
                FRAC(N) = 0.30
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
              itypexx= 1
          endif
!
!yliu DBL          if(rnum.eq.7) then
!yliu DBL              if(N.eq.2) then
!yliu DBL                FRAC(N) = 0.70
!yliu DBL                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
!yliu DBL              else if(N.eq.4) then
!yliu DBL                FRAC(N) = 0.30
!yliu DBL                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
!yliu DBL              else
!yliu DBL                FRAC(N)=1.0E-6
!yliu DBL                LAI(N)=1.0E-6
!yliu DBL              endif
!yliu DBL              itypexx=2
!yliu DBL          endif
!
          if(rnum.eq.9) then
              if(N.eq.2) then
                FRAC(N) = VCOVER0(N,MONTH,1)
                LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
              else
                FRAC(N)=1.0E-6
                LAI(N)=1.0E-6
              endif
              itypexx = 2
          endif
!zzq 2009Feb23: update vegetation type 8.
         datafs % xvtype=itypexx
         ITYPE = datafs % xvtype
          ENDDO
        endif
!yliu DBL
          if(ITYPE.eq.4) then
          itypexx = 4
          rnum=num_section(xlon,xlat,4)
          DO N=1,NPFT
            if(rnum.eq.41.or.rnum.eq.42) then
                if(N.eq.4) then
                   FRAC(N) = 0.70
                   LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
                else if(N.eq.5) then
                   FRAC(N) = 0.30
                   LAI(N)=ZLT0(N,MONTH,1)*GREEN0(N,MONTH,1)
                else
                   FRAC(N)=1.0E-6
                   LAI(N)=1.0E-6
                endif
                itypexx=4
            endif
            datafs % xvtype=itypexx
            ITYPE = datafs % xvtype
          ENDDO
          endif
!
!
         end subroutine
!=======================================================================
!
      SUBROUTINE vegcnvt(vg,xxlat)
!                                                          
!=======================================================================
        real::vg,xxlat,veg1
     
         veg1=vg
         if(vg.eq.1) veg1=1
!yliu DBL         if(vg.eq.2)then
!yliu DBL                 if(abs(xxlat).gt.45) then
!yliu DBL                      veg1=2
!yliu DBL                 else
!yliu DBL                      veg1=1
!yliu DBL                 endif
!yliu DBL         endif
         if(vg.eq.2) veg1=7
         if(vg.eq.3)then
                 if(abs(xxlat).gt.30) then
                      veg1=2
                 else
                      veg1=1
                 endif
         endif
         if(vg.eq.4) veg1=2
         if(vg.eq.5) veg1=2
         if(vg.eq.6) then
              if(abs(xxlat).gt.30) then
                     veg1=5    !c3
              else
                     veg1=4    !c4
              endif
          endif 
          if(vg.eq.7) veg1=3   ! C3
          if(vg.eq.8) veg1=4   ! C4
          if(vg.eq.9) veg1=5
!zzq2012dec06          if(vg.eq.10) veg1=6
          if(vg.eq.10) then
                 veg1=6
              if(abs(xxlat).lt.50) veg1=5 
          endif
!yliu DBL          if(vg.eq.11) veg1=7
!yliu DBL          if(vg.eq.12) veg1=8
!yliu DBL          if(vg.eq.13) veg1=9
          if(vg.eq.11) veg1=8
          if(vg.eq.12) veg1=9
          if(vg.eq.13) veg1=10
!
          vg=veg1

      END SUBROUTINE
!=======================================================================
      function num_section(xlo,xla,xtp) result(f_result)
      integer:: snum,f_result,xtp
      real::xlo,xla
       
!Austrilia region
      if((xlo.ge.60.and.xlo.lt.180).and.(xla.ge.-90.and.xla.lt.0)) then
        snum=1
      endif
!
!Southern America
      if((xlo.ge.180.and.xlo.lt.310).and.(xla.ge.-90.and.xla.lt.0)) then
        snum=2
      endif
!
!Northern America
      if((xlo.ge.180.and.xlo.lt.310).and.(xla.ge.0.and.xla.lt.90)) then
        snum=3
      endif
!
!India+middle Southern Peninsula
      if((xlo.ge.60.and.xlo.lt.100).and.(xla.ge.0.and.xla.lt.30)) then
        snum=4
      endif
!
!Middle Europe
      if((xlo.ge.50.and.xlo.lt.100).and.(xla.ge.30.and.xla.lt.90)) then
        snum=5
      endif
!
!near to north China
      if((xlo.ge.100.and.xlo.lt.180).and.(xla.ge.30.and.xla.lt.90)) then
        snum=6
      endif
!
!near to south China
      if((xlo.ge.100.and.xlo.lt.180).and.(xla.ge.0.and.xla.lt.30)) then
        snum=7
      endif
!
!Afrca
      if((xlo-360.ge.-50.and.xlo-360.lt.0.or.xlo.ge.0.and.xlo.lt.60)   &
           .and.(xla.ge.-90.and.xla.lt.30)) then
        snum=8
      endif
!
!Western Europe
      if((xlo-360.ge.-50.and.xlo-360.lt.0.or.xlo.ge.0.and.xlo.lt.60)   &
           .and.(xla.ge.30.and.xla.lt.90)) then
        snum=9
      endif
!ly2014.3.19
!India
      if(xtp.eq.4) then
      if((xlo.ge.60.and.xlo.lt.100).and.(xla.ge.0.and.xla.lt.30)) then
        snum=41
      endif
!Afrca
      if((xlo-360.ge.-20.and.xlo-360.lt.0.or.xlo.ge.0.and.xlo.lt.60) &
           .and.(xla.ge.10.and.xla.lt.15)) then
        snum=42
      endif
      endif
      f_result = snum

      END FUNCTION
END MODULE
