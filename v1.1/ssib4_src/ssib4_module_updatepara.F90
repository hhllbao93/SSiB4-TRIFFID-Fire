MODULE ssib4_module_updatepara
      implicit none
contains
!=======================================================================
!
      SUBROUTINE updatepara(ITYPE,HT_TRIF,GREEN_TRIF,LAI,SAI,  &
                       Z0,XDD,RBC,RDC,Z2,Z1,ZLT1,GREEN1,VCOVER1)
!                                        Zhengqiu Zhang, July 23,2008
!=======================================================================
!  
!-------------------------------------------------------------------
      integer:: ITYPE
      real:: ZLT1,GREEN1,VCOVER1
      real:: HT_TRIF,GREEN_TRIF,LAI,SAI
      real:: Z0,XDD,RBC,RDC,z1,z2
      real:: xdd1,rbc1,rdc1
!references for z2 and z1 when obtained the parameters
!yliu09May2017 addtype --
!      real,dimension(9):: z1x,z2x
!      data z2x/35.0,37, 0.6, 3.6, 5, 0.6,  0.1, 20.0,  0.1/
!      data z1x/1.0, 17, 0.1, 2,   2, 0.1, 0.01, 11.5, 0.01/
      real,dimension(10):: z1x,z2x
      data z2x/35.0,37, 0.6, 3.6, 5, 0.6, 20.0,  0.1, 20.0,  0.1/
      data z1x/1.0, 17, 0.1, 2,   2, 0.1,  1.2, 0.01, 11.5, 0.01/
!!!--
!     
       z2=HT_TRIF
       z1 = z2*z1x(itype)/z2x(itype)
       if(z1.lt.0.1) z1=0.1
       if(z2.le.z1) z2 =z1+0.05
!
!zzq this is the fitting 07/29/08
!zzq      z0 = 0.10*z2
!
!From Sellers et al. (1996),PP710.
       if(lai.gt.0.1) then
!zzq20081002    z0=z2*(1-0.91*exp(-0.0075*lai))
          z0=z2*(1-0.905*exp(-0.0035*lai))
       else
          z0=0.0906*z2
       endif

      z0 = max(z0,0.01122)

! for type 1
      if(itype.eq.1) then
        if((lai.ge.0.53).and.(z2.ge.1.40)) then
          xdd=0.3703*z2*(1-exp(-0.2*lai))+0.3436*z2-0.0389*lai-0.3132
          rbc=21.7890/lai-0.0038*z2-0.1072*lai+2.3485
          rdc=1.9194*lai*z2*log(1.+10./z2)+1.9681*z2+6.6121*lai+7.5986
        else
          xdd=0.3730*lai+0.2562*z2+0.00060
          rbc=-69.1119*lai-0.0026*z2+80.0
          rdc=-10.5262*lai+1.8547*z2+22.5
        end if
! for type 2
      else if(itype.eq.2) then
        if((lai.ge.0.54).and.(z2.ge.1.48)) then
          xdd=0.2732*z2*(1-exp(-0.2*lai))+0.5309*z2-0.0500*lai+0.1993
          rbc=2.9273/lai-0.0003*z2-0.0130*lai+0.2792
          rdc=4.7607*lai*z2*log(1.+20./z2)+4.0913*z2-7.5162*lai+61.2243
        else
          xdd=1.8497*lai+0.3669*z2+0.00060
          rbc=-17.2360*lai-0.0002*z2+15.0
          rdc=94.2568*lai+4.7720*z2+22.5
        end if
! for type 3
      else if(itype.eq.3) then
        if((lai.ge.0.29).and.(z2.ge.0.02)) then
          xdd=0.4339*z2*(1-exp(-0.4*lai))+0.1591*z2-0.0026*lai+0.0412
          rbc=9.7871/lai+2.2529*z2-0.0727*lai+0.2241
          rdc=29.9079*lai*z2*log(1.+1/z2)+1.1139*z2-6.0005*lai+22.4101
        else
          xdd=0.1562*lai+0.1833*z2+0.00060
          rbc=-159.0758*lai+2.0040*z2+80.0
          rdc=-3.5268*lai+12.2101*z2+22.5
        end if
! for type 4
      else if(itype.eq.4) then
        if((lai.ge.0.35).and.(z2.ge.0.14)) then
          xdd=0.3868*z2*(1-exp(-0.4*lai))+0.4039*z2-0.0123*lai+0.0291
          rbc=26.1456/lai-0.0172*z2-0.0891*lai+1.7106
          rdc=2.6074*lai*z2*log(1.+20./z2)+18.4615*z2-0.0353*lai+22.3821
        else
          xdd=0.2585*lai+0.1889*z2+0.00060
          rbc=-210.2962*lai-0.0071*z2+150.0
          rdc=9.1484*lai+9.0603*z2+22.5
        end if
! for type 5
      else if(itype.eq.5) then
        if((lai.ge.0.10).and.(z2.ge.0.20)) then
          xdd=0.8867*z2*(1-exp(-0.4*lai))+0.2650*z2-0.0485*lai+0.0294
          rbc=5.0966/lai-0.0571*z2-0.0942*lai+0.8352
          rdc=5.5541*lai*z2*log(1.+15./z2)+13.1841*z2-1.7100*lai+24.7814
        else
          xdd=0.8390*lai+0.0213*z2+0.00060
          rbc=-482.2008*lai-0.0040*z2+100.0
          rdc=52.2827*lai+1.0660*z2+22.5
        end if
! for type 6
      else if(itype.eq.6) then
        if((lai.ge.0.13).and.(z2.ge.0.02)) then
          xdd=0.9716*z2*(1-exp(-0.2*lai))+0.1325*z2+0.0031*lai+0.0300
          rbc=6.1682/lai+2.3327*z2-0.0818*lai-0.1133
          rdc=7.0649*lai*z2*log(1.+10/z2)+6.7953*z2-0.7850*lai+20.5034
        else
          xdd=0.2652*lai+0.1418*z2+0.00060
          rbc=-405.5650*lai+2.1106*z2+100.0
          rdc=-14.2372*lai+8.9501*z2+22.5
        end if
!yliu09May2017 --
      else if(itype.eq.7) then
        if((lai.ge.1.50).and.(z2.ge.1.48)) then
        xdd=-0.0554*z2*(1-exp(-1.5316*lai))+0.7847*z2-0.0530*lai+0.0127
        rbc=26.0688/lai-0.2045*z2+0.0301*lai+7.9867
        rdc=2.0310*lai*z2*log(1.+56.8035/z2)+3.0895*z2-20.8635*lai+105.7302
      else
        xdd=1.3316*lai+0.7164*z2-1.6684
        rbc=-33.4086*lai-0.1842*z2+70.1764
        rdc=22.4314*lai+5.6418*z2+70.0215
      end if
      rdc =amin1(rdc,400.)

      end if
!
      rbc =amax1(rbc,0.2)
      rdc =amax1(rdc,22.5)
!
!zzq..
!note if itype is not within the 6, take the original parameters.
!
!yliu09May2017 addtype       if (ITYPE.le.6)THEN
       if (ITYPE.le.7)THEN
           zlt1=LAI+SAI
           vcover1=1.0d0
           green1=green_trif
       end if

       END SUBROUTINE
END MODULE
