!=========================================================================
! This module is for storing data in memory
! designed by Zhengqiu Zhang, Sept26, 2013
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
MODULE module_records
     use ssib4_module_comssib4sfc
     use module_comsdef
     use module_string
     use module_netcdf
!
     implicit none
     integer,parameter::nvdesnum=171,regnum=2                   !Huilin add 6+1+12+3+2+7 variables
     integer,dimension(REGNUM):: nvdes_tag,nvdes,vartag
     character(len=80),dimension(REGNUM)::stitle
     character(len=80),dimension(REGNUM,nvdesnum)::gvarlongname
     character(len=20),dimension(REGNUM,nvdesnum)::gvarname
     character(len=20),dimension(REGNUM,nvdesnum)::gvarunit
!        
     integer,parameter:: MYNUM=171,MAXHR=8                      !Huilin add 6+12+3+2+7 variables
     real,dimension(lonf2s:lonf2e,latg2s:latg2e,MYNUM):: vouts,vcount
     real,dimension(lonf2s:lonf2e,latg2s:latg2e) :: landmsk
     real,dimension(lonf2s:lonf2e,latg2s:latg2e,MYNUM,MAXHR)::voutsh,counth
!
     real, dimension(lonf2s:lonf2e)                        :: xlongit
     real, dimension(latg2s:latg2e)                        :: xlatit
     type(ssibdataios),dimension(lonf2s:lonf2e,latg2s:latg2e) :: ssibios
     type(ssibdataout),dimension(lonf2s:lonf2e,latg2s:latg2e) :: dataout
!
contains
      subroutine initial
       integer::i,j,k,ihr
!
       do i=lonf2s,lonf2e
       do j=latg2s,latg2e
          landmsk(i,j)= 0.0
       enddo
       enddo
!
       do k=1,MYNUM
       do i=lonf2s,lonf2e
       do j=latg2s,latg2e
          vouts(i,j,k)= 0.0
          vcount(i,j,k)= 0.0
          do ihr=1,MAXHR
            voutsh(i,j,k,ihr)= 0.0
            counth(i,j,k,ihr)= 0.0
          enddo
       enddo
       enddo
       enddo
!
      end subroutine
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine recsh(nrg,dtt,nhr,n,x,s1,s2,s3)
!
      character(len=*):: s1,s2,s3
      integer:: i,j,k,n,nhr,ihr,ihour,nrg,kk
      real,dimension(LONF2S:LONF2E,LATG2S:LATG2E,n)::x
      real:: hourstep,dtt
!
       hourstep= dtt/3600.
!
       if(hourstep.eq.0) then
          write(*,*) "error for hourstep=",hourstep
       endif
!
       ihour=mod(nhr,24)
       ihr = int(ihour/hourstep)+1

       do k=1,n
!
       call regs(nrg,s1,s2,s3,n,k)
       kk = nvdes(nrg)
!
        do i=lonf2s,lonf2e
        do j=latg2s,latg2e
          if(x(i,j,k).ne.-999.0.and.landmsk(i,j).eq.1) then
            voutsh(i,j,kk,ihr)=voutsh(i,j,kk,ihr)+x(i,j,k)
            counth(i,j,kk,ihr)= counth(i,j,kk,ihr)+1
          endif
         enddo
       enddo
       enddo
!
      end subroutine
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine records(nrg,n,x,s1,s2,s3)
!
      character(len=*):: s1,s2,s3
      integer,intent(in)::nrg
      integer:: i,j,k,n,kk
      real,dimension(LONF2S:LONF2E,LATG2S:LATG2E,n)::x
!
       do k=1,n
!
       call regs(nrg,s1,s2,s3,n,k)
       kk = nvdes(nrg)
!
       do i=LONF2S,LONF2E
        do j=LATG2S,LATG2E
          if(x(i,j,k).ne.-999.0.and.landmsk(i,j).eq.1) then
            vouts(i,j,kk)=vouts(i,j,kk)+x(i,j,k)
            vcount(i,j,kk)= vcount(i,j,kk)+1
          endif
         enddo
       enddo
       enddo
!       
       end subroutine
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine seton_reg_global(nrg,stit)
      character(Len=*):: stit
      integer:: nrg
!
      stitle(nrg) =stit
      vartag(nrg) = 99
!
      end subroutine
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine setoff_reg_global(nrg) 
!
      integer:: nrg
!
      vartag(nrg) = -1
!
      end  subroutine                                            
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
      subroutine seton_reg_records(nrg)
!
      integer:: nrg
!
      nvdes(nrg) = 0
      nvdes_tag(nrg) = 99
!
      end subroutine
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
      subroutine setoff_reg_records(nrg)
!
      integer:: nrg
!
      nvdes_tag(nrg) = -1
!
      end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      subroutine regs(nrg,varname,varlongname,varunit,n,k)
!
      character(len=*):: varname,varlongname,varunit
      character(len=1):: str
      integer:: nrg,n,kk,k,nlen1,nlen2,nlen3
      integer:: len_trim

      if(nvdes_tag(nrg).ne.-1) then
          nvdes(nrg) = nvdes(nrg)+1
       endif

      if(vartag(nrg).eq.-1) return
      kk = nvdes(nrg)
      nlen1=len_trim(varname);
      nlen2=len_trim(varlongname);
      nlen3=len_trim(varunit);

      if(nlen1.gt.20.or.nlen2.gt.70.or.nlen3.gt.20) then
         write(*,*) 'recorded string is too long!' 
         stop
      endif   
      gvarunit(nrg,kk) = varunit(1:nlen3)

      if(n.gt.1) then
        call sfromi(str,1,k)
        gvarname(nrg,kk) = varname(1:nlen1)//str
        gvarlongname(nrg,kk) = varlongname(1:nlen2)//str
      else
        gvarname(nrg,kk) = varname
        gvarlongname(nrg,kk) = varlongname
      endif
!
      end subroutine
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
      subroutine timestring(ibgnyr,ibgnmn,ibgndy,fbgnhr,s)
!
      integer:: ibgnyr,ibgnmn,ibgndy
      real:: fbgnhr
      integer:: ierr
! * strings for creating time string in *ctl file.
      character(len=*):: s
      character(len=4):: sbgnyr
      character(len=3),dimension(12):: smonth
      character(len=2):: sbgndy,sbgnhr
! * months in words from digital   
      data smonth/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/
    
      call sfromi(sbgnyr,4,ibgnyr)
      call sfromi(sbgndy,2,ibgndy)
      call sfromi(sbgnhr,2,int(fbgnhr))
      s=sbgnhr//'Z'//sbgndy//smonth(ibgnmn)//sbgnyr

      end subroutine
END MODULE
