!
!This module is designed mainly for outputing model result in NetCdf.
!  zzq, sept 28,2013
!
MODULE module_writenc
      use ssib4_module_comssib4sfc
      use module_records
!
      implicit none
      integer:: monthly_ncid
      integer:: monthlyh_ncid(48)
!
contains
!=========================================================================
       subroutine netcdf_openfile(         &
          title,filename,                  &
          longitude,latitude,nrg,num,      &
          ibgnyr,ibgnmn,ibgndy,fbgnhr,ncid,irestart)
!
      integer:: nrg, num,ncid 
      integer:: ibgnyr,ibgnmn,ibgndy
      real   :: fbgnhr
      character(len=*) :: title,filename
      character(len=20):: var_names(num)
      character(len=80):: var_longnames(num)
      character(len=20):: var_units(num)
      character(len=30):: time_units
      character(len=80):: strname
      character(len=30):: strunit
      real,dimension(num):: value_mis
      real,dimension(LONF2S:LONF2E):: longitude
      real,dimension(LATG2S:LATG2E)::latitude
      real,dimension(1):: z
      integer:: k,len1,len2,len_trim
      integer:: irestart

!       print*,'outfile name:',filename
       z(1)=1
!
       do k=1,num
          value_mis(k) = -999
          var_names(k) = gvarname(nrg,k)
          strname = gvarlongname(nrg,k)
          strunit = gvarunit(nrg,k)
          len1=len_trim(strname)
          len2=len_trim(strunit)
          var_longnames(k) = strname(1:len1)//'('//strunit(1:len2)//')'
          var_units(k) = strunit
        enddo
!
        call  netcdf_time_units(                      &
           ibgnyr,ibgnmn,ibgndy,fbgnhr,time_units)
!        print*,"time unit:",time_units
!
        call netcdf_open_wf(                          &
         filename,title,                              &
         LONF2S,LONF2E,LATG2S,LATG2E,                 &
         1,1,num,                                     &
         ncid,                                        &
         var_names, var_longnames,                    &
         var_units,time_units,                        &
         value_mis,                                   &
         longitude, latitude, z, irestart)
!
       end subroutine
!=========================================================================
!
      subroutine netcdf_time_units(ibgnyr,ibgnmn,ibgndy,fbgnhr,s)
!
      integer:: ibgnyr,ibgnmn,ibgndy
      real   :: fbgnhr
      integer:: ierr
! * strings for creating time string in *ctl file.
      character(len=*):: s
      character(len=4):: sbgnyr
      character(len=2):: sbgnmn,sbgndy,sbgnhr
!
      call sfromi(sbgnyr,4,ibgnyr)
      call sfromi(sbgnmn,2,ibgnmn)
      call sfromi(sbgndy,2,ibgndy)
      call sfromi(sbgnhr,2,int(fbgnhr))
      s = "months since "//trim(sbgnyr)//"-"//trim(sbgnmn)//   &
          "-"//trim(sbgndy)//" "//trim(sbgnhr)//":00:00"
!
      end subroutine
!=======================================================================
       subroutine netcdf_inital_monthly(                       &
         expnum,outdir,nrg, ibgnyr,ibgnmn,ibgndy,fbgnhr,irestart)
!
      integer:: ibgnyr,ibgnmn,ibgndy
      real   :: fbgnhr
      character(len=100):: title,longfname,fname
      character(len=*):: outdir,expnum
      integer:: nrg,num,len_trim
      integer:: irestart
!
      if(vartag(nrg).eq.-1) return
!
       title="monthly averages"
       fname='/mon'//expnum(1:len_trim(expnum))//'.nc'
       longfname=outdir(1:len_trim(outdir))//fname(1:len_trim(fname))
       num = nvdes(nrg)
!
       call netcdf_openfile(title,longfname,   &
          xlongit,xlatit,nrg,num,              &
          ibgnyr,ibgnmn,ibgndy,fbgnhr,monthly_ncid,irestart)
!
       end subroutine
!=========================================================================
      subroutine netcdf_writeoutm(nrg,it)
!
      integer:: nrg,it,i,j,k
      real(kind=4),dimension(LONF2S:LONF2E,LATG2S:LATG2E):: r4vouts
      character(len=20):: varname
      integer:: len_trim
!
       do k=1,nvdes(nrg)
       do i=LONF2S,LONF2E
        do j=LATG2S,LATG2E
            if(vcount(i,j,k).ne.0.and.landmsk(i,j).eq.1) then
               r4vouts(i,j)=vouts(i,j,k)/vcount(i,j,k)
            else
               r4vouts(i,j)=-999
            endif
         enddo
       enddo
!
      varname = gvarname(nrg,k);
      call netcdf_write1(LONF2S,LONF2E,LATG2S,LATG2E,it,   &
         monthly_ncid,varname(1:len_trim(varname)),r4vouts)
      enddo
      call netcdf_write_timestep(monthly_ncid,it,1.0)
!
      end subroutine
!=========================================================================
      subroutine netcdf_writeoutm_close(nrg)
      integer:: nrg
!
      call netcdf_close(monthly_ncid)
!
      end  subroutine 
!=========================================================================
       subroutine netcdf_inital_monthlyh(              &
         expnum,outdir,nrg, dtt, ibgnyr,ibgnmn,ibgndy,fbgnhr,irestart)
!
      character(len=*):: outdir,expnum
      character(len=100):: title,longfname,fname
      character(len=2):: shr
      integer:: ibgnyr,ibgnmn,ibgndy
      real   :: fbgnhr
      integer:: nrg,num,ihr,mymaxhr,nhr,len_trim
      real:: dtt,hourstep
      integer:: irestart

      if(vartag(nrg).eq.-1) return
      title="monthly averages for steptimes"
      hourstep = dtt/3600.
      mymaxhr = int(24/hourstep)
!
      num = nvdes(nrg)
      do ihr=1,mymaxhr
        nhr = int((ihr-1)*hourstep)
        call sfromi(shr,2,nhr)            
        fname='/mon'//expnum(1:len_trim(expnum))//'_h'//shr//'.nc'
        longfname=outdir(1:len_trim(outdir))//fname(1:len_trim(fname))

        call netcdf_openfile(title,longfname,        &
          xlongit,xlatit,nrg,num,                    &
          ibgnyr,ibgnmn,ibgndy,fbgnhr,monthlyh_ncid(ihr),irestart)
      enddo
!
      end subroutine
!
!=========================================================================
      subroutine netcdf_writeouth(nrg,dtt,it)
!
      integer:: nrg,ihr, i,j,k,it,mymaxhr,nhr
      real,dimension(LONF2S:LONF2E,LATG2S:LATG2E):: r4vouts
      character(len=20):: varname;
      integer:: len_trim
      real:: dtt,hourstep
!
      hourstep = dtt/3600.
      mymaxhr = int(24/hourstep)
!
! ********* output set is like *********
!  00 03  06  09  11  15  17  21
!***************************************
!
      do 1000 ihr=1,mymaxhr
!
       nhr = int((ihr-1)*hourstep)
       do k=1,nvdes(nrg)
       do i=LONF2S,LONF2E
        do j=LATG2S,LATG2E
            if(counth(i,j,k,ihr).ne.0.and.landmsk(i,j).eq.1) then
               r4vouts(i,j)=voutsh(i,j,k,ihr)/counth(i,j,k,ihr)
            else
               r4vouts(i,j)=-999
            endif
         enddo
       enddo
!
      varname=gvarname(nrg,k); 
      call netcdf_write1(LONF2S,LONF2E,LATG2S,LATG2E,it,      &
        monthlyh_ncid(ihr),varname(1:len_trim(varname)),r4vouts)
      enddo
      call netcdf_write_timestep(monthlyh_ncid(ihr),it,1.0)
 1000  continue
!
      end subroutine
!=========================================================================
      subroutine netcdf_writeouth_close(nrg,dtt)
!
      integer:: nrg,ihr,mymaxhr
      real:: dtt
!
      mymaxhr = int(24*3600/dtt)
!
      do ihr=1,mymaxhr
          call netcdf_close(monthlyh_ncid(ihr))
      enddo
!
      end subroutine
END MODULE
