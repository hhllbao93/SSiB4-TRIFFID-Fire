!***************************************************************
! The programme is used to read and write a NetCDF file.
! Zhengqiu Zhang. Jun01,2009  
!***************************************************************
!
!***************************************************************
MODULE module_netcdf
      implicit none
#include <netcdf.inc>

contains
      subroutine netcdf_close(ncid)
      integer ncid
!      
      call hdlerr(nf_close(ncid))
!
      end subroutine
!***************************************************************
      subroutine netcdf_read0( &
            longmin,longmax,latmin,latmax,it,ncid,name,var)
      integer ncid,varid
      integer longmin,longmax,latmin,latmax,zmin,zmax,it
      real var(longmin:longmax,latmin:latmax)
      integer startB(3),countB(3)
      character(len=*) name
!
       zmin=1
       zmax=1
!
       startB(1)=1
       startB(2)=1
       startB(3)=it
       countB(1)=longmax-longmin+1
       countB(2)=latmax-latmin+1
       countB(3)=1
!
      call hdlerr(NF_INQ_VARID(ncid,name,varid))
      call hdlerr(NF_GET_VARA_REAL(ncid,varid,startB,countB,var))
!
      end subroutine
!***************************************************************
      subroutine netcdf_read1( &
            longmin,longmax,latmin,latmax,it,ncid,name,var)
      integer ncid,varid
      integer longmin,longmax,latmin,latmax,zmin,zmax,it
      real var(longmin:longmax,latmin:latmax)
      integer startB(4),countB(4)
      character(len=*) name
!
       zmin=1
       zmax=1
!
       startB(1)=1
       startB(2)=1
       startB(3)=1
       startB(4)=it
       countB(1)=longmax-longmin+1
       countB(2)=latmax-latmin+1
       countB(3)=zmax-zmin+1
       countB(4)=1
!
      call hdlerr(NF_INQ_VARID(ncid,name,varid))
      call hdlerr(NF_GET_VARA_REAL(ncid,varid,startB,countB,var))
!
      end subroutine
!***************************************************************
      subroutine netcdf_write1(longmin,longmax,latmin,latmax,it,ncid,name,var)
      integer ncid,varid
      integer longmin,longmax,latmin,latmax,zmin,zmax,it
      real,dimension(longmin:longmax,latmin:latmax) :: var
      integer,dimension(4):: startB,countB
      character(len=*) name
!
       zmin=1
       zmax=1
!
       startB(1)=1
       startB(2)=1
       startB(3)=1
       startB(4)=it
       countB(1)=longmax-longmin+1
       countB(2)=latmax-latmin+1
       countB(3)=zmax-zmin+1
       countB(4)=1
!
      call hdlerr(NF_INQ_VARID(ncid,name,varid))
      call hdlerr(NF_PUT_VARA_REAL(ncid,varid,startB,countB,var))
!
      end subroutine
!***************************************************************
      subroutine netcdf_write_timestep(ncid,nt,timestep) 
      integer ncid,timeid,nta,nt
      real time(nt)
      real timestep
      integer,dimension(1):: startA,countA
      integer it

      call hdlerr( NF_INQ_VARID(ncid,'time',timeid))
      call hdlerr( NF_INQ_DIMLEN(ncid,timeid,nta)) 

      startA(1)=1
      do it=1,nt
        countA(1)=it
        time(it)=(it-1)*timestep
      enddo

      call hdlerr(nf_put_vara_real(ncid,timeid,startA,countA,time))

      end subroutine
!***************************************************************
      subroutine trim_nf_put_att_text(ncid,varid,units,text)
      integer ncid,varid
      character(len=*) units
      character(len=*) text
      integer nlen,len_trim
!
      nlen=len_trim(text) 
      call hdlerr( nf_put_att_text(ncid, varid, &
                   units, nlen,text(1:nlen)))
      end subroutine
!***************************************************************
      subroutine netcdf_open_rf(filename, ncid) 
      character(len=*) filename
      integer ncid
      call hdlerr( NF_OPEN(filename, NF_NOCLOBBER, ncid))

      end subroutine
!***************************************************************
      subroutine netcdf_open_wf(                &
         filename,title,                        &
         longmin,longmax,latmin,latmax,         &
         zmin,zmax,num,                         &
         ncid,                                  &
         var_names, var_longnames,              &
         var_units,time_units,                  &
         value_mis,                             &
         longitude, latitude, altitude,irestart)
      integer longmin,longmax,latmin,latmax,zmin,zmax,num
      integer ncid
      integer, dimension(num):: varid

      integer::               &
        longitude_dim,longid, &
        latitude_dim,latid,   &
        altitude_dim,zid,     &
        time_dim,timeid
!
! variable shapes
      integer,dimension(1)::  &
        longitude_dims,       &
        latitude_dims,        &
        altitude_dims,        &
        time_dims
      integer,dimension(num) :: var_dims

      integer startA(1),countA(1),startB(4),countB(4)
      real longitude(longmin:longmax)
      real latitude(latmin:latmax)
      real altitude(zmin:zmax)

      character(len=*) filename,title
      character(len=*) time_units
      character(len=*) var_names(num)
      character(len=*) var_longnames(num)
      character(len=*) var_units(num)
      real value_mis(num)
      integer i,status
      integer longcnt,latcnt,zcnt
      integer irestart
!
!enter define mode
!zzq      status= nf_create(filename, NF_64BIT_OFFSET, ncid)
      if (irestart.eq.2) then
      status= nf_create(  &
!        filename, IOR(NF_NOCLOBBER,NF_64BIT_OFFSET), ncid)
        filename,IOR(NF_NOCLOBBER,NF_SHARE), ncid)
      else
      status= nf_create(  &
!        filename, IOR(NF_NOCLOBBER,NF_64BIT_OFFSET), ncid)
        filename,IOR(NF_CLOBBER,NF_SHARE), ncid)
      end if

!check NetCDF existence and open the file
      if(status.eq.NF_EEXIST) then
        call hdlerr( NF_OPEN(filename,IOR(NF_WRITE,NF_SHARE), ncid))
        do i=1,num
         call hdlerr( NF_INQ_VARID(ncid,var_names(i),varid(i)))
        enddo
        return
      endif

      longcnt=longmax-longmin+1
      latcnt=latmax-latmin+1
      zcnt=zmax-zmin+1
! define dimensions
      call hdlerr( nf_def_dim(ncid,'longitude',longcnt,longitude_dim))
      call hdlerr( nf_def_dim(ncid,'latitude',latcnt,latitude_dim))
      call hdlerr( nf_def_dim(ncid,'altitude',zcnt, altitude_dim))
      call hdlerr( nf_def_dim(ncid,'time',NF_UNLIMITED,time_dim))

! define variables
      longitude_dims(1) = longitude_dim
      latitude_dims(1)  = latitude_dim
      altitude_dims(1)  = altitude_dim
      time_dims(1)      = time_dim

      call hdlerr( nf_def_var(ncid, 'longitude',           &
                  NF_REAL, 1, longitude_dims, longid))
      call hdlerr( nf_def_var(ncid, 'latitude',            &
                  NF_REAL, 1, latitude_dims, latid))
      call hdlerr( nf_def_var(ncid,'altitude',             &
                  NF_REAL, 1, altitude_dims,zid))
      call hdlerr( nf_def_var(ncid,'time',                 &
                  NF_REAL,1,time_dim,timeid))

      var_dims(4) = time_dim
      var_dims(3) = altitude_dim
      var_dims(2) = latitude_dim
      var_dims(1) = longitude_dim

      do i=1,num
          call hdlerr(nf_def_var(ncid,                    &
              var_names(i),NF_REAL,4,var_dims,varid(i)))
      enddo

! assign attributes
      call trim_nf_put_att_text(ncid, longid,             &
                   'units', 'degrees_east')
      call trim_nf_put_att_text(ncid, longid,             &
                  'point_spacing', 'even')
      call trim_nf_put_att_text(ncid, latid,              &
                  'units', 'degrees_north')
      call trim_nf_put_att_text(ncid, latid,              &
                  'point_spacing', 'even')
      call trim_nf_put_att_text(ncid, timeid,             &
                  'units', time_units)
      call trim_nf_put_att_text(ncid, zid,                &
                  'units', 'level')

      do i=1,num
      call trim_nf_put_att_text(ncid, varid(i),           &
               'name', var_names(i))
      call trim_nf_put_att_text(ncid, varid(i),           &
               'long_name', var_longnames(i))
      call trim_nf_put_att_text(ncid, varid(i),           &
               'units', var_units(i))
      call hdlerr( nf_put_att_real(ncid, varid(i),        &
               'missing_value',NF_REAL, 1, value_mis(i)))
      enddo

      call trim_nf_put_att_text(ncid, NF_GLOBAL,'title', title)
      call trim_nf_put_att_text(ncid, NF_GLOBAL,          &
                'author', 'Zhengqiu Zhang')

! leave define mode
      call hdlerr( nf_enddef(ncid))
!
      startA(1)=1
      countA(1)=longcnt
      call hdlerr( nf_put_vara_real(ncid,longid,          &
              startA,countA,longitude))
!
      startA(1)=1
      countA(1)=latcnt
      call hdlerr( nf_put_vara_real(ncid, latid,          &
              startA,countA,latitude))
!
      startA(1)=1
      countA(1)=zcnt
      call hdlerr( nf_put_vara_real(ncid, zid,            &
              startA,countA,altitude))
!
      end subroutine
!***************************************************************
      subroutine hdlerr (status)
      integer status
!
!Add handling code here.
!
      IF (STATUS .NE. NF_NOERR) THEN
        PRINT *, NF_STRERROR(STATUS)
        STOP 'Stopped!!!!!!!!!'
      ENDIF
!
      end subroutine

END MODULE
