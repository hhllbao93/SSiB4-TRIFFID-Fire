!Reading data for 1 degree.
MODULE module_read1deg_clim
       use module_comsdef
       use ssib4_module_comsconst
       use ssib4_module_vegcnvt
       use ssib4_module_comssib4sfc
       use module_netcdf
       use module_string
!
      implicit none
      integer ncid_dswrf
      integer ncid_dlwrf
      integer ncid_pres
      integer ncid_tas
      integer ncid_prcp
      integer ncid_shum
      integer ncid_wind
      integer ncid_light
      integer ncid_dp
      integer ncid_agr
      integer ncid_gdp
      integer ncid_peatf
!Huilin add ncid_light ncid_Dp ncid_agr ncid_gdp ncid_peatf
contains
!
      subroutine readforcing(dtt,iyear1,krec,ks,ke,dswrf,dlwrf,  &
           pres,pmr,ppl,ppc,tm,qm,umm,vmm,zwind,                 &
           light,dp,agr,gdp,peatf)
!
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)::  &
          dswrf,dlwrf,pres,pmr,ppl,ppc,tm,qm,umm,vmm,&
          zwind,light,dp,agr,gdp,peatf
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)::  & 
          prcp,wind,tprec,sig,delsig
      real,dimension(360,180)::                      &
          dswrf1,dlwrf1,pres1,prcp1,tm1,qm1,wind1,   &
          light1,dp1,agr1,gdp1,peatf1
!Huilin add ncid_light ncid_Dp ncid_agr ncid_gdp ncid_peatf
!
      integer::iyear,iyear1
      real::dtt
      integer:: krec,ks,ke,indirlen,indirfirelen
      integer:: I,J
      character(len=100):: indir
      character(len=100):: indirfire
      character(len=4):: syear
      character(len=100):: file_dswrf
      character(len=100):: file_dlwrf
      character(len=100):: file_pres
      character(len=100):: file_tas
      character(len=100):: file_prcp
      character(len=100):: file_shum
      character(len=100):: file_wind
      character(len=100):: file_light
      character(len=100):: file_dp
      character(len=100):: file_agr
      character(len=100):: file_gdp
      character(len=100):: file_peatf
!       
      iyear = 1948
!=========Reading netcdf data ===================
!
!     READ IN THE ATMOSPHERIC data
!
      if(krec.eq.ks) then

      call sfromi(syear,4,iyear)
!
!      indir='/home7/xue/jzhang/ncdata_4806clim/'
      indir='/gpfs/fs1/work/hlhuang/data/ncdata_4806clim/'
!!!!      indir='/glade/fs1/work/zhangzq/ncdata_4806clim/'
!     indir='/gpfs/fs1/scratch/yeliu/data/ncdata/clim/'
      indirfire='/glade/p/univ/ucla0014/hlhuang/ncdata_4808/'
      indirlen=len_trim(indir)
      file_dswrf= indir(1:indirlen)//"dswrf_3hourly_4872clim.nc"
      file_dlwrf= indir(1:indirlen)//"dlwrf_3hourly_4872clim.nc"
      file_pres = indir(1:indirlen)//"pres_3hourly_4872clim.nc"
      file_tas  = indir(1:indirlen)//"tas_3hourly_4872clim.nc"
      file_prcp = indir(1:indirlen)//"prcp_3hourly_4872clim.nc"
      file_shum = indir(1:indirlen)//"shum_3hourly_4872clim.nc"
      file_wind = indir(1:indirlen)//"wind_3hourly_4872clim.nc"
      indirfirelen=len_trim(indirfire)
!     file_light= indir(1:indirlen)//"light_3hourly_4806clim.nc"
!     file_dp   = indir(1:indirlen)//"dp_3hourly_4806clim.nc"
!     file_agr  = indir(1:indirlen)//"agr_3hourly_4806clim.nc"
!     file_gdp  = indir(1:indirlen)//"gdp_3hourly_4806clim.nc"
!     file_peatf= indir(1:indirlen)//"peatf_3hourly_4806clim.nc"
      file_light= indirfire(1:indirfirelen)//"light_3hourly_1948-1948.nc"
      file_dp   = indirfire(1:indirfirelen)//"dp_3hourly_1948-1948.nc"
      file_agr  = indirfire(1:indirfirelen)//"agr_3hourly_1948-1948.nc"
      file_gdp  = indirfire(1:indirfirelen)//"gdp_3hourly_1948-1948.nc"
      file_peatf= indirfire(1:indirfirelen)//"peatf_3hourly_1948-1948.nc"
!Huilin add light Dp agr gdp peatf

      call netcdf_open_rf(file_dswrf,NCID_DSWRF)
      call netcdf_open_rf(file_dlwrf, NCID_DLWRF)
      call netcdf_open_rf(file_pres, NCID_PRES)
      call netcdf_open_rf(file_tas, NCID_TAS)
      call netcdf_open_rf(file_prcp, NCID_PRCP)
      call netcdf_open_rf(file_shum,NCID_SHUM)
      call netcdf_open_rf(file_wind, NCID_WIND)
      call netcdf_open_rf(file_light,NCID_LIGHT)
      call netcdf_open_rf(file_dp, NCID_DP)
      call netcdf_open_rf(file_agr,NCID_AGR)
      call netcdf_open_rf(file_gdp, NCID_GDP)
      call netcdf_open_rf(file_peatf, NCID_PEATF)
!Huilin add light Dp agr gdp peatf

       endif

      call netcdf_read0(1,360,1,180,krec, ncid_dswrf,"dswrf",dswrf1)
      call netcdf_read0(1,360,1,180,krec, ncid_dlwrf,"dlwrf",dlwrf1)
      call netcdf_read0(1,360,1,180,krec, ncid_pres,"pres",pres1)
      call netcdf_read0(1,360,1,180,krec, ncid_tas,"tas",tm1)
      call netcdf_read0(1,360,1,180,krec, ncid_prcp,"prcp",prcp1)
      call netcdf_read0(1,360,1,180,krec, ncid_shum,"shum",qm1)
      call netcdf_read0(1,360,1,180,krec, ncid_wind,"wind",wind1)
      call netcdf_read1(1,360,1,180,krec, ncid_light,"light",light1)
      call netcdf_read1(1,360,1,180,krec, ncid_dp,"dp",dp1)
      call netcdf_read1(1,360,1,180,krec, ncid_agr,"agr",agr1)
      call netcdf_read1(1,360,1,180,krec, ncid_gdp,"gdp",gdp1)
      call netcdf_read1(1,360,1,180,krec, ncid_peatf,"peatf",peatf1)
!Huilin add light Dp agr gdp peatf

      do i=lonf2s,lonf2e
      do j=latg2s,latg2e
        dswrf(i,j) = dswrf1(i,j)
        dlwrf(i,j) = dlwrf1(i,j)
        pres(i,j)  = pres1(i,j)
        tm(i,j)    = tm1(i,j)
        prcp(i,j)  = prcp1(i,j)
        qm(i,j)    = qm1(i,j)
        wind(i,j)  = wind1(i,j)
        light(i,j) = light1(i,j)
        dp(i,j)    = dp1(i,j)
        agr(i,j)   = agr1(i,j)
        gdp(i,j)   = gdp1(i,j)
        peatf(i,j) = peatf1(i,j)
      enddo
      enddo
!
      do i=lonf2s,lonf2e
      do j=latg2s,latg2e
!=========================================================================
! assume all precipitation is large scale.
! the input precipitation is in m/s.
!=========================================================================
         TPREC(i,j) = prcp(i,j)*DTT
         PPL(i,j) = TPREC(i,j)
         PPC(i,j) = TPREC(i,j)-PPL(i,j)

!=========================================================================
! reduced the wind from 10m to 40m
!=========================================================================
         umm(i,j)  = wind(i,j)+1.3   !reduced the wind from 10m to 40m
         vmm(i,j) = 0.0
!
         PMR(i,j)=pres(i,j)
         sig(i,j) = (EXP(AKAPPA*LOG(PMR(i,j)/pres(i,j))))**(1/AKAPPA)
         delsig(i,j) = (pres(i,j)-pmr(i,j))/pres(i,j)
!zzq         ZWIND(i,j)=  GASR/GRAV*DELSIG(i,j)*0.5*TM(i,j)

! assume reference height is 40m.
         ZWIND(i,j)= 40.0
      enddo
      enddo
!
      if(krec.eq.ke) then
         call netcdf_close(NCID_DSWRF)
         call netcdf_close(NCID_DLWRF)
         call netcdf_close(NCID_PRES)
         call netcdf_close(NCID_TAS)
         call netcdf_close(NCID_PRCP)
         call netcdf_close(NCID_SHUM)
         call netcdf_close(NCID_WIND)
         call netcdf_close(NCID_LIGHT)
         call netcdf_close(NCID_DP)
         call netcdf_close(NCID_AGR)
         call netcdf_close(NCID_GDP)
         call netcdf_close(NCID_PEATF)
      endif
!
      END SUBROUTINE
!
!----------------------------------------------------------------------
!
      subroutine param_ini(ssibv,dtt,   &
         ibgnyr,ibgnmn,ibgndy,fbgnhr,xlon,xlat)
!
!----------------------------------------------------------------------
      type(ssibdataios),dimension(lonf2s:lonf2e,latg2s:latg2e):: ssibv
      real,dimension(lonf2s:lonf2e):: xlon
      real,dimension(latg2s:latg2e):: xlat
      real,dimension(lonf2s:lonf2e,latg2s:latg2e)::  &
        swdown,rlwdown,pres,pmr,ppl,ppc,tm,qm,umm,vmm, &
        zwind,light,dp,agr,gdp,peatf 

      real,dimension(360,180):: ssibveg1,ssibvegx
      real,dimension(360,180):: fracveg1,fracvegx
      real,dimension(lonf2s:lonf2e,latg2s:latg2e,NTYPES):: fracveg
      real,dimension(360)::xlon1
      real,dimension(180)::xlat1
      integer :: ibgnyr,ibgnmn,ibgndy
      real    :: fbgnhr,dtt
      integer :: i,j,k,il,ntyp,it
      character(len=2):: ss
!
! starting date: January 1, 1948
      DTT   = 3600*3
      ibgnyr = 1948
      ibgnmn = 1
      ibgndy = 1
      fbgnhr = 0.0

!yliu      open(66,file='data_input/vegmap_180x360.asc',form='formatted')
!yliu      rewind (66)
!yliu      read (66,888) ssibveg1
!yliu      close (66)
!yliu  888 format(16f4.0)

      open(66,file='data_input/newmap/dom_type.asc',form='formatted')
      rewind (66)
      read (66,*) ssibveg1
      close (66)

      do k=1,NTYPES
         if ( k .lt. 10 ) write(ss,"(I1)") k
         if ( k .ge. 10 ) write(ss,"(I2)") k
         open(70,file='data_input/newmap/frac'//trim(ss)//'.asc', &
         form='formatted')
         read (70,*) fracveg1
         close (70)
         do i=1,360
         do j=1,180
             fracvegx(i,j)= fracveg1(i,180-j+1)
         enddo
         enddo
         do i=lonf2s,lonf2e
         do j=latg2s,latg2e
             fracveg(i,j,k)= fracvegx(i,j)
         enddo
         enddo
      enddo

      do k=1,NTYPES
         do i=lonf2s,lonf2e
         do j=latg2s,latg2e
            ssibv(i,j) % xveg_fracn(k) = fracveg(i,j,k)
         enddo
         enddo
      enddo
      do i=lonf2s,lonf2e
      do j=latg2s,latg2e
         ssibv(i,j) %xveg_fracn(8) = fracveg(i,j,8) + fracveg(i,j,10)
      enddo
      enddo

      do i=1,360
      do j=1,180
          xlon1(i)=i-1
          xlat1(j)=-89.5+j-1
          ssibvegx(i,j)= ssibveg1(i,180-j+1)
      enddo
      enddo
!
      do i=lonf2s,lonf2e
      do j=latg2s,latg2e
          xlon(i)=xlon1(i)
          xlat(j)= xlat1(j)
          ssibv(i,j) % xvtype= ssibvegx(i,j)
!yliu          call vegcnvt(ssibv(i,j)%xvtype, xlat(j))
      enddo
      enddo
!
! the following will be overwritten if running is the restart.
!
      call readforcing(dtt,ibgnyr,1,1,1,swdown,rlwdown,   &
           pres,pmr,ppl,ppc,tm,qm,umm,vmm,zwind,          &
           light,dp,agr,gdp,peatf)

      do j=latg2s,latg2e
      do i=lonf2s,lonf2e
        if(ssibv(i,j) % xvtype.gt.0.0) then
            ntyp=int(ssibv(i,j) %xvtype+0.5)
!yliu09May2017          if(ntyp.le.7) then
          if(ntyp.le.8) then
            ssibv(i,j) % xschsel = 1.0
          else
            ssibv(i,j) % xschsel = 2.0
          endif
!
!zzq 2009Feb23
!yliu09May2017          if(ntyp.eq.8) ssibv(i,j) % xschsel = 1.0
          if(ntyp.eq.9) ssibv(i,j) % xschsel = 1.0

          do it=1,NTYPES
!
            ssibv(i,j) % xtA(it)=tm(i,j)
            ssibv(i,j) % tc0(it)=tm(i,j)
            ssibv(i,j) % tg0(it)=tm(i,j)
            ssibv(i,j) % td0(it)=tm(i,j)
            ssibv(i,j) % tcm(it)=tm(i,j)
            ssibv(i,j) % tgm(it)=tm(i,j)
            ssibv(i,j) % tdm(it)=tm(i,j)
!
            ssibv(i,j) % w01(it) = 0.5
            ssibv(i,j) % w02(it) = 0.5
            ssibv(i,j) % w03(it) = 0.5
            ssibv(i,j) % wm1(it) = ssibv(i,j) % w01(it)
            ssibv(i,j) % wm2(it) = ssibv(i,j) % w02(it)
            ssibv(i,j) % wm3(it) = ssibv(i,j) % w03(it)
!
            ssibv(i,j) % cap01(it)=0.0
            ssibv(i,j) % cap02(it)=0.0
            ssibv(i,j) % capm1(it)=0.0
            ssibv(i,j) % capm2(it)=0.0
!
           enddo
         end if
        enddo
        enddo
!  *************************************************************
       print*, 'ssibveg is read!'
       end subroutine
END MODULE
