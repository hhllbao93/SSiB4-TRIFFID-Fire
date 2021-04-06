!
!This module is designed by zzq, mainly for model output in binary file.
!  zzq, sept 28,2013
!
MODULE module_writectl
      use module_records
contains
!
!====================================================
      subroutine writeouth(nrg,outdir,dtt,iyr,imn)
!
      character(len=*) :: outdir
      character(len=80):: outfile
      character(len=4) :: syear
      character(len=2) :: smon,shr
      integer:: iyr,imn,nhr,nrg
      integer:: nid,i,j,k,outdirlen,mymaxhr,latxlon
      real,dimension(LONF2S:LONF2E,LATG2S:LATG2E)::r4vouts
      real:: hourstep,dtt
!
      hourstep= dtt/3600.
      mymaxhr = int(24/hourstep)
!
! ********* output set is like *********
!  00 03  06  09  11  15  17  21
!***************************************
      nid=22
      call sfromi(syear,4,iyr)
      call sfromi(smon,2,imn)

      do 1000 ihr=1,mymaxhr
!
      nhr = int((ihr-1)*hourstep)
      call sfromi(shr,2,nhr)
      outdirlen=len_trim(outdir)
      outfile=outdir(1:outdirlen)//'/f'//syear//smon//'h'//shr
      latxlon= (lonf2e-lonf2s+1)*(latg2e-latg2s+1)
!
      open(nid,file=outfile,                &
          access='direct',recl=latxlon*4,   &
          form='unformatted',status='unknown')
!
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
      WRITE(nid,rec=k)((r4vouts(i,j),i=LONF2S,LONF2E), j=LATG2S,LATG2E)

      enddo
      close(nid)
1000  continue

      end subroutine
!====================================================
      subroutine writeout(nrg,outdir,iyr,imn)
!
      character(len=*) :: outdir
      character(len=80):: outfile
      character(len=4) :: syear
      character(len=2) :: smon
      integer:: iyr,imn
      integer:: nid,i,j,k,outdirlen
      real(kind=4),dimension(LONF2S:LONF2E,LATG2S:LATG2E)::r4vouts
!
      nid=22
      call sfromi(syear,4,iyr)
      call sfromi(smon,2,imn)
!
      outdirlen=len_trim(outdir)
      outfile=outdir(1:outdirlen)//'/f'//syear//smon

      open(nid,file=outfile,                                          & 
          access='direct',recl=(lonf2e-lonf2s+1)*(latg2e-latg2s+1)*4, & 
          form='unformatted',status='unknown')
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
      WRITE(nid,rec=k)((r4vouts(i,j),i=LONF2S,LONF2E), j=LATG2S,LATG2E)
      enddo
!
      close(nid)

      end subroutine
!
!=========================================================================
!
      subroutine write_script1deg(                            &
         nrg,dataf,nrec,ibgnyr,ibgnmn,ibgndy,fbgnhr,strstep,  &
         lons,lone,lats,late)
!
!=========================================================================
      character(len=*):: dataf,strstep
      character(len=80):: tstr
      character(len=80),dimesnion(3):: s
      integer:: lons,lone,lats,late
      integer:: NFID,nrec,nx,ny
      integer:: k,nrg
      integer:: ibgnyr,ibgnmn,ibgndy
      real:: fbgnhr

      nx = lone-lons +1
      ny = late-lats +1
 
      NFID=777
      open(NFID,file=dataf, status='unknown')

      write(NFID,'("dset &0")')
      write(NFID,'("undef -999")')
      write(NFID,5) stitle(nrg)
5     format("Title",2x,A<LEN_TRIM(stitle(nrg))>)
      write(NFID,'("options template")')
      write(NFID,10) nx,lons-1.0,1.0
      write(NFID,20) ny,-89.5+lats-1,1.0
10    format("xdef",1x,I4,1x,"linear",1x,f8.1,1x,f4.1)
20    format("ydef",1x,I4,1x,"linear",1x,f8.1,1x,f4.1)
      write(NFID,'("zdef 1 linear 1 1")')
      call timestring(ibgnyr,ibgnmn,ibgndy,fbgnhr,tstr)
      write(NFID,30) nrec,tstr,strstep
30    format("TDEF",2x,I<len_digit(nrec)>,2x,   &
           "LINEAR",2x,A<LEN_TRIM(tstr)>,2x,A<LEN_TRIM(strstep)>)
      write(NFID,40) nvdes(nrg)
40    format("VARS ", I<len_digit(nvdes(nrg))>)
      do k=1,nvdes(nrg)
       s(1)= gvarname(nrg,k)
       s(2)= gvarlongname(nrg,k)
       s(3)= gvarunit(nrg,k)
       write(NFID,50) s(1),k,s(2),s(3)
      enddo
50    format(A13,"0 99",2x,i<len_digit(k)>,     &
         2x,A<LEN_TRIM(s(2))>,'(',A<LEN_TRIM(s(3))>,')')
      write(NFID,'("ENDVARS")')
      close(NFID)

      end subroutine
!=========================================================================
!
      subroutine writescript(   &
          nrg,dataf,nrec,ibgnyr,ibgnmn,ibgndy,fbgnhr,strstep)
!
!=========================================================================
      character(len=*):: dataf,strstep
      character(len=80):: tstr
      character(len=80),dimension(3):: s
      integer:: NFID,nrec,n
      integer:: k,nrg
      integer:: nrec,ibgnyr,ibgnmn,ibgndy
      real   :: fbgnhr

      NFID=777
      open(NFID,file=dataf, status='unknown')
      write(NFID,'("dset &0")')
      write(NFID,'("undef -999")')
      write(NFID,5) stitle(nrg)
5     format("Title",2x,A<LEN_TRIM(stitle(nrg))>)
      write(NFID,'("options template")')
      write(NFID,'("xdef 1 linear 1 1")')
      write(NFID,'("ydef 1 linear 1 1")')
      write(NFID,'("zdef 1 linear 1 1")')
      call timestring(ibgnyr,ibgnmn,ibgndy,fbgnhr,tstr)
      write(NFID,30) nrec,tstr,strstep
30    format("TDEF",2x,I<len_digit(nrec)>,2x, &
           "LINEAR",2x,A<LEN_TRIM(tstr)>,2x,A<LEN_TRIM(strstep)>)
      write(NFID,40) nvdes(nrg)
40    format("VARS ", I<len_digit(nvdes(nrg))>)
      do k=1,nvdes(nrg)
       s(1)= gvarname(nrg,k)
       s(2)= gvarlongname(nrg,k)
       s(3)= gvarunit(nrg,k)
       write(NFID,50) s(1),k,s(2),s(3)
      enddo
50    format(A13,"0 99",2x,i<len_digit(k)>,   &
         2x,A<LEN_TRIM(s(2))>,'(',A<LEn_TRIM(s(3))>,')')
      write(NFID,'("ENDVARS")')
      close(NFID)

      end subroutine
!=========================================================================
!
      subroutine writestub  &
           (nrg,dataf,nrec,ibgnyr,ibgnmn,ibgndy,fbgnhr,strstep)
!
!=========================================================================
      character(len=*):: dataf,strstep
      character(len=80):: tstr
      character(len=80),dimension(3):: s
      integer:: NFID,nrec,nrg
      integer:: k
      integer:: nrec,ibgnyr,ibgndy
      real   :: fbgnhr

      NFID=777
      open(NFID,file=dataf, status='unknown')
      call timestring(ibgnyr,ibgnmn,ibgndy,fbgnhr,tstr)
      write(NFID,10) nrec,tstr,strstep
10    format("TDEF",2x,I<len_digit(nrec)>,2x,   & 
            "LINEAR",2x,A<LEN_TRIM(tstr)>,2x,A<LEN_TRIM(strstep)>)
      write(NFID,15) nvdes(nrg)
15    format("VARS ", I<len_digit(nvdes(nrg))>)
      do k=1,nvdes(nrg)
       s(1)= gvarname(nrg,k)
       s(2)= gvarlongname(nrg,k)
       s(3)= gvarunit(nrg,k)
       write(NFID,20) s(1),k,s(2),s(3)
      enddo
20    format(A13,"0 99",2x,i<len_digit(k)>,      &
        2x,A<LEN_TRIM(s(2))>,'(',A<LEN_TRIM(s(3))>,')')
      write(NFID,'("ENDVARS")')
      close(NFID)

      end subroutine
!====================================================
      subroutine duplicate(outdir,sfile1,sfile2,iyr,imn)
!
      character(len=*):: outdir,sfile1,sfile2
      integer:: iyr,imn
      character(len=120):: strcmd
      character(len=4):: syear
      character(len=2):: smon
!
      call sfromi(syear,4,iyr)
      call sfromi(smon,2,imn)
      strcmd=outdir(1:len_trim(outdir))          & 
        //'/'//sfile1(1:len_trim(sfile1))        &
        //' '//outdir(1:len_trim(outdir))        &
        //'/'//sfile2(1:len_trim(sfile2))        & 
        //syear//smon
      call system('cp -p '//strcmd);
      
      end subroutine
END MDOULE
