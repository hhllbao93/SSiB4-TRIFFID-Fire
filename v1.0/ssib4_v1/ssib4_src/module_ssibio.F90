MODULE  module_ssibio
       use module_records
       use module_time
!
       implicit none
!
contains
     subroutine ssibin(ssibv,irestart,outdir,dtt,    &
             iyear,imonth,iday,tdayhr,itime,fhour)
!
      type(ssibdataios),dimension(lonf2s:lonf2e,latg2s:latg2e) :: ssibv
      character(len=*),intent(in):: outdir
      real,intent(in)            :: dtt
      integer,intent(out)        :: iyear,imonth,iday   
      real,intent(out)           :: tdayhr
      real,intent(out)           :: fhour
      integer:: ibgnyear,ibgnmonth,ibgnday
      real   :: bgndayhr
!      
      integer:: nfsibi,itime,irestart,ierr
      character(len=200):: strfile
!
! irestart = 0 for starting the model.
! irestart = 1 for starting the model and read ssib.inp.
! irestart = 2 for restart
!
!     ***************************************************
!     initialize ssib variables
!     ***************************************************
!
       itime = 1
       fhour = 0.0
!
       if(irestart.eq.0) return
!trying to read the file  
       nfsibi = 86   ! for restart       
       fhour = 0.0   ! it is important to set fhour = 0 here 
       strfile = outdir(1:len_trim(outdir))//'/ssib.inp'
       open(unit=nfsibi,file=strfile, form='unformatted',status='old',err=600)
!
      rewind nfsibi
      read(nfsibi) ibgnyear,ibgnmonth,ibgnday,bgndayhr,fhour
      write(6,*) ' file ssib.inp opened. unit=',nfsibi
      write(6,*) ' ibgnyr,ibgnmn,ibgndy,bgndayhr,fhour=', &
              ibgnyear,ibgnmonth,ibgnday,bgndayhr,fhour
      go to 601
  600 close(nfsibi)
      write(6,*) ' file ssib.inp not existed, started from beginning.'
      return  
  601 continue
!
      itime = 2
!
      if(irestart.eq.1) then 
         fhour = 0    ! in this case time is the same as initialized
      else
         iyear  = ibgnyear
         imonth = ibgnmonth
         iday   = ibgnday
         tdayhr = bgndayhr 
         FHOUR =FHOUR + DTT/3600.   ! go ahead one time step
         call timeforward(fhour,iyear,imonth,iday,tdayhr,ierr)
      endif
!
      read(nfsibi) ssibv
      close(nfsibi)

!  *************************************************************
    END SUBROUTINE
!---------------------------------------------------------------------------------------------
    subroutine ssibout(ssibv,outdir,ibgnyear,ibgnmonth,ibgnday,bgndayhr,fhour)
      type(ssibdataios),dimension(lonf2s:lonf2e,latg2s:latg2e) :: ssibv
      character(len=*):: outdir
      integer,intent(in):: ibgnyear,ibgnmonth,ibgnday
      real,intent(in)   :: bgndayhr,fhour
       integer::nfsibo
       character(len=200) strfile
!
       nfsibo = 87
       strfile = outdir(1:len_trim(outdir))//'/ssib.inp'
       open(unit=nfsibo,file=strfile,form='unformatted',err=900)
        go to 901
  900   continue
        write(*,*) ' error in opening file sib.inp'
        call abort
  901   continue

#ifdef DEBUG_PRINT
        write(*,*) ' file ssib.inp opened. unit=',nfsibo
        print 99, ibgnyear,ibgnmonth,ibgnday,bgndayhr,fhour
99    format(1h ,'output ssib ibgnyr,ibgnmn,ibgndy,bgndayhr,fhour=', &
                4(1x,f4.1),2x,f10.1)
#endif
        rewind nfsibo
        write(nfsibo) ibgnyear,ibgnmonth,ibgnday,bgndayhr,fhour
        write(nfsibo) ssibv
        close(nfsibo)
!
!  *************************************************************
       end subroutine
END MODULE

