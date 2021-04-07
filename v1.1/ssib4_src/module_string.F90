!
!This module is designed mainly for Fortran string processing.
!  zzq, sept 28,2013
!
MODULE module_string
    implicit none
contains
!=========================================================================
!
      function len_trim1(string) result(f_result)
!
!=========================================================================
      character(len=*) string
      integer k,f_result

      f_result=0
      do k=len(string),1,-1
         f_result=k
         if(isalphabet(string(k:k)).eq.1) goto 40
      enddo
   40    continue

      end function 
!=========================================================================
!
      function len_digit(x) result(f_result)
!
!=========================================================================
      integer f_result,k,x,y

      f_result=1
      y=x
      do k=1, 100
         y=y/10
         f_result=k
         if(y.lt.1) goto 40
      enddo
   40    continue

      end function 
!=========================================================================
!
      function strchr(string,c) result(f_result)
!
!=========================================================================
      character(len=*):: string
      character(len=1):: c
      integer::k,f_result

      f_result=-1 
      do k=1,len(string)
          if(c.eq.string(k:k)) then
             f_result = k
             return
          endif
      enddo
!
      end function 
!=========================================================================
!
      function isdigit(c) result(f_result)
!
!=========================================================================
      character(len=1):: c
      integer::f_result

      if(ichar(c).ge.48.and.ichar(c).le.57) then
          f_result=1
      else
          f_result =0
      endif

      end function 
!=========================================================================
!
      function isalphabet(c) result(f_result)
!
!=========================================================================
      character(len=1) c
      integer n,f_result

      n=ichar(c)
      if(n.ge.33.and.n.le.127) then
          f_result=1
      else
          f_result =0
      endif

      end function 
!=========================================================================
!
      function atoi(string) result(f_result)
! 
!=========================================================================
      character(len=*):: string
      character(len=1):: c
      integer:: i,n,f_result
!
      f_result=0
      do i=1,len(string)
         c= string(i:i)
         n=ichar(c)-48
         if(isdigit(c).eq.1) then
             f_result = f_result*10+n
         endif
      enddo
!     
      end function 
!=========================================================================
!
      function atof(string) result(f_result)
! 
!=========================================================================

      character(len=*):: string
      character(len=255),dimension(2):: strs
      character(len=1):: c
      integer:: n,negative,ks,i
      real:: f_result
!
        negative = -1
        ks =1 
        n=1
!      
      do i=1,len(string)
        c= string(i:i)
       if(c.eq.char(45)) then
           negative =i
           ks =i+1
       else if(c.eq.char(46)) then
           n = 2
           ks =i+1
       else if(isdigit(c).eq.1) then
           strs(n) = string(ks:i)
        endif
      enddo
!
      if(n.eq.1) then
          f_result = atoi(strs(1))
      else
          call compact(strs(2),n)
          f_result=atoi(strs(1))+atoi(strs(2))/10**n
      endif
!
      if(negative.ne.-1) f_result = -f_result
!    
      end function 
!=========================================================================
!
      subroutine sformat(string,sfrmt,x)
!
!=========================================================================
      character(len=*):: sfrmt,string  
      real,intent(in)::x
!
!sp controls whether the plus sign is included for positive numbers, 
! i.e.,spf12.6 and spi4 for float and integer.
!
      write(string,sfrmt) x
      string ='('//string(1:len_trim(string))//'/'
!
      end subroutine
!
!=========================================================================
!
      subroutine sfromi(s,num,x)
!
!=========================================================================
      character(len=*) s
      integer x,y,num,a,i

      s=char(0)
      y=x
      do i=1,num
        a=mod(y,10)
        y=y/10
        if(a.ge.0.and.a.le.9)  s=char(48+a)//s
      enddo
      end subroutine
!
!=========================================================================
!
      subroutine compact(string,lens)
!
!=========================================================================
      character(len=*) string
      character(len=255) str
      character c
      integer i,k,flag,lens
!    
      flag  =0 
      k=0

      do i=1,len(string)
          c=string(i:i)
          if(IsAlphabet(c).eq.1) then
             if(flag.eq.0) then
               str=c
               flag =1
           else
               str=str(1:k)//c
           endif
           k = k+1
         endif
        enddo
        lens = k
        string = str
!	
      end subroutine
!
!=========================================================================
!
        subroutine uppercase(string) 
!
!=========================================================================
        character(len=*) string
        integer k,n
!
         do k=1,len(string) 
           n=ichar(string(k:k))
          if(n.ge.97.and.n.le.122) then
            string(k:k) = char(n-32) 
           endif
         enddo
!
      end subroutine
!
!=========================================================================
!
        subroutine lowercase(string) 
!
!=========================================================================
        character(len=*) string
        integer k,n
!
        do k=1,len(string) 
          n=ichar(string(k:k))
          if(n.ge.65.and.n.le.90) then
          string(k:k) = char(n+32) 
          endif
        enddo
!
      end subroutine
!=========================================================================
!
        subroutine breaks(string,cmarker,strs,n)
!
!=========================================================================
        character(len=*) string
        character(len=*) strs(*)
        character c,cmarker
        integer ismarker,i1,ks,k,n

        do k=1,len(string)
         c=string(k:k)
         i1=k
         if(c.ne.cmarker) goto 40
         enddo
   40    continue
!	
        ismarker =0
        n=1
        ks=i1
!
        do k=i1,len(string)
         c=string(k:k)
         if(c.ne.cmarker) then
           if(ismarker.eq.1) then
                n = n+1
               ks = k
             endif
             strs(n)=string(ks:k)
             ismarker = 0
         else
             ismarker = 1
         endif    
      enddo
!
      end subroutine
!=========================================================================
!
      subroutine strtimestep(s,nstepsec)
!
!=========================================================================
      character(len=*) s
      character(len=10) s1
      integer nstepsec,n
      integer nmn,nhr,ndy,nmo,ny


      nmn = nstepsec/60
      n = len_digit(nmn)
      call sfromi(s1,n,nmn)
      s = s1(1:len_trim(s1))//'mn'

      if(nmn.gt.60) then
        nhr = nmn/60
        n = len_digit(nhr)
        call sfromi(s1,n,nhr)
        s = s1(1:len_trim(s1))//'hr'
      endif
     
      if(nhr.gt.24) then
        ndy = nstepsec/86400
        n = len_digit(ndy)
        call sfromi(s1,n,ndy)
        s = s1(1:len_trim(s1))//'dy'
      endif

      end subroutine

END MODULE
