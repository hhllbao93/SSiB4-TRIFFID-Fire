      program vegmap
      real ssibveg1(360,180),ssibvegx(360,180) 
       
      open(66,file='vegmap_180x360.asc',form='formatted')
      open(77,file='vegmap_73x37.asc',form='formatted')
      rewind (66)
      read (66,888) ssibveg1
      close (66)
888   format(16f4.0)

      do i=1,360
      do j=1,180
          ssibvegx(i,j)= ssibveg1(i,180-j+1)
      enddo
      enddo

      write(77,888) ((ssibvegx(i,j),i=1,360,5),j=1,180,5)
      close(77)
      stop
      end

