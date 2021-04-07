      program vegmap
      real ssibveg1(360,180), ssibvegx(360,180),ssiblai(360,180,12) 
      real zlt1(12),zlt2(12),zlt3(12),zlt4(12),zlt5(12),zlt6(12)
      real zlt7(12),zlt8(12),zlt9(12),zlt10(12),zlt11(12),zlt12(12)
      real zlt(9,12),vc(9)

data zlt1/4.53781,4.53781,4.53781,4.53781,4.53781,4.53781,4.53781,4.53781,4.53781,4.53781,4.53781,4.53781/
data zlt2/5.99999,6.30000,6.80000,6.99998,7.19999,7.29996,7.10000,7.00003,6.79998,6.30002,5.99998,5.99999/
data zlt3/0.40000,0.50000,0.60000,0.70000,1.29999,3.00000,3.50001,1.49999,0.80000,0.60000,0.50000,0.40000/
data zlt4/1.11595,1.01520,0.73292,0.79373,0.77415,1.02609,2.27048,4.14710,3.48744,2.34535,1.42727,1.47393/
data zlt5/0.46117,0.30744,0.16232,0.16364,0.14659,0.14659,0.14659,0.70139,0.98104,0.59175,0.47143,0.46117/
data zlt6/0.17040,0.17040,0.17040,0.17040,0.17040,0.31567,1.59873,1.35057,0.11181,0.17040,0.17040,0.17040/
data zlt7/0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000/
data zlt8/0.01000,0.01000,0.03000,0.12000,0.30000,0.47000,0.45000,0.34000,0.13000,0.03000,0.01000,0.01000/
data zlt9/0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000/

      vc(1)=0.98
      vc(2)=0.75
      vc(3)=0.9
      vc(4)=0.3
      vc(5)=0.1
      vc(6)=0.3
      vc(7)=0.0001
      vc(8)=0.075
      vc(9)=0.0001

      do im=1,12
        zlt(1,im)=zlt1(im)
        zlt(2,im)=zlt2(im)
        zlt(3,im)=zlt3(im)
        zlt(4,im)=zlt4(im)
        zlt(5,im)=zlt5(im)
        zlt(6,im)=zlt6(im)
        zlt(7,im)=zlt7(im)
        zlt(8,im)=zlt8(im)
        zlt(9,im)=zlt9(im)
       enddo

      open(66,file='vegmap_180x360.asc',form='formatted')
      open(77,file='gradslai.gra',  &
          access='direct',recl=180*360*4,  &
          form='unformatted',status='unknown')

      rewind (66)
      read (66,888) ssibveg1
      close (66)
888   format(16f4.0)

      do i=1,360
      do j=1,180
          ssibvegx(i,j)= ssibveg1(i,180-j+1)
          xlat=-89.5+j-1
          call vegcnvt(ssibvegx(i,j), xlat)
          ity=ssibvegx(i,j)
          do im=1,12
              if(ity.ne.0.and.ityp.ne.8) then
                ssiblai(i,j,im)=zlt(ity,im)
               else
                ssiblai(i,j,im)=99
               endif
           enddo
      enddo
      enddo

      do im=1,12
         WRITE(77,rec=im) ((ssiblai(i,j,im),i=1,360),j=1,180)
      enddo

      close(77)
      end program
