!=======================================================================
!
      SUBROUTINE vegcnvt(vg,xxlat)
!                                                          
!=======================================================================
     
         veg1=vg
         if(vg.eq.1) veg1=1
         if(vg.eq.2)then
                 if(abs(xxlat).gt.45) then
                      veg1=2
                 else
                      veg1=1
                 endif
         endif
         if(vg.eq.3)then
                 if(abs(xxlat).gt.30) then
                      veg1=2
                 else
                      veg1=1
                 endif
         endif
         if(vg.eq.4) veg1=2
         if(vg.eq.5) veg1=2
         if(vg.eq.6) then
              if(abs(xxlat).gt.30) then
                     veg1=5    !c3
              else
                     veg1=4    !c4
              endif
          endif 
          if(vg.eq.7) veg1=3   ! C3
          if(vg.eq.8) veg1=4   ! C4
          if(vg.eq.9) veg1=5
          if(vg.eq.10) veg1=6
          if(vg.eq.11) veg1=7
          if(vg.eq.12) veg1=8
          if(vg.eq.13) veg1=9
!
          vg=veg1

       return
       END
