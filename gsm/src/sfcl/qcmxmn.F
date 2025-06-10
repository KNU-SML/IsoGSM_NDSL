      subroutine qcmxmn(fld,ijdim,loismsk,vmaxmin,svar)
c
      implicit none
c
c  loismsk=0  ... ocean without sea ice
c          1  ... land without snow
c          2  ... sea ice without snow
c          3  ... land with snow
c          4  ... sea ice with snow
c
      integer ijdim
      real    fld(ijdim)
      integer loismsk(ijdim)
      character*8 svar
c
c  first index for max and min
c  second index for loismsk
c
      real vmaxmin(2,5)
c
      integer kount,ij,n
      real perr
c
c  if criteria max .lt. min, do not check
c
      if(vmaxmin(1,1).lt.vmaxmin(2,1)) return
c
c  check against land-sea mask and ice cover mask
c
      kount=0
      do ij=1,ijdim
        do n=1,5
           if(loismsk(ij).eq.n-1) then
             if(fld(ij).gt.vmaxmin(1,n)) then
               fld(ij)=vmaxmin(1,n)
               kount=kount+1
             elseif(fld(ij).lt.vmaxmin(2,n)) then
               fld(ij)=vmaxmin(2,n)
               kount=kount+1
             endif
           endif
        enddo
      enddo
c 
c  print results
c
#ifdef DBG
      perr=float(kount)/float(ijdim)*100.
      write(6,100) svar,perr
  100 format(a10,":% qc points=",f4.1)
#endif
c
      return
      end
