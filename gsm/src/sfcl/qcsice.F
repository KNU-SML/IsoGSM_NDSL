      subroutine qcsice(ais,glacir,amxice,slmask,idim,jdim)
c
      implicit none
c
      integer idim,jdim
      real ais(idim,jdim),glacir(idim,jdim),
     1     amxice(idim,jdim),slmask(idim,jdim)
c
      integer i,j,ij,kount,ip,im,jp,jm
      real    perr
c
c  check sea-ice cover mask against land-sea mask
c
      kount=1
      do j=1,jdim
      do i=1,idim
        if(slmask(i,j).eq.0..and.glacir(i,j).eq.1..and.
     1     ais(i,j).ne.1.) then
          ais(i,j)=1.
          kount=kount+1
        endif
        if(slmask(i,j).eq.1..and.ais(i,j).eq.1.) then
          ais(i,j)=0.
          kount=kount+1
        endif
        if(slmask(i,j).eq.0..and.amxice(i,j).eq.0..and.
     1     ais(i,j).eq.1.) then
          ais(i,j)=0.
          kount=kount+1
        endif
      enddo
      enddo
c
c  remove isolated open ocean surrounded by sea ice and/or land
c
c     ij=0
c     do j=1,jdim
c       do i=1,idim
c         ij=ij+1
c         ip=i+1
c         im=i-1
c         jp=j+1
c         jm=j-1
c         if(jp.gt.jdim) jp=jdim-1
c         if(jm.lt.1) jm=2
c         if(ip.gt.idim) ip=idim-1
c         if(im.lt.1) im=2
c         if(slmask(i,j).eq.0..and.ais(i,j).eq.0.) then
c           if((slmask(ip,jp).eq.1..or.ais(ip,jp).eq.1.).and.
c    1         (slmask(i ,jp).eq.1..or.ais(i ,jp).eq.1.).and.
c    2         (slmask(im,jp).eq.1..or.ais(im,jp).eq.1.).and.
c    3         (slmask(ip,j ).eq.1..or.ais(ip,j ).eq.1.).and.
c    4         (slmask(im,j ).eq.1..or.ais(im,j ).eq.1.).and.
c    5         (slmask(ip,jm).eq.1..or.ais(ip,jm).eq.1.).and.
c    6         (slmask(i ,jm).eq.1..or.ais(i ,jm).eq.1.).and.
c    7         (slmask(im,jm).eq.1..or.ais(im,jm).eq.1.)) then
c               ais(i,j)=1.
cifdef DBG
c             write(6,*) ' isolated open sea point surrounded by',
c    1                   ' sea ice or land modified to sea ice'
cendif
c           endif
c         endif
c       enddo
c     enddo
c
c  print results
c
#ifdef DBG
      perr=float(kount)/float(idim*jdim)*100.
      write(6,100) perr
#endif
  100 format("qcsice:% qc points=",f4.1)
      return
      end
