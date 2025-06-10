      subroutine nfdr_init(idim,jdim,
     1                    cumppt,avethbndy,ythmois,kdi,fm100)
c
      implicit none
c
      integer idim,jdim,i,j
      real cumppt(idim,jdim),avethbndy(idim,jdim),ythmois(idim,jdim),
     *     fm100 (idim,jdim)
      integer kdi(idim,jdim)
c
c  set dead start history fields
c
      do j=1,jdim
       do i=1,idim
         cumppt   (i,j)=0.0
         avethbndy(i,j)=15.0
         ythmois  (i,j)=15.0
         kdi      (i,j)=100
         fm100    (i,j)=10.
       enddo
      enddo
c
      return
      end
