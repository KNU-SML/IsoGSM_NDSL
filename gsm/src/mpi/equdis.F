      subroutine equdis(ind,len,ncut,lenarr,lendef)
c$$$  subprogram documentation block
c
c subprogram: 	equdis
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: cut len into ncut pieces with load balancing by 
c           symmetric distribution
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	equdiv(len,ncut,lenarr)
c
c    input argument lists:
c	ind	- integer spread direction: 1 for regular,
c                                          -1 for reverse
c	len	- integer total length 
c	ncut	- integer number of subgroup
c
c    output argument list:
c	lenarr	- integer (ncut) length of each subgroup
c	lendef	- integer (len) redefine the index 
c 
c subprograms called: none
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
      integer ind,len,ncut,lenarr,lendef
     *       ,nn,n,i,lens,lene
      dimension lenarr(ncut),lendef(len)
      integer, allocatable :: lentmp(:)
      allocate(lentmp(len))
c
      do i=1,ncut
        lenarr(i)=0
      enddo
      if( ind.eq.1 ) then
        lens=1
        lene=len
      else
        lens=len
        lene=1
      endif
      i=1
      n=1
      do nn=lens,lene,ind
        lenarr(n)=lenarr(n)+1
        lentmp(nn)=n
        n=n+i
        if(n.eq.ncut+1) then
          i=-1
          n=n+i
        endif
        if(n.eq.0) then
          i=1
          n=n+i
        endif
      enddo
c
      n=0
      do i=1,ncut
        do nn=1,len
          if(lentmp(nn).eq.i) then
            n=n+1
            lendef(n)=nn
          endif
        enddo
      enddo
c
      deallocate(lentmp)
c
      return
      end
