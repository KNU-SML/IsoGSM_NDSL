      subroutine equdiv(len,ncut,lenarr)
c$$$  subprogram documentation block
c
c subprogram: 	equdiv
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: cut len into ncut pieces with load balancing
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	equdiv(len,ncut,lenarr)
c
c    input argument lists:
c	len	- integer total length 
c	ncut	- integer number of subgroup
c
c    output argument list:
c	lenarr	- integer (ncut) length of each subgroup
c 
c subprograms called: none
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
      integer len,ncut,lenarr
     *       ,n0,n1,n
      dimension lenarr(ncut)
      n0=len/ncut
      n1=mod(len,ncut)
      do n=1,n1
        lenarr(n)=n0+1
      enddo
      do n=n1+1,ncut
        lenarr(n)=n0
      enddo
      return
      end

