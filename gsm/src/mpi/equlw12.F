      subroutine equlw12(jcap,nrow,lwvdis,lnpdis,jcaprm)
c$$$  subprogram documentation block
c
c subprogram: 	equlw12
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: about equal division of all triangular truncation waves
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	equlw12(jcap,nrow,lwvdis,lnpdis,jcaprm)
c
c    input argument lists:
c	jcap	- integer spectral wavenumber
c	nrow	- integer number of nrow
c
c    output argument list:
c	lwvdis	- integer (nrow) length of wave in l
c	lnpdis	- integer (nrow) length of wave for l and n
c	jcaprm	- integer remain of unbalance in l
c 
c subprograms called:
c   equdiv	- to compute about equal number of subgroup 
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
      integer jcap,nrow,lwvdis,lnpdis,jcaprm
     *       ,jcaphf,nr,n,l,m,lw1pnt,lw2pnt
      dimension lwvdis(nrow),lnpdis(nrow)
      integer, allocatable :: lw1str(:),lw1end(:)
      integer, allocatable :: lw2str(:),lw2end(:)
c
      allocate (lw1str(nrow))
      allocate (lw1end(nrow))
      allocate (lw2str(nrow))
      allocate (lw2end(nrow))
c
      jcaphf=(jcap+1)/2
      jcaprm=mod(jcap+1,2)
      call equdiv(jcaphf,nrow,lwvdis)
c
      lw1pnt=1
      lw2pnt=jcap+2
c
      do n=1,nrow
        lw2end(n)=lw2pnt-1
        lw2pnt=lw2pnt-lwvdis(n)
        lw2str(n)=lw2pnt
        lw1str(n)=lw1pnt
        lw1pnt=lw1pnt+lwvdis(n)
        lw1end(n)=lw1pnt-1
      enddo
      lw1end(nrow)=lw1end(nrow)+jcaprm
c
      do nr=1,nrow
c
        m=0
        do l=lw1str(nr),lw1end(nr)
          do n=l,jcap+1
            m=m+1
          enddo
        enddo
        do l=lw2str(nr),lw2end(nr)
          do n=l,jcap+1
            m=m+1
          enddo
        enddo
        lnpdis(nr)=m
c
      enddo
      deallocate (lw1str)
      deallocate (lw1end)
      deallocate (lw2str)
      deallocate (lw2end)
c
      return
      end
