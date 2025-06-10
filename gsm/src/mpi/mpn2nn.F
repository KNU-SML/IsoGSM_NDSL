#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpn2nn(a,lntp,b,lntpp,nvar)
c$$$  subprogram documentation block
c
c subprogram: 	mpn2nn
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  transpose (lntp,nvar) to (lntpp,nvar)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpn2nn(a,lntp,b,lntpp,nvar)
c
c    input argument lists:
c	a	- real (lntp,nvar) partial field 
c	lntp	- integer partial spectral grid 
c	lntpp	- integer sub partial spectral grid 
c	nvar	- integer total set of fields
c
c    output argument list:
c	b	- real (lntpp,nvar) sub partial field 
c 
c subprograms called:
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer lntpp,lntp,nvar,n,m,ipe0,lntstr0
      real    a,b
      dimension a(lntp,nvar),b(lntpp,nvar)
c
c option for 1-d decomposition
      if( ncol.eq.1 ) then
        do n=1,nvar
          do m=1,lntp
            b(m,n) = a(m,n)
          enddo
        enddo
        return
      endif
c
c otherwise
c cut the part and through away the rest.
      ipe0=int(mype/ncol)*ncol
      lntstr0=lntstr(mype)*2-lnpstr(ipe0)*2
      do n=1,nvar
        do m=1,lntlen(mype)*2
          b(m,n)=a(m+lntstr0,n)
        enddo
      enddo
c
      return
      end
