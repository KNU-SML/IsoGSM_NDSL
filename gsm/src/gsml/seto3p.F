#define JCAP jcap
#define LATG latg
#define NPES npes
c
c prgmmr:  yifeng cui added RSM MP     org:sdsc     date:02-08-05
c
      subroutine seto3p (idate, fhour, ko3p, ko3l)
      
#include <paramodel.h>
#ifdef MP
#ifdef RSM
#include <comrmpi.h>
#else
#include <commpi.h>
#endif
#endif
c
      parameter (jo3=18,ko3=46,blat=-85.0,dphi=10.0)
c
      common /o3pcom/prdin, disin
c
      dimension prdin(jo3,ko3), disin(jo3,ko3)
      dimension idate(4)
      data nskipo/-1/
      save nskipo
c
cfpp$ noconcur r
c
      call incdte(idate(4),idate(2),idate(3),idate(1),
     1     jyy,jmm,jdd,jhh,nint(fhour))
      mon=jmm
      iday=jdd
c
      nskip = (mon-1)*3 + (iday-1) / 10
      if (iday .eq. 31) nskip = nskip -1
c
      if (nskip .ne. nskipo) then
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef RSM
         open(unit=ko3p,file='r_o3prod',form='formatted',status='old')
         open(unit=ko3l,file='r_o3loss',form='formatted',status='old')
#else
         open(unit=ko3p,file='o3prod',form='formatted',status='old')
         open(unit=ko3l,file='o3loss',form='formatted',status='old')
#endif
         rewind (ko3p)
         rewind (ko3l)
         do i=1,nskip
            read (ko3p,*) prdin
            read (ko3l,*) disin
         enddo
         read (ko3p,*) prdin
         read (ko3l,*) disin
#ifdef MP
      endif
#endif
#ifdef MP
#ifdef RSM
      call rmpbcastr(prdin,jo3*ko3)
      call rmpbcastr(disin,jo3*ko3)
#else
      call mpbcastr(prdin,jo3*ko3)
      call mpbcastr(disin,jo3*ko3)
#endif
#endif
c
         nskipo = nskip
      endif
c
      return
      end
