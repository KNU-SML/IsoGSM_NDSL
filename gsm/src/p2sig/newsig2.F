#define CP cp
#define RD rd
      subroutine newsig2(si,sl,del)
#include<paramodel.h>
c
c this routine sets coordinates for levels and use phillips method
c
#include <modlsigs.h>
c
      parameter(kdim=levmdl,kdimp=kdim+1,kdimm=kdim-1)
      dimension si(kdimp),sl(kdim),del(kdim)
      real rk,rk1,rkinv
      real, allocatable :: ci(:),cl(:),rpi(:)
c
      allocate (ci(kdimp),cl(kdim),rpi(kdimm))
c
      ci(1) = 0.
      do k=1,kdim
        del(k)=delmdl(k)
        ci(k+1)=ci(k)+delmdl(k)
      enddo
      ci(kdimp)=1.
c
      rk  = rd_/cp_
      rk1 = rk + 1.
      rkinv=1./rk
c
      levs=kdim
c
      do li=1,kdimp
        si(li) = 1. - ci(li)
      enddo
c
      do le=1,kdim
        dif = si(le)**rk1 - si(le+1)**rk1
        dif = dif / (rk1*(si(le)-si(le+1)))
        sl(le) = dif**rkinv
        cl(le) = 1. - sl(le)
      enddo
c
c     compute pi ratios for temp. matrix.
c
      do le=1,kdimm
        rpi(le) = (sl(le+1)/sl(le))
      enddo
      do le=1,kdimm
        rpi(le) = rpi(le)**rk
      enddo
c
      do 6 le=1,kdimp
      print 100, le, ci(le), si(le)
100   format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)
    6 continue
c
      print 200
  200 format (1h0)
      do 7 le=1,kdim
      print 300, le, cl(le), sl(le), del(le)
  300 format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,
     1 'del=', f6.3)
    7 continue
c
      print 400, (rpi(le), le=1,kdimm)
  400 format (1h0, 'rpi=', (18(1x,f6.3)) )
c
      deallocate (ci,cl,rpi)
      return
      end
