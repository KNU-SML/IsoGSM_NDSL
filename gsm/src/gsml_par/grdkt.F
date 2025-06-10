#define NSTYPE nstype
      subroutine grdkt
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgdft/ dfkt(ngrid,ntype)
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
      do k = 1, ntype
        if(satpsi(k).gt.0.) then
        dynw = tsat(k) * .05
        f1 = log10(satpsi(k)) + b(k) * log10(tsat(k)) + 2.
        do i = 1, ngrid
          theta = float(i-1) * dynw
          theta = min(tsat(k),theta)
          if(theta.gt.0.) then
           
            pf = f1 - b(k) * log10(theta)
          else
            pf = 5.2
          endif
          if(pf.le.5.1) then
            dfkt(i,k) = exp(-(2.7+pf)) * 420.
          else
            dfkt(i,k) = .1744
          endif
        enddo
        endif
      enddo
      return
      end
