#define NSTYPE nstype
      function funcdf(theta,ktype)
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgdf/ dfk(ngrid,ntype)
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
      w = (theta / tsat(ktype)) * 20. + 1.
      kw = w
      kw = min(kw,21)
      kw = max(kw,1)
      funcdf = dfk(kw,ktype)
     &         + (w - kw) * (dfk(kw+1,ktype) - dfk(kw,ktype))
      return
      end
