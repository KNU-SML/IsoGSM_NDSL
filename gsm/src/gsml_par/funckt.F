#define NSTYPE nstype
      function funckt(theta,ktype)
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgkt/ ktk(ngrid,ntype)
      real ktk
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
      w = (theta / tsat(ktype)) * 20. + 1.
      kw = w
      kw = min(kw,21)
      kw = max(kw,1)
      funckt = ktk(kw,ktype)
     &         + (w - kw) * (ktk(kw+1,ktype) - ktk(kw,ktype))
      return
      end
