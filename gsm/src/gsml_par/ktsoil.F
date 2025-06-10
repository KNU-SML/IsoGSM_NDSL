#define NSTYPE nstype
      real function ktsoil(theta,ktype)
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgdft/ dfkt(ngrid,ntype)
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
      w = (theta / tsat(ktype)) * 20. + 1.
      kw = w
      kw = min(kw,21)
      kw = max(kw,1)
      ktsoil = dfkt(kw,ktype)
     &         + (w - kw) * (dfkt(kw+1,ktype) - dfkt(kw,ktype))
      return
      end
