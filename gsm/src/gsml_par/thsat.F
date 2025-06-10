#define NSTYPE nstype
      function thsat(ktype)
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
      thsat = tsat(ktype)
      return
      end
