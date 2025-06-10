#define NSTYPE nstype
      function twlt(ktype)
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
      twlt = .1
      return
      end
