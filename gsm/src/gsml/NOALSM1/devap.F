#include <define.h>
      function devap (etp1,smc,zsoil,shdfac,smcmax,bexp,
     &                dksat,dwsat,smcdry,smcref,smcwlt,fxexp)
c
      implicit none
c
c ----------------------------------------------------------------------
c function devap
c ----------------------------------------------------------------------
c calculate direct soil evaporation
c ----------------------------------------------------------------------
      real bexp
      real devap
      real dksat
      real dwsat
      real etp1
      real fx
      real fxexp
      real shdfac
      real smc
      real smcdry
      real smcmax
      real zsoil
      real smcref
      real smcwlt
      real sratio
c ----------------------------------------------------------------------
c direct evap a function of relative soil moisture availability, linear
c when fxexp=1.
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c fx > 1 represents demand control
c fx < 1 represents flux control
c ----------------------------------------------------------------------
#ifndef NOAH_JUN05
      fx = ( (smc - smcdry) / (smcmax - smcdry) )**fxexp
      fx = max ( min ( fx, 1. ) ,0. )
#else
      sratio = (smc - smcdry) / (smcmax - smcdry)
      if(sratio.gt.0.) then
       fx=sratio**fxexp
       fx = max ( min ( fx, 1. ) ,0. )
      else
       fx=0.
      endif
#endif
c ----------------------------------------------------------------------
c allow for the direct-evap-reducing effect of shade
c ----------------------------------------------------------------------
      devap = fx * ( 1.0 - shdfac ) * etp1
c ----------------------------------------------------------------------
c end function devap
c ----------------------------------------------------------------------
      return
      end
