c
c   the following common blocks contain pretabulated co2 transmission
c       functions, evaluated using the methods of fels and
c       schwarzkopf (1981) and schwarzkopf and fels (1985),
c***common co2bd3 contains co2 transmission functions and temperature
c   and pressure derivatives for the 560-800 cm-1 band. also included
c   are the standard temperatures and the weighting function. these
c   data are in block data bd3:
c         co251    =  transmission fctns for t0 (std. profile)
c                       with p(sfc)=1013.25 mb
c         co258    =  transmission fctns. for t0 (std. profile)
c                       with p(sfc)= ^810 mb
c         cdt51    =  first temperature derivative of co251
c         cdt58    =  first temperature derivative of co258
c         c2d51    =  second temperature derivative of co251
c         c2d58    =  second temperature derivative of co251
c         co2m51   =  transmission fctns for t0 for adjacent pressure
c                        levels, with no pressure quadrature. used for
c                        nearby layer computations. p(sfc)=1013.25 mb
c         co2m58   =  same as co2m51,with p(sfc)= ^810 mb
c         cdtm51   =  first temperature derivative of co2m51
c         cdtm58   =  first temperature derivative of co2m58
c         c2dm51   =  second temperature derivative of co2m51
c         c2dm58   =  second temperature derivative of co2m58
c         stemp    =  standard temperatures for model pressure level
c                        structure with p(sfc)=1013.25 mb
c         gtemp    =  weighting function for model pressure level
c                        structure with p(sfc)=1013.25 mb.
c         b0       =  temp. coefficient used for co2 trans. fctn.
c                        correction for t(k). (see ref. 4 and bd3)
c         b1       =  temp. coefficient, used along with b0
c         b2       =  temp. coefficient, used along with b0
c         b3       =  temp. coefficient, used along with b0
c
      common/co2bd3/co251(lp1,lp1),co258(lp1,lp1),cdt51(lp1,lp1),
     1 cdt58(lp1,lp1),c2d51(lp1,lp1),c2d58(lp1,lp1),co2m51(l),
     2 co2m58(l),cdtm51(l),cdtm58(l),c2dm51(l),c2dm58(l),
     3 stemp(lp1),gtemp(lp1),b0,b1,b2,b3
c
c   common co2bd2 contains co2 transmission functions and temperature
c   and pressure derivatives for the 560-670 cm-1 part of the 15 um
c   co2 band.  these data are in block data bd2.
c         co231    =  transmission fctns for t0 (std. profile)
c                       with p(sfc)=1013.25 mb
c         co238    =  transmission fctns. for t0 (std. profile)
c                       with p(sfc)= ^810 mb
c         cdt31    =  first temperature derivative of co231
c         cdt38    =  first temperature derivative of co238
c         c2d31    =  second temperature derivative of co231
c         c2d38    =  second temperature derivative of co231
c
      common / co2bd2 / co231(lp1),co238(lp1),cdt31(lp1),
     1 cdt38(lp1),c2d31(lp1),c2d38(lp1)
c
c   common co2bd4 contains co2 transmission functions and temperature
c   and pressure derivatives for the 670-800 cm-1 part of the 15 um
c   co2 band.  these data are in block data bd4.
c         co271    =  transmission fctns for t0 (std. profile)
c                       with p(sfc)=1013.25 mb
c         co278    =  transmission fctns. for t0 (std. profile)
c                       with p(sfc)= ^810 mb
c         cdt71    =  first temperature derivative of co271
c         cdt78    =  first temperature derivative of co278
c         c2d71    =  second temperature derivative of co271
c         c2d78    =  second temperature derivative of co271
c
      common / co2bd4 / co271(lp1),co278(lp1),cdt71(lp1),
     1 cdt78(lp1),c2d71(lp1),c2d78(lp1)
c
c   common co2bd5 contains co2 transmission functions for the 2270-
c   2380 part of the 4.3 um co2 band. these data are in block data bd5.
c         co211    =  transmission fctns for t0 (std. profile)
c                       with p(sfc)=1013.25 mb
c         co218    =  transmission fctns. for t0 (std. profile)
c                       with p(sfc)= ^810 mb
c
      common / co2bd5 / co211(lp1),co218(lp1)
