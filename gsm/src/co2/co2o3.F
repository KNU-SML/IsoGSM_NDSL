#define LEVS levs
      program co2o3
c$$$  main program documentation block
c
c main program:  co2o3
c   prgmmr: cpc          org: w/np51     date: 01-03-31
c
c abstract:
cccc  program co2o3 = consolidation of a number of dan schwarzkopf,gfdl
c                     codes to produce a file of co2 hgt data
c                     for any vertical coordinate (read by subroutine
c                     conrad in the gfdl radiation codes)-k.a.c. jun89.
cnov89--updated (nov 89) for latest gfdl lw radiation.....k.a.c.
c program history log:
c
c namelists:
c
c input files:
c   unit
c
c output files:
c   unit   66
c
c subprograms called:
c     co2ptz
c     co2int
c     co2in1
c     co2ins
c
c attributes:
c   language: fortran
c
c$$$
#include <paramodel.h>
      parameter (l=levs_,lp1=l+1,lp2=l+2)
      dimension sgtemp(lp1,2),co2d1d(l,6),co2d2d(lp1,lp1,6)
cnov89
      dimension co2iq2(lp1,lp1,6),co2iq3(lp1,lp1,6),co2iq5(lp1,lp1,6)
cnov89
      dimension t41(lp2,2),t42(lp1),
     1          t43(lp2,2),t44(lp1)
      dimension t20(lp1,lp1,3),t21(lp1,lp1,3)
      dimension t22(lp1,lp1,3),t23(lp1,lp1,3)
      dimension sglvnu(lp1),siglnu(l)
      dimension stemp(lp1),gtemp(lp1)
      dimension cdtm51(l),co2m51(l),c2dm51(l)
      dimension cdtm58(l),co2m58(l),c2dm58(l)
      dimension cdt51(lp1,lp1),co251(lp1,lp1),c2d51(lp1,lp1)
      dimension cdt58(lp1,lp1),co258(lp1,lp1),c2d58(lp1,lp1)
cnov89
      dimension cdt31(lp1),co231(lp1),c2d31(lp1)
      dimension cdt38(lp1),co238(lp1),c2d38(lp1)
      dimension cdt71(lp1),co271(lp1),c2d71(lp1)
      dimension cdt78(lp1),co278(lp1),c2d78(lp1)
      dimension co211(lp1),co218(lp1)
      equivalence (cdt31(1),co2iq2(1,1,1)),(co231(1),co2iq2(1,1,2))
      equivalence (c2d31(1),co2iq2(1,1,3)),(cdt38(1),co2iq2(1,1,4))
      equivalence (co238(1),co2iq2(1,1,5)),(c2d38(1),co2iq2(1,1,6))
      equivalence (cdt71(1),co2iq3(1,1,1)),(co271(1),co2iq3(1,1,2))
      equivalence (c2d71(1),co2iq3(1,1,3)),(cdt78(1),co2iq3(1,1,4))
      equivalence (co278(1),co2iq3(1,1,5)),(c2d78(1),co2iq3(1,1,6))
      equivalence (co211(1),co2iq5(1,1,2)),(co218(1),co2iq5(1,1,5))
cnov89
      equivalence (stemp(1),sgtemp(1,1)),(gtemp(1),sgtemp(1,2))
      equivalence (cdtm51(1),co2d1d(1,1)),(co2m51(1),co2d1d(1,2))
      equivalence (c2dm51(1),co2d1d(1,3)),(cdtm58(1),co2d1d(1,4))
      equivalence (co2m58(1),co2d1d(1,5)),(c2dm58(1),co2d1d(1,6))
      equivalence (cdt51(1,1),co2d2d(1,1,1)),(co251(1,1),co2d2d(1,1,2))
      equivalence (c2d51(1,1),co2d2d(1,1,3)),(cdt58(1,1),co2d2d(1,1,4))
      equivalence (co258(1,1),co2d2d(1,1,5)),(c2d58(1,1),co2d2d(1,1,6))
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_co2',2001,0000,0000,'np51')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c===>  get sgtemp and output which used to be on units 41,42,43,44....
      lread = 0
      call co2ptz(sgtemp,t41,t42,t43,t44,sglvnu,siglnu,lread)
c===>  interpolate desired co2 data from the detailed(109,109) grid..
c         ir=1,iq=1 is for common /co2bd3/ in radiation code...
c           for the consolidated 490-850 cm-1 band...
cnov89
      ico2tp=61
cnov89
      ir = 1
c===>  ratio is the actual co2 mixing ratio in units of 330 ppmv.....
c        e.g. sinv in coeint has been obtained using rco2=330 ppmv
      ratio = 1.0
      read(5,*,end=599) co2ppm
      ratio = co2ppm / 330.
599   continue
      nmethd = 2
      call co2int(ico2tp,t41,t42,t22,ratio,ir,nmethd)
      ir = 1
      nmethd = 1
      call co2int(ico2tp,t41,t42,t20,ratio,ir,nmethd)
      ir = 1
      nmethd = 2
      call co2int(ico2tp,t43,t44,t23,ratio,ir,nmethd)
      ir = 1
      nmethd = 1
      call co2int(ico2tp,t43,t44,t21,ratio,ir,nmethd)
c===>    fill up the co2d1d array
c       the following gets co2 transmission functions and
c         their derivatives for tau(i,i+1),i=1,levs,
c         where the values are not obtained by quadrature but are the
c         actual transmissivities,etc,between a pair of pressures. these
c         are used only for nearby layer calculations including h2o..
c
      iq = 1
      call co2in1(t20,t21,co2d1d,iq)
c
c===>    fill up the co2d2d array
c    the following gets co2 transmission functions and their derivatives
c        from 109-level line-by-line calculations made using the 1982
c        mcclatchy tape (12511 lines),consolidated,interpolated
c        to the mrf vertical coordinate,and re-consolidated to a
c        200 cm-1 bandwidth. the interpolation method is described in
c        schwarzkopf and fels (j.g.r.,1985).
c
      call co2ins(t22,t23,co2d2d,iq)
c
cnov89
c===>  interpolate desired co2 data from the detailed(109,109) grid..
c         ir=2,iq=2 is for common /co2bd2/ in radiation code...
c           for the consolidated 490-670 cm-1 band...
      ico2tp=62
      ir = 2
      nmethd = 2
      call co2int(ico2tp,t41,t42,t22,ratio,ir,nmethd)
      call co2int(ico2tp,t43,t44,t23,ratio,ir,nmethd)
      iq = 2
      call co2ins(t22,t23,co2iq2,iq)
c===>  interpolate desired co2 data from the detailed(109,109) grid..
c         ir=3,iq=3 is for common /co2bd4/ in radiation code...
c           for the consolidated 670-850 cm-1 band...
      ico2tp=63
      ir = 3
      nmethd = 2
      call co2int(ico2tp,t41,t42,t22,ratio,ir,nmethd)
      call co2int(ico2tp,t43,t44,t23,ratio,ir,nmethd)
      iq = 3
      call co2ins(t22,t23,co2iq3,iq)
c---      following code not working and not needed yet
c===>  interpolate desired co2 data from the detailed(109,109) grid..
c         ir=4,iq=5 is for common /co2bd5/ in radiation code...
c           for the 4.3 micron band....this for rco2=300 ppmv...
c not used yet      ico2tp=65
c not used yet      ir = 4
c not used yet      ratio = 0.909091
c not used yet      nmethd = 2
c not used yet      call co2int(ico2tp,t41,t42,t22,ratio,ir,nmethd)
c not used yet      call co2int(ico2tp,t43,t44,t23,ratio,ir,nmethd)
c not used yet      iq = 5
c not used yet      call co2ins(t22,t23,co2iq5,iq)
cnov89
c...     write data to disk..
c            ...since these codes are compiled with autodbl,the co2 data
c               is converted to single precision in a later job step..
      rewind 66
      write(66) stemp
      write(66) gtemp
      write(66) cdtm51
      write(66) co2m51
      write(66) c2dm51
      write(66) cdtm58
      write(66) co2m58
      write(66) c2dm58
      write(66) cdt51
      write(66) co251
      write(66) c2d51
      write(66) cdt58
      write(66) co258
      write(66) c2d58
cnov89
      write(66) cdt31
      write(66) co231
      write(66) c2d31
      write(66) cdt38
      write(66) co238
      write(66) c2d38
      write(66) cdt71
      write(66) co271
      write(66) c2d71
      write(66) cdt78
      write(66) co278
      write(66) c2d78
cjun93
      write(66) co2ppm*1.e-6
c not used yet      write(66) co211
c not used yet      write(66) co218
cnov89
      rewind 66
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tage('clim_co2')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop
      end
