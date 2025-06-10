#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAPP lcapp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define MLVARK mlvark
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTRAC ntrac
#define SLVARK slvark
#define TWOJ1 twoj1
#define NTOTAL ntotal
        program gsm
c$$$  main program documentation block
c
c main program: gsm         make global forecast with spectral model
c   prgmmr: sela             org: nmc23       date: 81-01-01
c
c abstract: make global forecast with spectral model.
c
c program history log:
c   81-01-01  sela
c   99-05-21  hann-ming henry juang   add mpi in one loop version
c   99-06-22  hann-ming henry juang   add mpi in all loop version
c   00-06-01  hann-ming henry juang   change gsm into gsmini, gsmstep, and gsmend
c   00-07-01  song-you hong  add tracers and prognostic clouds
c   00-09-01  song-you hong  add icecloud for radiation
c   01-03-31  hann-ming juang   add w3tag for nco implementation
Clu_rev: Noah LSM implementation
c   02-16-03  change i/o to include fixed fields and new state variables
c   02-19-03  merge progtm/progt2
c
c input files:
c   unit   11    sigma file (analysis or at time t-dt)
c   unit   12    sigma file (at time t if not analysis)
c   unit   14    surface file
c   unit   15    co2 constants (dependent on vertical resolution)
c   unit   24    mountain variance (dependent on horizontal resolution)
c   unit   43    cloud tuning
c
c output files:
c   unit   51    sigma file (at time t-dt)
c   unit   52    sigma file (at time t)
c   unit   53    surface file
c   unit   61    initial zonal diagnostics
c   unit   63    flux diagnostics
c   unit   64    final zonal diagnostics
c   unit   67    grid point diagnostics
c
c attributes:
c   language: cray fortran
c
c$$$
c
#include <paramodel.h>
#include <comfibm.h>
c
#include <comcon.h>
#include <comgpd.h>
c
#ifdef LFM
#include <comlfm.h>
#endif
c
#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
#endif
c
#ifdef MP
#include <commpi.h>
#endif
#ifdef CLM_CWF
      parameter (npmx=1000)
      common /comfcwrk/cgs(LONF2S,levs_,latg2_),
     +                 cwgt(npmx,lonf2_,latg2_)
      common /comfcwrki/inx(npmx,lonf2_,latg2_),jnx(npmx,lonf2_,latg2_),
     +                 npnt(lonf2_,latg2_), lwconv
      logical lwconv
#endif
#ifdef SAS_DIAG
c ...............................................................
      common/comsdiag/dcu(LONF2S,levs_,latg2_),dcv(LONF2S,levs_,latg2_),
     +                dct(LONF2S,levs_,latg2_),dcq(LONF2S,levs_,latg2_),
     +                dch(LONF2S,levs_,latg2_),
     +                fcu(LONF2S,levs_,latg2_),fcd(LONF2S,levs_,latg2_),
     +                deltb(LONF2S,latg2_),delqb(LONF2S,latg2_),
     +                delhb(LONF2S,latg2_),cbmf(LONF2S,latg2_),
     +                dlt(LONF2S,levs_,latg2_),dlq(LONF2S,levs_,latg2_),
     +                dlh(LONF2S,levs_,latg2_)
#endif
#ifdef CONDENS
c ...............................................................
      common/comcond/
     +       qcicnv(LONF2S,levs_,latg2_), qrscnv(LONF2S,levs_,latg2_),
     +       qcilrg(LONF2S,levs_,latg2_), qrslrg(LONF2S,levs_,latg2_),
     +       taucld(LONF2S,levs_,latg2_),
     +       cldwp(LONF2S,levs_,latg2_),  cldip(LONF2S,levs_,latg2_)
#endif

#ifdef MP
      call mpinit(strwtime)
#else
#ifdef T3E
      strwtime=secondr()
#endif
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c initial call to global spectral model for getcon etc.
c
      call gsmini
c
c
c     start time loop
c
      do 10000 jdt=limlow,maxstp
c
#ifdef TIMEF
        tstep0=timef()*0.001
#endif
c
c call one time step for gsm related works (can be forward-step or centered-step)
c
        call gsmstep
c
#ifdef TIMEF
        tstep=timef()*0.001-tstep0
#ifndef NOPRINT
        if(mype.eq.master) 
     1   print*,' time step, time per step: ',jdt,tstep
#endif
#endif
c
10000 continue
c
#ifdef MP
#ifdef COUPLE_GSM_AND_ROMS
      call send_flx_recv_sst
#endif
#endif
c
#ifdef DBG
      if(mype.eq.master) print *,' time loop ended '
#endif
c
c  end  time  loop
c
      call gsmend
c
#ifdef MP
      call mpfine(endwtime)
      print *,' MPI wall time is ',endwtime-strwtime
#else
#ifdef T3E
      endwtime=secondr()
      print *,' NONMPI wall time is ',endwtime-strwtime
#endif
#endif
      stop
      end
