#define BGF bgf
#define BORDER border
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1 jgrd1
#define JGRD1P jgrd1p
#define JGRD jgrd
#define JWAV1 jwav1
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLWAVP llwavp
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNT2 lnt2
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac

      subroutine gradini
#include <paramodel.h>
c
#ifndef RSM
#include <comfgrid.h>
#include <comfphys.h>
#include <comfver.h>
#include <comio.h>
#include <comgrad.h>
#include <radiag.h>
#ifdef MP
#include <commpi.h>
#endif
#else
#include <rdparm.h>
#include <rscomf.rerun.h>
#include <rscomltb.h>
#include <rscommap.h>
#include <rscomgrad.h>
#ifdef MP
#include <comrmpi.h>
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#else
#define IGRD12S igrd12_
#define JGRD12S jgrd12_
#endif
#endif
c
      save rad1st
      logical rad1st
      data rad1st/.true./
c
      if(rad1st) then
#ifdef RSM
        rrs2=1.0
#endif
        runrad = .true.
#ifndef SW99UPDATE
        nalaer=49
        kalb=1
        istrat=1
        ko3       = 1 ! nasa o3
        ibnd      = 1 ! =1:use one nir band, =2:use three nir bands
        iswsrc(1) = 1 ! aerosol
        iswsrc(2) = 0 ! o2
        iswsrc(3) = 0 ! co2
        iswsrc(4) = 1 ! water vapor
        iswsrc(5) = 1 ! o3
#ifdef ICECLOUD
        icfc   = 1
        icwp   = 0
        if (ncldg_.gt. 0) icwp = 1
#else
        icfc = 0
        icwp = 0
#endif
#else
        nfalb =49
        nfaer =47
        kalb=1
        istrat=1
        ko3       = 1 ! nasa o3
        ibnd      = 2 ! =1:use one nir band, =2:use three nir bands
        iswsrc(1) = 1 ! aerosol
        iswsrc(2) = 1 ! o2
        iswsrc(3) = 1 ! co2
        iswsrc(4) = 1 ! water vapor
        iswsrc(5) = 1 ! o3
        icfc   = 1
        icwp   = 0
        if (ncldg_.gt. 0) icwp = 1
#endif
#ifdef RSM
        jsno=0
#endif
c
c   get cld-rh relations in tabular form for day 1
c
#ifdef MP
        if( mype.eq.master ) then
#endif
          call crhtab(rhcl,ier)
#ifdef MP
        endif
#ifndef RSM
        call mpbcasti(ier,1)
#else
        call rmpbcasti(ier,1)
#endif
#endif
        if (ier.lt.0) then
          istrat = 0
#ifdef MP
          if( mype.eq.master ) then
#endif
          print 1113
 1113     format(1h ,'===>tuning tables not available..abort')
#ifdef MP
          endif
#endif
#ifdef MP
#ifdef RSM
          call rmpabort
#else
          call mpabort
#endif
#else
          call abort
#endif
        end if
#ifdef MP
        nall=nbin*nlon*nlat*mcld*nseal
#ifndef RSM
        call mpbcastr(rhcl,nall)
#else
        call rmpbcastr(rhcl,nall)
#endif
#endif
c
#ifdef DBG
        if(mype.eq.master) print 1114,istrat
 1114   format(1h ,'..for diagnosed clds....istrat = ',i4)
#endif
#ifdef RSM
c
        do j=1,JGRD12S
          do i=1,IGRD12S
            sinlar(i,j)=sin(rlat(i,j))
            coslar(i,j)=sqrt(1. e 0 - sinlar(i,j)*sinlar(i,j))
          enddo
        enddo
#endif
c
        iunco2 = 15
        call gcljms(si)
c
cyh95.. add new data initialization routines for rad and aerosols
c
        call grlwsw(sl,iunco2)
#ifndef SW99UPDATE
#ifndef ICECLOUD
        call gaeros(si,sl)
#else
        call gaeros8(si,sl)
#endif
#else
        call gaeros99(si,sl)
#endif
#ifndef RSM
c
c   specify the latitude where permanent snow resides poleward
c     jsno=latitude closest to 70 deg n-indicating extent of perm
c     snow cover. will change as a fcn. of latitude structure
c
        jsno=(latg_+1)/9
#endif
c
c  get interval (hrs) between short-wave radiation calls.....
c
        raddt = 3600. * dtswav
        dtlw = 3600. * dtlwav
c
c  paerf when SWRMDC =.false.
c  alvsf,alnsf,alvwf,alnwf from the albaer.snl will be discarded
c   when NEWALB = .true.
c
#ifndef RSM
        call insur2(nalaer,alvsf,alnsf,alvwf,alnwf,
#ifndef SW99UPDATE
     1            facsf,facwf,paerf)
#else
     1            facsf,facwf)
#endif
#else
        call rinsur2(IGRD12S,JGRD12S,nalaer,alvsf,alnsf,alvwf,alnwf,
#ifndef SW99UPDATE
     1            facsf,facwf,paerf)
#else
     1            facsf,facwf)
#endif
#endif
c
        rad1st=.false.
      endif
c
      return
      end
