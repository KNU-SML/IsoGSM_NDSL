#define NPES npes
#define NCOL ncol
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define NROW nrow
#define LEVSP levsp
#define LONFP lonfp
#define LATGP latgp
#define JCAPP jcapp
#define LNTP lntp
#define LLNP llnp
      subroutine mpdimset(jcap,levs,lonf,latg)
c$$$  subprogram documentation block
c
c subprogram: 	mpdimset
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: preset all starting point and length for 
c           all pe for global spectral model.
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpdimset(jcap,levs,lonf,latg)
c
c    input argument lists:
c	jcap	- integer spectral wavenumber
c	levs	- integer vertical layer number
c	lonf	- integer gaussian grid for longitude
c	latg	- integer gaussian grid for latitude
c
c    output argument list:
c 
c subprograms called:
c   dimset	- to compute all dimension, starting point for all pe
c
c attributes:
c    language: fortran 90
c$$$
c
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer jcap,levs,lonf,latg
     *       ,levmax,lonmax,latmax,lntmax,lnpmax,lwvmax
c
      call dimset(jcap,levs,lonf,latg,
     *            npes,ncol,nrow,
     *            levstr,levlen,levmax,
     *            lonstr,lonlen,lonmax,
     *            latstr,latlen,latmax,
     *            lwvstr,lwvlen,lwvmax,
     *            lntstr,lntlen,lntmax,
     *            lnpstr,lnplen,lnpmax,
     *            lwvdef,latdef)
c
      if( levmax .ne. levsp_ ) then
        print *,' levmax levsp_ ',levmax,levsp_
        print *,' Error in mpdimset for levmax '
        call mpabort
      endif
c
      if( lonmax .ne. lonfp_ ) then
        print *,' lonmax lonfp ',lonmax,lonfp_
        print *,' Error in mpdimset for lonmax '
        call mpabort
      endif
c
      if( latmax .ne. latgp_ ) then
        print *,' latmax latgp_ ',latmax,latgp_
        print *,' Error in mpdimset for latmax '
        call mpabort
      endif
c
      if( lwvmax .ne. jcapp_ ) then
        print *,' lwvmax jcapp_ ',lwvmax,jcapp_
        print *,' Error in mpdimset for lwvmax '
        call mpabort
      endif
c
      if( lntmax .ne. lntp_ ) then
        print *,' lntmax lntp_ ',lntmax,lntp_
        print *,' Error in mpdimset for lntmax '
        call mpabort
      endif
c
      if( lnpmax .ne. llnp_ ) then
        print *,' lnpmax llnp_ ',lnpmax,llnp_
        print *,' Error in mpdimset for lnpmax '
        call mpabort
      endif
c
      return
      end
