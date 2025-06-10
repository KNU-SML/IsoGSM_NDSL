#define IOSGB iosgb
      program sgb
c$$$  main program documentation block
c
c main program:  sgb         transform sigma to sigma grib
c   prgmmr: iredell          org: w/nmc23     date: 95-01-19
c
c abstract: program transforms sigma input to sigma grib output.
c   by default, temperature, specific humidity, diverence and vorticity,
c   zonal and meridional winds, surface pressure and its gradient,
c   and surface orography and its gradient are all output on the
c   model physical gaussian grid and model sigma levels.
c   parameters and grid output can be controlled by namelist input.
c
c program history log:
c   01-03-31  hann-ming juang   add w3tag for nco implementation
c
c namelists:
c   namsgb:      parameters determining output format
c     io         number of longitude points (default: set by rdsgh)
c     jo         number of latitude points (default: set by rdsgh)
c     ncpus      number of parallel processes (default: environment)
c     mxbit      maximum number of bits to pack data (default: 16)
c     ids(255)   decimal scaling of packed data
c                (set to less than -128 to skip parameter altogether)
c                (default: set by subprogram idsdef)
c     icen       forecast center identifier (default: 7)
c     icen2      forecast sub-center identifier (default: 0)
c     igen       model generating code (default: from sigma file)
c
c input files:
c   unit   11    sigma file(s)
c
c output files:
c   unit   51    sigma grib1 file(s)
c
c subprograms called:
c   gncpus       get environment number of parallel processes
c   idsdef       set defaults for decimal scaling
c   rdsgh        read a sigma file header
c   sgb1         transform one sigma file to sigma grib
c
c attributes:
c   language: fortran
c
c$$$
#include <paramodel.h>
      parameter(levmax=100)
c
      dimension idate(4)
      dimension si(levmax+1),sl(levmax)
      dimension ids(255)
      logical lhgrad
      data lhgrad/.false./
c
      namelist/namsgb/ ids,icen,icen2,igen,fhs,fhe,fhinc,lhgrad
      data ids/255*0/
      data icen/7/,icen2/0/,igen/0/
      data fhs/0./,fhe/24./,fhinc/24./
c
#ifdef NCO+TAG
c
      call w3tagb('clim_sgb',2001,0000,0000,'np51')
c
#endif
c
      call gncpus(ncpus)
      print *,'ncpus=',ncpus
      call idsdef(2,ids)
c
c  idrt=4 .. gaussian, idrt=0 lat/lon
c
      if(io_.eq.36.or.io_.eq.72.or.io_.eq.120.or.
     1   io_.eq.144.or.io_.eq.180.or.io_.eq.360) then
         idrt=0
         write(6,*) 'lat/lon grid assumed'
      else
         idrt=4
         write(6,*) 'gaussian grid assumed'
      endif
      read(*,namsgb,end=5)
   5  continue
c
      nsig=11
      nsgb=51
c
      fh=fhs
      dowhile (fh.le.fhe)
        call rdsgh(nsig,fh,fhour,idate,si,sl,iret,ntrace,ncloud)
        if(iret.ne.0) then
          print *,'read sigma file error, iret=',iret
          call abort
        endif
        if(fh.ne.fhour) then
           print *,'warning!!! fh on filename inconsistent'
           print *,'fh=',fh,' fhour=',fhour
        endif
c
        call sgb1(idrt,fhour,idate,nsig,si,sl,
     $            ntrace,ncloud,
     &            nsgb,ncpus,ids,icen,icen2,igen,lhgrad)
        close(nsig) 
        close(nsgb) 
        fh=fh+nint(fhinc)
      enddo
#ifdef NCO_TAG
c
      call w3tage('clim_sgb')
c
#endif
      stop
      end
