      program mtn
c$$$  main program documentation block
c
c main program:  mtn  terrain maker for global spectral model
c   prgmmr: iredell       org: w/nmc23       date: 92-04-16
c
c abstract: this program creates 4 terrain-related files
c   computed from the navy 5-minute terrain dataset.
c   the model physics grid parameters and spectral truncation
c   are read by this program as input.
c   the 4 files produced are respectively:
c     1) sea-land mask on model physics grid used by asfc program
c     2) gridded orography on model physics grid used by asfc program
c     3) mountain variance on model physics grid used by fcst program
c     4) spectral orography in spectral domain used by chgr program
c
c program history log:
c   92-04-16  iredell
c   96-01-29  kanamitsu  modified to read from ascii file
c   98-03-10  kanamitsu  modified to do smoothing (hoskins)
c   98-13-10  hong       modified for subgrid orography used in kim's scheme
c   98-13-10  hong       modified for high-resolution gtopo new orography 
c   01-03-31  juang	 add w3tag for nco implementation
c
c usage:
c
c   input files:
c     unit5      - physics longitudes (im), physics latitudes (jm),
c                  spectral truncation (nm) and rhomboidal flag (nr)
c                  respectively read in free format.
c     unit11     - packed navy 10-minute terrain dataset.
c
c   output files:
c     unit51     - sea-land mask (im,jm), lon(im,jm), lat(im,jm)
c     unit52     - gridded orography (im,jm)
c     unit53     - mountain variance (im,jm)
c     unit54     - spectral orography ((nm+1)*((nr+1)*nm+2))
c     unit55     - unfiltered gridded orography (im,jm)
c
c   subprograms called:
c     unique:
c     tersub     - main subprogram
c     sphert     - spherical transform
c     gausslat   - compute gaussian latitudes
c     equallat   - compute equally-spaced latitudes
c     legendre   - compute legendre polynomials
c     fftfax     - fft (library call can be substituted)
c     rfftmlt    - fft (library call can be substituted)
c     library:
c     gbytes     - unpack bits
c     qpass      - fft utility
c     rpass      - fft utility
c
c   remarks: fortran 9x extensions are used.
c
c attributes:
c   language fortran
c
c$$$
#ifdef CRAY_THREAD
cfpp$ noconcur f
#endif
C
      parameter(nv=1)
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_mtn',2001,0000,0000,'np51   ')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c
      read(5,*) im,jm,nm,nr
      print*, im,jm,nm,nr,nv
      if(nr.ne.0) then
        print *,'rohmboidal truncation not allowed'
        stop
      endif
      nw=(nm+1)*((nr+1)*nm+2)
      call tersub(im,jm,nm,nr,nw,nv)
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tage('clim_mtn')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop
      end
