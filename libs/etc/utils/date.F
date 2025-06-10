#include <machine.h>
      program date
c$$$  main program documentation block
c
c main program:  date    get date from file
c   prgmmr: kanamitsu          org: w/np51     date: 01-03-31
c
c abstract: get date from model files
c
c program history log:
c   01-03-31  hann-ming juang  add w3tag calls for nco implementation
c
c namelists:
c   namin:      parameters determining new date
c
c input files:
c   unit   11  sigma file(s)
c
c output files:
c   unit   51  sigma file
c
c subprograms called:
c
c attributes:
c   language: fortran
c
c$$$
      dimension idate(4),lab(4)
c
      character*120 fnam
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_date',2001,0000,0000,'np51')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c   
      read(5,200) fnam
  200 format(A120)
c
      open(unit=11,file=fnam,form='unformatted',err=900)
c
      read(11) 
c
      read(11) fhour,idate
c
      ihour=nint(fhour)
c
#ifdef NCO_TAG
      write(51,100) idate(1),idate(2),idate(3),idate(4),ihour
#endif
      write(6,100) idate(1),idate(2),idate(3),idate(4),ihour
100   format(4i6,i10)
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL w3tage('clim_date')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop
c
  900 continue
      write(6,*) 'file open error in program date'
      call abort
      end
