#include <machine.h>
      program inchour
c$$$  main program documentation block
c
c main program:  inchour    get increment hour
c   prgmmr: kanamitsu          org: w/np51     date: 01-03-31
c
c abstract:  compute increment hour
c
c program history log:
c   01-03-31  hann-ming juang  add w3tag calls for nco implementation
c
c namelists:
c   namin:      parameters determining new date
c
c input files:
c
c output files:
c
c subprograms called:
c
c attributes:
c   language: fortran
c
c$$$
c
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_inchour',2001,0000,0000,'np51')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c
c  given two dates, compute hour increment
c
      read(5,*) iys,ims,ids,ihs,iye,ime,ide,ihe
      call compjd(iye,ime,ide,ihe,0,jde,fjde)
      call compjd(iys,ims,ids,ihs,0,jds,fjds)
      inc=(float(jde-jds)+fjde-fjds)*24
#ifdef NCO_TAG
      write(51,*) inc
#endif
      print *,inc
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL w3tage('clim_inchour')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop
      end
c
      subroutine compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)
c
      dimension ndm(12)
      data jdor/2415019/,jyr19/1900/
      data ndm/0,31,59,90,120,151,181,212,243,273,304,334/
c     
      jd=jdor
      jyrm9=jyr-jyr19
      lp=jyrm9/4
c
c     For year 1900 and after...
c
      if(jyrm9.ge.0) then
         if(lp.gt.0) then
            jd=jd+1461*lp-1
         endif
         ny=jyrm9-4*lp
c
      else
c
c     For years before 1900 years
c
         jd=jd+1461*(lp-1)+366
         ny=jyrm9-4*(lp-1)-1
         jd=jd+abs((jyrm9+1)/100)-abs((jyr-2000+1)/400)
         if(mod(jyr,100).eq.0) then
            if(mod(jyr,400).ne.0) then
               jd=jd+1
            endif
         endif
c
      endif
c
      ic=0
c
c     For year 1900 and after..
c
      if(jyrm9.ge.0) then
c        
c        For years that are not mulitples of 4
         if(ny.gt.0) then
            jd=jd+365*ny+1
            if(lp.eq.0) then
c              Between 1901-1903
               jd=jd-1
            else
c              2100,2200,2300,2500 are not leap, but 2000, 2400 are
               if(mod(jyr-ny,400).ne.0.and.mod(jyr-ny,100).eq.0) then
                  jd=jd-1
               endif
            endif
c
         else
c
c        For years that are multiples of 4
            if(jmnth.gt.2) ic=1
c           For year 1900
            if(lp.eq.0) ic=0
c           2100,2200,2300,2500 are not leap, but 2000, 2400 are           
            if(mod(jyr,400).ne.0.and.mod(jyr,100).eq.0) then
               ic=0
            endif
         endif
      else
c
c     For years before 1900 years
c
         if(ny.lt.3) then
            jd=jd+365*ny
         else
            jd=jd+365*ny
            if(jmnth.gt.2) ic=1
            if(mod(jyr,400).ne.0.and.mod(jyr,100).eq.0) then
               if(jmnth.gt.2) then
                  ic=0
               endif               
            endif
         endif
      endif
c
      jd=jd+ndm(jmnth)+jday+ic
      if(jhr.ge.12) then
        fjd=.041666667e0*float(jhr-12)+.00069444444e0*float(jmn)
      else
        jd=jd-1
        fjd=.5e0+.041666667e0*float(jhr)+.00069444444e0*float(jmn)
      endif
      return
      end
