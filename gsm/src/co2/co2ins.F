#define LEVS levs
      subroutine co2ins(t22,t23,t66,iq)
c     *********************************************************
c       save data on permanent data set denoted by co222 ******
c          ..... k.campana   march 1988,october 1988...
c          ..... k.campana   december 1988-cleaned up for launcher
c          ..... k.campana   november 1989-altered for new radiation
c     *********************************************************
#include <paramodel.h>
      parameter (l=levs_,lp1=l+1)
      dimension t22(lp1,lp1,3),t23(lp1,lp1,3),t66(lp1,lp1,6)
      dimension dcdt8(lp1,lp1),dcdt10(lp1,lp1),co2po(lp1,lp1),
     * co2800(lp1,lp1),co2po1(lp1,lp1),co2801(lp1,lp1),co2po2(lp1,lp1),
     * co2802(lp1,lp1),n(lp1),d2ct8(lp1,lp1),d2ct10(lp1,lp1)
ccc   itin=22
ccc   itin1=23
co222  latest code had  iq=1
ccc      iq=4
1011  format (4f20.14)
ccc      read (itin,1011) ((co2po(i,j),i=1,lp1),j=1,lp1)
ccc      read (itin1,1011) ((co2800(i,j),i=1,lp1),j=1,lp1)
ccc      read (itin,1011) ((co2po1(i,j),i=1,lp1),j=1,lp1)
ccc      read (itin1,1011) ((co2801(i,j),i=1,lp1),j=1,lp1)
ccc      read (itin,1011) ((co2po2(i,j),i=1,lp1),j=1,lp1)
ccc      read (itin1,1011) ((co2802(i,j),i=1,lp1),j=1,lp1)
      do 300 j=1,lp1
        do 300 i=1,lp1
          co2po(i,j) = t22(i,j,1)
cnov89
          if (iq.eq.5) go to 300
cnov89
          co2po1(i,j) = t22(i,j,2)
          co2po2(i,j) = t22(i,j,3)
  300 continue
      do 301 j=1,lp1
        do 301 i=1,lp1
          co2800(i,j) = t23(i,j,1)
cnov89
          if (iq.eq.5) go to 301
cnov89
          co2801(i,j) = t23(i,j,2)
          co2802(i,j) = t23(i,j,3)
  301 continue
c***the following code is rewritten so that the radiative bands
c   are:
c        iq=1    560-800     (consol.=490-850)
c        iq=2    560-670     (consol.=490-670)
c        iq=3    670-800     (consol.=670-850)
c        iq=4    560-760 (original code)   (consol.=490-850)
cnov89
c        iq=5   2270-2380    (consol.=2270-2380)
cnov89
c  the following loop obtains transmission functions for bands
c  used in radiative model calculations,with the equivalent
c  widths kept from the original consolidated co2 tf s.
cnov89
c      note: although the band transmission functions are
c  computed for all radiative bands, as of 9/28/88, they
c  are written out in full only for the full 15 um band cases
c  (iq=1,4).  in other cases, the transmissivities (1,k) are
c  written out, as these are the only ones needed for cts
c  calculations.  also, for the 4.3 um band (iq=5) the temp.
c  derivative terms are not written out, as they are unused.
cnov89
      if (iq.eq.1) then
         c1=1.5
         c2=0.5
      endif
      if (iq.eq.2) then
        c1=18./11.
        c2=7./11.
      endif
      if (iq.eq.3) then
        c1=18./13.
        c2=5./13.
      endif
      if (iq.eq.4) then
        c1=1.8
        c2=0.8
      endif
cnov89
      if (iq.eq.5) then
        c1=1.0
        c2=0.0
      endif
cnov89
      do 1021 i=1,lp1
      do 1021 j=1,lp1
      co2po(j,i)=c1*co2po(j,i)-c2
      co2800(j,i)=c1*co2800(j,i)-c2
cnov89
      if (iq.eq.5) go to 1021
cnov89
      co2po1(j,i)=c1*co2po1(j,i)-c2
      co2801(j,i)=c1*co2801(j,i)-c2
      co2po2(j,i)=c1*co2po2(j,i)-c2
      co2802(j,i)=c1*co2802(j,i)-c2
1021  continue
cnov89
      if (iq.ge.1.and.iq.le.4) then
cnov89
      do 1 j=1,lp1
      do 1 i=1,lp1
      dcdt8(i,j)=.02*(co2801(i,j)-co2802(i,j))*100.
      dcdt10(i,j)=.02*(co2po1(i,j)-co2po2(i,j))*100.
      d2ct8(i,j)=.0016*(co2801(i,j)+co2802(i,j)-2.*co2800(i,j))*1000.
      d2ct10(i,j)=.0016*(co2po1(i,j)+co2po2(i,j)-2.*co2po(i,j))*1000.
1     continue
cnov89
      endif
cnov89
co222 *********************************************************
ccc       rewind 66
c        save cdt51,co251,c2d51,cdt58,co258,c2d58..on tempo file
ccc       write (66) dcdt10
ccc       write (66) co2po
ccc       write (66) d2ct10
ccc       write (66) dcdt8
ccc       write (66) co2800
ccc       write (66) d2ct8
ccc       rewind 66
cnov89
      if (iq.eq.1.or.iq.eq.4) then
cnov89
      do 400 j=1,lp1
       do 400 i=1,lp1
        t66(i,j,1) = dcdt10(i,j)
        t66(i,j,2) = co2po(i,j)
        t66(i,j,3) = d2ct10(i,j)
        t66(i,j,4) = dcdt8(i,j)
        t66(i,j,5) = co2800(i,j)
        t66(i,j,6) = d2ct8(i,j)
  400 continue
cnov89
      else
      do 409 i=1,lp1
        t66(i,1,2) = co2po(1,i)
        t66(i,1,5) = co2800(1,i)
        if (iq.eq.5) go to 409
        t66(i,1,1) = dcdt10(1,i)
        t66(i,1,3) = d2ct10(1,i)
        t66(i,1,4) = dcdt8(1,i)
        t66(i,1,6) = d2ct8(1,i)
  409 continue
      endif
cnov89
co222 *********************************************************
      return
      end
