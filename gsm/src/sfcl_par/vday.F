      subroutine vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
#include <paramodel.h>
c
c julian day of the middle of each month
c
      dimension dayhf(13)
      data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,
     1           196.5,227.5,258.0,288.5,319.0,349.5,380.5/
      save dayhf
c
c number of days in a month
c
      dimension mjday(12)
      data mjday/31,28,31,30,31,30,31,31,30,31,30,31/
      save mjday
c
c julian day of the first day of the month
c
      dimension fjday(12)
c
      logical ijordr
c
      iret=0
c
      monend=9999
c
      fjday(1)=1.
      do mon=2,12
        fjday(mon)=fjday(mon-1)+float(mjday(mon-1))
      enddo
c
c  get julian day of the iy/im/id/ih provided.
c
      rjday=0.
      imm=im-1
      if(imm.ge.1) then
        do mon=1,imm
          rjday=rjday+mjday(mon)
        enddo
      endif
      rjday=rjday+id+float(ih)/24.0
      rjday=rjday+fh/24.0
      rjday = mod(rjday,365.)
      if(rjday.eq.0.) rjday = 365.
      if(rjday.le.0..or.rjday.gt.365.) then
        print *,'wrong rjday',rjday
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      endif
      do mon=1,11
        if(rjday.lt.fjday(mon+1)) then
          mfmon=mon
          go to 10
        endif
      enddo
      mfmon=12
   10 continue
      if(rjday.lt.dayhf(1)) rjday=rjday+365.
c
c  compute jy,jm,jd,jh of forecast
c
      jy=iy
      jm=im
      jd=id
      incdy=int(fh/24.)
      jh=ih+mod(fh,24.)
      incdy=incdy+jh/24
      jh=mod(jh,24)
      if(incdy.ge.1) then
      do incd=1,incdy
        jd=jd+1
        if(jm.eq.4.or.jm.eq.6.or.jm.eq.9.or.jm.eq.11) then
          mondy=30
        elseif(jm.eq.2) then
          if(mod(jy,4).eq.0) then
            mondy=29
c
c           Check 100/400 year leap
            if(mod(jy,100).eq.0.and.mod(jy,400).ne.0) then
               mondy=28
            endif
c
          else
            mondy=28
          endif
        else
          mondy=31
        endif
        if(jd.gt.mondy) then
          jm=jm+1
          jd=1
          if(jm.gt.12) then
            jy=jy+1
            jm=1
          endif
        endif
      enddo
      endif
c     write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
c
      return
      end
