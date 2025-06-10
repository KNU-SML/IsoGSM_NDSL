#define PI pi
      subroutine albsno(imx2,latrue,jsno,
     1           albdoa,rlat,albedr,slmskr,ssnow,tgr,tar)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
#include <paramodel.h>
      dimension albdoa(imx2),albedr(imx2),slmskr(imx2)
      dimension ssnow(imx2),rlat(imx2)
c
c added by bob grumbine for sea ice albedo algorithm
c
      real tgr(imx2), tar(imx2)
c
c modified by hmh juang for simplicity  and for regional model use
c
      parameter(snodeg=70.* pi_ /180.)
c
c the following determines surface albedo (albdoa),where snow exists.
c
        snochk=pi_
        if(jsno.eq.0) snochk=snodeg  ! jsno=0 for rsm
c
        do 350 i=1, imx2
c
c limit background albedo (in case snow leaves greenland)
c
          albdoa(i)=min(albedr(i),0.6 e 0)
          if(slmskr(i).eq.1.0 e 0) then
c
c check lat<jsno for global and abs(rlat)>snodeg for regional
c
            if(latrue.lt.jsno .or. abs(rlat(i)).gt.snochk) then
              if(ssnow(i).gt.0. e 0) albdoa(i)=0.75 e 0
            else
              if(ssnow(i).ge.1. e 0) then
                albdoa(i)=0.6 e 0
              else if(ssnow(i).gt.0.0 e 0) then
                albdoa(i)=albdoa(i)+sqrt(ssnow(i))*(0.6 e 0-albdoa(i))
              endif
            endif
          else if(slmskr(i).eq.2.0 e 0) then
            if (ssnow(i) .gt. 0.0) then
              if (tgr(i) .lt. 273.16 - 5.) then
                albdoa(i) = 0.8
              else if (tgr(i) .le. 273.16) then
                albdoa(i) = 0.65 + 0.03*(273.16 - tgr(i))
              else
                albdoa(i) = 0.65
              endif
             else
              if (tgr(i) .lt. 271.2 .or. tar(i) .lt. 273.16) then
                albdoa(i) = 0.65
              else if (tar(i) .lt. 273.16+5.) then
                albdoa(i) = 0.65 - 0.04*(tar(i) -273.16)
              else
                albdoa(i) = 0.45
              endif
            endif
          endif
          albdoa(i) = max(albdoa(i),.06 e 0)
350     continue
      return
      end
