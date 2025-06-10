      subroutine nfdr_kbdi(ppt,annual,cumppt,tmax,kdi)
c
      implicit none
      integer kdi 
      real ppt,annual,cumppt,tmax
c
      integer net,idq
      real pptnet,xkdi
c
c     This subroutine calculates the daily values of the Keetch-
c     Byram Drought Index as described in For. Serv. Pub. SE-38.
c     patterned after Bill Mains sbr "Keetch" dated 6/27/74
c
c     ppt is today's rainfall amount in inches
c     annual is average annual ppt (rain + water equivalent of snow)
c     cumppt is cumulative ppt
c     tmax is today's maximum temp in deg f
c     kdi is the keetch-byram d.i.
c     daylit is the hours of daylight - see entry m100 in sbr mi100
c
      pptnet=0.
      net=0
c
c     If no precip today skip to stmt 3
c
      if(ppt.le.0.0) go to 3
c 
c     Some precip today - check yesterdays precip amt
c
      if(cumppt.gt.0.20) go to 1
c
c     Still not enough
c
      cumppt=cumppt+ppt
c 
c     If we still dont have more than .20 inches of ppt, skip to stmt 4
c
      if(cumppt.le.0.20) go to 4
c 
c     Now we have enough, get net amt and adjust
c
      pptnet=cumppt-0.20
      go to 2
c 
c     All of todays rain counts as net ppt to reduce the kb-di
c
    1 continue
      pptnet=ppt
c 
c     Adjust yesterdays kdi with todays net rain, holding kdi >= 0
c
    2 continue
      net=(100.0*pptnet)+0.0005
      kdi=kdi-net
      if(kdi.lt.0)kdi=0
      go to 4
c 
c     With no ppt, set yesterdays ppt to zero
c
    3 continue
      cumppt=0.0
c 
c     Drying only takes place on days warmer than 50 deg f
c
    4 continue
      if(tmax.lt.50.) return
c 
c     We are now ready to calculate the drying factor
c
      xkdi=real(kdi)
      idq=(800.-xkdi)*(0.9676*exp(0.0486*tmax)-8.299)*0.001/
     *    (1.0+10.88*exp(-0.04409*annual))+0.5
c 
c     Update the drought index
c
      kdi=kdi+idq
c 
      return
      end
