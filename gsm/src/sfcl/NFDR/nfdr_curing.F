      subroutine nfdr_curing(fm1,annual,fm1000,wherb,fherbc,
     &                       fmwood,w1dp,wherbc)
c
      implicit none
c
c  fm1000 is the 1000 hour fuel moisture in percent
c  w1d is the 1 hour fuel load from the fuel model
c  wherb is the herbaceous fuel load from the fuel model
c  fherbc is the live herbaceous fuel moisture in percent
c  fmwood is the live woody fuel moisture in percent
c  w1dp is fuel model 1 hour load plus cured live herbaceous load
c  wherbc is remaining live load
c  fcc is the fraction of the live herbaceous load that is cured
c  Live herb load varies from 30-250% as fm1000 varies from 15-25%
c  Live woody load varies from 80-200% as fm1000 varies from 15-25%
c  15% 1000 hr moisture is estimated to be about a minimum value for Hawaii
c  25% 1000 hr moisture is standard NFDRS for setting max live moistures
c
      real fm1,annual,fm1000,wherb,fherbc,fmwood,w1dp,wherbc
c
      real fcc
      integer ifccpct
c 
c  Set the herb and live woody fuel moisture in percent for...
c  Climate class 1 (mean annual ppt <= 16 inches
c
        if(annual .le. 16.0)then
          fherbc=(fm1000-7.0)*9.5652 + 30.0
          fmwood=(fm1000-7.0)*5.2174 + 80.0
        endif
c      
c  Climate class 2 (16.0 <mean annual ppt <= 29 inches      
c
        if(annual .gt. 16.0 .and. annual .le. 29.0)then
          fherbc=(fm1000-9.0)*10.4762 + 30.0
          fmwood=(fm1000-9.0)*5.7143 + 80.0
        endif
c    
c  Climate class 3 (29.0 <mean annual ppt <= 69 inches       
c
        if(annual .gt. 29.0 .and. annual .le. 69.0)then
          fherbc=(fm1000-12.0)*12.2222 + 30.0
          fmwood=(fm1000-12.0)*6.6667 + 80.0
        endif
c    
c  Climate class 4 (mean annual ppt > 69 inches)     
c
        if(annual .gt.69.0)then
          fherbc=(fm1000-15.0)*14.6667 + 30.0
          fmwood=(fm1000-15.0)*8.0 + 30.0
        endif
c             
      if(fherbc.lt.30. ) fherbc=fm1
      if(fherbc.gt.250.) fherbc=250.
      if(fmwood.gt.200.) fmwood=200.0
      if(fmwood.lt.80.0) fmwood=80.0
c           
c  Do the load transfer as herb moisture ranges from 30% to 120%
c
      fcc=-0.0111*fherbc+1.33
      if(fcc.gt.1.)fcc=1.
      if(fcc.lt.0.)fcc=0.
      ifccpct = nint(fcc*100.0)
      w1dp=w1dp+wherb*fcc
      wherbc=wherb*(1.-fcc)
c 
      return
      end
