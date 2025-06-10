      subroutine nfdr_m100(jdate,pptdur,tmax,tmin,rhmax,rhmin,ij,
     *                bndry,fm100hr)
c
      implicit none
c
c  This subroutine calculates the average boundry conditions for the
c  past 24-hours and the 100-hr-tl fuel moisture. the boundry con-
c  ditions are weighted averages of the eqmcs calculated from the
c  db and rh extreme values. Philat is used to calculate
c  daylength which is the basis of the weighting function.
c  Philat is from Thornthwait via Bill Main at East Lansing, MI.
c
c  input
c
      integer jdate,ij
      real pptdur,tmax,tmin,rhmax,rhmin
c
c output
c
      real bndry,fm100hr
c
      real fr100,xfact,ymc100,decl,phi,tla,daylit,emc,emc1,emc2
      real bndry1
      integer lat
c
c  fr100 equals (1.0-.87*exp(-24./100.)) = 0.3156
c
      fr100=0.3156
c
c  Calculate declination for the daylength calculation.
c
      xfact = float(jdate-80)*0.01745
c
c  At this point, fm100hr is yesterdays value for the ij
c
      ymc100 = fm100hr
c
c  Determine the approximate latitude for this ij.
c
      lat=48
      if(ij .gt. 230 .and. ij .le. 460) lat=44
      if(ij .gt. 460 .and. ij .le. 690) lat=40
      if(ij .gt. 690 .and. ij .le. 920) lat=36
      if(ij .gt. 920 .and. ij .le. 1150) lat=31
      if(ij .gt. 1150) lat=27
c 
c  Convert degrees to radians by mult. by 0.01745
c
      decl=23.5*sin(xfact)
      decl=decl*.01745
      phi=-tan(float(lat)*.01745)
c
      tla=phi*sin(decl)
      if(abs(tla).lt.0.01)tla=0.01
c
      tla=atan(sqrt(1.0-tla*tla)/tla)
      if(tla.lt.0.0)tla=tla+3.141593
c
c  Hours of daylight (daylit).
c
      daylit = tla*7.64
c 
c  Float percip. duration is carried as pptdur.
c  pptdur = float(pptdur)
c 
c  Calculate eqmcs for 24-hr extremes
c
      call nfdr_eqmc(tmax,rhmin,emc1)
      call nfdr_eqmc(tmin,rhmax,emc2)
c 
c  Weighted average, emc
c
      emc=(daylit*emc1+(24.0-daylit)*emc2)/24.0
c 
c  24-hour boundry conditions,bndry, and the 100-hr-tl fuel moist.,
c  fm100.
c
      bndry=((24.-pptdur)*emc+(2.7*pptdur+76.)*pptdur)/24.
      bndry1=((24.-pptdur)*emc+(0.5*pptdur+41.)*
     +        pptdur)/24.0
      fm100hr = (bndry1-ymc100)*fr100+ymc100
c
c  now fm100hr is todays value for this ij.
c
      return
      end
