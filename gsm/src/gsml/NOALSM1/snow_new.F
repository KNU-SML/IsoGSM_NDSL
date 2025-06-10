      subroutine snow_new (temp,newsn,snowh,sndens)
c
      implicit none
c      
c ----------------------------------------------------------------------
c subroutine snow_new
c ----------------------------------------------------------------------
c calculate snow depth and densitity to account for the new snowfall.
c new values of snow depth & density returned.
c
c temp    air temperature (k)
c newsn   new snowfall (m)
c snowh   snow depth (m)
c sndens  snow density (g/cm3=dimensionless fraction of h2o density)
c ----------------------------------------------------------------------
      real sndens
      real dsnew
      real snowhc
      real hnewc
      real snowh
      real newsn
      real newsnc
      real temp 
      real tempc
c ----------------------------------------------------------------------
c conversion into simulation units      
c ----------------------------------------------------------------------
      snowhc = snowh*100.
      newsnc = newsn*100.
      tempc = temp-273.15
c ----------------------------------------------------------------------
c calculating new snowfall density depending on temperature
c equation from gottlib l. 'a general runoff model for snowcovered
c and glacierized basin', 6th nordic hydrological conference,
c vemadolen, sweden, 1980, 172-177pp.
c-----------------------------------------------------------------------
      if (tempc .le. -15.) then
        dsnew = 0.05
      else                                                      
        dsnew = 0.05+0.0017*(tempc+15.)**1.5
      endif
c ----------------------------------------------------------------------
c adjustment of snow density depending on new snowfall      
c ----------------------------------------------------------------------
      hnewc = newsnc/dsnew
      sndens = (snowhc*sndens+hnewc*dsnew)/(snowhc+hnewc)
      snowhc = snowhc+hnewc
      snowh = snowhc*0.01
c ----------------------------------------------------------------------
c end subroutine snow_new
c ----------------------------------------------------------------------
      return
      end
