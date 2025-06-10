      subroutine fwindx(im,jm,u10m,v10m,t2m,q2m,ps,fwi)
c
c  compute equilibrium moisture content emc
c  and fire weather index fwi following fosberg (1978) 
c
c  variables and units used in computation of fwi (3 mar 1990)
c
c  w   = windspeed in miles per hour
c  tfahr = temperature in degrees fahrenheit
c  rh    = relative humidity in percent
c
c  units of input variables
c
c  u10m,v10m = wind speed meters per second
c  t2m     = degrees kelvin
c  q2m     = gm/gm
c  hm    = percent
c  psfc  = kPascal
c
      implicit none
      integer im,jm
      real u10m(im,jm),v10m(im,jm),t2m(im,jm),q2m(im,jm)
      real ps(im,jm)
      real fwi(im,jm)
c
      real w(im,jm),hm(im,jm),tfahr(im,jm),emc(im,jm)
c
      real rd,rv,eps,epsm1
      parameter(rd= 2.8705e+2 ,rv= 4.6150e+2 ,eps=rd/rv,epsm1=rd/rv-1.)
c
      real fpvs
      real es, qs, a, eta 
      integer i,j
c
c  get 10m wind speed
c
      do j=1,jm
        do i=1,im
           w(i,j)=sqrt(u10m(i,j)*u10m(i,j)+v10m(i,j)*v10m(i,j))
        enddo
      enddo
c
c  wind speed in miles per hour
c
      do j=1,jm
        do i=1,im
          w(i,j)=w(i,j)*3600./1609.3
        enddo
      enddo
c
c  get 2m relative humidity.  always use saturtion over ice if t2m<273.16
c
      do j=1,jm
        do i=1,im
          es=fpvs(t2m(i,j))
          qs=eps*es/(ps(i,j)+epsm1*es)
          hm(i,j)=min(max(q2m(i,j)/qs,0.),1.)*100.
        enddo
      enddo
c
c  temperature in fahrenheit
c
      do j=1,jm
        do i=1,im
          tfahr(i,j)=(t2m(i,j)-273.16)*1.80+32.0
        enddo
      enddo
c
c  compute emc
c
      do j=1,jm
        do i=1,im
          if(hm(i,j).lt.10.)then
            emc(i,j)=0.03229+(0.281073-0.000578*tfahr(i,j))*hm(i,j)
          else
            if(hm(i,j).le.50.)then
              emc(i,j)=2.22749+0.160107*hm(i,j)-0.014784*tfahr(i,j)
            else
              emc(i,j)=21.0606+(0.005565*hm(i,j)-0.00035*tfahr(i,j)-
     1                 0.483199)*hm(i,j)
            endif
          endif
        enddo
      enddo
c
c  compute fire weather index
c
      do j=1,jm
        do i=1,im
          a=emc(i,j)/30.0
          eta=1.+(-2.+(1.5-0.5*a)*a)*a
          fwi(i,j)=eta*sqrt(1.+w(i,j)*w(i,j))/0.3002
        enddo
      enddo
c
      return
      end
