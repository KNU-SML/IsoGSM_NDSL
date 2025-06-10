#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONFP lonfp
      subroutine flxco299(ipts,swc,u1,du,nu,swh,w1,dw,nw,tbl,daytm,dflx)                                 
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c********************************************************************           
c  compute the absorption due to co2. ref: chou (j. climate, 1990,              
c     209-217)                                                                  
c     updated sep. 1999 based on nasa/tm-1999-104606, vol 15.
c  the effect of co2 absorption below the cloud top is neglected.               
c  input variables:                                                             
c     swc,swh     : column amount of co2 and water vapor                        
c     u1,du,w1,dw : coefficients
c     tbl         : look up co2 absorption table
c     nu,nw       : table dimensions
c     daytm       : daytime flag                                                
c  output variables:                                                            
c     dflx        : flux reduction due to co2 for clear sky                     
c                                                                               
c********************************************************************           
c                                                                               
#include <paramodel.h>
#include <rdparm99.h>
                      d i m e n s i o n                                         
     1  swc(imbx,lp1),   swh(imbx,lp1)                               
     2, dflx(imbx,lp1),  tbl(nu,nw)                                
      logical daytm(imax)                                                       
c
c ... table look-up for the reduction of clear-sky solar
c
      x1 = u1 - 0.5*du
      y1 = w1 - 0.5*dw
      do k=2,lp1
      do i=1,ipts
       if (daytm(i)) then
        clog = swc(i,k)
        wlog = swh(i,k)
        ic = int( (clog - x1)/du + 1.0e0)
        iw = int( (wlog - y1)/dw + 1.0e0)
        ic = max(2, min(nu, ic))
        iw = max(2, min(nw, iw))
        ic1 = ic - 1
        iw1 = iw - 1
        dc = clog - float(ic-2)*du - u1
        dd = wlog - float(iw-2)*dw - w1
        x2 = tbl(ic1,iw1) + (tbl(ic1,iw)-tbl(ic1,iw1))/dw * dd
        y2 = x2 + (tbl(ic,iw1) - tbl(ic1,iw1))/du * dc

        dflx(i,k) = dflx(i,k) + y2
       end if
      end do
      end do
c                                                                               
      return                                                                    
      end                                                                       
