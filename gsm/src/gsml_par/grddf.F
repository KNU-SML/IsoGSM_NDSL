#define NSTYPE nstype
      subroutine grddf                                                          
#include <paramodel.h>
      parameter(ntype=nstype_,ngrid=22)
      common /comgdf/ dfk(ngrid,ntype)                                          
      common /comgkt/ ktk(ngrid,ntype)                                          
      real ktk                                                                  
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)         
      do k = 1, ntype                                                           
        dynw = tsat(k) * .05                                                    
        f1 = b(k) * satkt(k) * satpsi(k) / tsat(k) ** (b(k) + 3.)               
        f2 = satkt(k) / tsat(k) ** (b(k) * 2. + 3.)                             
c                                                                               
c  convert from m/s to kg m-2 s-1 unit                                          
c                                                                               
        f1 = f1 * 1000.                                                         
        f2 = f2 * 1000.                                                         
        do i = 1, ngrid                                                         
          theta = float(i-1) * dynw                                             
          theta = min(tsat(k),theta)                                            
          dfk(i,k) = f1 * theta ** (b(k) + 2.)                                  
          ktk(i,k) = f2 * theta ** (b(k) * 2. + 3.)                             
        enddo                                                                   
      enddo                                                                     
      return                                                                    
      end                                                                       
