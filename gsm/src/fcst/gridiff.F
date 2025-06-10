#include <define.h>
      subroutine gridiff(ims2,imx2,kmx,deltim,sl,pslap,gt0,gq0,ntotal)
c                                                                               
      dimension sl(imx2,kmx)                                                       
      dimension pslap(imx2) 
      dimension gt0(imx2,kmx)
      dimension gq0(imx2,kmx*ntotal)
c                                                                               
      do k=1,kmx                                                              
        kd=max(k-1,1)                                                           
        ku=min(k+1,kmx)                                                       
        do j=1,ims2                                                           
          gt0(j,k)=gt0(j,k)+pslap(j)*deltim*                                    
     1            (gt0(j,ku)-gt0(j,kd))*sl(j,k)/(sl(j,ku)-sl(j,kd))                   
        enddo                                                                   
      enddo                                                                     
c
c diffusion is available to q only now.
c if you want to use surface horizontal diffusion to cloud water, 
c use "do n = 1,ntotal"
c
#ifdef ISOTOPE
      do n = 1,ntotal
#else
      do n = 1,1
#endif
        ki = (n-1)*kmx + 1
        ke = n*kmx 
        do k=1,kmx                                                              
          kk = (n-1)*kmx + k
          kd=max(kk-1,ki)                                                           
          ku=min(kk+1,ke)                                                       
          kd1=kd-(n-1)*kmx
          ku1=ku-(n-1)*kmx
          do j=1,ims2                                                           
            gq0(j,kk)=gq0(j,kk)+pslap(j)*deltim*                                    
     1               (gq0(j,ku)-gq0(j,kd))*sl(j,k)/(sl(j,ku1)-sl(j,kd1))                   
          enddo                                                                   
        enddo                                                                     
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
