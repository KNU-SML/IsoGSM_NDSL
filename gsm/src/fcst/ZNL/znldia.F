#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define PI pi

      subroutine znldia(nznl,rhour,idate,kdt,latb2,levs,                        
     &                  wgb,colrab,del)                                         
c                                                                               
#include <paramodel.h>
      dimension idate(4)                                                        
      dimension wgb(latb2),colrab(latb2),del(levs)                              
c                                                                               
#include <comznl.h>
c                                                                               
#ifdef DYNAMIC_ALLOC
      dimension znlm(nlb,levs+1,nrm)                                            
#else
      parameter(nlevp1=levs_+1)                                                 
      dimension znlm(nlb,nlevp1,nrm)                                            
#endif
                                                                                
      dimension znls(nlb,nst,nrs)                                               
      dimension weim(nlb),weis(nlb,nst)                                         
c.......................................................................        
      do n=1,nlb                                                                
      weim(n)=0.                                                                
      do m=1,nst                                                                
      weis(n,m)=0.                                                              
      enddo                                                                     
      do l=1,levs+1                                                             
      do m=1,nrm                                                                
      znlm(n,l,m)=0.                                                            
      enddo                                                                     
      enddo                                                                     
      do m=1,nst                                                                
      do k=1,nrs                                                                
      znls(n,m,k)=0.                                                            
      enddo                                                                     
      enddo                                                                     
      enddo                                                                     
c                                                                               
      do lat=1,latb2                                                            
        nlbh=nlb/2                                                              
        wlb=pi_/(2*nlbh)                                                        
        jb=colrab(lat)/wlb+1                                                    
        jn=jb+1                                                                 
        js=nlb-jb+1                                                             
c                                                                               
        weim(jn)=weim(jn)+zwm(1,lat)*wgb(lat)                                   
        weim(js)=weim(js)+zwm(2,lat)*wgb(lat)                                   
        do n=1,nrm                                                              
          if(zhm(n).gt.0.) then                                                 
            w=wgb(lat)*zhm(n)                                                   
            do k=1,levs                                                         
              znlm(jn,k,n)=znlm(jn,k,n)+zdm(1,k,n,lat)*w                        
              znlm(js,k,n)=znlm(js,k,n)+zdm(2,k,n,lat)*w                        
            enddo                                                               
          endif                                                                 
        enddo                                                                   
c                                                                               
        do k=1,nst                                                              
          weis(jn,k)=weis(jn,k)+zws(1,k,lat)*wgb(lat)                           
          weis(js,k)=weis(js,k)+zws(2,k,lat)*wgb(lat)                           
        enddo                                                                   
        do n=1,nrs                                                              
          do k=1,nst                                                            
            znls(jn,k,n)=znls(jn,k,n)+zds(1,k,n,lat)*wgb(lat)                   
            znls(js,k,n)=znls(js,k,n)+zds(2,k,n,lat)*wgb(lat)                   
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
c                                                                               
      do j=2,nlb                                                                
        weim(1)=weim(1)+weim(j)                                                 
      enddo                                                                     
      do k=1,levs                                                               
        do j=2,nlb                                                              
          do n=1,nrm                                                            
            znlm(1,k,n)=znlm(1,k,n)+znlm(j,k,n)                                 
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
      do j=1,nlb                                                                
        do k=1,levs                                                             
          do n=1,nrm                                                            
            znlm(j,levs+1,n)=znlm(j,levs+1,n)+znlm(j,k,n)*del(k)                
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
      do k=1,levs+1                                                             
        do j=1,nlb                                                              
          do n=1,nrm                                                            
            if(weim(j).ne.0.) znlm(j,k,n)=znlm(j,k,n)/weim(j)                   
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
c                                                                               
      do j=2,nlb                                                                
        do k=1,nst                                                              
          weis(1,k)=weis(1,k)+weis(j,k)                                         
        enddo                                                                   
      enddo                                                                     
      do k=1,nst                                                                
        do j=2,nlb                                                              
          do n=1,nrs                                                            
            znls(1,k,n)=znls(1,k,n)+znls(j,k,n)                                 
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
      do k=1,nst                                                                
        do j=1,nlb                                                              
          do n=1,nrs                                                            
            if(weis(j,k).ne.0.) znls(j,k,n)=znls(j,k,n)/weis(j,k)               
          enddo                                                                 
        enddo                                                                   
      enddo                                                                     
c                                                                               
      do j=1,nlb                                                                
        if(weis(1,1).ne.0.) znls(j,1,nsslmsk)=100.*weis(j,1)/weis(1,1)          
      enddo                                                                     
      do k=2,nst                                                                
        do j=1,nlb                                                              
          if(weis(j,1).ne.0.) znls(j,k,nsslmsk)=100.*weis(j,k)/weis(j,1)        
        enddo                                                                   
      enddo                                                                     
c                                                                               
      call znlprt(nznl,rhour,idate,kdt,levs+1,znlm,znls)                        
c                                                                               
      return                                                                    
      end                                                                       
