#define CP cp
#define LEVS levs
#define RD rd
      subroutine bnmc(kmx,si,sl,cd,ci,cq,cql)                                   
#include <paramodel.h>
      parameter(rd=rd_,cp=cp_,rocp=rd/cp)                                       
      dimension si(kmx+1),sl(kmx)                                               
      dimension cd(kmx,kmx+1),ci(kmx+1,kmx),cq(kmx+1,kmx),cql(kmx,kmx)          
c local array                                                                   
      dimension del(levs_),slk(levs_),alfa(2:levs_),beta(levs_-1)               
      km=kmx                                                                    
      do k=1,km                                                                 
        del(k)=si(k+1)-si(k)                                                    
        slk(k)=sl(k)**rocp                                                      
      enddo                                                                     
      do k=2,km                                                                 
        alfa(k)=0.5*(1.-slk(k-1)/slk(k))/rocp                                   
      enddo                                                                     
      do k=1,km-1                                                               
        beta(k)=0.5*(slk(k+1)/slk(k)-1.)/rocp                                   
      enddo                                                                     
                                                                                
      do kd=1,km                                                                
        do ki=1,km+1                                                            
          cd(kd,ki)=0                                                           
        enddo                                                                   
      enddo                                                                     
      do kd=1,km                                                                
        cd(kd,kd)=-1/del(kd)                                                    
        cd(kd,kd+1)=1/del(kd)                                                   
      enddo                                                                     
                                                                                
      do ki=1,km+1                                                              
        do kd=1,km                                                              
          ci(ki,kd)=0                                                           
        enddo                                                                   
      enddo                                                                     
      do ki=2,km                                                                
        ci(ki,ki-1)=0.5                                                         
        ci(ki,ki)=0.5                                                           
      enddo                                                                     
                                                                                
      do ki=1,km+1                                                              
        do kd=ki,km                                                             
          cq(ki,kd)=0                                                           
        enddo                                                                   
        do kd=1,ki-1                                                            
          cq(ki,kd)=del(kd)                                                     
        enddo                                                                   
      enddo                                                                     
                                                                                
      cql(1,1)=del(1)-beta(1)*si(2)                                             
      do kt=2,km-1                                                              
        cql(1,kt)=del(kt)-alfa(kt)*si(kt)-beta(kt)*si(kt+1)                     
      enddo                                                                     
      cql(1,km)=del(km)-alfa(km)*si(km)                                         
      do kz=2,km                                                                
        cql(kz,1)=cql(1,1)+beta(1)                                              
        do kt=2,kz-1                                                            
          cql(kz,kt)=cql(1,kt)+alfa(kt)+beta(kt)                                
        enddo                                                                   
        cql(kz,kz)=cql(1,kz)+alfa(kz)                                           
        do kt=kz+1,km                                                           
          cql(kz,kt)=cql(1,kt)                                                  
        enddo                                                                   
      enddo                                                                     
                                                                                
      return                                                                    
      end                                                                       
