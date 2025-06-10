#define LATG latg
#define LEVS levs
#define LONF lonf
      subroutine omegas(dphi,dlam,ug,vg,dg,del,vvel,ps,sl)                      
#include <paramodel.h>
      save                                                                      
      parameter (idim=lonf_,jdim=latg_,                                         
     1           kdim=levs_,kdimp=kdim+1,kdimm=kdim-1)                          
c....   code lifted from post (mcp1840) jun 88--computes vvel (cb/sec)          
c....    input ps in cb,output vvel in cb/sec                                   
c....   do loops altered for better vectorization possibilities..k.a.c.         
      dimension dphi(idim),dlam(idim),                                          
     1               ug(idim,kdim),vg(idim,kdim),                               
     2 dg(idim,kdim),del(kdim),ps(idim),sl(kdim)                                
c...   vvel contains omega in layers on return from subroutine...               
      dimension vvel(idim,kdim)                                                 
      dimension cg(idim,kdim)                                                   
      dimension db(idim,kdim),cb(idim,kdim),dot(idim,kdimp)                     
      do 1 k=1,kdimp                                                            
        do 1 lo=1,idim                                                          
          dot(lo,k) = 0.                                                        
    1 continue                                                                  
c...  compute c=v(true)*del(ln(ps))                                             
      do 5 le=1,kdim                                                            
        do 4 lo=1,idim                                                          
          cg(lo,le)=ug(lo,le)*dlam(lo)+vg(lo,le)*dphi(lo)                       
    4   continue                                                                
    5 continue                                                                  
      do 10 lo=1,idim                                                           
        db(lo,1)=del(1)*dg(lo,1)                                                
        cb(lo,1)=del(1)*cg(lo,1)                                                
   10 continue                                                                  
      do 6 le=1,kdimm                                                           
        do 6 lo=1,idim                                                          
          db(lo,le+1)=db(lo,le)+del(le+1)*dg(lo,le+1)                           
          cb(lo,le+1)=cb(lo,le)+del(le+1)*cg(lo,le+1)                           
    6 continue                                                                  
c                                                                               
c...    sigma dot computed only at interior interfaces                          
c  this is negative sigma dot                                                   
c                                                                               
      do 7 k=1,kdimm                                                            
        do 7 lo=1,idim                                                          
          dot(lo,k+1)=dot(lo,k)+del(k)                                          
     1               *(db(lo,kdim)+cb(lo,kdim)-dg(lo,k)-cg(lo,k))               
    7 continue                                                                  
      do 8 k=1,kdim                                                             
        do 8 lo=1,idim                                                          
          vvel(lo,k)=  sl(k)*(cg(lo,k)-cb(lo,kdim)-db(lo,kdim))-                
     1                0.5*(dot(lo,k+1)+dot(lo,k))                               
          vvel(lo,k)=vvel(lo,k)*ps(lo)                                          
ccc       vvel(lo,k)=vvel(lo,k)*ps(lo)*10.                                      
    8 continue                                                                  
      return                                                                    
      end                                                                       
