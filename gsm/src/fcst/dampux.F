#define JCAP jcap
#define LEVH levh
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define NCLDG ncldg
#define NTOTAL ntotal
#define NTRAC ntrac
#define RERTH rerth
      subroutine dampux(div,vor,tem,rmx,deltim,fdamp,arn,spdmax,
     1                  snnp1,lnts2,lnoffset)                
#include <paramodel.h>
#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
      dimension fdamp(LNT22S),arn(LNT22S),spdmax(levs_), 
     1 rmx(LNT22S,levh_),snnp1(lnt2_),
     3 vor(LNT22S,levs_), div(LNT22S,levs_),tem(LNT22S,levs_)
c.................................................................              
      do j=1,lnts2                                                          
        jj=lnoffset+j
        arn(j)=snnp1(jj)+0.25e0                                                   
        arn(j)= sqrt(arn(j))                                                      
        arn(j)=arn(j)-0.5e0                                                      
      enddo
      alfa=2.5e0                                                               
      beta=rerth_*1.009e0/deltim                                               
      alfadt=alfa*deltim/rerth_                                                 
c.................................................                              
      do k=1,levs_                                                           
        rncrit=beta/spdmax(k)                                                     
        coef=alfadt*spdmax(k)                                                     
        do j=1,lnts2                                                          
          if (arn(j).gt.rncrit) then                                                
            div(j,k) =div(j,k)/(1.+(arn(j)-rncrit)*coef)                              
            vor(j,k) =vor(j,k)/(1.+(arn(j)-rncrit)*coef)                              
            tem(j,k) =tem(j,k)/(1.+(arn(j)-rncrit)*coef)                              
          end if                                                                    
        enddo
#ifndef NISLQ
        do it=1,ntotal_
          is=(it-1)*levs_                                            
          l=is+k                                                    
            do j=1,lnts2                                           
              if (arn(j).gt.rncrit) then                          
                rmx(j,l) =rmx(j,l)/(1.+(arn(j)-rncrit)*coef)     
              end if                                            
            enddo
        enddo
#endif
      enddo
      return                                                   
      end                                                     
