#define IO io
#define JO jo
      subroutine rowsep(a)                                                      
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
#include <paramodel.h>
      dimension a(2*io_,(jo_+1)/2)                                              
      dimension b(io_,jo_)                                                      
c                                                                               
      do j=1,(jo_+1)/2                                                          
        do i=1,io_                                                              
          b(i,j)=a(i,j)                                                         
          b(i,jo_-j+1)=a(io_+i,j)                                               
        enddo                                                                   
      enddo                                                                     
      do ij=1,io_*jo_                                                           
        a(ij,1)=b(ij,1)                                                         
      enddo                                                                     
      return                                                                    
      end                                                                       
