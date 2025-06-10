#include <define.h>
      subroutine gncpus(ncpus)                                                  
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram: gncpus         gets environment number of cpus                    
c   prgmmr: iredell          org: w/nmc23     date: 94-08-19                    
c                                                                               
c abstract: gets and returns the environment variable ncpus,                    
c   designating the number of processors over which to parallelize.             
c                                                                               
c program history log:                                                          
c   94-08-19  iredell                                                           
c                                                                               
c usage:    call gncpus(ncpus)                                                  
c   output arguments:                                                           
c     ncpus        integer number of cpus                                       
c                                                                               
c subprograms called:                                                           
c   getenv       get environment variable                                       
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      integer getenv                                                            
      character*8 cncpus                                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      ncpus=1                                                              
#define DEFAULT
#ifdef GETENV
#undef DEFAULT
      iret=getenv('ncpus',cncpus)                                        
      if(iret.eq.1) then                                                   
        read(cncpus,'(bn,i8)',iostat=ios) ncpus                            
        ncpus=max(ncpus,1)                                                 
        print *,'ncpus=',ncpus                                             
      endif                                                                
#endif
#ifdef PXFGETENV
      iret=0                                                             
      call pxfgetenv('ncpus',5,cncpus,linval,jret)                       
      if(jret.eq.0) iret=1                                               
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
