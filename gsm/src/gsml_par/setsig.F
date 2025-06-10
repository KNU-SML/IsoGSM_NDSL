#define CP cp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define RD rd
      subroutine setsig(ci, si, del, sl, cl, rpi, n1)                          
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    setsig      sets up model sigma structure.                     
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-05                     
c                                                                               
c abstract: sets up model sigma structure based on vertical                     
c   sigma spacing defined in the subroutine.                                    
c                                                                               
c program history log:                                                          
c   88-04-05  joseph sela                                                       
c                                                                               
c usage:    call setsig (ci, si, del, sl, cl, rpi)                              
c                                                                               
c   output argument list:                                                       
c     ci       - array of 1.0-si at each level.                                 
c     si       - array of sigma value at each level.                            
c     del      - array of sigma spacing at each layer.                          
c     sl       - array of sigma at midpoint of sigma layers.                    
c     cl       - array of 1.0-sl at each layer midpoint.                        
c     rpi      - array of pi ratios needed in thermodynamic equation.           
c                                                                               
c   output files:                                                               
c     output   - printout file.                                                 
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      real rk,rk1,rkr                                                           
      dimension ci(levp1_), si(levp1_),                                         
     1 del(levs_), sl(levs_), cl(levs_), rpi(levm1_)                            
      integer idate(4)                                                          
c                                                                               
#ifndef NOPRINT
      print 98, n1                                                             
98    format (1h0, 'begin setsig - getting sigs from unit',i4)                  
#endif
      close(n1)
      open (unit=n1,file='sigit ',form='unformatted',err=999)
      go to 998
  999 continue
      print *,'error opening sigit in getcon'
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
  998 continue
      rewind n1
      read(n1)                                                                 
      read(n1) fhour,idate,si,sl                                               
      rewind n1                                                                
      do 1 li=1,levp1_                                                          
1     ci(li) = 1.e0 - si(li)                                                   
      do 3 le=1,levs_                                                           
      cl(le) = 1.e0 - sl(le)                                                   
      del(le) = si(le) - si(le+1)                                               
3     continue                                                                  
c     compute pi ratios for temp. matrix.                                       
      rk = rd_/cp_                                                              
      do 4 le=1,levm1_                                                          
      base = sl(le+1)/sl(le)                                                    
4     rpi(le) = base**rk                                                        
#ifdef DBG
      do 5 le=1,levp1_                                                          
      print 100, le, ci(le), si(le)                                             
100   format (1h , 'level=', i2, 2x, 'ci=', f6.3, 2x, 'si=', f6.3)              
5     continue                                                                  
      print 97                                                                  
97    format (1h0)                                                              
      do 6 le=1,levs_                                                           
      print 101, le, cl(le), sl(le), del(le)                                    
101   format (1h , 'layer=', i2, 2x, 'cl=', f6.3, 2x, 'sl=', f6.3, 2x,          
     1 'del=', f6.3)                                                            
6     continue                                                                  
      print 102, (rpi(le), le=1,levm1_)                                         
102   format (1h0, 'rpi=', (18(1x,f6.3)) )                                      
#endif
      return                                                                    
      end                                                                       
