#define LATG latg
#define LONF lonf
      program cltune                                                            
c$$$  main program documentation block
c
c main program:  cltune      cloud tuning program
c   prgmmr: cpc          org: w/np51     date: 00-00-00
c
c abstract: 
c
c program history log:
c   01-03-31  
c
c namelists:
c   namtune:     
c	fnsig	sigma file unit
c	vvcldin	data index
c	ieee	ieee (1) or not (0)
c
c input files:
c   unit   10	
c
c output files:

c subprograms called:
c   getrh       get relative humidity
c   rhfreq      request relative humidity
c
c attributes:
c   language: fortran
c
c$$$
c                                                                               
#include <paramodel.h>
      save                                                                      
c                                                                               
      parameter(idim=lonf_,jdim=latg_,mcld=3)                                   
c                                                                               
      parameter(maxfile=500)                                                    
c                                                                               
      dimension slmask(idim,jdim)                                               
      dimension idate(4)                                                        
      dimension rhmdl(idim,jdim,mcld)                                           
c                                                                               
      character*120 fnsig(maxfile)                                              
c                                                                               
      namelist/namtune/fnsig,vvcldin,ieee                                       
c                                                                               
      data fnsig/maxfile*'          '/                                          
      data vvcldin/-9.99/                                                       
      data ieee/1/                                                              
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_cldtune',2001,0000,0000,'np51')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      read (5,namtune)                                                          
c                                                                               
      do n=1,maxfile                                                            
        if(fnsig(n)(1:10).eq.'          ') then                                 
          go to 10                                                              
        endif                                                                   
      enddo                                                                     
  10  continue                                                                  
      nfiles=n-1                                                                
      print *,'nfiles=',nfiles                                                  
      if(nfiles.le.0) then                                                      
        call abort                                                              
      endif                                                                     
c                                                                               
      rewind 10                                                                 
      read (10) slmask                                                          
c                                                                               
      do n=1,nfiles                                                             
        call getrh(fnsig(n),vvcldin,fhour,idate,rhmdl,ieee)                     
        call rhfreq(slmask,n,nfiles,fhour,idate,rhmdl)                          
      enddo                                                                     
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL w3tage('clim_cldtune')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop                                                                      
      end                                                                       
