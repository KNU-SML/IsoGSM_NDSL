#define IO io
#define JCAP jcap
#define JO jo
#define KO ko
#define LEVS levs
#include <define.h>
      program pgb                                                               
c$$$  main program documentation block                                          
c                                                                               
c main program:  pgb         transform sigma to pressure grib                   
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: program transforms sigma input to pressure grib output.             
c   the output consists of data on a regular lat/lon grid.                      
c   geopotential height, wind components, relative humidity,                    
c   temperature and vertical velocity are output on mandatory pressures.        
c   also output are sundry fields consisting of                                 
c   precipitable water, three lower level relative humidities,                  
c   lower level potential temperature and wind components,                      
c   surface temperature, pressure, omega and relative humidity,                 
c   and tropopause temperature, pressure, wind components and shear.            
c   first nampgb namelist is read to determine output format.                   
c   then a sigma (grid or spectral) file is read from unit 11 and               
c   the program produces and writes a pressure grib1 file to unit 51.           
c   then a sigma file is read from unit 12 and                                  
c   the program produces and writes a pressure grib1 file to unit 52.           
c   the program continues until an empty input file is encountered.             
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c   01-03-31  juang    add w3tag for nco implementation
c                                                                               
c namelists:                                                                    
c   nampgb:      parameters determining output format                           
c     io         number of longitude points                                     
c     jo         number of latitude points                                      
c     ko         number of pressure levels                                      
c     po(ko)     pressures in mb (default: mandatory levels)                    
c     ncpus      number of parallel processes (default: 1)                      
c     mxbit      maximum number of bits to pack data (default: 16)              
c     ids(255)   decimal scaling of packed data                                 
c                (default: set by subprogram idsdef)                            
c     pot(255)   highest pressure in mb to output data                          
c                as a function of parameter indicator                           
c                (default: 300 for rh, 100 for omega, 0 otherwise)              
c     icen       forecast center identifier (default: 7)                        
c     icen2      forecast sub-center identifier (default: 0)                    
c     igen       model generating code (default: from sigma file)               
c                                                                               
c input files:                                                                  
c   unit   11-?  sigma file(s)                                                  
c                                                                               
c output files:                                                                 
c   unit   51-?  pressure grib1 file(s)                                         
c                                                                               
c subprograms called:                                                           
c   idsdef       set defaults for decimal scaling                               
c   gpvs         compute saturated vapor pressure table                         
c   gtdp         compute dewpoint temerature table                              
c   gthe         compute equivalent potential temperature table                 
c   gtma         compute moist adiabat table                                    
c   rdsgh        read a sigma file header                                       
c   pgb1         transform one sigma file to pressure grib                      
c                                                                               
c attributes:                                                                   
c   language: fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(levmax=100,komax=100)                                           
      dimension idate(4)                                                        
      dimension si(levmax+1),sl(levmax)                                         
      dimension po(komax)                                                       
      dimension ids(255),pot(255)                                               
c                                                                               
      namelist/nampgb/ ids,pot,icen,icen2,igen,fhs,fhe,fhinc
      parameter(io=io_, jo=jo_, ko=ko_, mxbit=16)                               
      parameter(io2=2*io_, io22=2*io_+6,johf=(jo_+1)/2)                         
      parameter(jcap=jcap_,levs=levs_)                                          
      parameter(nc=(jcap_+1)*(jcap_+2)+1,nctop=(jcap_+1)*2)                     
c                                                                               
#include <postplevs.h>
c                                                                               
      data ncpus/1/                                                             
c                                                                               
      data ids/255*0/,pot/255*0./                                               
      data icen/7/,icen2/1/,igen/0/                                             
      data fhs/0./,fhe/120./,fhinc/24./
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_pgb',2001,0000,0000,'np51   ')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c                                                                               
c  set defaults and read namelist                                               
c                                                                               
      call idsdef(2,ids)                                                        
      pot(33)=po(ko)                                                            
      pot(34)=po(ko)                                                            
      pot( 7)=po(ko)                                                            
      pot(11)=po(ko)                                                            
      pot(52)=100.                                                              
      pot(39)=100.                                                              
      pot(152)=100.                                                              
      pot(153)=100.                                                              
      pot(154)=po(ko)                                                              
      read(*,nampgb,end=5)                                                      
5     continue                                                                  
      call gpvs                                                                 
      call gtdp                                                                 
      call gthe                                                                 
      call gtma                                                                 
      nsig=11                                                                   
      npgb=50                                                                   
c
      fh=fhs
      dowhile(fh.le.fhe)
c                                                                               
c  read sigma header record                                               
c                                                                               
        call rdsgh(nsig,fh,fhour,idate,si,sl,iret,ntrace,ncloud)
        if(iret.ne.0) then
          print *,'read sigma file error'
          call abort
        endif
c                                                                               
c  transform to pressure grib and attempt to read next sigma header             
c                                                                               
        iymdh=idate(4)*1000000+idate(2)*10000+idate(3)*100+idate(1)             
        print *,' posting date ',iymdh,'+',nint(fhour),                         
     &          '   sigma spectral t',jcap,' l',levs,                           
     &          '   pressure grid ',io,'x',jo,'x',ko                            
        call pgb1(fhour,idate,nsig,si,sl,po,
     &            npgb,ncpus,ids,pot,icen,icen2,igen,ntrace,ncloud)                           
        close(nsig)                                                             
        close(npgb)                                                             
        fh=fh+nint(fhinc)
      enddo                                                                     
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tage('clim_pgb')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop                                                                      
      end                                                                       
