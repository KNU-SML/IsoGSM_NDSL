#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs

      subroutine znlprt(nznl,rhour,idate,kdt,nlv,znlm,znls)                     
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    znlprt      print zonal diagnostics.                           
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: print zonal diagnostics.                                            
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell                                                      
c                                                                               
c usage:    call znlprt (rhour,idate,kdt,nlb,nst,nlv,                           
c    &            nrd,nrm,nrs,cnd,cnm,cns,isd,ism,iss,znld,znlm,znls)           
c   input argument list:                                                        
c     rhour    - current forecast hour                                          
c     idate    - initial hour,month,day,year                                    
c     kdt      - forecast step if positive or -100-diab.ini.step or             
c                -4 for initial-dt data in diab.ini                             
c                -3 for initial    data in diab.ini                             
c                -2 for initial-dt data                                         
c                -1 for initial    data                                         
c                0  for data after initialization                               
c     nlb      - number of latitude bands                                       
c     nst      - number of surface types                                        
c     nlv      - number of levels                                               
c     nrd      - number of 3d dynamics fields                                   
c     nrm      - number of 3d physics fields                                    
c     nrs      - number of 2d physics fields                                    
c     cnd      - character*8 identifications of 3d dynamics fields              
c     cnm      - character*8 identifications of 3d physics fields               
c     cns      - character*8 identifications of 2d physics fields               
c     isd      - status flags of 3d dynamics fields                             
c     ism      - status flags of 3d physics fields                              
c     iss      - status flags of 2d physics fields                              
c     znld     - zonal diagnostic array of 3d dynamics fields                   
c     znlm     - zonal diagnostic array of 3d physics fields                    
c     znls     - zonal diagnostic array of 2d physics fields                    
c                                                                               
c   subprograms called:                                                         
c     ipwrcn   - function to determine order of magnitude of output             
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
#include <paramodel.h>
#include <comznl.h>
      integer idate(4)                                                          
      character*8 cnm,cns                                                       
      dimension cnm(nrm),cns(nrs)                                               
      dimension znlm(nlb,nlv,nrm)                                               
      dimension znls(nlb,nst,nrs)                                               
#ifdef DYNAMIC_ALLOC
      dimension wrkm(nlb,nlv),wrks(nlb,nst)                                     
#else
      parameter(nnlv=levs_+1)                                                   
      dimension wrkm(nlb,nnlv),wrks(nlb,nst)                                    
#endif
      character*32 labz                                                         
      character*8 clb(6)                                                        
      character*4 cst(6)                                                        
      character*80 cfmt                                                         
      data cnm/'u       ','v       ','tv      ','q       ','vot**2  ',          
     &         'div**2  ','omega   ','t       ','rh      ','ke      ',          
     &         'dtconv  ','dtlarg  ','dtshal  ','dtvrdf  ','dqconv  ',          
     &         'dqshal  ','dqvrdf  ','duvrdf  ','dvvrdf  ','dthsw   ',          
     &         'dthlw   ','cloud   ','cvcld   '/                                
      data cns/'rain    ','rainc   ','dtsfc   ','dqsfc   ','dusfc   ',          
     &         'dvsfc   ','rcov    ','rcovc   ','ps      ','tskin   ',          
     &         'wetness ','snow    ','tg1     ','tg2     ','tg3     ',          
     &         'sfcsw   ','sfclw   ','rhs     ','tvs     ','ts      ',          
     &         'qs      ','zorl    ','slmsk   ','dugwd   ','dvgwd   ',          
     &         'duasfc  ','duagwd  ','duamtn  ','ua      ','uap     ',          
     &         'ep      ','cldwork '/                                           
      data clb/'90n-90s ','90n-60n ','60n-30n ','30n-30s ',                     
     &                    '30s-60s ','60s-90s '/                                
      data cst/'mean',' lnd','slnd',' ice','sice',' sea'/                       
      data ncol/12/                                                             
c                                                                               
      print 900,kdt                                                             
c                                                                               
      do 40 n=1,nrm                                                             
      ipwr=ipwrcn(cnm(n))                                                       
      print 910,cnm(n),ipwr,rhour,idate,kdt                                     
      nrow=(nlv-1)/ncol+1                                                       
      do 30 krow=1,nrow                                                         
      k1=(krow-1)*ncol+1                                                        
      if(krow.lt.nrow) then                                                     
        k2=k1-1+ncol                                                            
        write(cfmt,921) ncol                                                    
        print cfmt,(k,k=k1,k2)                                                  
      elseif(k1.lt.nlv) then                                                    
        k2=nlv                                                                  
        write(cfmt,922) k2-k1                                                   
        print cfmt,(k,k=k1,k2-1)                                                
      else                                                                      
        k2=nlv                                                                  
        write(cfmt,923)                                                         
        print cfmt                                                              
      endif                                                                     
      write(cfmt,930) -ipwr,k2-k1+1                                             
      print cfmt,(clb(j),(znlm(j,k,n),k=k1,k2),j=1,nlb)                         
30    continue                                                                  
40    continue                                                                  
c                                                                               
      do 50 n=1,nrs                                                             
      ipwr=ipwrcn(cns(n))                                                       
      print 910,cns(n),ipwr,rhour,idate,kdt                                     
      write(cfmt,920) nst                                                       
      print cfmt,(cst(k),k=1,nst)                                               
      write(cfmt,930) -ipwr,nst                                                 
      print cfmt,(clb(j),(znls(j,k,n),k=1,nst),j=1,nlb)                         
50    continue                                                                  
c                                                                               
      if(nznl.gt.0) then                                                        
        labz='znl 92/6'                                                         
        fstep=kdt                                                               
        do j=1,nlv                                                              
        do i=1,nlb                                                              
        wrkm(i,j)=0.                                                            
        enddo                                                                   
        enddo                                                                   
        do j=1,nst                                                              
        do i=1,nlb                                                              
        wrks(i,j)=0.                                                            
        enddo                                                                   
        enddo                                                                   
        write(nznl)labz                                                         
        write(nznl)rhour,fstep,idate,                                           
     &             (((znlm(j,k,n),j=1,nlb),k=1,nlv),n= 1,10),                   
     &             (((znlm(j,k,n),j=1,nlb),k=1,nlv),n=20,21),                   
     &             (  wrkm                         ,n=13,30),                   
     &             (((znls(j,k,n),j=1,nlb),k=1,nst),n= 9,15),                   
     &             (  wrks                         ,n= 8,11),                   
     &             (((znls(j,k,n),j=1,nlb),k=1,nst),n=16,22),                   
     &             (  wrks                         ,n=19,29),                   
     &             (((znls(j,k,n),j=1,nlb),k=1,nst),n=23,23),                   
     &             (((znlm(j,k,n),j=1,nlb),k=1,nlv),n=11,19),                   
     &             (  wrkm                         ,n=10,30),                   
     &             (((znls(j,k,n),j=1,nlb),k=1,nst),n= 1, 8),                   
     &             (((znls(j,k,n),j=1,nlb),k=1,nst),n=24,30),                   
     &             (  wrks                         ,n=16,30)                    
        close(nznl)                                                             
      endif                                                                     
      return                                                                    
900   format('0','zonally averaged diagnostics',2x,'kdt=',i4)                   
910   format(1x,a8,' (10**',i3,')',                                             
     &       ' fhour=',f6.1,' idate= (',4i4,')',' kdt=',i4)                     
920   format("(5x,'--lat-- ',",i2,"(5x,a4))")                                   
921   format("(5x,'--lat-- ',",i2,"(5x,'k=',i2))")                              
922   format("(5x,'--lat-- ',",i2,"(5x,'k=',i2),5x,' sum')")                    
923   format("(5x,'--lat-- ',                   5x,' sum')")                    
930   format("(5x,a8,",i3,"p",i2,"f9.2))")                                      
      end                                                                       
