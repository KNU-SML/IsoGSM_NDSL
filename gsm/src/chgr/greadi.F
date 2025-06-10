#define IJCAP ijcap
#define ILATG ilatg
#define ILEVS ilevs
#define ILONF ilonf
#define IGEN igen
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine greadi(n,fhour,idate,gzi,qi,tei,dii,zei,rqi
     &       ,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
     &       ,iret)
#include <paramodel.h>
      save                                                                      
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    greadi      read input sigma file for chgr.                    
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: reads a sigma file.                                                 
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell  docblock written (prehistorical program)            
c                                                                               
c usage:    call greadi(n,fhour,idate,gzi,qi,tei,dii,zei,rqi,iret)              
c   input argument list:                                                        
c     n        - logical unit number to read                                    
c                                                                               
c   output argument list:                                                       
c     fhour    - forecast hour                                                  
c     idate    - initial hour,month,day,year                                    
c     gzi      - orography                                                      
c     qi       - ln(psfc)                                                       
c     tei      - temperature                                                    
c     dii      - divergence                                                     
c     zei      - vorticity                                                      
c     rqi      - moisture                                                       
c     iret     - return code (0: normal, 1: end of file, 2: i/o error)          
c                                                                               
c   subprograms called:                                                         
c     setsig   - set sigma values                                               
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
c                                                                               
#include <parmchgr.h>
c                                                                               
      parameter(idimti=idimi*2,jdimhfi=jdimi/2,ijdimi=idimi*jdimi,              
     1          mdimi=(mwavei+1)*(mwavei+2),                                    
     3          mdimvi=(mwavei+1)*(mwavei+4),                                   
     4          mwvp2i=(mwavei+1)*2,                                            
     5          mdimhfi=mdimi/2,mdimvhi=mdimvi/2,                               
     6          kdimpi=kdimi+1,kdimmi=kdimi-1)                                  
c                                                                               
      dimension idate(4)                                                        
c hmhj 2                                                                        
      common /inpver/ siin(kdimpi),slin(kdimi),delin(kdimi),                    
     1                ciin(kdimpi),clin(kdimi),rpiin(kdimmi)                    
     2               ,ak(kdimpi),bk(kdimpi)

      common /label/ lab(8)                                                     
      character*4 lab                                                           
c                                                                               
      dimension gzi(mdimi),qi(mdimi),                                           
     1          tei(mdimi,kdimi),dii(mdimi,kdimi),zei(mdimi,kdimi),             
     2          rqi(mdimi,kdimqi)                                               
c
      parameter(kdum=201-ilevs_-1-ilevs_,kdum2=21,kens=2)
      dimension dummy(kdum),dummy2(kdum2),ensemble(kens)
c                                                                               
c     spectral input data file format                                           
c          lab                                                                  
c          hour,idate(4),siin(kdimpi),slin(kdimi)                               
c          zlni qi tei dii zei rqi                                              
c                                                                               
      rewind n                                                                  
      read(n,end=100,err=200) lab                                               
      print *,' end read lab.'                                                  
c     read(n) fhour,idate,siin,slin                                             
c
      read(n,err=201)fhour,idate,siin,slin
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
      if(tracers.eq.0.) then
        tracers=1.
        print *,'tracers reset to 1.'
      endif
      itrace=nint(tracers)
      icloud=nint(clouds)
      print *,'tread unit,fhour,idate=',n,fhour,idate
      print *,' number of tracers input = ',itrace
      print *,' number of cloud input   = ',icloud
      goto 202
201   continue
      rewind n
      read(n) lab
      read(n,err=201)fhour,idate,siin,slin
c
      do i=1,kdum
        dummy(i)=0.
      enddo
      waves=ijcap_
      xlayers=ilevs_
      trun=1.
      order=2.
      realform=1.
      gencode=igen_
      rlond=ilonf_
      rlatd=ilatg_
      rlonp=ilonf_
      rlatp=ilatg_
      rlonr=ilonf_
      rlatr=ilatg_
      tracers=1.
      clouds=0.
      pdryini=0.
      subcen=0.
      do i=1,kens
        ensemble(i)=0.
      enddo
      ppid=0.
      slid=0.
      vcid=0.
      vmid=0.
      vtid=0.
      do k=1,kdum2
        dummy2(k)=0.
      enddo
c
      itrace=nint(tracers)
      icloud=nint(clouds)
c
      print *,'tread old format unit,fhour,idate=',n,fhour,idate
      print *,' number of tracers input = ',itrace
      print *,' number of cloud input   = ',icloud
202   continue
c
      rmaxsiin=siin(1)
      do k=1,kdimpi
       rmaxsiin=max(rmaxsiin,siin(k))
      enddo
c
      if(rmaxsiin.gt.1. .or. rmaxsiin.eq.0.) then
       call setakbk(kdimi,siin,slin,ak,bk,delin)
      else
       call setsig(ciin,siin,delin,slin,clin,rpiin)
      endif

      read(n) gzi                                                               
      read(n) qi                                                                
c                                                                               
      do 10 k=1,kdimi                                                           
      read(n)(tei(i,k),i=1,mdimi)                                               
   10 continue                                                                  
      do 20 k=1,kdimi                                                           
      read(n)(dii(i,k),i=1,mdimi)                                               
      read(n)(zei(i,k),i=1,mdimi)                                               
   20 continue                                                                  
c                                                                               
      do 30 k=1,kdimi*min(itrace,ntotal_)                                    
      read(n)(rqi(i,k),i=1,mdimi)                                               
   30 continue                                                                  
*** set zero for tracers not included in input ***
#ifdef ISOTOPE
#ifndef NOFRAC
      do k=kdimi+1,kdimi*2
        do i=1,mdimi
          rqi(i,k)=0.
        enddo
      enddo
      do k=kdimi*2+1,kdimqi
        do i=1,mdimi
          rqi(i,k)=rqi(i,mod(k-1,kdimi)+1)
        enddo
      enddo
#else
      do k=kdimi*itrace+1,kdimqi
        do i=1,mdimi
          rqi(i,k)=0.
        enddo
      enddo
#endif
c
#else
      if (itrace.lt.ntrac_) then
        do k=kdimi*itrace+1,kdimqi
          do i=1,mdimi
            rqi(i,k)=0.
          enddo
        enddo
      endif
#endif
c                                                                               
      iret=0                                                                    
      return                                                                    
c                                                                               
100   continue                                                                  
      iret=1                                                                    
      return                                                                    
c                                                                               
200   continue                                                                  
      iret=2                                                                    
      return                                                                    
c                                                                               
      end                                                                       
