#define LEVP1 levp1
#define LEVS levs
      subroutine rdsgh(nsig,fh,fhour,idate,si,sl,iret,ntrace,ncloud)
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    rdsgh       read sigma file header record                      
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: reads the header record from the sigma file.                        
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:      call rdsgh(nsig,fhour,idate,si,sl,iret)                           
c                                                                               
c   input argument list:                                                        
c     nsig     - integer unit from which to read header                         
c                                                                               
c   output argument list:                                                       
c     fhour    - real forecast hour                                             
c     idate    - integer (4) date                                               
c     si       - real (levs+1) sigma interfaces                                 
c     sl       - real (levs) sigma levels                                       
c                                                                               
c   input files:                                                                
c     nsig     - sigma file                                                     
c                                                                               
c subprograms called:                                                           
c   maxfac       return maximum prime factor                                    
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      character*32 clabe                                                        
      dimension idate(4)                                                        
      dimension si(levs_+1),sl(levs_)                                           
      dimension dummy(201-levp1_-levs_)
      dimension ensemble(2),dummy2(21)
      character*3 fni
      parameter(fni='sig',nchi=3)
c
      character*80 fno
#ifdef ASSIGN
      character*80 asgnstr
#endif
c
      call fnam(fni,nchi,fh,fno,ncho)
#ifdef ASSIGN
      write(asgnstr,'(9hassign u:,I2,)') nsig
      call assign('assign -R')
      call assign(asgnstr)
#endif
      open(unit=nsig,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      write(6,*) ' error in opening file ',fno(1:ncho)
      call abort
  901 continue
      write(6,*) ' file ',fno(1:ncho),' opened. unit=',nsig

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  read and extract header record                                               
c  read sigma spectral file header and determine gaussian grid                  
c     print *,'reading lab'                                                     
c     call flush(6)                                                             
      read(nsig,end=91,err=92) clabe                                            
c     print *,'reading fhour,....'                                              
c     call flush(6)                                                             
      read(nsig,err=201)fhour,idate,si,sl
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
      ntrace=nint(tracers)
      ncloud=nint(clouds)
      if(fh.ne.fhour) then
        print *,'fh and fhour does not match'
        call abort
      endif
      goto 202
201   continue
      rewind nsig
      read(nsig,end=91,err=92) clabe                                            
      read(nsig,err=201)fhour,idate,si,sl
      ntrace=1
      ncloud=0
202   continue
#ifndef NOPRINT
      print *,'tread unit,fhour,idate=',n,fhour,idate
      print *,' number of tracers input = ',ntrace
      print *,' number of cloud input   = ',ncloud
#endif
      iret=0
      return                                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  end of file encountered                                                      
91    iret=1                                                                    
      return                                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  i/o error encountered                                                        
92    iret=2                                                                    
      return                                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      end                                                                       
