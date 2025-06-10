#define IGRD igrd
#define IGRD1 igrd1
#define IO2 io2
#define JGRD jgrd
#define JGRD1 jgrd1
#define JO2 jo2
#define LATG latg
#define LONF lonf
      subroutine gribit2(f,lbm,idrt,im,jm,mxbit,colat1,                          
     &                  ilpds,iptv,icen,igen,ibms,ipu,itl,il1,il2,              
     &                  iyr,imo,idy,ihr,iftu,ip1,ip2,itr,                       
     &                  ina,inm,icen2,ids,iens,                                 
     &                  xlat1,xlon1,xlat2,xlon2,delx,dely,oritru,proj,          
     &                  grib,lgrib,ierr)                                        
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    gribit      create grib message                                
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: create a grib message from a full field.                            
c   at present, only global latlon grids and gaussian grids                     
c   and regional polar projections are allowed.                                 
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c   94-05-04  juang (for gsm and rsm use)                                       
c   98-01-28  juang chnage y2k and add lambert
c                                                                               
c usage:    call gribit(f,lbm,idrt,im,jm,mxbit,colat1,                          
c    &                  ilpds,iptv,icen,igen,ibms,ipu,itl,il1,il2,              
c    &                  iyr,imo,idy,ihr,iftu,ip1,ip2,itr,                       
c    &                  ina,inm,icen2,ids,iens,                                 
c    &                  xlat1,xlon1,delx,dely,oritru,proj,                      
c    &                  grib,lgrib,ierr)                                        
c   input argument list:                                                        
c     f        - real (im*jm) field data to pack into grib message              
c     lbm      - logical (im*jm) bitmap to use if ibms=1                        
c     idrt     - integer data representation type                               
c                (0 for latlon or 4 for gaussian or 5 for polar)                
c     im       - integer longitudinal dimension                                 
c     jm       - integer latitudinal dimension                                  
c     mxbit    - integer maximum number of bits to use (0 for no limit)         
c     colat1   - real first colatitude of grid if idrt=4 (radians)              
c     ilpds    - integer length of the pds (usually 28)                         
c     iptv     - integer parameter table version (usually 1)                    
c     icen     - integer forecast center (usually 7)                            
c     igen     - integer model generating code                                  
c     ibms     - integer bitmap flag (0 for no bitmap)                          
c     ipu      - integer parameter and unit indicator                           
c     itl      - integer type of level indicator                                
c     il1      - integer first level value (0 for single level)                 
c     il2      - integer second level value                                     
c     iyr      - integer year                                                   
c     imo      - integer month                                                  
c     idy      - integer day                                                    
c     ihr      - integer hour                                                   
c     iftu     - integer forecast time unit (1 for hour)                        
c     ip1      - integer first time period                                      
c     ip2      - integer second time period (0 for single period)               
c     itr      - integer time range indicator (10 for single period)            
c     ina      - integer number included in average                             
c     inm      - integer number missing from average                            
c     icen2    - integer forecast subcenter                                     
c                (usually 0 but 1 for reanal or 2 for ensemble)                 
c     ids      - integer decimal scaling                                        
c     iens     - integer (5) ensemble extended pds values                       
c                (application,type,identification,product,smoothing)            
c                (used only if icen2=2 and ilpds>=45)                           
c     xlat1    - real first point of regional latitude (radians)                
c     xlon1    - real first point of regional longitude (radians)               
c     xlat2    - real last  point of regional latitude (radians)                
c     xlon2    - real last  point of regional longitude (radians)               
c     delx     - real dx on 60n for regional (m)                                
c     dely     - real dy on 60n for regional (m)                                
c     proj     - real polar projection flag 1 for north -1 for south            
c                     mercater projection 0                                     
c     oritru   - real orientation of regional polar projection or               
c                     truth for regional mercater projection                    
c                                                                               
c   output argument list:                                                       
c     grib     - character (lgrib) grib message                                 
c     lgrib    - integer length of grib message                                 
c                (no more than 100+ilpds+im*jm*(mxbit+1)/8)                     
c     ierr     - integer error code (0 for success)                             
c                                                                               
c subprograms called:                                                           
c   gtbits     - compute number of bits and round data appropriately            
c   w3fi72     - engrib data into a grib1 message                               
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      real f(im*jm)                                                             
      logical lbm(im*jm)                                                        
      character grib(*)                                                         
#ifdef DYNAMIC_ALLOC
      integer ibm(im*jm),ipds(100),igds(100),ibds(100)              
      real fr(im*jm)                                                            
      character pds(ilpds)                                                      
#else
#ifdef RSM
#if ( ( _igrd_ + 1 )*( _jgrd_ + 1 ) >= ( _io2_ * _jo2_ ) )
#define NIBM igrd1_*jgrd1_
#else
#define NIBM io2_*jo2_
#endif
#else
#if ( ( _lonf_ * _latg_ ) >= ( _io2_ * _jo2_ ) )
#define NIBM lonf_*latg_
#else
#define NIBM io2_*jo2_
#endif
#endif
      integer ibm(NIBM),ipds(100),igds(100),ibds(100)                           
      real fr(NIBM)                                                             
      character pds(100)                                                        
#endif
c                                                                               
      integer iens(5),kprob(2),kclust(16),kmembr(80)                            
      dimension xprob(2)                                                        
c                                                                               
#ifdef REAL4_W3LIB
      integer*4 iens4(5),kprob4(2),kclust4(16),kmembr4(80)                      
#ifdef DYNAMIC_ALLOC
      integer*4 ibm4(im*jm)
#else
      integer*4 ibm4(NIBM)
#endif
      integer*4 ipds4(100),igds4(100),ibds4(100)                                
      integer*4 nbit4,nf4,nfo4,lgrib4,ierr4                                     
#endif
c
      kprob(1)=0
      kprob(2)=0
c                                                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  determine grid parameters                                                    
      pi=acos(-1.)                                                              
      nf=im*jm                                                                  
      if(idrt.eq.0) then                                                        
        if(im.eq.144.and.jm.eq.73) then                                         
          igrid=2                                                               
        elseif(im.eq.360.and.jm.eq.181) then                                    
          igrid=3                                                               
        else                                                                    
          igrid=255                                                             
        endif                                                                   
        iresfl=128                                                              
        iscan=0                                                                 
        lat1=nint(90.e3)                                                        
        lon1=0                                                                  
        lati=nint(180.e3/(jm-1))                                                
        loni=nint(360.e3/im)                                                    
        igds09=-lat1                                                            
        igds10=-loni                                                            
        igds11=lati                                                             
        igds12=loni                                                             
        igds13=0                                                             
        igds14=0                                                             
        igds15=0                                                             
        igds16=0                                                             
        igds17=0                                                             
        igds18=0                                                             
      elseif(idrt.eq.4) then                                                    
        if(im.eq.192.and.jm.eq.94) then                                         
          igrid=98                                                              
        elseif(im.eq.384.and.jm.eq.190) then                                    
          igrid=126                                                             
        else                                                                    
          igrid=255                                                             
        endif                                                                   
        iresfl=128                                                              
        iscan=0                                                                 
        lat1=nint(90.e3-180.e3/pi*colat1)                                       
        lon1=0                                                                  
        lati=jm/2                                                               
        loni=nint(360.e3/im)                                                    
        igds09=-lat1                                                            
        igds10=-loni                                                            
        igds11=lati                                                             
        igds12=loni                                                             
        igds13=iscan                                                            
        igds14=0                                                             
        igds15=0                                                             
        igds16=0                                                             
        igds17=0                                                             
        igds18=0                                                             
      elseif(idrt.eq.5) then    ! polar projection                              
        igrid=255                                                               
        lat1=nint(180.e3/acos(-1.) * xlat1)                                     
        lon1=nint(180.e3/acos(-1.) * xlon1)                                     
        iresfl=0                                                                
        igds09=nint(oritru*1.e3)                                                
        igds10=delx                                                             
        igds11=dely                                                             
        if( nint(proj).eq.1  ) igds12=0         ! north polar proj              
        if( nint(proj).eq.-1 ) igds12=128       ! south polat proj              
        iscan=64                                                                
        igds13=iscan                                                            
        igds14=0                                                             
        igds15=0                                                             
        igds16=0                                                             
        igds17=0                                                             
        igds18=0                                                             
#ifdef RSM
      elseif(idrt.eq.3) then    ! lambert projection
        igrid=255
        lat1=nint(180.e3/acos(-1.) * xlat1)
        lon1=nint(180.e3/acos(-1.) * xlon1)
        iresfl=8
        igds09=nint(oritru*1.e3)
        igds10=delx
        igds11=dely
        if( nint(proj).eq.2  ) igds12=0         ! north proj
        if( nint(proj).eq.-2 ) igds12=128       ! south proj
        iscan=64
        igds13=iscan
        igds14=0
        igds15=nint(oritru*1.e3)
        igds16=igds15
        igds17=-90000
        igds18=0
#endif
      elseif(idrt.eq.1) then    ! mercater projection                           
        igrid=255                                                               
        lat1=nint(180.e3/acos(-1.) * xlat1)                                     
        lon1=nint(180.e3/acos(-1.) * xlon1)                                     
        iresfl=0                                                                
        igds09=nint(180.e3/acos(-1.) * xlat2)                                   
        igds10=nint(180.e3/acos(-1.) * xlon2)                                   
        igds11=delx                                                             
        igds12=dely                                                             
        igds13=nint(oritru*1.e3)                                                
        iscan=64                                                                
        igds14=iscan                                                            
        igds15=0                                                             
        igds16=0                                                             
        igds17=0                                                             
        igds18=0                                                             
      else                                                                      
        ierr=40                                                                 
        return                                                                  
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  reset time range parameter in case of overflow                               
      if(itr.ge.2.and.itr.le.5.and.ip2.ge.256) then                             
        jp1=ip2                                                                 
        jp2=0                                                                   
        jtr=10                                                                  
      else                                                                      
        jp1=ip1                                                                 
        jp2=ip2                                                                 
        jtr=itr                                                                 
      endif                                                                     
c for y2k
c     iyr4=iyr
c     if(iyr.le.100) iyr4=2050-mod(2050-iyr,100)
      iy =mod(iyr-1,100)+1
      ic =(iyr-1)/100+1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  fill pds parameters                                                          
      ipds(01)=ilpds    ! length of pds                                         
      ipds(02)=iptv     ! parameter table version id                            
      ipds(03)=icen     ! center id                                             
      ipds(04)=igen     ! generating model id                                   
      ipds(05)=igrid    ! grid id                                               
      ipds(06)=1        ! gds flag                                              
      ipds(07)=ibms     ! bms flag                                              
      ipds(08)=ipu      ! parameter unit id                                     
      ipds(09)=itl      ! type of level id                                      
      ipds(10)=il1      ! level 1 or 0                                          
      ipds(11)=il2      ! level 2                                               
      ipds(12)=iy       ! year                                                  
      ipds(13)=imo      ! month                                                 
      ipds(14)=idy      ! day                                                   
      ipds(15)=ihr      ! hour                                                  
      ipds(16)=0        ! minute                                                
      ipds(17)=iftu     ! forecast time unit id                                 
      ipds(18)=jp1      ! time period 1                                         
      ipds(19)=jp2      ! time period 2 or 0                                    
      ipds(20)=jtr      ! time range indicator                                  
      ipds(21)=ina      ! number in average                                     
      ipds(22)=inm      ! number missing                                        
      ipds(23)=ic       ! century                                               
      ipds(24)=icen2    ! forecast subcenter                                    
      ipds(25)=ids      ! decimal scaling                                       
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  fill gds and bds parameters                                                  
      igds(01)=0        ! number of vertical coords                             
      igds(02)=255      ! vertical coord flag                                   
      igds(03)=idrt     ! data representation type                              
      igds(04)=im       ! east-west points                                      
      igds(05)=jm       ! north-south points                                    
      igds(06)=lat1     ! latitude of origin                                    
      igds(07)=lon1     ! longitude of origin                                   
      igds(08)=iresfl   ! resolution flag                                       
      igds(09)=igds09   ! latitude of end or orientation                        
      igds(10)=igds10   ! longitude of end or dx in meter on 60n                
      igds(11)=igds11   ! lat increment or gaussian lats or dy in meter         
      igds(12)=igds12   ! longitude increment or projection                     
      igds(13)=igds13   ! scanning mode or lat of intercut on earth for mercater
      igds(14)=igds14   ! not used or scanning mode for mercater                
      igds(15)=igds15   ! not used or cut latitude near pole for lambert
      igds(16)=igds16   ! not used or cut latitude near equator for lambert
      igds(17)=igds17   ! not used or lat of south pole for lambert
      igds(18)=igds18   ! not used or lon of south pole for lambert
      ibds(1)=0       ! bds flags                                               
      ibds(2)=0       ! bds flags                                               
      ibds(3)=0       ! bds flags                                               
      ibds(4)=0       ! bds flags                                               
      ibds(5)=0       ! bds flags                                               
      ibds(6)=0       ! bds flags                                               
      ibds(7)=0       ! bds flags                                               
      ibds(8)=0       ! bds flags                                               
      ibds(9)=0       ! bds flags                                               
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  fill bitmap and count valid data.  reset bitmap flag if all valid.           
      nbm=nf                                                                    
      if(ibms.ne.0) then                                                        
        nbm=0                                                                   
        do i=1,nf                                                               
          if(lbm(i)) then                                                       
            ibm(i)=1                                                            
            nbm=nbm+1                                                           
          else                                                                  
            ibm(i)=0                                                            
          endif                                                                 
        enddo                                                                   
        if(nbm.eq.nf) ipds(7)=0                                                 
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  round data and determine number of bits                                      
      if(nbm.eq.0) then                                                         
        do i=1,nf                                                               
          fr(i)=0.                                                              
        enddo                                                                   
        nbit=0                                                                  
      else                                                                      
        call gtbits(ipds(7),ids,nf,ibm,f,fr,fmin,fmax,nbit)                     
c       write(0,'("gtbits:",4i4,4x,2i4,4x,2g16.6)')                             
c    &   ipu,itl,il1,il2,ids,nbit,fmin,fmax                                     
        if(mxbit.gt.0) nbit=min(nbit,mxbit)                                     
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  create product definition section                                            
#ifdef REAL4_W3LIB
      do i=1,100                                                                
        ipds4(i)=ipds(i)                                                        
      enddo                                                                     
      call w3fi68(ipds4,pds)                                                    
#else
      call w3fi68(ipds,pds)                                                     
#endif
      if(icen2.eq.2.and.ilpds.ge.45) then                                       
        ilast=45                                                                
#ifdef REAL4_W3LIB
        do i=1,5                                                                
          iens4(i)=iens(i)                                                      
        enddo                                                                   
        do i=1,2                                                                
          kprob4(i)=kprob(i)                                                    
        enddo                                                                   
        do i=1,16                                                               
          kclust4(i)=kclust(i)                                                  
        enddo                                                                   
        do i=1,80                                                               
          kmembr4(i)=kmembr(i)                                                  
        enddo                                                                   
        ilast4=ilast                                                            
        call pdsens(iens4,kprob4,xprob,kclust4,kmembr4,ilast4,pds)              
#else
        call pdsens(iens,kprob,xprob,kclust,kmembr,ilast,pds)                   
#endif
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  create grib message                                                          
#ifdef REAL4_W3LIB
      nbit4=nbit                                                                
      do i=1,100                                                                
        ipds4(i)=ipds(i)                                                        
      enddo                                                                     
      do i=1,100                                                                
        igds4(i)=igds(i)                                                        
      enddo                                                                     
      do i=1,nf                                                                 
        ibm4(i)=ibm(i)                                                          
      enddo                                                                     
      nf4=nf                                                                    
      do i=1,100                                                                
        ibds4(i)=ibds(i)                                                        
      enddo                                                                     
      call w3fi72(0,fr,0,nbit4,1,ipds4,pds,                                     
     &            1,255,igds4,0,0,ibm4,nf4,ibds4,                               
     &            nfo4,grib,lgrib4,ierr4)                                       
      nfo=nfo4                                                                  
      lgrib=lgrib4                                                              
      ierr=ierr4                                                                
#else
      call w3fi72(0,fr,0,nbit,1,ipds,pds,                                       
     &            1,255,igds,0,0,ibm,nf,ibds,                                   
     &            nfo,grib,lgrib,ierr)                                          
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
