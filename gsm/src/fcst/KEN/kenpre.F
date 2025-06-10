#define LATG2 latg2
#define LATG latg
#define LEVS levs
#define LONF lonf
#define LPNT lpnt
#define LTSTP ltstp
#define MLVARK mlvark
#define PI pi
#define SLVARK slvark
#include <define.h>
      subroutine kenpre(con,colrad,idim,jdim2,n1)                             
c....                                                                           
c  this routine computes the kendata grid point indices igrd,jgrd for           
c   fcst grid, for the                       
c   npoint points,given the longitude and latitude of same (con).               
c   colrad is the colatitude of the fcst grid (dimensioned jdim2,which          
c   is half of the total latitudinal grid points),                              
c....                                                                           
#include <paramodel.h>
#include <comgpd.h>
#ifdef KEN
      dimension colrad(jdim2),con(1700)                                         
      dimension alat(nptken),alon(nptken),kpoi(nptken)                          
      dimension slimsk(lonf_,latg_),blatf(latg_)
      dimension colrar(latg2_),wgr(latg2_),wgrcs(latg2_),rrs2(latg2_)           
      data pi/pi_/                                                              
c....   begin here..                                                            
      do 2 k=1,nstken                                                           
       do 2 j=1,nptken                                                          
        do 2 i=1,nvrken                                                         
         svdata(i,j,k) = 0.                                                     
    2 continue                                                                  
c--- get slmsk because we wish to get nearest point                             
c       of same sfc type...this array not available til after step1             
c       if we were not trying to cover ourselves for out-board radi8            
c       ,then this code could be called from step1...                           
      open (unit=n1,file='sfci',form='unformatted',err=999)
      go to 998
  999 continue
      print *,'error opening sfci in kenpre'
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
      read(n1) ghour,id1,id2,id3,id4                                          
#ifndef NOPRINT
   99 format(1h ,'fhour, idate=',f6.2,2x,4(1x,i4))                              
      print *,'in kenpre read slmsk from unit=',n1                            
      print 99,ghour, id1,id2,id3,id4                                           
#endif
      read(n1)                                                                
      read(n1)                                                                
      read(n1)                                                                
      read(n1)                                                                
      read(n1)                                                                
      read(n1)                                                                
c.....  skip cv, cvb, cvt, albedo                                               
      read(n1)                                                                
      read(n1)                                                                
      read(n1)                                                                
      read(n1)                                                                
      read(n1) slimsk                                                         
      rewind n1                                                               
cccc           call row1ns(slimsk)                                              
      call glats (latg2_, colrar, wgr, wgrcs, rrs2)                             
      dxf = 360. / lonf_                                                        
      dxr = 360. / lonf_                                                        
      ilonf = lonf_                                                             
      jlatg2 = latg2_                                                           
      jlatg = latg_                                                             
      jfp1 = jlatg2 + 1                                                         
c----    get latitude of gaussian grids                                         
      do 3 j=1,jlatg2                                                           
        blatf(j) =(pi /2. - colrad(j)) * 180. / pi                              
    3 continue                                                                  
      blatf(jfp1) = -blatf(jlatg2)                                              
c...    put lat/lon into useable arrays (max=200),where                         
c         npoint gt 0 implies npoint lat/lon s in con and                       
c           if abs(lat) between   0, 90 look for nearest point                  
c                       between 100,190 look for nearest land point             
c                       between 200,290 look for nearest sea point              
c         npoint lt 0 implies lat/lon of center of region ,only..               
c                    lat,lon=con(1301),con(1501)                                
c           let xy=abs(npoint) and always be 2 digits                           
c                              and do not differentiate land/sea,               
c             then x between 1,9 means create array of every x points           
c                   (i.e. x=1 means every point,x=3 means every 3rd,..          
c              and y between 0,9 means create (y+1,y+1) array..                 
c           thus xy can have values 10-99                                       
c....                                                                           
      npute = -1                                                                
      if (npoint.lt.0) then                                                     
       xy = abs(npoint)                                                         
       if (xy.lt.10..or.xy.gt.99.) then                                         
        npute = 0                                                               
#ifndef NOPRINT
        print 98,npoint                                                         
   98   format(1h ,' num(1300)=',i6,'out of -range, so set=1')                  
#endif
        npoint = 1                                                              
       else                                                                     
        npoint = 1                                                              
        iskp = xy/10                                                            
        iy   = xy - iskp*10 + 1                                                 
#ifndef NOPRINT
        print 97,iy,iy,iskp                                                     
   97   format(1h ,' prepare regional (',i2,',',i2,') array - every',           
     1             i2,' points')                                                
#endif
        npute = iy * iy                                                         
       end if                                                                   
      end if                                                                    
      do 5 k = 1, npoint                                                        
       ils = -1                                                                 
       ylat = abs(con(k+1300))                                                  
       if (ylat.ge.100.and.ylat.le.190.) then                                   
c...     land point is desired...                                               
        ils = 1                                                                 
        sgn = con(k+1300) / ylat                                                
        con(k+1300) = ylat-100.                                                 
        if (sgn.lt.0.) con(k+1300) = - (ylat-100.)                              
       end if                                                                   
       if (ylat.ge.200.and.ylat.le.290.) then                                   
        ils = 0                                                                 
        sgn = con(k+1300) / ylat                                                
        con(k+1300) = ylat-200.                                                 
        if (sgn.lt.0.) con(k+1300) = - (ylat-200.)                              
       end if                                                                   
       xlat = con(k+1300)                                                       
       xlon = con(k+1500)                                                       
#ifndef NOPRINT
       if (npute.lt.0.and.ils.eq.-1) print 197,k,xlat,xlon                      
       if (npute.lt.0.and.ils.eq.0) print 198,k,xlat,xlon                       
       if (npute.lt.0.and.ils.eq.1) print 199,k,xlat,xlon                       
  197  format(1h ,' ==== station ',i4,' at latlon=',2f8.2,                      
     1            ' desired as nearest point')                                  
  198  format(1h ,' ==== station ',i4,' at latlon=',2f8.2,                      
     1            ' desired as ocean pt')                                       
  199  format(1h ,' ==== station ',i4,' at latlon=',2f8.2,                      
     1            ' desired as land pt')                                        
#endif
       if (npute.gt.0.and.k.gt.1) go to 195                                     
        alat(k) = con(k+1300)                                                   
        alon(k) = con(k+1500)                                                   
  195  continue                                                                 
       if (xlon.lt.0) xlon = 360. + con(k+1500)                                 
       if (npute.lt.0) then                                                     
         ils = -1                                                               
         call getij (xlat,xlon,slimsk,blatf,dxf,                                
     1               ils,ilonf,jlatg,ki,kj)                                     
       else                                                                     
         ils = -1                                                               
         call getij (xlat,xlon,slimsk,blatf,dxf,                                
     1               ils,ilonf,jlatg,ki,kj)                                     
       end if                                                                   
       igrd(k) = ki                                                             
       jgrd(k) = kj                                                             
       if(npute.gt.0) go to 5                                                   
        if(xlat.lt.0.) then                                                     
         igrd(k) = ki + ilonf                                                   
         jgrd(k) = jlatg + 1 - kj                                               
        endif                                                                   
    5 continue                                                                  
c....    regional block , i,j still in single latitude structure..              
      if (npute.gt.0) then                                                      
       iback = iy/2                                                             
c....    if iy = 1 the all we want is 1 point                                   
       if (iback.le.0) then                                                     
         npoint = 1                                                             
         go to 59                                                               
       end if                                                                   
       istarf = igrd(1) - iback*iskp                                            
       jstarf = jgrd(1) - iback*iskp                                            
       npoint = 0                                                               
       do 30 kyj=1,iy                                                           
        do 30 kxi=1,iy                                                          
         npoint = npoint + 1                                                    
         igrd(npoint) = istarf + (kxi-1)*iskp                                   
         jgrd(npoint) = jstarf + (kyj-1)*iskp                                   
   30  continue                                                                 
       do 32 n=1,npoint                                                         
         kpoi(n) = 0                                                            
         if (jgrd(n).gt.jlatg.or.jgrd(n).lt.1) go to 32                         
         if (igrd(n).gt.ilonf) igrd(n) = igrd(n) - ilonf                        
         if (igrd(n).lt.1) igrd(n) = igrd(n) + ilonf                            
         kpoi(n) = n                                                            
   32  continue                                                                 
c...    squeeze out the out of bounds points(kpoi=0)                            
       npp = 0                                                                  
       do 33 n=1,npoint                                                         
        if (kpoi(n).le.0) go to 33                                              
         npp = npp + 1                                                          
         igrd(npp) = igrd(kpoi(n))                                              
         jgrd(npp) = jgrd(kpoi(n))                                              
         if (jgrd(npp).gt.jlatg2) then                                          
          igrd(npp) = igrd(npp) + ilonf                                         
          jgrd(npp) = jlatg+1-jgrd(npp)                                         
         end if                                                                 
   33  continue                                                                 
       npoint = npp                                                             
      end if                                                                    
c...................  debug print                                               
   59 continue                                                                  
      do 60 k=1,npoint                                                          
       ig=igrd(k)                                                               
       jg=jgrd(k)                                                               
       iclnd=ig                                                                 
       jclnd=jg                                                                 
       if(igrd(k).le.ilonf) then                                                
        blat=90.-colrad(jgrd(k))*180./pi                                        
        blon=(igrd(k)-1)*360./ilonf                                             
        if(blon.gt.180.) blon=blon-360.                                         
       else                                                                     
        blat=colrad(jgrd(k))*180./pi-90.                                        
        blon=(igrd(k)-1-ilonf)*360./ilonf                                       
        if(blon.gt.180.) blon=blon-360.                                         
        iclnd=ig-ilonf                                                          
        jclnd=jlatg+1-jg                                                        
       endif                                                                    
#ifndef NOPRINT
       write(6,61) k,alat(k),alon(k),blat,blon                                  
       write(6,62) jgrd(k),igrd(k),slimsk(iclnd,jclnd)                          
60    continue                                                                  
61    format(' kenpre: k,orig lat-lon,compt lat-lon=',i4,4f8.2)                 
62    format('          ....jgrd,igrd,slmsk=',2i6,f6.1)                         
#endif
#endif
      return                                                                    
      end                                                                       
