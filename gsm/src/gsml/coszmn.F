#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine coszmn(imx2,jmx2,
     1                  dtswav,solhr,sinlat,coslat,sdec,cdec,slag,
     2                  xlon,coszen,ldg,coszdg)                       
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c                                                                               
c===>  compute mean cos solar zen angl over dtswav hrs                          
c....   cosine of solar zen angl for both n. and s. hemispheres.                
c        solhr=time(hrs) after 00z (greenwich time)..                           
c        xlon is east long(radians)..                                           
c        sinlat, coslat are sin and cos of latitude (n. hemisphere)             
c        sdec, cdec = the sine and cosine of the solar declination.             
c        slag = equation of time                                                
#include <paramodel.h>
      dimension xlon(imx2,jmx2),coszen(imx2,jmx2)
      logical ldg                                           
      dimension coszdg(imx2,jmx2)                    
      dimension sinlat(imx2,jmx2),coslat(imx2,jmx2)
c local array   
c
c yifeng cui added rsm mp    org:sdsc     date:02-09-05
c                                                                
#ifdef MP
#ifdef RSM
#define ILOTS igrd12p_
#else
#define ILOTS lonf2p_
#endif
#else
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      dimension coszn(ILOTS)                                                    
      dimension istsun(ILOTS)                                                   
      nstp = 6                                                                  
      istp = nstp*dtswav                                                        
      pid12 = (2.e0 * asin(1.e0)) / 12.e0                                    
c----                                                                           
      do 2000 j=1,jmx2                                                          
        do 20 i=1,imx2                                                          
          coszen(i,j) = 0.e0                                                   
          istsun(i) = 0                                                         
   20   continue                                                                
      do 1000 it=1,istp                                                         
       cns = pid12 * (solhr-12.e0+(it-1)*1.e0/nstp) +slag                     
       do 40 i=1,imx2                                                           
         ss=sinlat(i,j)*sdec                                                    
         cc=coslat(i,j)*cdec                                                    
         coszn(i) = ss + cc * cos(cns + xlon(i,j))                              
         coszen(i,j) = coszen(i,j) + max (0.e0, coszn(i))                      
         if(coszn(i).gt.0.e0) istsun(i) = istsun(i) + 1                        
   40  continue                                                                 
 1000 continue                                                                  
       do 1500 i=1,imx2                                                         
        if(ldg) coszdg(i,j) = coszen(i,j) / istp                                
        if(istsun(i).gt.0) coszen(i,j) = coszen(i,j) / istsun(i)                
 1500  continue                                                                 
 2000 continue                                                                  
      return                                                                    
      end                                                                       
