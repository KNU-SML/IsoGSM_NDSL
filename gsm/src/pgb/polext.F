      subroutine polext(mp,im,fnx,fsx,fn,fs)                                    
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    polext      extrapolate a field to the poles.                  
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31                     
c                                                                               
c abstract: a global horizontal field is exterpolated to the poles.             
c   polar scalars are the average of the closest latitude circle values.        
c   polar vector components are taken from the wavenumber 1 component           
c   extracted from the values on the closest latitude circle.                   
c   polar flags are copied from the closest prime meridian value.               
c                                                                               
c program history log:                                                          
c   93-04-28  iredell                                                           
c                                                                               
c usage:    call polext(mp,im,fnx,fsx,fn,fs)                                    
c   input argument list:                                                        
c     mp       - integer field parameter identifier                             
c                (0 for scalar, 1 for vector, 2 for flag)                       
c     im       - integer number of longitudes                                   
c     fnx      - real (im) field values on the closest latitude circle          
c                to the north pole                                              
c     fsx      - real (im) field values on the closest latitude circle          
c                to the south pole                                              
c                                                                               
c   output argument list:                                                       
c     fn       - real (im) field values extrapolated to the north pole          
c     fs       - real (im) field values extrapolated to the south pole          
c                                                                               
c attributes:                                                                   
c   language: ansi fortran 77                                                   
c                                                                               
c$$$                                                                            
      real fnx(im),fsx(im),fn(im),fs(im)                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
c  get polar values for scalars or vectors                                      
      pi=acos(-1.)                                                              
      if(mp.eq.0) then                                                          
c  full scalar                                                                  
        fnp=0.                                                                  
        fsp=0.                                                                  
        do 1010 i=1,im                                                          
          fnp=fnp+fnx(i)                                                        
          fsp=fsp+fsx(i)                                                        
1010      continue                                                              
        fnp=fnp/im                                                              
        fsp=fsp/im                                                              
        do 1020 i=1,im                                                          
          fn(i)=fnp                                                             
          fs(i)=fsp                                                             
1020    continue                                                                
      elseif(mp.eq.1) then                                                      
c  full vector                                                                  
        fnpc=0.                                                                 
        fnps=0.                                                                 
        fspc=0.                                                                 
        fsps=0.                                                                 
        do 1030 i=1,im                                                          
          ci=cos(2*pi*(i-1)/im)                                                 
          si=sin(2*pi*(i-1)/im)                                                 
          fnpc=fnpc+ci*fnx(i)                                                   
          fnps=fnps+si*fnx(i)                                                   
          fspc=fspc+ci*fsx(i)                                                   
          fsps=fsps+si*fsx(i)                                                   
1030    continue                                                                
        fnpc=2*fnpc/im                                                          
        fnps=2*fnps/im                                                          
        fspc=2*fspc/im                                                          
        fsps=2*fsps/im                                                          
        do 1040 i=1,im                                                          
          ci=cos(2*pi*(i-1)/im)                                                 
          si=sin(2*pi*(i-1)/im)                                                 
          fn(i)=fnpc*ci+fnps*si                                                 
          fs(i)=fspc*ci+fsps*si                                                 
1040    continue                                                                
      elseif(mp.eq.2) then                                                      
c  full flag                                                                    
        do 1050 i=1,im                                                          
          fn(i)=fnx(1)                                                          
          fs(i)=fsx(1)                                                          
1050    continue                                                                
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
