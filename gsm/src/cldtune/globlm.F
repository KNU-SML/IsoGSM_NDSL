      subroutine globlm(y,globl,rlat,wgt,wgta,                                  
     1                               kd,ld,colrad,ltitl,fhour,ipts)             
      save                                                                      
c...   compute area weighted means - global,hemispheric,etcccc                  
c..   input:                                                                    
c          y      = latitudinal means - ld latitudes,kd layers..                
c          colrad = co-latiudes in radians (n.h.)                               
c          ltitl  = printing label (character*20)                               
c          fhour  = forecast hour                                               
c          ipts   = number of points on lat circle for average..                
c..   output:   printed values (ft06f001)                                       
c            globl(1,kd)=global mean                                            
c            globl(2,kd)=northern hemisphere mean                               
c            globl(3,kd)=northern hemisphere polar (90-60deg) mean              
c            globl(4,kd)=northern hemisphere midlat(60-30deg) mean              
c            globl(5,kd)=northern hemisphere tropic(30-0 deg) mean              
c            globl(6,kd)=southern hemisphere mean                               
c            globl(7,kd)=southern hemisphere polar (90-60deg) mean              
c            globl(8,kd)=southern hemisphere midlat(60-30deg) mean              
c            globl(9,kd)=southern hemisphere tropic(30-0 deg) mean              
      dimension y(kd,ld),colrad(ld)                                             
      character*20 ltitl                                                        
      dimension globl(9,kd),rlat(ld),wgt(ld),wtsum(9)                           
      dimension ipts(ld),wgta(ld)                                               
      dimension globx(9,100)                                                    
      data latlon/ 0 /                                                          
      jj = ld / 2                                                               
      kk = kd                                                                   
cc    if (latlon.gt.0) go to 30                                                 
c...   compute area wgts and latitudes just poleward of 60deg (j60)             
c         and just poleward of 30deg (j30)                                      
      rad60 = 3.14159265 / 3.                                                   
      rad30 = 3.14159265 / 6.                                                   
c...   latitude in radians                                                      
      do 5 j=1,jj                                                               
        rlat(j) =  3.14159265 / 2. - colrad(j)                                  
    5 continue                                                                  
      do 10 j=1,jj                                                              
        if(rlat(j).lt.rad60) go to 15                                           
        j60 = j                                                                 
   10 continue                                                                  
      go to 200                                                                 
   15 do 20 j=1,jj                                                              
        if(rlat(j).lt.rad30) go to 25                                           
        j30 = j                                                                 
   20 continue                                                                  
      go to 205                                                                 
   25 if (latlon.le.0) print 9,jj,j60,j30                                       
    9 format(1h ,' jj=',i4,'  j60=',i4,'  j30=',i4)                             
      latlon = 10                                                               
cc 30 continue                                                                  
c...                                                                            
      do 185 k=1,kk                                                             
        do 185 i=1,9                                                            
          globx(i,k) = 0.                                                       
          globl(i,k) = 0.                                                       
  185 continue                                                                  
      do 26 n=1,9                                                               
        wtsum(n) = 0.                                                           
   26 continue                                                                  
      do 27 la=1,jj                                                             
        rlatn = 3.14159265e0 / 2.e0                                             
        if (la.gt.1) rlatn = 0.5e0 * (rlat(la)+rlat(la-1))                      
        rlats = 0.e0                                                            
        if (la.lt.jj) rlats = 0.5e0 * (rlat(la)+rlat(la+1))                     
        wgt(la) = abs(sin(rlatn)-sin(rlats))                                    
        wgt(ld+1-la) = abs(sin(rlatn)-sin(rlats))                               
   27 continue                                                                  
c...    compute area for row (allows for changing num pts)                      
      do 29 la=1,ld                                                             
        wgta(la) = wgt(la) * ipts(la)                                           
   29 continue                                                                  
      do 28 la=1,jj                                                             
c..    hemispheric area                                                         
        wtsum(2) = wtsum(2) + wgta(la)                                          
        wtsum(6) = wtsum(6) + wgta(ld+1-la)                                     
c..    latitude strip areas                                                     
        if (la.gt.j60) go to 31                                                 
        wtsum(3) = wtsum(3) + wgta(la)                                          
        wtsum(7) = wtsum(7) + wgta(ld+1-la)                                     
        go to 28                                                                
   31   if (la.gt.j30) go to 32                                                 
        wtsum(4) = wtsum(4) + wgta(la)                                          
        wtsum(8) = wtsum(8) + wgta(ld+1-la)                                     
        go to 28                                                                
   32   wtsum(5) = wtsum(5) + wgta(la)                                          
        wtsum(9) = wtsum(9) + wgta(ld+1-la)                                     
   28 continue                                                                  
c..    global area                                                              
        wtsum(1) = wtsum(2) + wtsum(6)                                          
      do 35 k = 1,kk                                                            
        do 35 j = 1,9                                                           
          globl(j,k) = 0.e0                                                     
   35 continue                                                                  
      do 75 la=1,jj                                                             
       do 45 k=1,kk                                                             
        globl(2,k) = globl(2,k) + wgt(la)*y(la,k)                               
        globl(6,k) = globl(6,k) + wgt(ld+1-la)*y(ld+1-la,k)                     
        globl(1,k) = globl(2,k) + globl(6,k)                                    
   45  continue                                                                 
       if (la.gt.j60) go to 55                                                  
       do 50 k=1,kk                                                             
        globl(3,k) = globl(3,k) + wgt(la)*y(la,k)                               
        globl(7,k) = globl(7,k) + wgt(ld+1-la)*y(ld+1-la,k)                     
   50  continue                                                                 
       go to 75                                                                 
   55  if (la.gt.j30) go to 65                                                  
       do 60 k=1,kk                                                             
        globl(4,k) = globl(4,k) + wgt(la)*y(la,k)                               
        globl(8,k) = globl(8,k) + wgt(ld+1-la)*y(ld+1-la,k)                     
   60  continue                                                                 
       go to 75                                                                 
   65  do 70 k=1,kk                                                             
        globl(5,k) = globl(5,k) + wgt(la)*y(la,k)                               
        globl(9,k) = globl(9,k) + wgt(ld+1-la)*y(ld+1-la,k)                     
   70  continue                                                                 
   75 continue                                                                  
      do 85 k=1,kk                                                              
        do 80 i=1,9                                                             
ccc       globx(i,k) = globl(i,k) / wtsum(1)                                    
          if (wtsum(i).gt.0.) globl(i,k) = globl(i,k) / wtsum(i)                
   80   continue                                                                
   85 continue                                                                  
      print 93                                                                  
      print 92,ltitl,fhour                                                      
      print 97                                                                  
      print 95                                                                  
      do 90 ken=1,kk                                                            
        k = kk+1-ken                                                            
        print 100,k,(globl(i,k),i=1,9)                                          
   90 continue                                                                  
      print 93                                                                  
ccc   print 98                                                                  
ccc   print 95                                                                  
ccc   do 11 ken=1,kk                                                            
ccc     k = kk+1-ken                                                            
ccc     print 100,k,(globx(i,k),i=1,9)                                          
cc 11 continue                                                                  
ccc   print 93                                                                  
      return                                                                    
   92 format (1h ,10x,' area weighted means',2x,a20,' fcsthr=',f6.1)            
   93 format (1h ,11x,' -------------------------------------')                 
   95 format (1h ,4x,'k',5x,'global',4x,'n hemi',4x,'nhpolr',4x,                
     1        'nhmidl',4x,'nhtrpc',4x,'s hemi',4x,                              
     2        'shpolr',4x,'shmidl',4x,'shtrpc')                                 
   97 format (1h ,15x,'====>  weighted relative to specified region')           
cc 98 format (1h ,15x,'====>  weighted relative to entire globe')               
  100 format (1h ,i5,9f10.3)                                                    
  200 print 201                                                                 
      stop                                                                      
  205 print 206                                                                 
      stop                                                                      
  201 format (1h ,' stop -- can t find j60')                                    
  206 format (1h ,' stop -- can t find j30')                                    
      end                                                                       
