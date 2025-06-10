      subroutine cfreq(yin,iin,jin,kin,ilon,jlat,mlon,mlat,                     
     1               xin,cfrq,nbin,kpts,iunit,nsl,slmask,londia,latdia)         
c-mk 1               xin,cfrq,nbin,kpts,iunit,nsl,cland,londia,latdia)          
      save                                                                      
c-    *****************************************************************         
c     *     compute cumulative frequency distribution for input--     *         
c     *          yin = fraction(rh or cld) 0.0-1.0                    *         
c     *          iin,jin,kin = dimensions of yin                      *         
c     *          ilon = left and right lon of bins                    *         
c     *          jlat = top and bott lat of bins                      *         
c     *          mlon = lon dimension for regions                     *         
c     *          mlat = lat dimension for regions                     *         
c     *          nbin = number of bins for 0.0-1.0                    *         
c     *            example: mlon=20,mlat=10(70n-70s)..region=19deg    *         
c     *                       thus regions approx 19deg x 14deg       *         
c     *                     nbin=100 for binning every .01            *         
c     *          cland=land sea array..land=x, sea = ' '              *         
c     *         slmask=land sea array..land=1., sea=0.                *         
c     *      for output--                                             *         
c     *          cfrq = freq distrib (i,lo,la,k)                      *         
c     *            i=bin,lo=bin lon column,la=bin lat row,k=cld type  *         
c     *          iunit= disk to store raw distribution                *         
c     *      workspace--                                              *         
c     *          xin                                                  *         
c-    *          kpts(number of points in each region                 *         
c-    *****************************************************************         
      dimension yin(iin,jin,kin)                                                
      dimension xin(iin,jin,kin)                                                
c-mk  character*1 cland(iin,jin)                                                
      dimension   slmask(iin,jin)                                               
      dimension cfrq(nbin,mlon,mlat,kin,nsl)                                    
      dimension kpts(mlon,mlat,kin,nsl),ilon(mlon,2),jlat(mlat,2)               
      parameter (mcld=3,nlon=2,nlat=4,nseal=2)                                  
      dimension lmiss(nlon,nlat,mcld,nseal),kmiss(nlon,nlat,mcld,nseal)         
      character*1 land,sea                                                      
      data land/'x'/,sea/' '/                                                   
      do 5 n=1,nsl                                                              
      do 5 k=1,kin                                                              
       do 5 l=1,mlat                                                            
        do 5 j=1,mlon                                                           
         kpts(j,l,k,n) = 0                                                      
         kmiss(j,l,k,n) = 0                                                     
         lmiss(j,l,k,n) = 0                                                     
    5 continue                                                                  
      do 10 n=1,nsl                                                             
      do 10 k=1,kin                                                             
       do 10 l=1,mlat                                                           
        do 10 j=1,mlon                                                          
         do 10 i=1,nbin                                                         
          cfrq(i,j,l,k,n) = 0.                                                  
   10 continue                                                                  
c....  convert fraction to whole percent                                        
      do 15 k=1,kin                                                             
       do 15 j=1,jin                                                            
        do 15 i=1,iin                                                           
         yin(i,j,k) = yin(i,j,k)*100.                                           
   15 continue                                                                  
      xbscl = nbin/100.                                                         
      do 100 k=1,kin                                                            
       imiss = 0                                                                
       do 50 j=1,jin                                                            
        do 40 i=1,iin                                                           
c----   find out which region this point belongs                                
         ireg = 0                                                               
         jreg = 0                                                               
         do 35 jk=1,mlat                                                        
          if (j.ge.jlat(jk,1).and.j.lt.jlat(jk,2))then                          
           jreg = jk                                                            
           go to 36                                                             
          end if                                                                
   35    continue                                                               
   36    do 37 ik=1,mlon                                                        
          if (i.ge.ilon(ik,1).and.i.lt.ilon(ik,2))then                          
           ireg = ik                                                            
           go to 38                                                             
          end if                                                                
   37    continue                                                               
   38    if (ireg.eq.0.or.jreg.eq.0) then                                       
cccc       print 27,i,j,k                                                       
ccc27 format(1h ,'cannot region-place a gauss grid point, i,j,k=',3i4)          
           go to 40                                                             
         end if                                                                 
ccccc    if either data is missing (xin,yin lt 0.) do not use....               
ccc      if (yin(i,j,k).lt.0..or.xin(i,j,k).lt.0.) then                         
c....    if data is missing (yin lt 0.) do not use....                          
         isl=2                                                                  
         if (slmask(i,j).eq.1.) isl=1                                           
         if (yin(i,j,k).lt.0..and.yin(i,j,k).gt.-150.) then                     
          lmiss(ireg,jreg,k,isl)=lmiss(ireg,jreg,k,isl)+1                       
          imiss = imiss + 1                                                     
          go to 40                                                              
         elseif (yin(i,j,k).lt.-150.) then                                      
          kmiss(ireg,jreg,k,isl)=kmiss(ireg,jreg,k,isl)+1                       
          imiss = imiss + 1                                                     
          go to 40                                                              
         endif                                                                  
c.                                                                              
c----   find out which bin (i) for cfrq..take adv of xin=0,100                  
         nbnb = yin(i,j,k)*xbscl + 1                                            
         if(nbnb.gt.nbin)   nbnb=nbin                                           
         cfrq(nbnb,ireg,jreg,k,isl) = cfrq(nbnb,ireg,jreg,k,isl) + 1.           
c----   keep track of number of points during (k=1)                             
         kpts(ireg,jreg,k,isl) = kpts(ireg,jreg,k,isl) + 1                      
cc       if (ireg.eq.londia.and.jreg.eq.latdia.and.k.eq.1) then                 
cc        print 2020,ireg,jreg,k,isl,yin(i,j,k),xin(i,j,k)                      
c2020     format(1h ,' i,j,k=',3i4,'...lansea=',i2,' input=',f10.3,             
cc   1               '..other=',f10.4)                                          
cc       end if                                                                 
   40   continue                                                                
   50  continue                                                                 
       print 2050,imiss,k                                                       
 2050  format(1h ,' number pts classified as missing=',i8,'..lyr=',i3)          
  100 continue                                                                  
c                                                                               
      do k=1,kin                                                                
       do n=1,nsl                                                               
        print *,'cloud level=',k,' sealand=',n                                  
        do l=1,mlat                                                             
         do j=1,mlon                                                            
           print *,'  lat=',l,' lon=',j,                                        
     1             ' undmtn=',lmiss(j,l,k,n),                                   
     2             ' vvelfilt=',kmiss(j,l,k,n)                                  
         enddo                                                                  
        enddo                                                                   
       enddo                                                                    
      enddo                                                                     
c                                                                               
      write (iunit) cfrq                                                        
      write (iunit) kmiss                                                       
c                                                                               
      return                                                                    
      end                                                                       
