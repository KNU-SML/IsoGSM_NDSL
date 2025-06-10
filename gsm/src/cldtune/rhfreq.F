#define LATG latg
#define LONF lonf
#include <define.h>
      subroutine rhfreq(slmask,n,nfiles,fhrmrf,idate,rhmdl)                     
#include <paramodel.h>
      save                                                                      
c-    *****************************************************************         
c     *     driver code for computing cumulative frequency distrib    *         
c-----******* of model rh(max) only                                   *         
c-----*******(unlike opnl codes,which didnot use data if rtn missing) *         
c-----*******             from .radbench(gcfrqbl4)..kac..oct 96       *         
c     *          for mitchell/hahn method of tuning cloud             *         
c     *          parameterization with cloud observations..           *         
c     *     do for lon,lat array covring gauss grid..land/sea seperate*         
c     *          1.read model max rh in cloud domains  ,3 cld lyrs hml*         
c     *          2.compute frequency distribution                     *         
c     *          3.save raw rh and rtneph distribution + numpts       *         
c     *          4. nbn=100, jlon=2, jlat=4, nsl=2                    *         
c-    *****************************************************************         
      parameter (idim=lonf_,jdim=latg_,idim1=idim+1,jdim2=jdim+2,               
     1           mcld=3,                                                        
     2           nbin=100,nlon=2,nlat=4,nseal=2)                                
      dimension idate(4),kidate(4)                                              
      dimension rtcld (idim,jdim,mcld)                                          
      dimension rhmdl (idim,jdim,mcld)                                          
      dimension rhfd(nbin,nlon,nlat,mcld,nseal)                                 
      dimension rhcf(nbin,nlon,nlat,mcld,nseal)                                 
      dimension kpts(nlon,nlat,mcld,nseal)                                      
      dimension kmiss(nlon,nlat,mcld,nseal)                                     
      dimension ilon(nlon,2),jlat(nlat,2)                                       
      dimension jbi(nbin)                                                       
      dimension yave(jdim),globl(9),rlat(jdim),wgt(jdim)                        
      dimension jpts(jdim),wgta(jdim)                                           
      parameter(numdays=500)                                                    
      dimension icdays(numdays),kpttot(mcld)                                    
      dimension colrad(jdim)                                                    
      dimension colr(jdim),xlat(jdim)                                           
      dimension mdays(12)                                                       
      character*1 ptitle(20)                                                    
      dimension slmask(idim,jdim)                                               
      character*12 lbl,ltitle(mcld,3)                                           
      data mdays / 31,28,31,30,31,30,31,31,30,31,30,31 /                        
      data (ptitle(k),k=1,20) /' ',' ',' ',' ',' ',' ',' ',' ',                 
     1         ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/                 
      data ltitle/                                                              
     1  'rtnef lowcld','rtnef midcld','rtnef hi cld',                           
     1  'mrf hi  cld ','mrf mid cld ','mrf low cld ',                           
     1  'mrf low -rh-','mrf mid -rh-','mrf hi -rh- '/                           
      equivalence(lbl,ptitle(1))                                                
c                                                                               
      londia=1                                                                  
      latdia=3                                                                  
c                                                                               
      do j=1,jdim                                                               
        do i=1,idim                                                             
          slmask(i,j)=nint(slmask(i,j))                                         
        enddo                                                                   
      enddo                                                                     
c     call nntprt(slmask,idim,jdim,1.)                                          
c     rewind 10                                                                 
      khalf = jdim/2                                                            
      xhalf = khalf                                                             
      call glatss(khalf,colrad)                                                 
      do 152 l=1,khalf                                                          
        colr(l) = colrad(l)                                                     
        xlat(l) = (3.14159265/2. - colr(l))*180./3.14159265                     
  152 continue                                                                  
      ii = idim                                                                 
      jj = jdim                                                                 
c                                                                               
c     kk = kdim                                                                 
c find good estimate of 60deg lat row                                           
c                                                                               
      limlat = 60.                                                              
      xlatin = 90./xhalf                                                        
      latlim = (90.-limlat)/xlatin                                              
#ifdef DBG
      print 100,latlim                                                          
      print 101,(xlat(l),l=1,khalf)                                             
#endif
  100 format(1h1,'....row number of estimated 60deg =',i6)                      
  101 format(1h ,'lats....',15f7.1)                                             
c                                                                               
c compute left,right and top,bottom values of columns and rows                  
c    on gaussian grid for the (2,4) array covering 0-360,60n-60s                
c                                                                               
      klon = idim/nlon                                                          
      klat = (khalf-(latlim-1))*2/nlat                                          
      print 102,klon,klat                                                       
  102 format(1h1,'..column increment,klon=',i4,'..row incr,klat=',i4)           
      do 17 i=1,nlon                                                            
       ilon(i,1) = klon*(i-1)+1                                                 
       ilon(i,2) = ilon(i,1) + klon                                             
       if (i.eq.nlon) ilon(i,2) = idim                                          
   17 continue                                                                  
      do 18 j=1,nlat                                                            
       jlat(j,1) = latlim + (j-1)*klat                                          
       jlat(j,2) = jlat(j,1) + klat                                             
       if (j.eq.nlat) jlat(j,2) = jdim+1-latlim                                 
   18 continue                                                                  
      print 103,(ilon(l,1),l=1,nlon)                                            
      print 104,(ilon(l,2),l=1,nlon)                                            
      print 105,(jlat(l,1),l=1,nlat)                                            
      print 106,(jlat(l,2),l=1,nlat)                                            
  103 format(1h ,'...left column boundary=',2i5)                                
  104 format(1h ,'...rite column boundary=',2i5)                                
  105 format(1h ,'...top row boundary  =',4i5)                                  
  106 format(1h ,'...bot row boundary  =',4i5)                                  
c                                                                               
c     rhm = 76                                                                  
c     read (irhm) fhrmrf,idate                                                  
c                                                                               
      write(6,*) 'fhrmrf,idate from unit 76:',fhrmrf,idate                      
c                                                                               
      if (fhrmrf.lt.0.) stop                                                    
c                                                                               
c compute number of days since 1 jan 92, in order to see if                     
c verifying data is okay...assume initial hour = 00z                            
c                                                                               
      iy = idate(4) + 1900                                                      
      im = idate(2)                                                             
      id = idate(3)                                                             
      ifda = (fhrmrf+0.5)/24.-1                                                 
c                                                                               
c  read in h,m,l rhmaxx on gg                                                   
c   by prev step..store as l,m,h                                                
c                                                                               
c     do 2901 k=1,mcld                                                          
c                                                                               
      do 2901 kmc=mcld,1,-1                                                     
c                                                                               
c  k=1 is top, i.e., k=kmc is bottom                                            
c                                                                               
c      kmc = mcld+1-k                                                           
c      read (irhm) ((rhmdl(i,j,kmc),i=1,idim),j=1,jdim)                         
c                                                                               
       do 713 i=1,nbin                                                          
        jbi(i) = 0                                                              
  713  continue                                                                 
c                                                                               
       jmis  = 0                                                                
       do 723 j=1,jdim                                                          
        do 723 i=1,idim                                                         
         if (rhmdl(i,j,kmc).lt.0.) then                                         
          jmis = jmis + 1                                                       
         else                                                                   
          ibi = rhmdl(i,j,kmc)*100.+1                                           
          if (ibi.ge.nbin) ibi=nbin                                             
          jbi(ibi) = jbi(ibi) + 1                                               
         end if                                                                 
  723  continue                                                                 
       lbl = ltitle(kmc,3)                                                      
       call avlat(yave,rhmdl(1,1,kmc),idim,jdim,1,kmc,jpts)                     
       call globlm(yave,globl,rlat,wgt,wgta,                                    
     1              1,jdim,colr,ptitle,fhrmrf,jpts)                             
       print 731,kmc,jmis                                                       
  731  format(1h ,'kmc=',i1,                                                    
     1        ' global missing due to underground or vvel=',i10)                
#ifdef DBG
       print 732,jbi                                                            
#endif
  732  format(1h ,10i10)                                                        
 2901 continue                                                                  
c                                                                               
c... position the file..if already 15 days strip off oldest (1st), adj icdays.. 
c                                                                               
      icfrq = 66                                                                
      if(n.eq.1) rewind icfrq                                                   
c                                                                               
      icdays(1) = idate(3)*10000 + idate(2)*100 + idate(4)                      
      write (icfrq) nfiles,icdays                                               
      write (icfrq) fhrmrf,idate                                                
      write(6,*) 'fhrmrf=',fhrmrf,' idate=',idate                               
c                                                                               
c here we have bl,l,m,h rhmax in rhmdl...                                       
c     so get cumulative frequency distribution..rhcf                            
c                                                                               
      call cfreq(rhmdl,idim,jdim,mcld,ilon,jlat,nlon,nlat,                      
     1           rtcld,rhcf,nbin,kpts,icfrq,nseal,slmask,londia,latdia)         
c                                                                               
c                                                                               
c  raw distribution is written for                                              
c      rhcf during the call to cfreq                                            
c                                                                               
      write (icfrq) kpts                                                        
c                                                                               
      do 70 kcl=1,mcld                                                          
       kpttot(kcl) = 0                                                          
       do 70 l=1,nlat                                                           
        do 70 i=1,nlon                                                          
         kpttot(kcl) = kpttot(kcl)+kpts(i,l,kcl,1)+kpts(i,l,kcl,2)              
   70 continue                                                                  
      print 198,(kpttot(kcl),kcl=1,mcld)                                        
  198 format(1h ,' total num gg pts used globally..for type 1-mcld=',           
     1           3i8)                                                           
  200 format(1h ,' lon band=',i4,' latitude band=',i4,                          
     1           'cld type (low=1) =',i4,' sea-land=',i2)                       
  201 format(1h0,'num of pts used in box=',i8)                                  
  202 format(1h0,' model rh ',' obs rtcld')                                     
  203 format(2f10.2)                                                            
      do 255 nsl=1,2                                                            
      do 255 kcl=1,mcld                                                         
       print 264,kcl,nsl                                                        
       print 265,((kpts(i,l,kcl,nsl),i=1,nlon),l=1,nlat)                        
  255 continue                                                                  
  264 format(1h ,' number of gg points used in each area..by latitude',         
     1           '..for cloud type=',i4,'sealand=',i2)                          
  265 format(1h ,20i6)                                                          
  210 format(1h ,' no crit rh for lat=',i3,' and lon band=',i3,                 
     1           ' land(=1) sea=',i3)                                           
 1221 format(1h0,' critical rh for lon,lat arrays for cld type=',i3,            
     1           ' land(=1) sea=',i3)                                           
  211 format(1h ,15f6.2)                                                        
  412 format(1h ,'can t find rh for bin=',i4,' cld=',f10.2)                     
  419 format(1h1,'cloud typ= ',i4,' lon,lat bx=',2i4,'-landsea=',i2)            
  420 format(1h ,'bin=',i3,'-- camt=',f10.2,'-- rh=',f10.2)                     
      return                                                                    
 999  continue                                                                  
      write(6,*) 'error in reading unit 66'                                     
      call abort                                                                
  913 format(1h ,' -t-  k = ',i4,' max = ',e15.5,' min =',e15.5)                
  914 format(1h ,' -q-  k = ',i4,' max = ',e15.5,' min =',e15.5)                
  915 format(1h ,' -rh- k = ',i4,' max = ',e15.5,' min =',e15.5)                
 3001 format(5x,'input label= ',a8,2x,2z16,2x,a8)                               
 3002 format(5x,'input data dates time and forecast hour',/,10x,                
     1       4i15,f7.1)                                                         
      end                                                                       
