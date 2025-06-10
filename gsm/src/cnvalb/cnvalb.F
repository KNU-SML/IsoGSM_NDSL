#define LATG latg
#define LONF lonf
      program cnvalb                                                            
c$$$  main program documentation block
c
c main program:  cnvalb  reads vegitation and maps to gaussin grid
c   prgmmr: cpc          org: w/np51     date: 01-03-31
c
c abstract:
c***************************************************************                
c   this program reads modified breigleb vegitation type data                   
c     file in one degree lat/lon resolution, then converts to                   
c     surface albedo based on brieglebs scheme, and maps                       
c     the data to gaussin grid for forecast model use.                          
c                                                                               
c   input data:                                                                 
c     vegitation index file                  - unit 11                          
c                                                                               
c   output data:                             - unit 51                          
c     vis spectrum albedo (strong cosz dep.)                                    
c     vis spectrum albedo (weak cosz dep.)                                      
c     nir spectrum albedo (strong cosz dep.)                                    
c     nir spectrum albedo (weak cosz dep.)                                      
c     fraction of coverage for strong dep. in each grid                         
c     fraction of coverage for weak depen. in each grid                         
c                                                                               
c   input namelist parameters:                                                  
c     imax,jmax      - output data i,j (lon/lat) dimensions                     
c     igau           - =1 output on gaussing grid                               
c                      =0 output on regular lat/lon grid                        
c     kssn           - n.h. season index 1-4 (w,s,s,f); =5 annual mean.         
c                      =0 produce win,spr,sum,and fal all at once               
c     clab           - label for output title (a4)                              
c     ichk           - =1 print out check point values                          
c                      =0 do not print                                          
c     iout           - =1 output formatted data, scale 1e3                      
c                      =0 output unformatted data in fraction                   
c     outdr          - output data file name and path (a48)                     
c****************************************************************               
c                                                                               
c attributes:
c   language: fortran
c
c$$$
c    1 ( idm1=360,  jdm1=180,  j30=61,  idm2=384,  jdm2=196,                    
#include <paramodel.h>
                   p a r a m e t e r                                            
     1 ( idm1=360,  jdm1=180,  j30=61,  idm2=lonf_,  jdm2=latg_,                
     2   idm1p=idm1+1, jdm1p=jdm1+1,  jdm1h=jdm1/2,  nmax=13,                   
     3   idm2p=idm2+1, jdm2p=jdm2+1,  jdm2h=jdm2/2,                             
     4   ndr=48, nfl=60 )                                                       
                   d i m e n s i o n                                            
     1  albv1(nmax,5),albv2(idm2,jdm2,3),idx (idm1,jdm1)                        
     2, albn1(nmax,5),albn2(idm2,jdm2,3),frac(idm2,jdm2,3)                      
     3, iswx (nmax),  idata(idm2,jdm2),  issx(idm1,jdm1)                        
                   d i m e n s i o n                                            
     1  rlon(idm1p), rlat(jdm1p), alon(idm2p), alat(jdm2p)                      
     2, blat(jdm2), ir1(idm2), ir2(idm2), jr1(jdm2), jr2(jdm2)                  
                   c h a r a c t e r                                            
     1  cdat(idm1,jdm1)*1, ch1*1, title(2)*40, clab*4, cssn(5)*3                
     2, dssn(5)*3, outdr*48, cttv(3)*12, cttn(3)*12, cttf(3)*12                 
     3, cfl*60                                                                  
c                                                                               
      data imax/lonf_/,jmax/latg_/                                              
      data igau/1/,kssn/0/,clab/'    '/,ichk/1/,iout/1/                         
      data outdr/'./alb_out.asc'/                                               
c                                                                               
                   namelist /namein/                                            
     1  imax, jmax, igau, kssn, clab, ichk, iout, outdr                         
c                                                                               
      data  albv1                                                               
c    1  / .06,.09,.10,.10,.09,.35,.24,.06,.06,.07,.08,.08,.90                   
c    2,   .06,.10,.11,.12,.08,.35,.24,.06,.06,.06,.05,.05,.90                   
c    3,   .06,.09,.11,.14,.10,.35,.24,.06,.07,.07,.06,.06,.90                   
c    4,   .06,.09,.10,.13,.10,.35,.24,.06,.06,.05,.05,.05,.90                   
c    5,   .06,.09,.11,.14,.10,.35,.24,.06,.06,.06,.06,.06,.90 /                 
c2   1  / .06,.08,.09,.09,.15,.34,.21,.06,.06,.06,.05,.06,.75                   
c2   2,   .06,.08,.10,.11,.13,.35,.21,.06,.06,.06,.06,.06,.75                   
c2   3,   .06,.09,.10,.12,.09,.36,.22,.06,.06,.06,.06,.06,.75                   
c2   4,   .06,.08,.10,.11,.12,.34,.21,.06,.06,.06,.06,.06,.75                   
c2   5,   .06,.08,.10,.11,.12,.35,.21,.06,.06,.06,.06,.06,.75 /                 
     1  / .06,.08,.09,.12,.14,.35,.20,.06,.06,.07,.06,.06,.75                   
     2,   .06,.08,.10,.12,.12,.35,.20,.06,.06,.06,.06,.06,.75                   
     3,   .06,.08,.10,.13,.08,.35,.20,.06,.06,.05,.06,.06,.75                   
     4,   .06,.08,.10,.12,.12,.35,.20,.06,.06,.06,.06,.06,.75                   
     5,   .06,.08,.10,.12,.12,.35,.20,.06,.06,.06,.06,.06,.75 /                 
      data  albn1                                                               
c    1  / .06,.27,.29,.28,.22,.50,.40,.20,.20,.24,.32,.30,.65                   
c    2,   .06,.30,.35,.32,.24,.50,.40,.20,.21,.26,.24,.22,.65                   
c    3,   .06,.30,.35,.40,.31,.50,.40,.20,.25,.32,.31,.29,.65                   
c    4,   .06,.28,.32,.33,.27,.50,.40,.20,.19,.22,.24,.22,.65                   
c    5,   .06,.28,.33,.34,.27,.51,.40,.22,.20,.28,.31,.27,.65 /                 
c    1  / .06,.21,.25,.25,.34,.48,.33,.22,.22,.21,.26,.27,.60                   
c    2,   .06,.24,.28,.28,.28,.50,.33,.22,.23,.23,.26,.25,.60                   
c    3,   .06,.25,.30,.31,.21,.51,.35,.22,.24,.25,.31,.27,.60                   
c    4,   .06,.24,.28,.28,.28,.50,.33,.22,.23,.23,.27,.26,.60                   
c    5,   .06,.24,.28,.28,.28,.50,.33,.22,.23,.23,.28,.26,.60 /                 
     1  / .06,.19,.25,.26,.30,.45,.30,.24,.22,.21,.27,.23,.50                   
     2,   .06,.22,.28,.28,.25,.50,.32,.23,.22,.21,.29,.25,.50                   
     3,   .06,.24,.30,.30,.18,.50,.30,.22,.21,.21,.31,.26,.50                   
     4,   .06,.22,.28,.28,.25,.50,.32,.23,.22,.21,.29,.25,.50                   
     5,   .06,.22,.28,.28,.25,.50,.31,.23,.22,.21,.29,.25,.50 /                 
      data  iswx / 3,1,1,1,1,1,1,2,2,2,2,2,2 /                                  
      data  cssn / 'win','spr','sum','fal','anu' /                              
     1,     dssn / 'win','spr','sum','fal','anu' /                              
      data  cttv /'vis s albedo','vis w albedo','vis o albedo'/                 
     1,     cttn /'nir s albedo','nir w albedo','nir o albedo'/                 
     2,     cttf /'frac cover s','frac cover w','frac cover o'/                 
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_cnvalb',2001,0000,0000,'np51')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c --- preset default seasonal index = 5 (anu)                                   
c                                                                               
      do 10 j=1,jdm1                                                            
      do 10 i=1,idm1                                                            
        issx(i,j) = 5                                                           
  10  continue                                                                  
c                                                                               
c --- read in namelist                                                          
c                                                                               
      read (5,namein)                                                           
      write(6,namein)                                                           
      if (imax.gt.idm2 .or. jmax.gt.jdm2) then                                  
        write(6,12) imax,jmax,idm2,jdm2                                         
  12    format(2x,'defined output array size exceeds limits',                   
     1        /4x,'imax,jmax=',2i4,'  idm2,jdm2=',2i4)                          
        call abort                                                              
      end if                                                                    
      if (kssn.lt.0 .or. kssn.gt.5) then                                        
        write(6,14) kssn                                                        
  14    format(2x,'season index out off range, kssn =',i4)                      
        call abort                                                              
      end if                                                                    
c                                                                               
      imaxp = imax + 1                                                          
      jmaxp = jmax + 1                                                          
      pi = atan(1.0e0) * 4.0e0                                                  
      tpi = 2.0 * pi                                                            
      hpi = 0.5 * pi                                                            
      raddeg = 180.0 / pi                                                       
      degrad = pi / 180.0                                                       
c                                                                               
c --- read in vegtype index file, convert char to integer                       
c                                                                               
      read(11,16) title                                                         
      write(6,16) title                                                         
  16  format(2a40)                                                              
      read(11,20) ((cdat(i,j),i=1,idm1),j=1,jdm1)                               
  20  format(80a1)                                                              
      do 25 j=1,jdm1                                                            
      do 25 i=1,idm1                                                            
        ch1 = cdat(i,j)                                                         
        if      (ch1.ge.'0' .and. ch1.le.'9') then                              
          idx(i,j) = ichar(cdat(i,j)) - ichar('0') + 1                          
        else if (ch1.ge.'A' .and. ch1.le.'C') then                              
          idx(i,j) = ichar(cdat(i,j)) - ichar('A') + 11                         
        else                                                                    
          write(6,22) i,j, cdat(i,j)                                            
  22      format(3x,'error in input index file: i,j =',2i4,                     
     1           3x,'cdat(i,j) = ',a1)                                          
          call abort                                                            
        end if                                                                  
  25  continue                                                                  
c --- check output file name and path                                           
      nn = ndr                                                                  
  28  if (outdr(nn:nn) .ne. ' ') go to 32                                       
      if (nn .eq. 1) then                                                       
        write(6,30)                                                             
  30    format(3x,'empty outdr - supply one')                                   
        call abort                                                              
      end if                                                                    
      nn = nn - 1                                                               
      go to 28                                                                  
  32  continue                                                                  
c                                                                               
c --- set up input data lon, lat of each grid box edges in radians              
c     matthews data are at the center of grid boxes, so the edge               
c     of the first box is at the 0 deg                                          
        drlon = tpi / float(idm1)                                               
        drlat = pi / float(jdm1)                                                
c                                                                               
      do 35 i=1,idm1p                                                           
        rlon(i) = drlon * (i - 1)                                               
  35  continue                                                                  
c --- the edge of the first lat box is at 90 deg                                
      do 40 j=1,jdm1p                                                           
        rlat(j) = hpi - drlat * (j - 1)                                         
  40  continue                                                                  
c --- check print                                                               
      if (ichk .eq. 1) then                                                     
        write(6,42)                                                             
  42    format(/3x,'input lon values:')                                         
        write(6,44) (raddeg*rlon(i),i=1,idm1p)                                  
  44    format(16f5.1)                                                          
        write(6,46)                                                             
  46    format(/3x,'input lat values:')                                         
        write(6,44) (raddeg*rlat(j),j=1,jdm1p)                                  
      end if                                                                    
c --- rem: model lon starts at greenwich, so the edge of                        
c     the first lon box is at the -0.5 of the box width                         
      dalon = tpi / float(imax)                                                 
      do 50 i=1,imaxp                                                           
        alon(i) = dalon * (i - 1.5)                                             
  50  continue                                                                  
c --- blat contains model lat grid                                              
      if (igau .eq. 1) then                                                     
        call glats(jmax,blat)                                                   
        alat(1) = hpi                                                           
        do 60 j=1,jmax-1                                                        
          alat(j+1) = 0.5 * (blat(j) + blat(j+1))                               
  60    continue                                                                
        alat(jmaxp) = -hpi                                                      
      else                                                                      
        dalat = pi / float(jmax)                                                
        do 70 j=1,jmaxp                                                         
          alat(j) = hpi - dalat * (j - 1)                                       
          if (j.le.jmax) blat(j) = alat(j) - 0.5*dalat                          
  70    continue                                                                
      end if                                                                    
c --- check print                                                               
      if (ichk .eq. 1) then                                                     
        write(6,72)                                                             
  72    format(/3x,'output lon values (box edges):')                            
        write(6,44) (raddeg*alon(i),i=1,imaxp)                                  
        write(6,74)                                                             
  74    format(/3x,'output lat values (box centers):')                          
        write(6,44) (raddeg*blat(j),j=1,jmax)                                   
        write(6,76)                                                             
  76    format(/3x,'output lat values (box edges):')                            
        write(6,44) (raddeg*alat(j),j=1,jmaxp)                                  
      end if                                                                    
c                                                                               
c --- find starting and ending indeces of the original                          
c     1 deg data for each of the output grid boxes                              
      ir = 1                                                                    
      do 100 ia=1,imax                                                          
        ir1(ia) = 0                                                             
        ir2(ia) = 0                                                             
        if (alon(ia) .lt. 0.0) then                                             
          xlon = tpi + alon(ia)                                                 
          ir0 = idm1                                                            
  80      if (xlon.ge.rlon(ir0) .and.                                           
     1        xlon.le.rlon(ir0+1)) then                                         
            ir1(ia) = ir0                                                       
            go to 90                                                            
          else                                                                  
            ir0 = ir0 - 1                                                       
            go to 80                                                            
          end if                                                                
        end if                                                                  
c                                                                               
  90    if (ir1(ia).eq.0 .and. alon(ia).ge.rlon(ir)                             
     1                   .and. alon(ia).le.rlon(ir+1)) then                     
            ir1(ia) = ir                                                        
        end if                                                                  
        if (ir2(ia).eq.0 .and. alon(ia+1).ge.rlon(ir)                           
     1                   .and. alon(ia+1).le.rlon(ir+1)) then                   
          ir2(ia) = ir                                                          
          go to 100                                                             
        else                                                                    
          if (ir .ge. idm1) then                                                
            ir2(ia) =ir                                                         
            go to 120                                                           
          else                                                                  
            ir = ir + 1                                                         
            go to 90                                                            
          end if                                                                
        end if                                                                  
 100  continue                                                                  
c                                                                               
 120  if (ichk .eq. 1) then                                                     
        write(6,122)                                                            
 122    format(/3x,'mapping index ir1:')                                        
        write(6,124) (ir1(i),i=1,imax)                                          
 124    format(20i4)                                                            
        write(6,126)                                                            
 126    format(/3x,'mappint index ir2:')                                        
        write(6,124) (ir2(i),i=1,imax)                                          
      end if                                                                    
c                                                                               
      jr = 1                                                                    
      do 160 ja=1,jmax                                                          
        jr1(ja) = 0                                                             
        jr2(ja) = 0                                                             
 150    if (jr1(ja).eq.0 .and. alat(ja).le.rlat(jr)                             
     1                   .and. alat(ja).ge.rlat(jr+1)) then                     
          jr1(ja) = jr                                                          
        end if                                                                  
        if (jr2(ja).eq.0 .and. alat(ja+1).le.rlat(jr)                           
     1                   .and. alat(ja+1).ge.rlat(jr+1)) then                   
          jr2(ja) = jr                                                          
          go to 160                                                             
        else                                                                    
          if (jr .ge. jdm1) then                                                
            jr2(ja) = jr                                                        
            go to 180                                                           
          else                                                                  
            jr = jr + 1                                                         
            go to 150                                                           
          end if                                                                
        end if                                                                  
 160  continue                                                                  
c                                                                               
 180  if (ichk .eq. 1) then                                                     
        write(6,182)                                                            
 182    format(/3x,'mapping index jr1:')                                        
        write(6,124) (jr1(j),j=1,jmax)                                          
        write(6,184)                                                            
 184    format(/3x,'mapping index jr2:')                                        
        write(6,124) (jr2(j),j=1,jmax)                                          
      end if                                                                    
c                                                                               
c --- set season index, reverse nh and sh                                       
c                                                                               
      if (kssn .eq. 0) then                                                     
        iss1 = 1                                                                
        iss2 = 4                                                                
      else                                                                      
        iss1 = kssn                                                             
        iss2 = kssn                                                             
      end if                                                                    
c                                                                               
c --- loop over seasons                                                         
c                                                                               
      do 1000 issn=iss1,iss2                                                    
c                                                                               
        write(6,186) cssn(issn)                                                 
 186    format(/3x,'start to create data set : ',a3)                            
c                                                                               
        if (issn .eq. 1) then                                                   
          jssn = 3                                                              
        else if (issn .eq. 2) then                                              
          jssn = 4                                                              
        else if (issn .eq. 3) then                                              
          jssn = 1                                                              
        else if (issn .eq. 4) then                                              
          jssn = 2                                                              
        else                                                                    
          jssn = 5                                                              
        end if                                                                  
c --- reverse nh and sh beyond 30 degree h-s                                    
c     within n30 - s30 always set up as 5 (anu)                                 
c       do 190 j1=1,j30                                                         
c         j2 = jdm1p - j1                                                       
c       do 190 i=1,idm1                                                         
c         issx(i,j1) = issn                                                     
c         issx(i,j2) = jssn                                                     
c190    continue                                                                
c --- reverse nh and sh indeces                                                 
        do 190 j1=1,jdm1h                                                       
          j2 = jdm1p - j1                                                       
        do 190 i=1,idm1                                                         
          issx(i,j1) = issn                                                     
          issx(i,j2) = jssn                                                     
 190    continue                                                                
c                                                                               
        do 200 k=1,3                                                            
        do 200 j=1,jdm2                                                         
        do 200 i=1,idm2                                                         
          albv2(i,j,k) = 0.0                                                    
          albn2(i,j,k) = 0.0                                                    
          frac (i,j,k) = 0.0                                                    
 200    continue                                                                
c                                                                               
c--------------------------------------------------------------                 
c     map data to output grid                                                   
c--------------------------------------------------------------                 
c                                                                               
        do 500 ja=1,jmax                                                        
c --- fraction of coverage in j direction (cos weighted)                        
          jjr1 = jr1(ja)                                                        
          jjr2 = jr2(ja)                                                        
          if (jjr1 .eq. jjr2) then                                              
            fjr1 = 1.0                                                          
            fjr2 = 1.0                                                          
          else                                                                  
            fjr1 = (cos(alat(ja)) - cos(rlat(jjr1+1)))                          
     1           / (cos(rlat(jjr1)) - cos(rlat(jjr1+1)))                        
            fjr2 = (cos(rlat(jjr2)) - cos(alat(ja+1)))                          
     1           / (cos(rlat(jjr2)) - cos(rlat(jjr2+1)))                        
          end if                                                                
c                                                                               
          do 450 ia = 1,imax                                                    
c --- fraction of coverage in i direction                                       
            iir1 = ir1(ia)                                                      
            iir2 = ir2(ia)                                                      
            if (iir1 .eq. iir2) then                                            
              fir1 = 1.0                                                        
              fir2 = 1.0                                                        
            else                                                                
              fir1 = (rlon(iir1+1) - alon(ia))                                  
     1             / (rlon(iir1+1) - rlon(iir1))                                
              fir2 = (alon(ia+1) - rlon(iir2))                                  
     1             / (rlon(iir2+1) - rlon(iir2))                                
            end if                                                              
c                                                                               
c --- first latitude band                                                       
c                                                                               
            indx = idx(iir1,jjr1)                                               
            mssn = issx(iir1,jjr1)                                              
            isw  = iswx(indx)                                                   
            albv2(ia,ja,isw) = albv1(indx,mssn) * fir1*fjr1                     
            albn2(ia,ja,isw) = albn1(indx,mssn) * fir1*fjr1                     
            frac (ia,ja,isw) = fir1 * fjr1                                      
c                                                                               
            do 250 ir=iir1+1,iir2-1                                             
              indx = idx(ir,jjr1)                                               
              mssn = issx(ir,jjr1)                                              
              isw  = iswx(indx)                                                 
              albv2(ia,ja,isw) = albv2(ia,ja,isw)                               
     1                         + albv1(indx,mssn) * fjr1                        
              albn2(ia,ja,isw) = albn2(ia,ja,isw)                               
     1                         + albn1(indx,mssn) * fjr1                        
              frac (ia,ja,isw) = frac (ia,ja,isw) + fjr1                        
 250        continue                                                            
c                                                                               
            indx = idx(iir2,jjr1)                                               
            mssn = issx(iir2,jjr1)                                              
            isw  = iswx(indx)                                                   
            albv2(ia,ja,isw) = albv2(ia,ja,isw)                                 
     1                       + albv1(indx,mssn) * fir2*fjr1                     
            albn2(ia,ja,isw) = albn2(ia,ja,isw)                                 
     1                       + albn1(indx,mssn) * fir2*fjr1                     
            frac (ia,ja,isw) = frac (ia,ja,isw) + fir2*fjr1                     
c                                                                               
c --- consequitive latitudes                                                    
c                                                                               
            do 350 jr=jjr1+1,jjr2-1                                             
c                                                                               
              indx = idx(iir1,jr)                                               
              mssn = issx(iir1,jr)                                              
              isw  = iswx(indx)                                                 
              albv2(ia,ja,isw) = albv2(ia,ja,isw)                               
     1                         + albv1(indx,mssn) * fir1                        
              albn2(ia,ja,isw) = albn2(ia,ja,isw)                               
     1                         + albn1(indx,mssn) * fir1                        
              frac (ia,ja,isw) = frac (ia,ja,isw) + fir1                        
c                                                                               
              do 300 ir=iir1+1,iir2-1                                           
                indx = idx(ir,jr)                                               
                mssn = issx(ir,jr)                                              
                isw  = iswx(indx)                                               
                albv2(ia,ja,isw) = albv2(ia,ja,isw)                             
     1                           + albv1(indx,mssn)                             
                albn2(ia,ja,isw) = albn2(ia,ja,isw)                             
     1                           + albn1(indx,mssn)                             
                frac (ia,ja,isw) = frac (ia,ja,isw) + 1.0                       
 300          continue                                                          
c                                                                               
              indx = idx(iir2,jr)                                               
              mssn = issx(iir2,jr)                                              
              isw  = iswx(indx)                                                 
              albv2(ia,ja,isw) = albv2(ia,ja,isw)                               
     1                         + albv1(indx,mssn) * fir2                        
              albn2(ia,ja,isw) = albn2(ia,ja,isw)                               
     1                         + albn1(indx,mssn) * fir2                        
              frac (ia,ja,isw) = frac (ia,ja,isw) + fir2                        
 350        continue                                                            
c                                                                               
c --- last latitude                                                             
c                                                                               
            indx = idx(iir1,jjr2)                                               
            mssn = issx(iir1,jjr2)                                              
            isw  = iswx(indx)                                                   
            albv2(ia,ja,isw) = albv2(ia,ja,isw)                                 
     1                       + albv1(indx,mssn) * fir1*fjr2                     
            albn2(ia,ja,isw) = albn2(ia,ja,isw)                                 
     1                       + albn1(indx,mssn) * fir1*fjr2                     
            frac (ia,ja,isw) = frac (ia,ja,isw) + fir1*fjr2                     
c                                                                               
            do 400 ir=iir1+1,iir2-1                                             
              indx = idx(ir,jjr2)                                               
              mssn = issx(ir,jjr2)                                              
              isw  = iswx(indx)                                                 
              if (indx .gt. 7) isw = 2                                          
              albv2(ia,ja,isw) = albv2(ia,ja,isw)                               
     1                         + albv1(indx,mssn) * fjr2                        
              albn2(ia,ja,isw) = albn2(ia,ja,isw)                               
     1                         + albn1(indx,mssn) * fjr2                        
              frac (ia,ja,isw) = frac (ia,ja,isw) + fjr2                        
 400        continue                                                            
c                                                                               
            indx = idx(iir2,jjr2)                                               
            mssn = issx(iir2,jjr2)                                              
            isw  = iswx(indx)                                                   
            if (indx .gt. 7) isw = 2                                            
            albv2(ia,ja,isw) = albv2(ia,ja,isw)                                 
     1                       + albv1(indx,mssn) * fir2*fjr2                     
            albn2(ia,ja,isw) = albn2(ia,ja,isw)                                 
     1                       + albn1(indx,mssn) * fir2*fjr2                     
            frac (ia,ja,isw) = frac (ia,ja,isw) + fir2*fjr2                     
c                                                                               
 450      continue                                                              
c                                                                               
 500    continue                                                                
c                                                                               
c --- find mean values                                                          
c                                                                               
        do 600 ja=1,jmax                                                        
        do 600 ia=1,imax                                                        
c                                                                               
          fract = frac(ia,ja,1) + frac(ia,ja,2) + frac(ia,ja,3)                 
          if (fract .le. 0.0) then                                              
            write(6,520) ia,ja,fract                                            
 520        format(2x,'error in total fractions',                               
     1             ' ia,ja,fract =',2i4,f6.2)                                   
            call abort                                                          
          end if                                                                
c                                                                               
          do 540 ka=1,3                                                         
            if (frac(ia,ja,ka) .gt. 0.0) then                                   
              albv2(ia,ja,ka) = albv2(ia,ja,ka) / frac(ia,ja,ka)                
              albn2(ia,ja,ka) = albn2(ia,ja,ka) / frac(ia,ja,ka)                
            else                                                                
              albv2(ia,ja,ka) = 0.0                                             
              albn2(ia,ja,ka) = 0.0                                             
            end if                                                              
 540      continue                                                              
c                                                                               
          do 560 ka = 1,3                                                       
            frac(ia,ja,ka) = frac(ia,ja,ka) / fract                             
 560      continue                                                              
c                                                                               
 600    continue                                                                
c                                                                               
c       cfl = outdr(1:nn)//'.'//dssn(issn)                                      
c
        cfl(1:nn)=outdr(1:nn)
        cfl(nn+1:nn+1)='.'
        cfl(nn+1+1:nn+1+3)=dssn(issn)(1:3)
        do n1=nn+5,60
           cfl(n1:n1)=' '
        enddo
c                                                                               
        if (iout .eq. 1) then                                                   
          write(6,602) cfl                                                      
 602      format(3x,'output formatted data to file: ',a60)                      
c                                                                               
#ifdef ASSIGN
          call assign('assign -R')
#endif
          open(51,file=cfl,form='formatted',status='unknown')                   
c                                                                               
          do 620 k=1,2                                                          
            write(51,604) cttv(k),clab,imax,jmax,cssn(issn)                     
 604        format(a12,1x,a4,' idm,jdm:',2i4,1x,a3,3x,                          
     1             'note: data from n to s,  scale 1e3')                        
            do 610 j=1,jmax                                                     
            do 610 i=1,imax                                                     
              idata(i,j) = nint(1000.0 * albv2(i,j,k))                          
 610        continue                                                            
            write(51,612) ((idata(i,j),i=1,imax),j=1,jmax)                      
 612        format(20i4)                                                        
 620      continue                                                              
          do 640 k=1,2                                                          
            write(51,604) cttn(k),clab,imax,jmax,cssn(issn)                     
            do 630 j=1,jmax                                                     
            do 630 i=1,imax                                                     
              idata(i,j) = nint(1000.0 * albn2(i,j,k))                          
 630        continue                                                            
            write(51,612) ((idata(i,j),i=1,imax),j=1,jmax)                      
 640      continue                                                              
          do 660 k=1,2                                                          
            write(51,604) cttf(k),clab,imax,jmax,cssn(issn)                     
            do 650 j=1,jmax                                                     
            do 650 i=1,imax                                                     
              idata(i,j) = nint(1000.0 * frac(i,j,k))                           
 650        continue                                                            
            write(51,612) ((idata(i,j),i=1,imax),j=1,jmax)                      
 660      continue                                                              
        else                                                                    
          write(6,662) cfl                                                      
 662      format(3x,'output unformatted data to file: ',a60)                    
#ifdef ASSIGN
          call assign('assign -R')
#endif
          open(51,file=cfl,form='unformatted',status='unknown')                 
c                                                                               
          do 670 k=1,2                                                          
            write(51) ((albv2(i,j,k),i=1,imax),j=1,jmax)                        
 670      continue                                                              
          do 680 k=1,2                                                          
            write(51) ((albn2(i,j,k),i=1,imax),j=1,jmax)                        
 680      continue                                                              
          do 690 k=1,2                                                          
            write(51) ((frac(i,j,k),i=1,imax),j=1,jmax)                         
 690      continue                                                              
        end if                                                                  
c                                                                               
1000  continue                                                                  
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tage('clim_cnvalb')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop                                                                      
      end                                                                       
