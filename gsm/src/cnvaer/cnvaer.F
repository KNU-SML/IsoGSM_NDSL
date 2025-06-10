#define LATG latg
#define LONF lonf
      program cnvaer                                                            
c$$$  main program documentation block
c
c main program:  cnvaer      reads aerosol and maps to gaussin grid
c   prgmmr: cpc          org: w/np51     date: 01-03-31
c
c abstract:
c***************************************************************                
c   this program reads aerosol type index data file in one                      
c     degree lat/lon resolution, then converts to aerosol                       
c     profile distribution, and maps to gaussin grid for                        
c     forecast model use.                                                       
c                                                                               
c   input data:                                                                 
c     aerosol type index file (20a4)               - unit 11                    
c       data in ncmc form, where n,m=1-5 for the major and the                  
c       secondary profile:                                                      
c     1.ubr; 2.cont-1; 3.mar-1; 4.cont-2; 5.mar-2                               
c       and c is char expression of percentage in 5 increment                   
c     'a'=0,'b'=5,... 'u'=100.                                                  
c                                                                               
c   output data:                                                                
c     aerosol profile map in feaction              - unit 51                    
c                                                                               
c   input namelist parameters:                                                  
c     imax,jmax      - output data i,j (lon/lat) dimensions                     
c     igau           - =1 output on gaussing grid                               
c                      =0 output on regular lat/lon grid                        
c     issn           - season index 1-4 (w,s,s,f); =5 annual mean.              
c     ichk           - =1 print out check point values                          
c                      =0 do not print                                          
c     iout           - =1 output formatted data, scale 1e3                      
c                      =0 output unformatted data in fraction                   
c****************************************************************               
c                                                                               
c attributes:
c   language: fortran
c
c$$$
c
c----1 ( idm1=360, jdm1=180, kdm=5,   idm2=384, jdm2=196,                      
#include <paramodel.h>
                   p a r a m e t e r                                            
     1 ( idm1=360, jdm1=180, kdm=5,   idm2=lonf_, jdm2=latg_,                   
     2   idm1p=idm1+1, jdm1p=jdm1+1,  nmax=13,                                  
     3   idm2p=idm2+1, jdm2p=jdm2+1,  jdm2h=jdm2/2 )                            
                   d i m e n s i o n                                            
     1  faer1(idm1,jdm1,kdm), faer2(idm2,jdm2,kdm), iaer2(idm2,jdm2)            
     2, fsum2(idm2,jdm2,kdm)                                                    
                   d i m e n s i o n                                            
     1  rlon(idm1p), rlat(jdm1p), alon(idm2p), alat(jdm2p)                      
     2, blat(jdm2), ir1(idm2), ir2(idm2), jr1(jdm2), jr2(jdm2)                  
                   c h a r a c t e r                                            
     1  cdat(idm1,jdm1)*4,  ch1*4,  title(2)*40,  clab*4                        
     2, cssn(5)*3,  ctyp(kdm)*6                                                 
c
      real*4 frc4
c                                                                               
      data imax/lonf_/,jmax/latg_/                                              
      data igau/1/,issn/5/,clab/'    '/,ichk/1/,iout/1/                         
c                                                                               
                   namelist /namein/                                            
     1  imax, jmax, igau, issn, clab, ichk, iout                                
c                                                                               
      data cssn/'win','spr','sum','fal','anu' /                                 
      data ctyp/'urban ','cont-1','mar-1 ','cont-2','mar-2 '/                   
      data m1typ / 3 /                                                          
c
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_cnvaer',2001,0000,0000,'np51)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
c                                                                               
c --- read in namelist                                                          
c                                                                               
      read (5,namein)                                                           
      write(6,namein)                                                           
      if (imax.gt.idm2 .or. jmax.gt.jdm2) then                                  
        write(6,2) imax,jmax,idm2,jdm2                                          
   2    format(2x,'defined output array size exceeds limits',                   
     1        /4x,'imax,jmax=',2i4,'  idm2,jdm2=',2i4)                          
        call abort                                                              
      end if                                                                    
      if (issn.lt.1 .or. issn.gt.5) then                                        
        write(6,4) issn                                                         
   4    format(2x,'season index out off range, issn =',i4)                      
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
      do 10 k=1,kdm                                                             
      do 10 j=1,jdm1                                                            
      do 10 i=1,idm1                                                            
        faer1(i,j,k) = 0.0                                                      
  10  continue                                                                  
c                                                                               
c --- read in aerosol type index file, convert char to integer                  
c                                                                               
      read(11,12) title                                                         
      write(6,12) title                                                         
  12  format(2a40)                                                              
      read(11,14) ((cdat(i,j),i=1,idm1),j=1,jdm1)                               
  14  format(20a4)                                                              
      do 20 j=1,jdm1                                                            
      do 20 i=1,idm1                                                            
        ch1 = cdat(i,j)                                                         
        if (ch1(1:1).lt.'1' .or. ch1(1:1).gt.'5' .or.                           
     1      ch1(3:3).lt.'1' .or. ch1(3:3).gt.'5') then                          
          write(6,16) i,j,cdat(i,j)                                             
  16      format(3x,'error in input index file: i,j =',2i4,                     
     1           3x,'cdat(i,j) = ',a4)                                          
          call abort                                                            
        else if (ch1(2:2).lt.'A' .or. ch1(2:2).gt.'U' .or.                      
     1           ch1(4:4).lt.'A' .or. ch1(4:4).gt.'U') then                     
          write(6,16) i,j, cdat(i,j)                                            
          call abort                                                            
        else                                                                    
          ity1 = ichar(ch1(1:1)) - ichar('0')                                   
          ity2 = ichar(ch1(3:3)) - ichar('0')                                   
          npc1 = 5 * (ichar(ch1(2:2)) - ichar('A'))                             
          npc2 = 5 * (ichar(ch1(4:4)) - ichar('A'))                             
          faer1(i,j,ity1) = 0.01 * npc1                                         
          faer1(i,j,ity2) = 0.01 * npc2                                         
        end if                                                                  
  20  continue                                                                  
c                                                                               
c --- check input data for errors                                               
c                                                                               
      do 40 j=1,jdm1                                                            
      do 40 i=1,idm1                                                            
        frc = 0.0                                                               
        do 30 k=1,kdm                                                           
          frc = frc + faer1(i,j,k)                                              
  30    continue                                                                
        frc4=frc
        frc=frc4
        if (frc .gt. 1.0) then                                                  
c --- total fraction can not exceed 1.0                                         
          write(6,32) i,j,(faer1(i,j,k),k=1,kdm)                                
  32      format(2x,'error in input data, i,j=',2i4,                            
     1           '  faer1(k),k=1,5 :'/5f8.3)                                    
          call abort                                                            
        else if (frc .lt. 1.0) then                                             
c --- fill mar-1 type to make total frc = 1.0                                   
          faer1(i,j,m1typ) = faer1(i,j,m1typ) + 1.0 - frc                       
        end if                                                                  
  40  continue                                                                  
c                                                                               
c --- set up input data lon, lat of each grid box edges in radians              
c     matthews data are at the center of grid boxes, so the edge               
c     of the first box is at the 0 deg                                          
      drlon = tpi / float(idm1)                                                 
      drlat = pi / float(jdm1)                                                  
c                                                                               
      do 45 i=1,idm1p                                                           
        rlon(i) = drlon * (i - 1)                                               
  45  continue                                                                  
c --- the edge of the first lat box is at 90 deg                                
      do 50 j=1,jdm1p                                                           
        rlat(j) = hpi - drlat * (j - 1)                                         
  50  continue                                                                  
c --- check print                                                               
      if (ichk .eq. 1) then                                                     
        write(6,52)                                                             
  52    format(/3x,'input lon values (box edges):')                             
        write(6,54) (raddeg*rlon(i),i=1,idm1p)                                  
  54    format(16f5.1)                                                          
        write(6,56)                                                             
  56    format(/3x,'input lat values (box edges):')                             
        write(6,54) (raddeg*rlat(j),j=1,jdm1p)                                  
      end if                                                                    
c --- rem: model lon starts at greenwich, so the edge of                        
c     the first lon box is at the -0.5 of the box width                         
      dalon = tpi / float(imax)                                                 
      do 60 i=1,imaxp                                                           
        alon(i) = dalon * (i - 1.5)                                             
  60  continue                                                                  
c --- blat contains model lat grid                                              
      if (igau .eq. 1) then                                                     
        call glats(jmax,blat)                                                   
        alat(1) = hpi                                                           
        do 65 j=1,jmax-1                                                        
          alat(j+1) = 0.5 * (blat(j) + blat(j+1))                               
  65    continue                                                                
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
        write(6,54) (raddeg*alon(i),i=1,imaxp)                                  
        write(6,74)                                                             
  74    format(/3x,'output lat values (box centers):')                          
        write(6,54) (raddeg*blat(j),j=1,jmax)                                   
        write(6,76)                                                             
  76    format(/3x,'output lat values (box edges):')                            
        write(6,54) (raddeg*alat(j),j=1,jmaxp)                                  
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
      do 200 k=1,kdm                                                            
      do 200 j=1,jdm2                                                           
      do 200 i=1,idm2                                                           
        faer2(i,j,k) = 0.0                                                      
 200  continue                                                                  
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
     1         / (cos(rlat(jjr1)) - cos(rlat(jjr1+1)))                          
          fjr2 = (cos(rlat(jjr2)) - cos(alat(ja+1)))                            
     1         / (cos(rlat(jjr2)) - cos(rlat(jjr2+1)))                          
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
     1           / (rlon(iir1+1) - rlon(iir1))                                  
            fir2 = (alon(ia+1) - rlon(iir2))                                    
     1           / (rlon(iir2+1) - rlon(iir2))                                  
          end if                                                                
c                                                                               
          do 430 k=1,kdm                                                        
c                                                                               
c --- first latitude band                                                       
c                                                                               
            faer2(ia,ja,k) = faer1(iir1,jjr1,k) * fir1*fjr1                     
            fsum2(ia,ja,k) = fir1*fjr1                                          
            do 250 ir=iir1+1,iir2-1                                             
              faer2(ia,ja,k) = faer2(ia,ja,k)                                   
     1                       + faer1(ir,jjr1,k) * fjr1                          
              fsum2(ia,ja,k) = fsum2(ia,ja,k) + fjr1                            
 250        continue                                                            
            faer2(ia,ja,k) = faer2(ia,ja,k)                                     
     1                     + faer1(iir2,jjr1,k) * fir2*fjr1                     
            fsum2(ia,ja,k) = fsum2(ia,ja,k) + fir2*fjr1                         
c                                                                               
c --- consequitive latitudes                                                    
c                                                                               
            do 350 jr=jjr1+1,jjr2-1                                             
c                                                                               
              faer2(ia,ja,k) = faer2(ia,ja,k)                                   
     1                       + faer1(iir1,jr,k) * fir1                          
              fsum2(ia,ja,k) = fsum2(ia,ja,k) + fir1                            
              do 300 ir=iir1+1,iir2-1                                           
                faer2(ia,ja,k) = faer2(ia,ja,k)                                 
     1                         + faer1(ir,jr,k)                                 
                fsum2(ia,ja,k) = fsum2(ia,ja,k) + 1.0                           
 300          continue                                                          
              faer2(ia,ja,k) = faer2(ia,ja,k)                                   
     1                       + faer1(iir2,jr,k) * fir2                          
              fsum2(ia,ja,k) = fsum2(ia,ja,k) + fir2                            
 350        continue                                                            
c                                                                               
c --- last latitude                                                             
c                                                                               
            faer2(ia,ja,k) = faer2(ia,ja,k)                                     
     1                     + faer1(iir1,jjr2,k) * fir1*fjr2                     
            fsum2(ia,ja,k) = fsum2(ia,ja,k) + fir1*fjr2                         
            do 400 ir=iir1+1,iir2-1                                             
              faer2(ia,ja,k) = faer2(ia,ja,k)                                   
     1                       + faer1(ir,jjr2,k) * fjr2                          
              fsum2(ia,ja,k) = fsum2(ia,ja,k) + fjr2                            
 400        continue                                                            
            faer2(ia,ja,k) = faer2(ia,ja,k)                                     
     1                     + faer1(iir2,jjr2,k) * fir2*fjr2                     
            fsum2(ia,ja,k) = fsum2(ia,ja,k) + fir2*fjr2                         
c                                                                               
 430      continue                                                              
c                                                                               
 450    continue                                                                
c                                                                               
 500  continue                                                                  
c                                                                               
c --- compute average                                                           
c                                                                               
      do 510 k=1,kdm                                                            
      do 510 j=1,jmax                                                           
      do 510 i=1,imax                                                           
        if (fsum2(i,j,k) .gt. 0.0) then                                         
          faer2(i,j,k) = faer2(i,j,k) / fsum2(i,j,k)                            
        else                                                                    
          write(6,502) i,j,k,fsum2(i,j,k)                                       
 502      format(2x,'error in fsum, i,j,k=',3i4,' fsum2 =',f10.3)               
          call abort                                                            
        end if                                                                  
 510  continue                                                                  
c                                                                               
c --- final check                                                               
c                                                                               
      do 550 j=1,jmax                                                           
      do 550 i=1,imax                                                           
c                                                                               
        do 520 k=1,kdm                                                          
          if (faer2(i,j,k) .lt. 0.01) faer2(i,j,k) = 0.0                        
          faer2(i,j,k) = min(1.0, faer2(i,j,k))                                 
 520    continue                                                                
c                                                                               
        fsea = faer2(i,j,3) + faer2(i,j,5)                                      
        flnd = faer2(i,j,1) + faer2(i,j,2) + faer2(i,j,4)                       
        if (flnd .gt. fsea) then                                                
          fac = 1.0 - fsea - flnd + faer2(i,j,2)                                
          if (fac .lt. 0.0) then                                                
            write(6,532) i,j,(faer2(i,j,k),k=1,kdm)                             
 532        format(2x,'error in output array, i,j=',2i3,                        
     1             '  faer2(k),k=1,5 :'/5f8.3)                                  
            call abort                                                          
          end if                                                                
          faer2(i,j,2) = fac                                                    
        else                                                                    
          fac = 1.0 - flnd - fsea + faer2(i,j,3)                                
          if (fac .lt. 0.0) then                                                
            write(6,532) i,j,(faer2(i,j,k),k=1,kdm)                             
            call abort                                                          
          end if                                                                
          faer2(i,j,3) = fac                                                    
        end if                                                                  
c                                                                               
 550  continue                                                                  
c --- write to output file                                                      
c                                                                               
      if (iout .eq. 1) then                                                     
c --- output formatted file, scaled to 1e3                                      
        do 650 k=1,kdm                                                          
          do 600 j=1,jmax                                                       
          do 600 i=1,imax                                                       
            iaer2(i,j) = nint(1000.0 * faer2(i,j,k))                            
 600      continue                                                              
          write(51,622) ctyp(k),clab,imax,jmax,cssn(issn)                       
 622      format('aerosol ',a6,1x,a4,' idm,jdm:',2i4,1x,a3,2x,                  
     1           'note: data from n to s,  scale 1e3')                          
          write(51,624) ((iaer2(i,j),i=1,imax),j=1,jmax)                        
 624      format(20i4)                                                          
 650    continue                                                                
      else                                                                      
c --- output unformatted data file                                              
        do 700 k=1,kdm                                                          
          write(51) ((faer2(i,j,k),i=1,imax),j=1,jmax)                          
 700    continue                                                                
      end if                                                                    
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tage('clim_cnvaer')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif

      stop                                                                      
      end                                                                       
