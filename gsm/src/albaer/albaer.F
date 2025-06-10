#define LATG latg
#define LONF lonf
      program albaer                                                            
c$$$  main program documentation block                                          
c                                                                               
c main program:  albaer      create one dataset from albedo and aerosol
c   prgmmr: cpc          org: w/np51     date: 01-03-31                    
c                                                                               
c abstract: the programm read in separate albedo data and aerosol data                
c       output from programs 'cnvalb.f' and 'cnvaer.f' and combines             
c       them into one data set for forcast model use.                           
c                                                                               
c namelists:                                                                    
c       imax,jmax   - i,j dimension of input data                               
c       infl        - =1 input data is formatted integer, scaled by 1e3         
c                     =0 input data is unformatted in fraction                  
c       ialb,iaer   - num of input albedo/aerosols data sets                    
c                       monthly=12, seasonal=4, annual mean=1                   
c
c input files:                                                                  
c       six surface 2-d data sets           from unit 11,12,...                 
c       five aerosol 2-d data sets          from unit 21,22,...                 
c                                                                               
c output files:                                                                 
c       eleven 2-d data set                 to unit 51                          
c                                                                               
c subprograms called:                                                           
c   abort       system abort 
c                                                                               
c attributes:                                                                   
c   language: fortran                                                      
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter (idim=lonf_, jdim=latg_, kd1=6, kd2=5)                          
      dimension idata(idim,jdim), fdata(idim,jdim)                              
      character title*80                                                        
      data imax/lonf_/,jmax/latg_/                                              
      data infl/1/,ialb/4/,iaer/1/                                              
      namelist /namein/ imax,jmax,infl,ialb,iaer                                
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tagb('clim_albaer',2001,0000,0000,'np51   ')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      read(5,namein)                                                            
      write(6,namein)                                                           
c                                                                               
      if (imax.gt.idim .or. jmax.gt.jdim) then                                  
        write(6,12) idim,jdim                                                   
  12    format(3x,'error in dimension specifications:'/                         
     1         3x,'maximum limits: idim,jdim =',2i4)                            
        call abort                                                              
      end if                                                                    
      if (infl .eq. 1) then                                                     
c --- loop through each of the six records in a file                            
        do 40 k=1,kd1                                                           
c --- loop through input files                                                  
        do 40 m=1,ialb                                                          
          nf = 10 + m                                                           
          read(nf,14) title                                                     
  14      format(a80)                                                           
          write(6,16) title                                                     
  16      format(3x,'input data file:'/a80)                                     
          write(6,22) k                                                         
  22      format(3x,'read input file group-1 k =',i2)                           
          read(nf,24) ((idata(i,j),i=1,imax),j=1,jmax)                          
  24      format(20i4)                                                          
          do 30 j=1,jmax                                                        
          do 30 i=1,imax                                                        
            fdata(i,j) = 0.001 * idata(i,j)                                     
            if (fdata(i,j).lt.0.0 .or. fdata(i,j).gt.1.0) then                  
              write(6,26) i,j,fdata(i,j)                                        
  26          format(3x,'error in data: i,j =',2i4,                             
     1               '  fdata(i,j) =',f8.4)                                     
              call abort                                                        
            end if                                                              
  30      continue                                                              
          if (k.le.4 .or. m.eq.ialb) then                                       
            write(51) ((fdata(i,j),i=1,imax),j=1,jmax)                          
            write(6,32) k, m                                                    
  32        format(3x,'write record num =',i3,' for inpur file',                
     1             i3)                                                          
          end if                                                                
  40    continue                                                                
c --- loop through each of the five records in a file                           
        do 60 k=1,kd2                                                           
c --- loop through input files                                                  
        do 60 m=1,iaer                                                          
          nf = 20 + m                                                           
          read(nf,14) title                                                     
          write(6,16) title                                                     
          write(6,42) k                                                         
  42      format(3x,'read input file group-2 k=',i2)                            
          read(nf,24) ((idata(i,j),i=1,imax),j=1,jmax)                          
          do 50 j=1,jmax                                                        
          do 50 i=1,imax                                                        
            fdata(i,j) = 0.001 * idata(i,j)                                     
            if (fdata(i,j).lt.0.0 .or. fdata(i,j).gt.1.0) then                  
              write(6,26) i,j,fdata(i,j)                                        
              call abort                                                        
            end if                                                              
  50      continue                                                              
          write(51) ((fdata(i,j),i=1,imax),j=1,jmax)                            
          write(6,32) k, m                                                      
  60    continue                                                                
      else                                                                      
        do 80 k=1,kd1                                                           
        do 80 m=1,ialb                                                          
          nf = 10 + m                                                           
          write(6,22) k                                                         
          read(nf) ((fdata(i,j),i=1,imax),j=1,jmax)                             
          do 70 j=1,jmax                                                        
          do 70 i=1,imax                                                        
            if (fdata(i,j).lt.0.0 .or. fdata(i,j).gt.1.0) then                  
              write(6,26) i,j,fdata(i,j)                                        
              call abort                                                        
            end if                                                              
  70      continue                                                              
          if (k.le.4 .or. m.eq.ialb) then                                       
            write(51) ((fdata(i,j),i=1,imax),j=1,jmax)                          
            write(6,32) k, m                                                    
          end if                                                                
  80    continue                                                                
        do 100 k=1,kd2                                                          
        do 100 m=1,iaer                                                         
          nf = 20 + m                                                           
          write(6,42) k                                                         
          read(nf) ((fdata(i,j),i=1,imax),j=1,jmax)                             
          do 90 j=1,jmax                                                        
          do 90 i=1,imax                                                        
            if (fdata(i,j).lt.0.0 .or. fdata(i,j).gt.1.0) then                  
              write(6,26) i,j,fdata(i,j)                                        
              call abort                                                        
            end if                                                              
  90    continue                                                                
          write(51) ((fdata(i,j),i=1,imax),j=1,jmax)                            
          write(6,32) k, m                                                      
 100    continue                                                                
      end if                                                                    
c                                                                               
#ifdef NCO_TAG
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call w3tage('clim_albaer')
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif
      stop                                                                      
      end                                                                       
