#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LATG latg
#define LEVS levs
#define LONFP lonfp
#define NPES npes
      subroutine conrad(nfile,rco2)                                             
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    conrad      initializes arrays for -lw- radiation              
c   prgmmr: kenneth campana  org: w/nmc23    date: 89-07-07                     
c                                                                               
c abstract: reads co2 transmission function data(from external file),           
c   which has been pre-computed for current vertical coordinate on              
c   the front-end machine. word conversion between front-end and c205           
c   occurs here. also call tabl86 to set up tables for lw calculation           
c   this code (conrad) is only called once...                                   
c                                                                               
c program history log:                                                          
c   84-01-01  fels and schwarzkopf,gfdl.                                        
c   89-07-07  kenneth campana - removed unnecessary code and added              
c                            reading and word conversion of co2 data.           
c   89-11-29  kenneth campana - commented co2 reads because they                
c                            are not yet ready for the new gfdl lw.             
c
c   02-09-09  yifeng cui      - rsm mp    org:sdsc  
c
c                                                                               
c usage:    call conrad(nfile)                                                  
c   input argument list:                                                        
c     nfile    - integer name of external co2 file.                             
c                                                                               
c attributes:                                                                   
c   language: fortran 200.                                                      
c   machine:  cyber 205.                                                        
c                                                                               
c$$$                                                                            
c    *******************************************************************        
c    *                           c o n r a d                           *        
c    *    read co2 transmission data from unit(nfile)for new vertical  *        
c    *      coordinate tests      ...                                  *        
c    *    these arrays used to be in block data    ...k.campana-mar 90 *        
c    *******************************************************************        
c                                                                               
#include <paramodel.h>
#include <rdparm.h>
#ifdef MP
#ifdef RSM
#include <comrmpi.h>
#else
#include <commpi.h>
#endif
#endif
      dimension sgtmp(lp1,2),co21d(l,6),co22d(lp1,lp1,6)                        
      dimension co21d3(lp1,6),co21d7(lp1,6)                                     
c                                                                               
c                 co2 data tables for user''s vertical coordinate               
c                                                                               
c   the following common blocks contain pretabulated co2 transmission           
c       functions, evaluated using the methods of fels and                      
c       schwarzkopf (1981) and schwarzkopf and fels (1985),                     
c-----  the 2-dimensional arrays are                                            
c                    co2 transmission functions and their derivatives           
c        from 109-level line-by-line calculations made using the 1982           
c        mcclatchy tape (12511 lines),consolidated,interpolated                 
c        to the nmc mrf vertical coordinatte,and re-consolidated to a           
c        200 cm-1 bandwidth. the interpolation method is described in           
c        schwarzkopf and fels (j.g.r.,1985).                                    
c-----  the 1-dim arrays are                                                    
c                  co2 transmission functions and their derivatives             
c          for tau(i,i+1),i=1,l,                                                
c            where the values are not obtained by quadrature,but are the        
c            actual transmissivities,etc,between a pair of pressures.           
c          these used only for nearby layer calculations including qh2o.        
c-----  the weighting function gtemp=p(k)**0.2*(1.+p(k)/30000.)**0.8/           
c         1013250.,where p(k)=pressure,nmc mrf(new)  l18 data levels for        
c         pstar=1013250.                                                        
c-----  stemp is us standard atmospheres,1976,at data pressure levels           
c        using nmc mrf sigmas,where pstar=1013.25 mb (ptz program)              
#include <co2dta.h>
c====>   begin here to get constants for radiation package                      
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef RSM
      open(unit=nfile,file='r_co2con',form='unformatted',status='old')
#else
      open(unit=nfile,file='co2con',form='unformatted',status='old')
#endif
      rewind nfile                                                              
c       read in pre-computed co2 transmission data....                          
      do 10 kk=1,2                                                              
        read(nfile) (sgtmp(i,kk),i=1,lp1)                                       
   10 continue                                                                  
      do 15 kk=1,6                                                              
        read(nfile) (co21d(i,kk),i=1,l)                                         
   15 continue                                                                  
      do 20 kk=1,6                                                              
        read(nfile) ((co22d(i,j,kk),i=1,lp1),j=1,lp1)                           
   20 continue                                                                  
      do 25 kk=1,6                                                              
        read(nfile) (co21d3(i,kk),i=1,lp1)                                      
   25 continue                                                                  
      do 30 kk=1,6                                                              
        read(nfile) (co21d7(i,kk),i=1,lp1)                                      
   30 continue                                                                  
c  read co2 concentration in ppm (defaulted in gradfs if missing)               
      read(nfile,end=31) rco2                                                   
   31 continue                                                                  
#ifndef NOPRINT
      print *,'co2 concentration is ',rco2                                      
#endif
      rewind nfile                                                              
#ifdef MP
      endif
#ifdef RSM
      call rmpbcastr(sgtmp,lp1*2)
      call rmpbcastr(co21d,l*6)
      call rmpbcastr(co22d,lp1*lp1*6)
      call rmpbcastr(co21d3,lp1*6)
      call rmpbcastr(co21d7,lp1*6)
      call rmpbcastr(rco2,1)
#else
      call mpbcastr(sgtmp,lp1*2)
      call mpbcastr(co21d,l*6)
      call mpbcastr(co22d,lp1*lp1*6)
      call mpbcastr(co21d3,lp1*6)
      call mpbcastr(co21d7,lp1*6)
      call mpbcastr(rco2,1)
#endif
#endif
      do 35 k=1,lp1                                                             
        stemp(k) = sgtmp(k,1)                                                   
        gtemp(k) = sgtmp(k,2)                                                   
   35 continue                                                                  
      do 40 k=1,l                                                               
        cdtm51(k) = co21d(k,1)                                                  
        co2m51(k) = co21d(k,2)                                                  
        c2dm51(k) = co21d(k,3)                                                  
        cdtm58(k) = co21d(k,4)                                                  
        co2m58(k) = co21d(k,5)                                                  
        c2dm58(k) = co21d(k,6)                                                  
   40 continue                                                                  
      do 45 j=1,lp1                                                             
       do 45 i=1,lp1                                                            
        cdt51(i,j) = co22d(i,j,1)                                               
        co251(i,j) = co22d(i,j,2)                                               
        c2d51(i,j) = co22d(i,j,3)                                               
        cdt58(i,j) = co22d(i,j,4)                                               
        co258(i,j) = co22d(i,j,5)                                               
        c2d58(i,j) = co22d(i,j,6)                                               
   45 continue                                                                  
      do 50 k=1,lp1                                                             
        cdt31(k) = co21d3(k,1)                                                  
        co231(k) = co21d3(k,2)                                                  
        c2d31(k) = co21d3(k,3)                                                  
        cdt38(k) = co21d3(k,4)                                                  
        co238(k) = co21d3(k,5)                                                  
        c2d38(k) = co21d3(k,6)                                                  
   50 continue                                                                  
      do 55 k=1,lp1                                                             
        cdt71(k) = co21d7(k,1)                                                  
        co271(k) = co21d7(k,2)                                                  
        c2d71(k) = co21d7(k,3)                                                  
        cdt78(k) = co21d7(k,4)                                                  
        co278(k) = co21d7(k,5)                                                  
        c2d78(k) = co21d7(k,6)                                                  
   55 continue                                                                  
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef DBG
      print 66,nfile                                                            
   66 format(1h ,'----read co2 transmission functions from unit ',i2)           
#endif
#ifdef MP
      endif
#endif
c                                                                               
c......    define tables for lw radiation                                       
      call table                                                                
c                                                                               
      return                                                                    
      end                                                                       
