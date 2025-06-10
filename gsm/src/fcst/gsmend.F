#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAPP lcapp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define MLVARK mlvark
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTRAC ntrac
#define SLVARK slvark
#define TWOJ1 twoj1


        subroutine gsmend
c$$$  subroutine program documentation block                                          
c                                                                               
c subroutine: gsmend         finish global forecast 
c   prgmmr: juang             org: w/np5       date: 00-06-01                    
c                                                                               
c abstract: finish global forecast with spectral model.                           
c                                                                               
c program history log:                                                          
c   81-01-01  sela                                                              
c   99-05-21  hann-ming henry juang   add mpi in one loop version
c   99-06-22  hann-ming henry juang   add mpi in all loop version
c   00-06-01  hann-ming henry juang   change gsm into gsmini, gsmstep and gsmend
c                                                                               
c input files:                                                                  
c   unit   11    sigma file (analysis or at time t-dt)                          
c   unit   12    sigma file (at time t if not analysis)                         
c   unit   14    surface file                                                   
c   unit   15    co2 constants (dependent on vertical resolution)               
c   unit   24    mountain variance (dependent on horizontal resolution)         
c   unit   43    cloud tuning                                                   
c                                                                               
c output files:                                                                 
c   unit   51    sigma file (at time t-dt)                                      
c   unit   52    sigma file (at time t)                                         
c   unit   53    surface file                                                   
c   unit   61    initial zonal diagnostics                                      
c   unit   63    flux diagnostics                                               
c   unit   64    final zonal diagnostics                                        
c   unit   67    grid point diagnostics                                         
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
c....                                                                           
#include <paramodel.h>
#include <comfibm.h>
#include <comcon.h>
#include <comgpd.h>
#include <comfgsm.h>
#ifdef LFM
#include <comlfm.h>
#endif
c
#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
#endif
c
#ifdef MP
#include <commpi.h>
#endif
c
c write sigitdt and sigit at the end of inchour
c
#ifdef SCM
      call wrisig(n1,thour,idate,qm,tem,dim,zem,rm ,sl,si,gz,z00,vvel,1)
      call wrisig(n1,thour,idate,q ,te ,di ,ze ,rq ,sl,si,gz,z00,vvel,2)
#else
#ifndef HYBRID
      call wrisig(n1,thour,idate,qm,tem,dim,zem,rm ,sl,si,gz,z00,1)
      call wrisig(n1,thour,idate,q ,te ,di ,ze ,rq ,sl,si,gz,z00,2)
#else
      call wrisig(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,1)
      call wrisig(n1,thour,idate,q ,te ,di ,ze ,rq,ak5,bk5,gz,z00,2)
#endif
#endif
c
c
      return
      end
