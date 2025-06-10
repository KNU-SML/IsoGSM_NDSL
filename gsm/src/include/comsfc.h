c
c  common for variables in the surface file
c
#ifndef RSM
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#else
#ifdef MP
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif
c
c  important:  do not change the order of the variables
c              do not insert any variable, simply add them
c              at the end
c   in addition this file must be consistent with comsfcmrg.h
c
#ifdef OSULSM1
      real tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,albedo,slmsk
      real plantr,canopy,f10m
#endif
c
#ifdef OSULSM2
      real tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,albedo,slmsk
      real vfrac,canopy,f10m,vtype,stype,facalf,uustar,ffmm,ffhh
#endif
c
#ifdef NOALSM1
      real tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,albedo,slmsk
      real vfrac,canopy,f10m,vtype,stype,facalf,uustar,ffmm,ffhh
      real prcp,srflag,snwdph,slc,shdmin,shdmax,slope,snoalb
#endif
c
#ifdef VICLSM1
      real tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,albedo,slmsk
      real vfrac,canopy,f10m,vtype,vroot,facalf,uustar,ffmm,ffhh
      real prcp,srflag
      real binf,ds,dsm,ws,cef,expt,kest,dph,bub,qrt,bkd,sld,wcr,wpw,smr
      real smx,dphn,smxn,expn,bubn,alpn,betn,gamn,flai,silz,snwz,sic
      real csno,rsno,tsf,tpk,sfw,pkw,lstsn
#endif
c
#ifdef VICLSM2
      real tsea,smc,sheleg,stc,tg3,zorl,cv,cvb,cvt,albedo,slmsk
      real vfrac,canopy,f10m,vtype,vroot,facalf,uustar,ffmm,ffhh
      real prcp,srflag
      real binf,ds,dsm,ws,cef,expt,kest,dph,bub,qrt,bkd,sld,wcr,wpw,smr
      real smx,dphn,smxn,expn,bubn,alpn,betn,gamn,silz,snwz
      real nveg,mvfr,mcnp,mvty,flai,msno,msmc,msic,mstc
      real csno,rsno,tsf,tpk,sfw,pkw,lstsn
#endif
c
#ifdef OSULSM1
      common /comsfc/
     &          tsea  (LONF2S,LATG2S),
     &          smc   (LONF2S,LATG2S,lsoil_),
     &          sheleg(LONF2S,LATG2S),
     &          stc   (LONF2S,LATG2S,lsoil_),
     &          tg3   (LONF2S,LATG2S),
     &          zorl  (LONF2S,LATG2S),
     &          cv    (LONF2S,LATG2S),
     &          cvb   (LONF2S,LATG2S),
     &          cvt   (LONF2S,LATG2S),
     &          albedo(LONF2S,LATG2S,lalbd_),
     &          slmsk (LONF2S,LATG2S),
     &          plantr(LONF2S,LATG2S),
     &          canopy(LONF2S,LATG2S),
     &          f10m  (LONF2S,LATG2S) 
#endif
c
#ifdef OSULSM2
      common /comsfc/
     &          tsea  (LONF2S,LATG2S),
     &          smc   (LONF2S,LATG2S,lsoil_),
     &          sheleg(LONF2S,LATG2S),
     &          stc   (LONF2S,LATG2S,lsoil_),
     &          tg3   (LONF2S,LATG2S),
     &          zorl  (LONF2S,LATG2S),
     &          cv    (LONF2S,LATG2S),
     &          cvb   (LONF2S,LATG2S),
     &          cvt   (LONF2S,LATG2S),
     &          albedo(LONF2S,LATG2S,lalbd_),
     &          slmsk (LONF2S,LATG2S),
     &          vfrac (LONF2S,LATG2S),
     &          canopy(LONF2S,LATG2S),
     &          f10m  (LONF2S,LATG2S),
     &          vtype (LONF2S,LATG2S),
     &          stype (LONF2S,LATG2S),
     &          facalf(LONF2S,LATG2S,2),
     &          uustar(LONF2S,LATG2S),
     &          ffmm  (LONF2S,LATG2S),
     &          ffhh  (LONF2S,LATG2S) 
#endif
c
#ifdef NOALSM1
      common /comsfc/
     &          tsea  (LONF2S,LATG2S),
     &          smc   (LONF2S,LATG2S,lsoil_),
     &          sheleg(LONF2S,LATG2S),
     &          stc   (LONF2S,LATG2S,lsoil_),
     &          tg3   (LONF2S,LATG2S),
     &          zorl  (LONF2S,LATG2S),
     &          cv    (LONF2S,LATG2S),
     &          cvb   (LONF2S,LATG2S),
     &          cvt   (LONF2S,LATG2S),
     &          albedo(LONF2S,LATG2S,lalbd_),
     &          slmsk (LONF2S,LATG2S),
     &          vfrac (LONF2S,LATG2S),
     &          canopy(LONF2S,LATG2S),
     &          f10m  (LONF2S,LATG2S), 
     &          vtype (LONF2S,LATG2S),
     &          stype (LONF2S,LATG2S),
     &          facalf(LONF2S,LATG2S,2),
     &          uustar(LONF2S,LATG2S),
     &          ffmm  (LONF2S,LATG2S),
     &          ffhh  (LONF2S,LATG2S),
     &          prcp  (LONF2S,LATG2S),
     &          srflag(LONF2S,LATG2S),
     &          snwdph(LONF2S,LATG2S),
     &          slc   (LONF2S,LATG2S,lsoil_),
     &          shdmin(LONF2S,LATG2S),
     &          shdmax(LONF2S,LATG2S),
     &          slope (LONF2S,LATG2S),
     &          snoalb(LONF2S,LATG2S) 
#endif
c
#ifdef VICLSM1
      common /comsfc/
     &          tsea  (LONF2S,LATG2S),
     &          smc   (LONF2S,LATG2S,lsoil_),
     &          sheleg(LONF2S,LATG2S),
     &          stc   (LONF2S,LATG2S,nsoil_),
     &          tg3   (LONF2S,LATG2S),
     &          zorl  (LONF2S,LATG2S),
     &          cv    (LONF2S,LATG2S),
     &          cvb   (LONF2S,LATG2S),
     &          cvt   (LONF2S,LATG2S),
     &          albedo(LONF2S,LATG2S,lalbd_),
     &          slmsk (LONF2S,LATG2S),
     &          vfrac (LONF2S,LATG2S),
     &          canopy(LONF2S,LATG2S),
     &          f10m  (LONF2S,LATG2S),
     &          vtype (LONF2S,LATG2S),
     &          vroot (LONF2S,LATG2S,lsoil_),
     &          facalf(LONF2S,LATG2S,2),
     &          uustar(LONF2S,LATG2S),
     &          ffmm  (LONF2S,LATG2S),
     &          ffhh  (LONF2S,LATG2S),
     &          prcp  (LONF2S,LATG2S),
     &          srflag(LONF2S,LATG2S),
     &          binf  (LONF2S,LATG2S),
     &          ds    (LONF2S,LATG2S),
     &          dsm   (LONF2S,LATG2S),
     &          ws    (LONF2S,LATG2S),
     &          cef   (LONF2S,LATG2S),
     &          expt  (LONF2S,LATG2S,lsoil_),
     &          kest  (LONF2S,LATG2S,lsoil_),
     &          dph   (LONF2S,LATG2S,lsoil_),
     &          bub   (LONF2S,LATG2S,lsoil_),
     &          qrt   (LONF2S,LATG2S,lsoil_),
     &          bkd   (LONF2S,LATG2S,lsoil_),
     &          sld   (LONF2S,LATG2S,lsoil_),
     &          wcr   (LONF2S,LATG2S,lsoil_),
     &          wpw   (LONF2S,LATG2S,lsoil_),
     &          smr   (LONF2S,LATG2S,lsoil_),
     &          smx   (LONF2S,LATG2S,lsoil_),
     &          dphn  (LONF2S,LATG2S,nsoil_),
     &          smxn  (LONF2S,LATG2S,nsoil_),
     &          expn  (LONF2S,LATG2S,nsoil_),
     &          bubn  (LONF2S,LATG2S,nsoil_),
     &          alpn  (LONF2S,LATG2S,nsoil_),
     &          betn  (LONF2S,LATG2S,nsoil_),
     &          gamn  (LONF2S,LATG2S,nsoil_),
     &          flai  (LONF2S,LATG2S),
     &          silz  (LONF2S,LATG2S),
     &          snwz  (LONF2S,LATG2S),
     &          sic   (LONF2S,LATG2S,lsoil_),
     &          csno  (LONF2S,LATG2S),
     &          rsno  (LONF2S,LATG2S),
     &          tsf   (LONF2S,LATG2S),
     &          tpk   (LONF2S,LATG2S),
     &          sfw   (LONF2S,LATG2S),
     &          pkw   (LONF2S,LATG2S),
     &          lstsn (LONF2S,LATG2S)
#endif
c
#ifdef VICLSM2
      common /comsfc/
     &          tsea  (LONF2S,LATG2S),
     &          smc   (LONF2S,LATG2S,lsoil_),
     &          sheleg(LONF2S,LATG2S),
     &          stc   (LONF2S,LATG2S,nsoil_),
     &          tg3   (LONF2S,LATG2S),
     &          zorl  (LONF2S,LATG2S),
     &          cv    (LONF2S,LATG2S),
     &          cvb   (LONF2S,LATG2S),
     &          cvt   (LONF2S,LATG2S),
     &          albedo(LONF2S,LATG2S,lalbd_),
     &          slmsk (LONF2S,LATG2S),
     &          vfrac (LONF2S,LATG2S),
     &          canopy(LONF2S,LATG2S),
     &          f10m  (LONF2S,LATG2S),
     &          vtype (LONF2S,LATG2S),
     &          vroot (LONF2S,LATG2S,kslmb),
     &          facalf(LONF2S,LATG2S,2),
     &          uustar(LONF2S,LATG2S),
     &          ffmm  (LONF2S,LATG2S),
     &          ffhh  (LONF2S,LATG2S),
     &          prcp  (LONF2S,LATG2S),
     &          srflag(LONF2S,LATG2S),
     &          binf  (LONF2S,LATG2S),
     &          ds    (LONF2S,LATG2S),
     &          dsm   (LONF2S,LATG2S),
     &          ws    (LONF2S,LATG2S),
     &          cef   (LONF2S,LATG2S),
     &          expt  (LONF2S,LATG2S,lsoil_),
     &          kest  (LONF2S,LATG2S,lsoil_),
     &          dph   (LONF2S,LATG2S,lsoil_),
     &          bub   (LONF2S,LATG2S,lsoil_),
     &          qrt   (LONF2S,LATG2S,lsoil_),
     &          bkd   (LONF2S,LATG2S,lsoil_),
     &          sld   (LONF2S,LATG2S,lsoil_),
     &          wcr   (LONF2S,LATG2S,lsoil_),
     &          wpw   (LONF2S,LATG2S,lsoil_),
     &          smr   (LONF2S,LATG2S,lsoil_),
     &          smx   (LONF2S,LATG2S,lsoil_),
     &          dphn  (LONF2S,LATG2S,nsoil_),
     &          smxn  (LONF2S,LATG2S,nsoil_),
     &          expn  (LONF2S,LATG2S,nsoil_),
     &          bubn  (LONF2S,LATG2S,nsoil_),
     &          alpn  (LONF2S,LATG2S,nsoil_),
     &          betn  (LONF2S,LATG2S,nsoil_),
     &          gamn  (LONF2S,LATG2S,nsoil_),
     &          silz  (LONF2S,LATG2S),
     &          snwz  (LONF2S,LATG2S),
     &          nveg  (LONF2S,LATG2S),
     &          mvfr  (LONF2S,LATG2S,msub_),
     &          mcnp  (LONF2S,LATG2S,msub_),
     &          mvty  (LONF2S,LATG2S,msub_),
     &          flai  (LONF2S,LATG2S,msub_),
     &          msno  (LONF2S,LATG2S,msub_),
     &          msmc  (LONF2S,LATG2S,kslmb),
     &          msic  (LONF2S,LATG2S,kslmb),
     &          mstc  (LONF2S,LATG2S,nslmb),
     &          csno  (LONF2S,LATG2S,msub_),
     &          rsno  (LONF2S,LATG2S,msub_),
     &          tsf   (LONF2S,LATG2S,msub_),
     &          tpk   (LONF2S,LATG2S,msub_),
     &          sfw   (LONF2S,LATG2S,msub_),
     &          pkw   (LONF2S,LATG2S,msub_),
     &          lstsn (LONF2S,LATG2S,msub_)
#endif
