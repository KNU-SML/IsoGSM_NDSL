#define CP cp
#define G g
#define HVAP hvap
#define RD rd
#define RV rv
      subroutine t93clds_driver
     * (im,ix,km,deltim,si,sl,del,ps,q,lw,t,a,
     *  dtsur,dqsur,fqs,fss,hpbl,slmsk,
     *  vvel,dlf,swhr,lwhr,
     *  cmf,du,solin,
     *  rn, lat)
c*****************************************
#include <paramodel.h>

c
#undef  MADJSS   ! moisture adjustment for supersaturation
#define MADJSS2  ! moisture adjustment for supersaturation
#undef  POUT
#define A0EVAP
c

      parameter (cb2pa=1.e3)
      parameter (hvap=hvap_, g=g_, cp=cp_, rd=rd_)
      parameter (rcp=1.d0/cp, rr=1.d0/rd, gor=g/rd, rkap=rd/cp)
      parameter (elocp=hvap/cp)
      parameter (pref=1.d5)  ! reference pressure = 1000 hPa = 1.d5 Pa
c
c input
c
      integer   im, ix, km
      real      deltim
      real      si(ix,km+1), sl(ix,km), del(ix,km), ps(ix)
      real      slk(km)
      real      dlf(ix,km), du(ix,km)
      real      dtsur(ix,km),dqsur(ix,km)
      real      solin(ix,km)
      real      vvel(ix,km), cmf(ix,km)
      real      swhr(ix,km), lwhr(ix,km)
      real      hpbl(ix), fqs(ix), fss(ix)
      real      slmsk(ix)  ! land/sea/ice mask (0:sea.1:land,2:ice)
      integer   lat
c
c input/output
c
      real      q(ix,km), lw(ix,km), t(ix,km), a(ix,km) 
c
c output
c
      real      rn(ix)
c
c local
c
      real      p(km+1),   pk(km) , ti(km),  th(km),  qi(km)
      real      lwi(km), clw(km), ciw(km), rho(km), dlfi(km), dui(km)
      real      z(km+1),   dtsuri(km), dqsuri(km), solini(km)
      real      wpa(km), cmfi(km), hrtr(km), hrsl(km)
      real      ai(km), pref
      real      prec, precl, precf
      real      dtliw(km), dqliw(km)
      real      sacv(km), sac(km), sabl(km), disp(km)

      real      dthor(km)  ! not used
      integer   iflag      ! not used

      real      coez(km)   ! for calculating z

#ifdef T93PR2
#define CLDT93V2
#endif

#ifdef CLDT93V2
      real      lfrc
#endif
c
#ifdef DBG
#ifdef MP
c     from commpi.h i should use include commpi.h
      integer npes,ncol,nrow,mype,master,msgtag
      common /commpi/ npes,ncol,nrow,mype,master,msgtag
c     from commpi.h
#endif
      integer   nnn(12)
#endif
c
c  moisture adjustment for supersaturation
c
#ifdef MADJSS
      parameter(rv=rv_, rg=1.d0/g)
      parameter(el2orc=hvap*hvap/(rv*cp),
     &          eps=rd/rv,epsm1=rd/rv-1.)
      parameter(epsq=2.e-12, qsmx=1.e4)
      parameter(rhcrit=0.8)
#endif
c---
#ifdef MADJSS2
      parameter(rv=rv_, rg=1.d0/g)
      parameter(el2orc=hvap*hvap/(rv*cp),
     &          eps=rd/rv,epsm1=rd/rv-1.)
      parameter(epsq=2.e-12, qsmx=1.e4)
      parameter(rhcrit=0.8)
#endif

c
c  output values at a point (ipos,jpos)
c
#ifdef POUT
      parameter (ipos=7,jpos=1)
      save iout 
      data iout /0/   ! after output iout=1
      real rh(km)
#ifndef MADJSS
#ifndef MADJSS2
      parameter(rv=rv_, rg=1.d0/g)
      parameter(el2orc=hvap*hvap/(rv*cp),
     &          eps=rd/rv,epsm1=rd/rv-1.)
      parameter(epsq=2.e-12, qsmx=1.e4)
#endif
#else
#endif
#endif

c
c     start
c
      dt=2*deltim
      rdt=1.0/dt
      dtmmd2m = dt/86400.d0*1.d-3  ! [mm/day] -> [m]
c
c     initialization
      rn(1:im)=0.0

c
c     main
c

      do i=1,im

!       calculate z : translate from moninp
      do k = 1, km-1
        coez(k) = log(si(i,k+1) / si(i,k)) / gor
      enddo
      do k=1,km
        slk(k)=sl(i,k)
      enddo

#ifdef MADJSS2
c
c     moisture adjustment for supersaturation  
c     -> lrgscl.F
c
        prec=0.0
        do k=km,1,-1
#ifdef ICE
          es=fpvs(t(i,k))
#else
          es=fpvs0(t(i,k))
#endif
          qs=eps*es/(sl(i,k)*ps(i)+epsm1*es)
          if(qs.lt.0.0) qs=qsmx
          qs = max(qs,epsq)
          qcond=(q(i,k)-qs)/(1.+el2orc*qs/(t(i,k)*t(i,k)))
          if(qcond.gt.0.) then
            dpovg=del(i,k)/g*ps(i)
            q(i,k)=q(i,k)-qcond
            t(i,k)=t(i,k)+qcond*elocp
            prec=prec+qcond*dpovg
          elseif(prec.gt.0.) then
c            rh=q(i,k)/qs
c            if (rh.lt.rhcrit) then
              dpovg=del(i,k)/g*ps(i)
              qevap=-qcond*(1.-exp(-0.32*sqrt(dt*prec)))
              rnevap=min(qevap*dpovg,prec)
              q(i,k)=q(i,k)+rnevap/dpovg
              t(i,k)=t(i,k)-rnevap/dpovg*elocp
              prec=prec-rnevap
c            endif
          endif
        end do
        rn(i) = rn(i) + prec  ! [m]
c
c     end: moisture adjustment for supersaturation  
c
#endif

        ti(1:km)     = t(i,1:km)            
        qi(1:km)     = max(q(i,1:km),0.0)           ! must be >=0                    
        lwi(1:km)    = max(lw(i,1:km),0.0)         
        ai(1:km)     = min(max(a(i,1:km),0.0),1.0)  ! cloud amount must be 0-1          
        dlfi(1:km)   = dlf(i,1:km)*rdt     ! [kg/kg] -> [kg/kg/s]
        dui(1:km)    = du(i,1:km)*rdt      ! [-] -> [/s]                 
        dtsuri(1:km) = dtsur(i,1:km)           
        dqsuri(1:km) = dqsur(i,1:km)           
        solini(1:km) = solin(i,1:km)           
        wpa(1:km)    = vvel(i,1:km) * cb2pa        ! [cb/s]->[Pa/s]
        cmfi(1:km)   = cmf(i,1:km)               
        hrtr(1:km)   = lwhr(i,1:km)              
        hrsl(1:km)   = swhr(i,1:km)              

        pk(1:km)  = ps(i)*sl(i,1:km)   * cb2pa       ! [cb]->[Pa]
        p(1:km+1) = ps(i)*si(i,1:km+1) * cb2pa       ! [cb]->[Pa] 

        rho(1:km) = rr * pk(1:km) / ti(1:km)      ! rho=p/(r*t)
        th(1:km)  = ti(1:km) * (pref/pk(1:km))**rkap ! theta [K]

!       divide lw between clw and ciw
        do k = 1, km
          tmp    = fliq(ti(k))        ! rate of liquid    
          clw(k) = lwi(k) * tmp
          ciw(k) = lwi(k) - clw(k)
        end do

!       calculate z : translate from moninp
        z(1) = 0.d0
        do k = 1, km-1
          z(k+1) = z(k) - ti(k) * coez(k)
        enddo
        z(km+1) = 5.e4  ! for instance

#ifdef CLDT93V2
        lfrc = min(slmsk(i),1.0)  ! =0:sea, =1:land&ice
#endif

#ifdef POUT
        if(i.eq.ipos.and.lat.eq.jpos) then
          if (iout.eq.0) then

          nf=98
          open (unit=nf,file='inputcld.txt',form='formatted',
     &          status='unknown')
          call pointout(nf,'sgk       ',km  ,sl)
          call pointout(nf,'p         ',km+1,p)
          call pointout(nf,'pk        ',km  ,pk)
          call pointout(nf,'t         ',km  ,ti)
          call pointout(nf,'th        ',km  ,th)
          call pointout(nf,'q         ',km  ,qi)
          call pointout(nf,'lw        ',km  ,lwi)
          call pointout(nf,'clw       ',km  ,clw)
          call pointout(nf,'ciw       ',km  ,ciw)
          call pointout(nf,'rho       ',km  ,rho)
          call pointout(nf,'dlf       ',km  ,dlfi)
          call pointout(nf,'du        ',km  ,dui)
          call pointout(nf,'z         ',km+1,z)
          call pointout(nf,'dtsur     ',km  ,dtsuri)
          call pointout(nf,'dqsur     ',km  ,dqsuri)
          call pointout(nf,'wpa       ',km  ,wpa)
          call pointout(nf,'cmf       ',km  ,cmfi)
          call pointout(nf,'hrtr      ',km  ,hrtr)
          call pointout(nf,'hrsl      ',km  ,hrsl)
          call pointout(nf,'hpbl      ',1   ,hpbl(i))
          call pointout(nf,'fqs       ',1   ,fqs(i))
          call pointout(nf,'fss       ',1   ,fss(i))
          call pointout(nf,'pref      ',1   ,pref)
          call pointout(nf,'delt      ',1   ,dt)
          call pointout(nf,'a         ',km  ,ai)
          do k=1,km
#ifdef ICE
            es=fpvs(ti(k))
#else
            es=fpvs0(ti(k))
#endif
            qs=eps*es/(sl(i,k)*ps(i)+epsm1*es)
            if(qs.lt.0.0) qs=qsmx
            qs = max(qs,epsq)
            rh(k)=qi(k)/qs
          end do
          call pointout(nf,'rh        ',km  ,rh)

c          do k=1,km-1
c            rh(k)=(th(k+1)-th(k))/(pk(k+1)-pk(k))
c            rh(k)=rh(k)*100.  ! [K/Pa]->[K/hPa]
c          end do
c          rh(km)=0.0
c          call pointout(nf,'dthdp     ',km  ,rh)
c          do k=1,km-1
c            rh(k)=(ti(k+1)-ti(k))/(z(k+1)-z(k))
c            rh(k)=rh(k)*1000.  ! [K/m]->[K/km]
c          end do
c          rh(km)=0.0
c          call pointout(nf,'dtdz      ',km  ,rh)
c
          
          close(nf)
          iout=1

          endif
        endif
#endif

#ifdef CLDT93V2
#define T93CLDS t93cldsv2
#else
#define T93CLDS t93clds
#endif

        call T93CLDS(       
     1               slk   ,p      ,pk     ,ti     ,th    ,  
#ifdef CLDT93V2
     2               qi    ,        clw    ,ciw    ,rho   ,dlfi   ,dui ,
#else
     2               qi    ,lwi    ,clw    ,ciw    ,rho   ,dlfi   ,dui ,
#endif
     3               z     ,dtsuri ,dqsuri ,solini ,wpa   , 
     4               cmfi  ,hrtr   ,hrsl   ,dthor  ,hpbl(i), 
     5               fqs(i),fss(i) ,pref   ,dt     ,ai    , 
     6               dtliw ,dqliw  ,prec   ,precl  ,precf ,iflag,
#ifdef CLDT93V2
     7               sacv  ,sac    ,sabl   ,disp   ,lfrc)
#else
     7               sacv  ,sac    ,sabl   ,disp )
#endif

#ifdef CLDT93V2
        lwi(1:km)=clw(1:km)+ciw(1:km)  ! total cloud water (liquid+ice)
#endif

#ifdef POUT
        if(i.eq.ipos.and.lat.eq.jpos) then
          if (iout.eq.1) then

          nf=99
          open (unit=nf,file='outputcld.txt',form='formatted',
     &          status='unknown')
          call pointout(nf,'a         ',km  ,ai)
          call pointout(nf,'lw        ',km  ,lwi)
          call pointout(nf,'clw       ',km  ,clw)
          call pointout(nf,'ciw       ',km  ,ciw)
          call pointout(nf,'prec      ',1   ,prec)
          call pointout(nf,'precl     ',1   ,precl)
          call pointout(nf,'precf     ',1   ,precf)
          call pointout(nf,'dtliw     ',km  ,dtliw)
          call pointout(nf,'dqliw     ',km  ,dqliw)
          call pointout(nf,'sacv      ',km  ,sacv)
          call pointout(nf,'sac       ',km  ,sac)
          call pointout(nf,'sabl      ',km  ,sabl)
          call pointout(nf,'disp      ',km  ,disp)

          close(nf)
          iout=2

          stop
          endif

        endif
#endif

c
c       output
c
        do k=1,km
          t(i,k)  = t(i,k) + dtliw(k)
          q(i,k)  = q(i,k) + dqliw(k)
          lw(i,k) = lwi(k)
          a(i,k)  = ai(k)
        end do
        rn(i)= rn(i) + prec * dtmmd2m      ! [m] prec=[mm/day]


#ifdef MADJSS
c
c     moisture adjustment for supersaturation  
c     -> lrgscl.F
c
        prec=0.0
        do k=km,1,-1
#ifdef ICE
          es=fpvs(t(i,k))
#else
          es=fpvs0(t(i,k))
#endif
          qs=eps*es/(sl(i,k)*ps(i)+epsm1*es)
          if(qs.lt.0.0) qs=qsmx
          qs = max(qs,epsq)
          qcond=(q(i,k)-qs)/(1.+el2orc*qs/(t(i,k)*t(i,k)))
          if(qcond.gt.0.) then
            dpovg=del(i,k)/g*ps(i)
            q(i,k)=q(i,k)-qcond
            t(i,k)=t(i,k)+qcond*elocp
            prec=prec+qcond*dpovg
          elseif(prec.gt.0.) then
c            rh=q(i,k)/qs
c            if (rh.lt.rhcrit) then
              dpovg=del(i,k)/g*ps(i)
              qevap=-qcond*(1.-exp(-0.32*sqrt(dt*prec)))
              rnevap=min(qevap*dpovg,prec)
              q(i,k)=q(i,k)+rnevap/dpovg
              t(i,k)=t(i,k)-rnevap/dpovg*elocp
              prec=prec-rnevap
c            endif
          endif
        end do
        rn(i) = rn(i) + prec  ! [m]
c
c     end: moisture adjustment for supersaturation  
c
#endif

#ifdef A0EVAP
c
c  a=0 -> lw=0 (evap)
c
        do k=1,km
          if (a(i,k).le.0.0.and.lw(i,k).gt.0.0) then
            q(i,k) = q(i,k) + lw(i,k)
            t(i,k) = t(i,k) - elocp * lw(i,k) 
            lw(i,k) = 0.0
          end if  
        end do
c
c  end
c
#endif
c
c  loop i end
c
      end do

#ifdef DBG
#ifdef MP
      if (mype.eq.master) then
#endif
      nnn(:)=0
      tmp1=0.0
      tmp2=0.0

      do k=1,km
      do i=1,im
        tmp1=max(tmp1,a(i,k))
        tmp2=min(tmp2,a(i,k))

        if (a(i,k).le.0.0) then
          nn=1
        else if (a(i,k).ge.1.0) then
          nn=12
        else
          nn=nint(a(i,k)*10.)+2
          nn=min(max(1,nn),12)
          nnn(nn)=nnn(nn)+1 
        endif
      end do
      end do
      print *,'cl1:',tmp1,tmp2
      print *,'cl2:',(nnn(n),n=1,12)

#ifdef MP
      endif
#endif
#endif

c  move water from vapor to liquid should the liquid amount be negative 
c         from [precpd.F]
c
cx      do k=1,km
cx      do i=1,im
cx          if (lw(i,k).lt.0.d0) then
cx            q(i,k)   = q(i,k) + lw(i,k)
cx            t(i,k)   = t(i,k) - elocp * lw(i,k)
cx            lw(i,k)  = 0.d0
cx          endif
cx      end do
cx      end do
c
      return
      end

      subroutine pointout
     &    (nf,vname,km,data)
c*********************************************************
      integer nf
      character*10 vname
      integer km
      real    data(km)
c
      write(nf,*) vname,km
      do k=1,km
        write(nf,*) k,data(k)
      end do
c
      return
      end
