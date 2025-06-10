#include <define.h>
      subroutine kf_eta_para (i, j,                           
     1                u0,v0,t0,qv0,p0,dzq,w0avg1d,         
     2                dt,dx,dxsq,rhoe,                     
     3                xlv0,xlv1,xls0,xls1,cp,r,g,          
     4                ep2,svp1,svp2,svp3,svpt0,            
     5                dqdt,dqidt,dqcdt,dqrdt,dqsdt,dtdt,   
     6                raincv,nca,ntst,                     
     7                p_qi,p_qs,p_first_scalar,warm_rain,  
     8                ids,ide, jds,jde, kds,kde,           
     9                ims,ime, jms,jme, kms,kme,           
     1                its,ite, jts,jte, kts,kte,
     2                ltop, klcl)

      implicit none

      integer, intent(in   ) :: ids,ide, jds,jde, kds,kde, 
     1                          ims,ime, jms,jme, kms,kme, 
     2                          its,ite, jts,jte, kts,kte, 
     3                          i,j,ntst,p_qi,p_qs,p_first_scalar

      logical, intent(in   ) :: warm_rain

      real, dimension( kts:kte ),                          
     1      intent(in   ) ::                           u0, 
     2                                                 v0, 
     3                                                 t0, 
     4                                                qv0, 
     5                                                 p0, 
     6                                               rhoe, 
     7                                                dzq, 
     8                                            w0avg1d

      real,  intent(in   ) :: dt,dx,dxsq
      real,  intent(in   ) :: xlv0,xlv1,xls0,xls1,cp,r,g
      real,  intent(in   ) :: ep2,svp1,svp2,svp3,svpt0

      real, dimension( kts:kte ), intent(inout) ::         
     1                                               dqdt, 
     2                                              dqidt, 
     3                                              dqcdt, 
     4                                              dqrdt, 
     5                                              dqsdt, 
     6                                               dtdt

      real,    dimension( ims:ime , jms:jme ),             
     1      intent(inout) ::                          nca

      real, dimension( ims:ime , jms:jme ),                
     1      intent(inout) ::                       raincv

!... common variables
      integer, parameter :: kfnt=250,kfnp=220
      real, dimension(kfnt,kfnp) :: ttab,qstab
      real, dimension(kfnp) :: the0k
      real, dimension(200) :: alu
      real :: rdpr,rdthk,plutop
      common/kflut/ ttab,qstab,the0k,
     *              alu,rdpr,rdthk,plutop

!...define local variables...

      real, dimension( kts:kte ) ::                        
     1      q0,z0,tv0,tu,tvu,qu,tz,tvd,                    
     2      qd,qes,thtes,tg,tvg,qg,wu,wd,w0,ems,emsd,      
     3      umf,uer,udr,dmf,der,ddr,umf2,uer2,             
     4      udr2,dmf2,der2,ddr2,dza,thta0,thetee,         
     5      thtau,theteu,thtad,theted,qliq,qice,          
     6      qlqout,qicout,pptliq,pptice,detlq,detic,      
     7      detlq2,detic2,ratio,ratio2


      real, dimension( kts:kte ) ::                        
     1      domgdp,exn,tvqu,dp,rh,eqfrc,wspd,              
     2      qdt,fxm,thtag,thpa,thfxout,                    
     3      thfxin,qpa,qfxout,qfxin,qlpa,qlfxin,           
     4      qlfxout,qipa,qifxin,qifxout,qrpa,              
     5      qrfxin,qrfxout,qspa,qsfxin,qsfxout,            
     6      ql0,qlg,qi0,qig,qr0,qrg,qs0,qsg


      real, dimension( kts:kte+1 ) :: omg
      real, dimension( kts:kte ) :: rainfb,snowfb
      real, dimension( kts:kte ) ::                        
     1      cldhgt,qsd,dilfrc,ddilfrc,tke,tgu,qgu,thteeg

! local vars

      real    :: p00,t00,rlf,rhic,rhbc,pie,         
     1           ttfrz,tbfrz,c5,rate
      real    :: gdry,rocp,aliq,bliq,                      
     1           cliq,dliq
      real    :: fbfrc,p300,dpthmx,thmix,qmix,zmix,pmix,   
     1           rocpq,tmix,emix,tlog,tdpt,tlcl,tvlcl,     
     2           cporq,plcl,es,dlp,tenv,qenv,tven,tvbar,   
     3           zlcl,wkl,wabs,trppt,wsigne,dtlcl,gdt,wlcl,
     4           tvavg,qese,wtw,rholcl,au0,vmflcl,upold,   
     5           upnew,abe,wklcl,ttemp,frc1,   
     6           qnewic,rl,r1,qnwfrz,effq,be,boterm,enterm,
     7           dzz,udlbe,rei,ee2,ud2,ttmp,f1,f2,         
     8           thttmp,qtmp,tmpliq,tmpice,tu95,tu10,ee1,  
     9           ud1,dptt,qnewlq,dumfdp,ee,tsat,           
     1           thta,vconv,timec,shsign,vws,pef, 
     2           cbh,rcbh,pefcbh,peff,peff2,tder,thtmin,   
     3           dtmltd,qs,tadvec,dpdd,frc,dpt,rdd,a1,     
     4           dssdt,dtmp,t1rh,qsrh,pptflx,cpr,cndtnf,   
     5           updinc,aincm2,devdmf,ppr,rced,dpptdf,     
     6           dmflfs,dmflfs2,rced2,ddinc,aincmx,aincm1, 
     7           ainc,tder2,pptfl2,fabe,stab,dtt,dtt1,     
     8           dtime,tma,tmb,tmm,bcoeff,acoeff,qvdiff,   
     9           topomg,cpm,dq,abeg,dabe,dfda,frc2,dr,     
     1           udfrc,tuc,qgs,rh0,rhg,qinit,qfnl,err2,    
     2           relerr,rlc,rls,rnc,fabeold,aincold,uefrc, 
     3           ddfrc,tdc,defrc,rhbar,dmffrc,dpmin,dilbe
      real    ::    astrt,tp,value,aintrp,tkemax,qfrz,
     1           qss,pptmlt,dtmelt,rhh,evac,binc

      integer :: indlu,nu,nuchm,nnn,klfs
      real    :: chmin,pm15,chmax,dtrh,rad,dppp
      real    :: tvdiff,dttot,absomg,absomgtc,frdp

      integer :: kx,k,kl

      integer :: ncheck
      integer, dimension (kts:kte) :: kcheck

      integer :: istop,ml,l5,kmix,low,                     
     1           lc,mxlayr,llfc,nlayrs,nk,                 
     2           kpbl,klcl,lcl,let,iflag,                  
     3           nk1,ltop,nj,ltop1,                        
     4           ltopm1,lvf,kstart,kmin,lfs,               
     5           nd,nic,ldb,ldt,nd1,ndk,                   
     6           nm,lmax,ncount,noitr,                     
     7           nstep,ntc,nchm,ishall,nshall
      logical :: iprnt

      data p00,t00/1.e5,273.16/
      data rlf/3.339e5/
      data rhic,rhbc/1.,0.90/
      data pie,ttfrz,tbfrz,c5/3.141592654,268.16,248.16,1.0723e-3/
      data rate/0.03/
!-----------------------------------------------------------
      iprnt=.false.
      gdry=-g/cp
      rocp=r/cp
      nshall = 0
      kl=kte
      kx=kte

!     aliq = 613.3
!     bliq = 17.502
!     cliq = 4780.8
!     dliq = 32.19
      aliq = svp1*1000.
      bliq = svp2
      cliq = svp2*svpt0
      dliq = svp3
!

!
!**********************************************************************
!                                                      ! ppt fb mods
!...option to feed convectively generated rainwater    ! ppt fb mods
!...into grid-resolved rainwater (or snow/graupel)     ! ppt fb mods
!...field.  "fbfrc" is the fraction of available       ! ppt fb mods
!...precipitation to be fed back (0.0 - 1.0)...        ! ppt fb mods
      fbfrc=0.0                                        ! ppt fb mods
!...mods to allow shallow convection...
      nchm = 0
      ishall = 0
      dpmin = 5.e3
!...
      p300=p0(1)-30000.
!
!...pressure perturbation term is only defined at mid-point of
!...vertical layers...since total pressure is needed at the top and
!...bottom of layers below, do an interpolation...
!
!...input a vertical sounding ... note that model layers are numbered
!...from bottom-up in the kf scheme...
!
      ml=0 
!sue  tmprpsb=1./psb(i,j)
!sue  cell=ptop*tmprpsb
!
!
! initialize the level index
!
      l5=1
      llfc=1
      ml=1
!
      do k=1,kx
!
!...if q0 is above saturation value, reduce it to saturation level...
!
         es=aliq*exp((bliq*t0(k)-cliq)/(t0(k)-dliq))
         qes(k)=0.622*es/(p0(k)-es)
         q0(k)=amin1(qes(k),qv0(k))
         q0(k)=amax1(0.000001,q0(k))
         ql0(k)=0.
         qi0(k)=0.
         qr0(k)=0.
         qs0(k)=0.
         rh(k) = q0(k)/qes(k)
         dilfrc(k) = 1.
         tv0(k)=t0(k)*(1.+0.608*q0(k))
!        rhoe(k)=p0(k)/(r*tv0(k))
!   dp is the pressure interval between full sigma levels...
         dp(k)=rhoe(k)*g*dzq(k)
! if turbulent kinetic energy (tke) is available from turbulent mixing scheme
! use it for shallow convection..for now, assume it is not available....
!         tke(k) = q2(i,j,nk)
         tke(k) = 0.
         cldhgt(k) = 0.
         if(p0(k).ge.500e2)l5=k
         if(p0(k).ge.p300)llfc=k
         if(t0(k).gt.t00)ml=k
      enddo

!
!...dzq is dz between sigma surfaces, dza is dz between model half level
        z0(1)=.5*dzq(1)
        do k=2,kl
          z0(k)=z0(k-1)+.5*(dzq(k)+dzq(k-1))
          dza(k-1)=z0(k)-z0(k-1)
        enddo   
        dza(kl)=0.

!
!
!  to save time, specify a pressure interval to move up in sequential
!  check of different ~50 mb deep groups of adjacent model layers in
!  the process of identifying updraft source layer (usl).  note that 
!  this search is terminated as soon as a buoyant parcel is found and 
!  this parcel can produce a cloud greater than specifed minimum depth
!  (chmin)...for now, set interval at 15 mb...
!
       ncheck = 1
       kcheck(ncheck)=1
       pm15 = p0(1)-15.e2
       do k=2,llfc
         if(p0(k).lt.pm15)then
           ncheck = ncheck+1
           kcheck(ncheck) = k
           pm15 = pm15-15.e2
         endif
       enddo
!
       nu=0
       nuchm=0
       usl:   do
           nu = nu+1
           if(nu.gt.ncheck)then 
             if(ishall.eq.1)then
               chmax = 0.
               nchm = 0
               do nk = 1,ncheck
                 nnn=kcheck(nk)
                 if(cldhgt(nnn).gt.chmax)then
                   nchm = nnn
                   nuchm = nk
                   chmax = cldhgt(nnn)
                 endif
               enddo
               nu = nuchm-1
               fbfrc=1.
               cycle usl
             else
               return
             endif
           endif      
           kmix = kcheck(nu)
           low=kmix

!...
           lc = low
!
!...assume that in order to support a deep updraft you need a layer of
!...unstable air at least 50 mb deep...to approximate this, isolate a
!...group of adjacent individual model layers, with the base at level
!...lc, such that the combined depth of these layers is at least 50 mb..
!   
           nlayrs=0
           dpthmx=0.
           nk=lc-1
           do 
             nk=nk+1   
             dpthmx=dpthmx+dp(nk)
             nlayrs=nlayrs+1
             if(dpthmx.gt.dpmin)then
               exit 
             endif
           end do    
           if(dpthmx.lt.dpmin)then 
             return
           endif
           kpbl=lc+nlayrs-1   

!
!...********************************************************
!...for computational simplicity without much loss in accuracy,
!...mix temperature instead of theta for evaluating convective
!...initiation (triggering) potential...
!          thmix=0.
           tmix=0.
           qmix=0.
           zmix=0.
           pmix=0.
!
!...find the thermodynamic characteristics of the layer by
!...mass-weighting the characteristics of the individual model
!...layers...
!
           do nk=lc,kpbl
             tmix=tmix+dp(nk)*t0(nk)
             qmix=qmix+dp(nk)*q0(nk)
             zmix=zmix+dp(nk)*z0(nk)
             pmix=pmix+dp(nk)*p0(nk)
           enddo   
!         thmix=thmix/dpthmx
          tmix=tmix/dpthmx
          qmix=qmix/dpthmx
          zmix=zmix/dpthmx
          pmix=pmix/dpthmx
          emix=qmix*pmix/(0.622+qmix)

!
!...find the temperature of the mixture at its lcl...
!
!        tlog=alog(emix/aliq)
! ...calculate dewpoint using lookup table...
!
          astrt=1.e-3
          ainc=0.075
          a1=emix/aliq
          tp=(a1-astrt)/ainc
          indlu=int(tp)+1
          value=(indlu-1)*ainc+astrt
          aintrp=(a1-value)/ainc
          tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
          tdpt=(cliq-dliq*tlog)/(bliq-tlog)
          tlcl=tdpt-(.212+1.571e-3*(tdpt-t00)-4.36e-4*(tmix-t00))
     1        *(tmix-tdpt)
          tlcl=amin1(tlcl,tmix)
          tvlcl=tlcl*(1.+0.608*qmix)
          zlcl = zmix+(tlcl-tmix)/gdry
          nk = lc-1

          do 
            nk = nk+1
            klcl=nk
            if(zlcl.le.z0(nk) .or. nk.gt.kl)then
              exit
            endif 
          enddo   
          if(nk.gt.kl)then
            return  
          endif
          k=klcl-1
          dlp=(zlcl-z0(k))/(z0(klcl)-z0(k))

!     
!...estimate environmental temperature and mixing ratio at the lcl...
!     
          tenv=t0(k)+(t0(klcl)-t0(k))*dlp
          qenv=q0(k)+(q0(klcl)-q0(k))*dlp
          tven=tenv*(1.+0.608*qenv)

!     
!...check to see if cloud is buoyant using fritsch-chappell trigger
!...function described in kain and fritsch (1992)...w0 is an
!...aproximate value for the running-mean grid-scale vertical
!...velocity, which gives smoother fields of convective initiation
!...than the instantaneous value...formula relating temperature
!...perturbation to vertical velocity has been used with the most
!...success at grid lengths near 25 km.  for different grid-lengths,
!...adjust vertical velocity to equivalent value for 25 km grid
!...length, assuming linear dependence of w on grid length...
!     
          if(zlcl.lt.2.e3)then
            wklcl=0.02*zlcl/2.e3
          else
            wklcl=0.02
          endif
          wkl=(w0avg1d(k)+(w0avg1d(klcl)-w0avg1d(k))*dlp)*dx/25.e3-wklcl
          if(wkl.lt.0.0001)then
            dtlcl=0.
          else 
            dtlcl=4.64*wkl**0.33
          endif
!
!...for eta model, give parcel an extra temperature perturbation based
!...the threshold rh for condensation (u00)...
!
!...for now, just assume u00=0.75...
!...!!!!!! for mm5, set dtrh = 0. !!!!!!!!
!         u00 = 0.75
!         if(u00.lt.1.)then
!           qslcl=qes(k)+(qes(klcl)-qes(k))*dlp
!           rhlcl = qenv/qslcl
!           dqssdt = qmix*(cliq-bliq*dliq)/((tlcl-dliq)*(tlcl-dliq))
!           if(rhlcl.ge.0.75 .and. rhlcl.le.0.95)then
!             dtrh = 0.25*(rhlcl-0.75)*qmix/dqssdt
!           elseif(rhlcl.gt.0.95)then
!             dtrh = (1./rhlcl-1.)*qmix/dqssdt
!           else
               dtrh = 0.
!           endif
!         endif   
!         if(ishall.eq.1)iprnt=.true.
!         iprnt=.true.
!         if(tlcl+dtlcl.gt.tenv)goto 45
!
      trigger:  if(tlcl+dtlcl+dtrh.lt.tenv)then   
!
! parcel not buoyant, cycle back to start of trigger and evaluate next potential usl...
!
            cycle usl
!
          else                   ! parcel is buoyant, determine updraft
!     
!...convective triggering criteria has been satisfied...compute
!...equivalent potential temperature
!...(theteu) and vertical velocity of the rising parcel at the lcl...
!     
          call envirtht(pmix,tmix,qmix,theteu(k),aliq,bliq,cliq,dliq)
!
!..modify calculation of initial parcel vertical velocity...jsk 11/26/97
!
            dttot = dtlcl+dtrh
            if(dttot.gt.1.e-4)then
              gdt=2.*g*dttot*500./tven
              wlcl=1.+0.5*sqrt(gdt)
              wlcl = amin1(wlcl,3.)
            else
              wlcl=1.
            endif
            plcl=p0(k)+(p0(klcl)-p0(k))*dlp
            wtw=wlcl*wlcl
!
            tvlcl=tlcl*(1.+0.608*qmix)
            rholcl=plcl/(r*tvlcl)
!        
            lcl=klcl
            let=lcl
! make rad a function of background vertical velocity...
            if(wkl.lt.0.)then
              rad = 1000.
            elseif(wkl.gt.0.1)then
              rad = 2000.
            else
              rad = 1000.+1000*wkl/0.1
            endif
!     
!*******************************************************************
!                                                                  *
!                 compute updraft properties                       *
!                                                                  *
!*******************************************************************
!     
!     
!...
!...estimate initial updraft mass flux (umf(k))...
!     
            wu(k)=wlcl
            au0=0.01*dxsq
            umf(k)=rholcl*au0
            vmflcl=umf(k)
            upold=vmflcl
            upnew=upold
!     
!...ratio2 is the degree of glaciation in the cloud (0 to 1),
!...uer is the envir entrainment rate, abe is available
!...buoyant energy, trppt is the total rate of precipitation
!...production...
!     
            ratio2(k)=0.
            uer(k)=0.
            abe=0.
            trppt=0.
            tu(k)=tlcl
            tvu(k)=tvlcl
            qu(k)=qmix
            eqfrc(k)=1.
            qliq(k)=0.
            qice(k)=0.
            qlqout(k)=0.
            qicout(k)=0.
            detlq(k)=0.
            detic(k)=0.
            pptliq(k)=0.
            pptice(k)=0.
            iflag=0
!     
!...ttemp is used during calculation of the linear glaciation
!...process; it is initially set to the temperature at which
!...freezing is specified to begin.  within the glaciation
!...interval, it is set equal to the updraft temp at the
!...previous model level...
!     
            ttemp=ttfrz
!     
!...enter the loop for updraft calculations...calculate updraft temp,
!...mixing ratio, vertical mass flux, lateral detrainment of mass and
!...moisture, precipitation rates at each model level...
!     
!     
            ee1=1.
            ud1=0.
            rei = 0.
            dilbe = 0.
      updraft:    do nk=k,kl-1
              nk1=nk+1
              ratio2(nk1)=ratio2(nk)
              frc1=0.
              tu(nk1)=t0(nk1)
              theteu(nk1)=theteu(nk)
              qu(nk1)=qu(nk)
              qliq(nk1)=qliq(nk)
              qice(nk1)=qice(nk)
              call tpmix2(p0(nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1),
     1               qice(nk1),qnewlq,qnewic,xlv1,xlv0)
!
!
!...check to see if updraft temp is above the temperature at which
!...glaciation is assumed to initiate; if it is, calculate the
!...fraction of remaining liquid water to freeze...ttfrz is the
!...temp at which freezing begins, tbfrz the temp below which all
!...liquid water is frozen at each level...
!
              if(tu(nk1).le.ttfrz)then
                if(tu(nk1).gt.tbfrz)then
                  if(ttemp.gt.ttfrz)ttemp=ttfrz
                  frc1=(ttemp-tu(nk1))/(ttemp-tbfrz)
                else
                  frc1=1.
                  iflag=1
                endif
                ttemp=tu(nk1)
!
!  determine the effects of liquid water freezing when temperature
!...is below ttfrz...
!
                qfrz = (qliq(nk1)+qnewlq)*frc1
                qnewic=qnewic+qnewlq*frc1
                qnewlq=qnewlq-qnewlq*frc1
                qice(nk1) = qice(nk1)+qliq(nk1)*frc1
                qliq(nk1) = qliq(nk1)-qliq(nk1)*frc1
                call dtfrznew(tu(nk1),p0(nk1),theteu(nk1),qu(nk1),qfrz, 
     1                    qice(nk1),aliq,bliq,cliq,dliq)
              endif
              tvu(nk1)=tu(nk1)*(1.+0.608*qu(nk1))
!
!  calculate updraft vertical velocity and precipitation fallout...
!
              if(nk.eq.k)then
                be=(tvlcl+tvu(nk1))/(tven+tv0(nk1))-1.
                boterm=2.*(z0(nk1)-zlcl)*g*be/1.5
                dzz=z0(nk1)-zlcl
              else
                be=(tvu(nk)+tvu(nk1))/(tv0(nk)+tv0(nk1))-1.
                boterm=2.*dza(nk)*g*be/1.5
                dzz=dza(nk)
              endif
              enterm=2.*rei*wtw/upold

              call condload(qliq(nk1),qice(nk1),wtw,dzz,boterm,enterm,
     1                  rate,qnewlq,qnewic,qlqout(nk1),qicout(nk1),g)
!
!...if vert velocity is less than zero, exit the updraft loop and,
!...if cloud is tall enough, finalize updraft calculations...
!
              if(wtw.lt.1.e-3)then
                exit
              else
                wu(nk1)=sqrt(wtw)
              endif
!...calculate value of theta-e in environment to entrain into updraft...
!
              call envirtht(p0(nk1),t0(nk1),q0(nk1),thetee(nk1),aliq,
     1                 bliq,cliq,dliq)
!
!...rei is the rate of environmental inflow...
!
              rei=vmflcl*dp(nk1)*0.03/rad
              tvqu(nk1)=tu(nk1)*(1.+0.608*qu(nk1)-qliq(nk1)-qice(nk1))
              if(nk.eq.k)then
                dilbe=((tvlcl+tvqu(nk1))/(tven+tv0(nk1))-1.)*dzz
              else
                dilbe=((tvqu(nk)+tvqu(nk1))/(tv0(nk)+tv0(nk1))-1.)*dzz
              endif
              if(dilbe.gt.0.)abe=abe+dilbe*g
!
!...if cloud parcels are virtually colder than the environment, minimal 
!...entrainment (0.5*rei) is imposed...
!
              if(tvqu(nk1).le.tv0(nk1))then   ! entrain/detrain if block
                ee2=0.5
                ud2=1.
                eqfrc(nk1)=0.
              else
                let=nk1
                ttmp=tvqu(nk1)
!
!...determine the critical mixed fraction of updraft and environmental air...
!
                f1=0.95
                f2=1.-f1
                thttmp=f1*thetee(nk1)+f2*theteu(nk1)
                qtmp=f1*q0(nk1)+f2*qu(nk1)
                tmpliq=f2*qliq(nk1)
                tmpice=f2*qice(nk1)
                call tpmix2(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,
     1                     qnewlq,qnewic,xlv1,xlv0)
                tu95=ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
                if(tu95.gt.tv0(nk1))then
                  ee2=1.
                  ud2=0.
                  eqfrc(nk1)=1.0
                else
                  f1=0.10
                  f2=1.-f1
                  thttmp=f1*thetee(nk1)+f2*theteu(nk1)
                  qtmp=f1*q0(nk1)+f2*qu(nk1)
                  tmpliq=f2*qliq(nk1)
                  tmpice=f2*qice(nk1)
                  call tpmix2(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,
     1                         qnewlq,qnewic,xlv1,xlv0)
                  tu10=ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
                  tvdiff = abs(tu10-tvqu(nk1))
                  if(tvdiff.lt.1.e-3)then
                    ee2=1.
                    ud2=0.
                    eqfrc(nk1)=1.0
                  else
                    eqfrc(nk1)=(tv0(nk1)-tvqu(nk1))*f1/(tu10-tvqu(nk1))
                    eqfrc(nk1)=amax1(0.0,eqfrc(nk1))
                    eqfrc(nk1)=amin1(1.0,eqfrc(nk1))
                    if(eqfrc(nk1).eq.1)then
                      ee2=1.
                      ud2=0.
                    elseif(eqfrc(nk1).eq.0.)then
                      ee2=0.
                      ud2=1.
                    else
!
!...subroutine prof5 integrates over the gaussian dist to determine the
!   fractional entrainment and detrainment rates...
!
                      call prof5(eqfrc(nk1),ee2,ud2)
                    endif
                  endif
                endif
              endif                  ! end of entrain/detrain if block
!
!
!...net entrainment and detrainment rates are given by the average fractional
!   values in the layer...
!
              ee2 = amax1(ee2,0.5)
              ud2 = 1.5*ud2
              uer(nk1)=0.5*rei*(ee1+ee2)
              udr(nk1)=0.5*rei*(ud1+ud2)
!
!...if the calculated updraft detrainment rate is greater than the total
!   updraft mass flux, all cloud mass detrains, exit updraft calculations...
!
              if(umf(nk)-udr(nk1).lt.10.)then
!
!...if the calculated detrained mass flux is greater than the total upd mass
!   flux, impose total detrainment of updraft mass at the previous model lvl..
!   first, correct abe calculation if needed...
!
                if(dilbe.gt.0.)then
                  abe=abe-dilbe*g
                endif
                let=nk
!               write(98,1015)p0(nk1)/100.
                exit 
              else
                ee1=ee2
                ud1=ud2
                upold=umf(nk)-udr(nk1)
                upnew=upold+uer(nk1)
                umf(nk1)=upnew
                dilfrc(nk1) = upnew/upold
!
!...detlq and detic are the rates of detrainment of liquid and
!...ice in the detraining updraft mass...
!
                detlq(nk1)=qliq(nk1)*udr(nk1)
                detic(nk1)=qice(nk1)*udr(nk1)
                qdt(nk1)=qu(nk1)
                qu(nk1)=(upold*qu(nk1)+uer(nk1)*q0(nk1))/upnew
              theteu(nk1)=(theteu(nk1)*upold+thetee(nk1)*uer(nk1))/upnew
                qliq(nk1)=qliq(nk1)*upold/upnew
                qice(nk1)=qice(nk1)*upold/upnew
!
!...pptliq is the rate of generation (fallout) of
!...liquid precip at a given model lvl, pptice the same for ice,
!...trppt is the total rate of production of precip up to the
!...current model level...
!
                pptliq(nk1)=qlqout(nk1)*umf(nk)
                pptice(nk1)=qicout(nk1)*umf(nk)
!
                trppt=trppt+pptliq(nk1)+pptice(nk1)
                if(nk1.le.kpbl)uer(nk1)=uer(nk1)+vmflcl*dp(nk1)/dpthmx
              endif
!
            end do updraft
!
!...check cloud depth...if cloud is tall enough, estimate the equilibriu
!   temperature level (let) and adjust mass flux profile at cloud top so
!   that mass flux decreases to zero as a linear function of pressure be
!   the let and cloud top...
!     
!...ltop is the model level just below the level at which vertical veloc
!   first becomes negative...
!     
            ltop=nk
            cldhgt(lc)=z0(ltop)-zlcl 
!
!...instead of using the same minimum cloud height (for deep convection)
!...everywhere, try specifying minimum cloud depth as a function of tlcl...
!
!
!
            if(tlcl.gt.293.)then
              chmin = 4.e3
            elseif(tlcl.le.293. .and. tlcl.ge.273)then
              chmin = 2.e3 + 100.*(tlcl-273.)
            elseif(tlcl.lt.273.)then
              chmin = 2.e3
            endif

!     
!...if cloud top height is less than the specified minimum for deep 
!...convection, save value to consider this level as source for 
!...shallow convection, go back up to check next level...
!     
!...try specifying minimum cloud depth as a function of tlcl...
!
!
!...do not allow any cloud from this layer if:
!
!...            1.) if there is no cape, or 
!...            2.) cloud top is at model level just above lcl, or
!...            3.) cloud top is within updraft source layer, or
!...            4.) cloud-top detrainment layer begins within 
!...                updraft source layer.
!
            if(ltop.le.klcl .or. ltop.le.kpbl .or. let+1.le.kpbl)then  
                                              ! no convection allowed
              cldhgt(lc)=0.
              do nk=k,ltop
                umf(nk)=0.
                udr(nk)=0.
                uer(nk)=0.
                detlq(nk)=0.
                detic(nk)=0.
                pptliq(nk)=0.
                pptice(nk)=0.
              enddo
!        
            elseif(cldhgt(lc).gt.chmin .and. abe.gt.1)then      
                                              ! deep convection allowed
              ishall=0
              exit usl
            else
!
!...to disallow shallow convection, comment out next line !!!!!!!!
              ishall = 1
              if(nu.eq.nuchm)then
                exit usl               
                                   ! shallow convection from this layer
              else
! remember this layer (by virtue of non-zero cldhgt) 
! as potential shallow-cloud layer
                do nk=k,ltop
                  umf(nk)=0.
                  udr(nk)=0.
                  uer(nk)=0.
                  detlq(nk)=0.
                  detic(nk)=0.
                  pptliq(nk)=0.
                  pptice(nk)=0.
                enddo
              endif
            endif
          endif trigger
        end do usl
      if(ishall.eq.1)then
        kstart=max0(kpbl,klcl)
        let=kstart
      endif
!     
!...if the let and ltop are the same, detrain all of the updraft mass fl
!   this level...
!     
      if(let.eq.ltop)then
        udr(ltop)=umf(ltop)+udr(ltop)-uer(ltop)
        detlq(ltop)=qliq(ltop)*udr(ltop)*upnew/upold
        detic(ltop)=qice(ltop)*udr(ltop)*upnew/upold
        uer(ltop)=0.
        umf(ltop)=0.
      else 
!     
!   begin total detrainment at the level above the let...
!     
      dptt=0.
      do nj=let+1,ltop
        dptt=dptt+dp(nj)
      enddo
      dumfdp=umf(let)/dptt
!     
!...adjust mass flux profiles, detrainment rates, and precipitation fall
!   rates to reflect the linear decrease in mass flx between the let and
!     
      do nk=let+1,ltop
!
!...entrainment is allowed at every level except for ltop, so disallow
!..entrainment at ltop and adjust entrainment rates between let and ltop
!...so the the dilution factor due to entyrianment is not changed but 
!...the actual entrainment rate will change due due forced total 
!...detrainment in this layer...
!
        if(nk.eq.ltop)then
          udr(nk) = umf(nk-1)
          uer(nk) = 0.
          detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)
          detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)
        else
          umf(nk)=umf(nk-1)-dp(nk)*dumfdp
          uer(nk)=umf(nk)*(1.-1./dilfrc(nk))
          udr(nk)=umf(nk-1)-umf(nk)+uer(nk)
          detlq(nk)=udr(nk)*qliq(nk)*dilfrc(nk)
          detic(nk)=udr(nk)*qice(nk)*dilfrc(nk)
        endif
        if(nk.ge.let+2)then
          trppt=trppt-pptliq(nk)-pptice(nk)
          pptliq(nk)=umf(nk-1)*qlqout(nk)
          pptice(nk)=umf(nk-1)*qicout(nk)
          trppt=trppt+pptliq(nk)+pptice(nk)
        endif
      enddo
      endif
!     
! initialize some arrays below cloud base and above cloud top...
!
      do nk=1,k
        if(nk.ge.lc)then
          if(nk.eq.lc)then
            umf(nk)=vmflcl*dp(nk)/dpthmx
            uer(nk)=vmflcl*dp(nk)/dpthmx
          elseif(nk.le.kpbl)then
            uer(nk)=vmflcl*dp(nk)/dpthmx
            umf(nk)=umf(nk-1)+uer(nk)
          else
            umf(nk)=vmflcl
            uer(nk)=0.
          endif
          tu(nk)=tmix+(z0(nk)-zmix)*gdry
          qu(nk)=qmix
          wu(nk)=wlcl
        else
          tu(nk)=0.
          qu(nk)=0.
          umf(nk)=0.
          wu(nk)=0.
          uer(nk)=0.
        endif
        udr(nk)=0.
        qdt(nk)=0.
        qliq(nk)=0.
        qice(nk)=0.
        qlqout(nk)=0.
        qicout(nk)=0.
        pptliq(nk)=0.
        pptice(nk)=0.
        detlq(nk)=0.
        detic(nk)=0.
        ratio2(nk)=0.
      call envirtht(p0(nk),t0(nk),q0(nk),thetee(nk),aliq,bliq,cliq,dliq)
         eqfrc(nk)=1.0
      enddo
!     
      ltop1=ltop+1
      ltopm1=ltop-1
!     
!...define variables above cloud top...
!     
      do nk=ltop1,kx
        umf(nk)=0.
        udr(nk)=0.
        uer(nk)=0.
        qdt(nk)=0.
        qliq(nk)=0.
        qice(nk)=0.
        qlqout(nk)=0.
        qicout(nk)=0.
        detlq(nk)=0.
        detic(nk)=0.
        pptliq(nk)=0.
        pptice(nk)=0.
        if(nk.gt.ltop1)then
          tu(nk)=0.
          qu(nk)=0.
          wu(nk)=0.
        endif
        thta0(nk)=0.
        thtau(nk)=0.
        ems(nk)=0.
        emsd(nk)=0.
        tg(nk)=t0(nk)
        qg(nk)=q0(nk)
        qlg(nk)=0.
        qig(nk)=0.
        qrg(nk)=0.
        qsg(nk)=0.
        omg(nk)=0.
      enddo
        omg(kx+1)=0.
        do nk=1,ltop
          ems(nk)=dp(nk)*dxsq/g
          emsd(nk)=1./ems(nk)
!     
!...initialize some variables to be used later in the vert advection sch
!     
          exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*qdt(nk)))
          thtau(nk)=tu(nk)*exn(nk)
          exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*q0(nk)))
          thta0(nk)=t0(nk)*exn(nk)
          ddilfrc(nk) = 1./dilfrc(nk)
          omg(nk)=0.
        enddo
!     if (xtime.lt.10.)then
!      write(98,1025)klcl,zlcl,dtlcl,ltop,p0(ltop),iflag,
!    * tmix-t00,pmix,qmix,abe
!      write(98,1030)p0(let)/100.,p0(ltop)/100.,vmflcl,plcl/100.,
!    * wlcl,cldhgt
!     endif
!     
!...compute convective time scale(timec). the mean wind at the lcl
!...and midtroposphere is used.
!     
        wspd(klcl)=sqrt(u0(klcl)*u0(klcl)+v0(klcl)*v0(klcl))
        wspd(l5)=sqrt(u0(l5)*u0(l5)+v0(l5)*v0(l5))
        wspd(ltop)=sqrt(u0(ltop)*u0(ltop)+v0(ltop)*v0(ltop))
        vconv=.5*(wspd(klcl)+wspd(l5))
!...for eta model, dx is a function of location...
!       timec=dx(i,j)/vconv
        timec=dx/vconv
        tadvec=timec
        timec=amax1(1800.,timec)
        timec=amin1(3600.,timec)
        if(ishall.eq.1)timec=2400.
        nic=nint(timec/dt)
        timec=float(nic)*dt
!     
!...compute wind shear and precipitation efficiency.
!     
        if(wspd(ltop).gt.wspd(klcl))then
          shsign=1.
        else
          shsign=-1.
        endif
        vws=(u0(ltop)-u0(klcl))*(u0(ltop)-u0(klcl))+(v0(ltop)-v0(klcl))*
     1      (v0(ltop)-v0(klcl))
        vws=1.e3*shsign*sqrt(vws)/(z0(ltop)-z0(lcl))
        pef=1.591+vws*(-.639+vws*(9.53e-2-vws*4.96e-3))
        pef=amax1(pef,.2)
        pef=amin1(pef,.9)
!     
!...precipitation efficiency is a function of the height of cloud base.
!     
        cbh=(zlcl-z0(1))*3.281e-3
        if(cbh.lt.3.)then
          rcbh=.02
        else
          rcbh=.96729352+cbh*(-.70034167+cbh*(.162179896+cbh*(- 
     1         1.2569798e-2+cbh*(4.2772e-4-cbh*5.44e-6))))
        endif
        if(cbh.gt.25)rcbh=2.4
        pefcbh=1./(1.+rcbh)
        pefcbh=amin1(pefcbh,.9)
!     
!... mean pef. is used to compute rainfall.
!     
        peff=.5*(pef+pefcbh)
        peff2 = peff                                ! jsk mods
       if(iprnt)then  
         write(98,1035)pef,pefcbh,lc,let,wkl,vws
!       call flush(98)   
       endif     
!        write(98,1035)pef,pefcbh,lc,let,wkl,vws
!*****************************************************************
!                                                                *
!                  compute downdraft properties                  *
!                                                                *
!*****************************************************************
!     
!     
       tder=0.
      devap:if(ishall.eq.1)then
         lfs = 1
       else
!
!...start downdraft about 150 mb above cloud base...
!
!        kstart=max0(kpbl,klcl)
!        kstart=kpbl                                  ! changed 7/23/99
         kstart=kpbl+1                                ! changed 7/23/99
         do nk = kstart+1,kl
           dppp = p0(kstart)-p0(nk)
!          if(dppp.gt.200.e2)then
           if(dppp.gt.150.e2)then
             klfs = nk
             exit 
           endif
         enddo
         klfs = min0(klfs,let-1)
         lfs = klfs
!
!...if lfs is not at least 50 mb above cloud base (implying that the 
!...level of equil temp, let, is just above cloud base) do not allow a
!...downdraft...
!
        if((p0(kstart)-p0(lfs)).gt.50.e2)then
          theted(lfs) = thetee(lfs)
          qd(lfs) = q0(lfs)
!
!...call tpmix2dd to find wet-bulb temp, qv...
!
          call tpmix2dd(p0(lfs),theted(lfs),tz(lfs),qss)
          thtad(lfs)=tz(lfs)*(p00/p0(lfs))**(0.2854*(1.-0.28*qss))
!     
!...take a first guess at the initial downdraft mass flux...
!     
          tvd(lfs)=tz(lfs)*(1.+0.608*qss)
          rdd=p0(lfs)/(r*tvd(lfs))
          a1=(1.-peff)*au0
          dmf(lfs)=-a1*rdd
          der(lfs)=dmf(lfs)
          ddr(lfs)=0.
          rhbar = rh(lfs)*dp(lfs)
          dptt = dp(lfs)
          do nd = lfs-1,kstart,-1
            nd1 = nd+1
            der(nd)=der(lfs)*ems(nd)/ems(lfs)
            ddr(nd)=0.
            dmf(nd)=dmf(nd1)+der(nd)
            theted(nd)=(theted(nd1)*dmf(nd1)+thetee(nd)*der(nd))/dmf(nd)
            qd(nd)=(qd(nd1)*dmf(nd1)+q0(nd)*der(nd))/dmf(nd)    
            dptt = dptt+dp(nd)
            rhbar = rhbar+rh(nd)*dp(nd)
          enddo
          rhbar = rhbar/dptt
          dmffrc = 2.*(1.-rhbar)
          dpdd = 0.
!...calculate melting effect
!... first, compute total frozen precipitation generated...
!
          pptmlt = 0.
          do nk = klcl,ltop
            pptmlt = pptmlt+pptice(nk)
          enddo
          if(lc.lt.ml)then
!...for now, calculate melting effect as if dmf = -umf at klcl, i.e., as
!...if dmffrc=1.  otherwise, for small dmffrc, dtmelt gets too large!
!...12/14/98 jsk...
            dtmelt = rlf*pptmlt/(cp*umf(klcl))
          else
            dtmelt = 0.
          endif
          ldt = min0(lfs-1,kstart-1)
!
          call tpmix2dd(p0(kstart),theted(kstart),tz(kstart),qss)
!
          tz(kstart) = tz(kstart)-dtmelt
          es=aliq*exp((bliq*tz(kstart)-cliq)/(tz(kstart)-dliq))
          qss=0.622*es/(p0(kstart)-es)
          theted(kstart)=tz(kstart)*(1.e5/p0(kstart))**(0.2854*
     1          (1.-0.28*qss))*
     2          exp((3374.6525/tz(kstart)-2.5403)*qss*(1.+0.81*qss))
!....  
          ldt = min0(lfs-1,kstart-1)
          do nd = ldt,1,-1
            dpdd = dpdd+dp(nd)
            theted(nd) = theted(kstart)
            qd(nd)     = qd(kstart)       
!
!...call tpmix2dd to find wet bulb temp, saturation mixing ratio...
!
            call tpmix2dd(p0(nd),theted(nd),tz(nd),qss)
            qsd(nd) = qss
!
!...specify rh decrease of 20%/km in downdraft...
!
            rhh = 1.-0.2/1000.*(z0(kstart)-z0(nd))
!
!...adjust downdraft temp, q to specified rh:
!
            if(rhh.lt.1.)then
              dssdt=(cliq-bliq*dliq)/((tz(nd)-dliq)*(tz(nd)-dliq))              
              rl=xlv0-xlv1*tz(nd)
              dtmp=rl*qss*(1.-rhh)/(cp+rl*rhh*qss*dssdt)
              t1rh=tz(nd)+dtmp
              es=rhh*aliq*exp((bliq*t1rh-cliq)/(t1rh-dliq))
              qsrh=0.622*es/(p0(nd)-es)
!
!...check to see if mixing ratio at specified rh is less than actual
!...mixing ratio...if so, adjust to give zero evaporation...
!
              if(qsrh.lt.qd(nd))then
                qsrh=qd(nd)
                t1rh=tz(nd)+(qss-qsrh)*rl/cp
              endif
              tz(nd)=t1rh
              qss=qsrh
              qsd(nd) = qss
            endif         
            tvd(nd) = tz(nd)*(1.+0.608*qsd(nd))
            if(tvd(nd).gt.tv0(nd).or.nd.eq.1)then
              ldb=nd
              exit
            endif
          enddo
          if((p0(ldb)-p0(lfs)) .gt. 50.e2)then   
                                         ! minimum downdraft depth! 
            do nd=ldt,ldb,-1
              nd1 = nd+1
              ddr(nd) = -dmf(kstart)*dp(nd)/dpdd
              der(nd) = 0.
              dmf(nd) = dmf(nd1)+ddr(nd)
              tder=tder+(qsd(nd)-qd(nd))*ddr(nd)
              qd(nd)=qsd(nd)
              thtad(nd)=tz(nd)*(p00/p0(nd))**(0.2854*(1.-0.28*qd(nd)))
            enddo
          endif
        endif
      endif devap
!
!...if downdraft does not evaporate any water for specified relative
!...humidity, no downdraft is allowed...
!
      d_mf:   if(tder.lt.1.)then
!           write(98,3004)i,j 
!3004       format(' ','no downdraft!;  i=',i3,2x,'j=',i3,'ishall =',i2)
          pptflx=trppt
          cpr=trppt
          tder=0.
          cndtnf=0.
          updinc=1.
          ldb=lfs
          do ndk=1,ltop
            dmf(ndk)=0.
            der(ndk)=0.
            ddr(ndk)=0.
            thtad(ndk)=0.
            wd(ndk)=0.
            tz(ndk)=0.
            qd(ndk)=0.
          enddo
          aincm2=100.
        else 
          ddinc = -dmffrc*umf(klcl)/dmf(kstart)
          updinc=1.
          if(tder*ddinc.gt.trppt)then
            ddinc = trppt/tder
          endif
          tder = tder*ddinc
          do nk=ldb,lfs
            dmf(nk)=dmf(nk)*ddinc
            der(nk)=der(nk)*ddinc
            ddr(nk)=ddr(nk)*ddinc
          enddo
         cpr=trppt
         pptflx = trppt-tder
         peff=pptflx/trppt
         if(iprnt)then
           write(98,*)'precip efficiency =',peff
!          call flush(98)   
         endif
!
!
!...adjust updraft mass flux, mass detrainment rate, and liquid water an
!   detrainment rates to be consistent with the transfer of the estimate
!   from the updraft to the downdraft at the lfs...
!     
!         do nk=lc,lfs
!           umf(nk)=umf(nk)*updinc
!           udr(nk)=udr(nk)*updinc
!           uer(nk)=uer(nk)*updinc
!           pptliq(nk)=pptliq(nk)*updinc
!           pptice(nk)=pptice(nk)*updinc
!           detlq(nk)=detlq(nk)*updinc
!           detic(nk)=detic(nk)*updinc
!         enddo
!     
!...zero out the arrays for downdraft data at levels above and below the
!...downdraft...
!     
         if(ldb.gt.1)then
           do nk=1,ldb-1
             dmf(nk)=0.
             der(nk)=0.
             ddr(nk)=0.
             wd(nk)=0.
             tz(nk)=0.
             qd(nk)=0.
             thtad(nk)=0.
           enddo
         endif
         do nk=lfs+1,kx
           dmf(nk)=0.
           der(nk)=0.
           ddr(nk)=0.
           wd(nk)=0.
           tz(nk)=0.
           qd(nk)=0.
           thtad(nk)=0.
         enddo
         do nk=ldt+1,lfs-1
           tz(nk)=0.
           qd(nk)=0.
           thtad(nk)=0.
         enddo
       endif d_mf
!
!...set limits on the updraft and downdraft mass fluxes so that the infl
!   into convective drafts from a given layer is no more than is availab
!   in that layer initially...
!     
       aincmx=1000.
       lmax=max0(klcl,lfs)
       do nk=lc,lmax
         if((uer(nk)-der(nk)).gt.1.e-3)then
           aincm1=ems(nk)/((uer(nk)-der(nk))*timec)
           aincmx=amin1(aincmx,aincm1)
         endif
       enddo
       ainc=1.
       if(aincmx.lt.ainc)ainc=aincmx
!     
!...save the relevent variables for a unit updraft and downdraft...they will 
!...be iteratively adjusted by the factor ainc to satisfy the stabilization
!...closure...
!     
       tder2=tder
       pptfl2=pptflx
       do nk=1,ltop
         detlq2(nk)=detlq(nk)
         detic2(nk)=detic(nk)
         udr2(nk)=udr(nk)
         uer2(nk)=uer(nk)
         ddr2(nk)=ddr(nk)
         der2(nk)=der(nk)
         umf2(nk)=umf(nk)
         dmf2(nk)=dmf(nk)
       enddo
       fabe=1.
       stab=0.95
       noitr=0
       istop=0
!
        if(ishall.eq.1)then                
                                ! first for shallow convection
!
! no iteration for shallow convection; if turbulent kinetic energy (tke) is available
! from a turbulence parameterization, scale cloud-base updraft mass flux as a function
! of tke, but for now, just specify shallow-cloud mass flux using tkemax = 5...
!
!...find the maximum tke value between lc and klcl...
!         tkemax = 0.
          tkemax = 5.
!          do 173 k = lc,klcl
!            nk = kx-k+1
!            tkemax = amax1(tkemax,q2(i,j,nk))
! 173      continue
!          tkemax = amin1(tkemax,10.)
!          tkemax = amax1(tkemax,5.)
!c         tkemax = 10.
!c...3_24_99...dpmin was changed for shallow convection so that it is the
!c...          the same as for deep convection (5.e3).  since this doubles
!c...          (roughly) the value of dpthmx, add a factor of 0.5 to calcu-
!c...          lation of evac...
!c         evac  = tkemax*0.1
          evac  = 0.5*tkemax*0.1
!         ainc = 0.1*dpthmx*dxij*dxij/(vmflcl*g*timec)
!          ainc = evac*dpthmx*dx(i,j)*dx(i,j)/(vmflcl*g*timec)
          ainc = evac*dpthmx*dxsq/(vmflcl*g*timec)
          tder=tder2*ainc
          pptflx=pptfl2*ainc
          do nk=1,ltop
            umf(nk)=umf2(nk)*ainc
            dmf(nk)=dmf2(nk)*ainc
            detlq(nk)=detlq2(nk)*ainc
            detic(nk)=detic2(nk)*ainc
            udr(nk)=udr2(nk)*ainc
            uer(nk)=uer2(nk)*ainc
            der(nk)=der2(nk)*ainc
            ddr(nk)=ddr2(nk)*ainc
          enddo
        endif                    ! otherwise for deep convection
! use iterative procedure to find mass fluxes...
      iter:     do ncount=1,10
!     
!*****************************************************************
!                                                                *
!           compute properties for compensational subsidence     *
!                                                                *
!*****************************************************************
!     
!...determine omega value necessary at top and bottom of each layer to
!...satisfy mass continuity...
!     
            dtt=timec
            do nk=1,ltop
              domgdp(nk)=-(uer(nk)-der(nk)-udr(nk)-ddr(nk))*emsd(nk)
              if(nk.gt.1)then
                omg(nk)=omg(nk-1)-dp(nk-1)*domgdp(nk-1)
                absomg = abs(omg(nk))
                absomgtc = absomg*timec
                frdp = 0.75*dp(nk-1)
                if(absomgtc.gt.frdp)then
                  dtt1 = frdp/absomg
                  dtt=amin1(dtt,dtt1)
                endif
              endif
            enddo
            do nk=1,ltop
              thpa(nk)=thta0(nk)
              qpa(nk)=q0(nk)
              nstep=nint(timec/dtt+1)
              dtime=timec/float(nstep)
              fxm(nk)=omg(nk)*dxsq/g
            enddo
!     
!...do an upstream/forward-in-time advection of theta, qv...
!     
        do ntc=1,nstep
!     
!...assign theta and q values at the top and bottom of each layer based
!...sign of omega...
!     
            do  nk=1,ltop
              thfxin(nk)=0.
              thfxout(nk)=0.
              qfxin(nk)=0.
              qfxout(nk)=0.
            enddo
            do nk=2,ltop
              if(omg(nk).le.0.)then
                thfxin(nk)=-fxm(nk)*thpa(nk-1)
                qfxin(nk)=-fxm(nk)*qpa(nk-1)
                thfxout(nk-1)=thfxout(nk-1)+thfxin(nk)
                qfxout(nk-1)=qfxout(nk-1)+qfxin(nk)
              else
                thfxout(nk)=fxm(nk)*thpa(nk)
                qfxout(nk)=fxm(nk)*qpa(nk)
                thfxin(nk-1)=thfxin(nk-1)+thfxout(nk)
                qfxin(nk-1)=qfxin(nk-1)+qfxout(nk)
              endif
            enddo
!     
!...update the theta and qv values at each level...
!     
            do nk=1,ltop
              thpa(nk)=thpa(nk)+(thfxin(nk)+udr(nk)*thtau(nk)+ddr(nk)* 
     1               thtad(nk)-thfxout(nk)-(uer(nk)-der(nk))*thta0(nk))*
     2               dtime*emsd(nk)
              qpa(nk)=qpa(nk)+(qfxin(nk)+udr(nk)*qdt(nk)+ddr(nk)*qd(nk)-
     1               qfxout(nk)-(uer(nk)-der(nk))*q0(nk))*dtime*emsd(nk)
            enddo   
          enddo   
          do nk=1,ltop
            thtag(nk)=thpa(nk)
            qg(nk)=qpa(nk)
          enddo
!     
!...check to see if mixing ratio dips below zero anywhere;  if so, borro
!...moisture from adjacent layers to bring it back up above zero...
!     
        do nk=1,ltop
          if(qg(nk).lt.0.)then
            if(nk.eq.1)then                             ! jsk mods
!              print *,' problem with kf scheme:  ' ! jsk mods
!              print *,'qg = 0 at the surface!!!!!!!'    ! jsk mods
            endif                                       ! jsk mods
            nk1=nk+1
            if(nk.eq.ltop)then
              nk1=klcl
            endif
            tma=qg(nk1)*ems(nk1)
            tmb=qg(nk-1)*ems(nk-1)
            tmm=(qg(nk)-1.e-9)*ems(nk  )
            bcoeff=-tmm/((tma*tma)/tmb+tmb)
            acoeff=bcoeff*tma/tmb
            tmb=tmb*(1.-bcoeff)
            tma=tma*(1.-acoeff)
            if(nk.eq.ltop)then
              qvdiff=(qg(nk1)-tma*emsd(nk1))*100./qg(nk1)
!              if(abs(qvdiff).gt.1.)then
!         print *,'!!!warning!!! cloud base water vapor changes by ',
!     1          qvdiff, 
!     2          '% when moisture is borrowed to prevent negative ',
!     3               'values in kain-fritsch'
!              endif
            endif
            qg(nk)=1.e-9
            qg(nk1)=tma*emsd(nk1)
            qg(nk-1)=tmb*emsd(nk-1)
          endif
        enddo
        topomg=(udr(ltop)-uer(ltop))*dp(ltop)*emsd(ltop)
        if(abs(topomg-omg(ltop)).gt.1.e-3)then
!       write(99,*)'error:  mass does not balance in kf scheme;        
!     1 topomg, omg =',topomg,omg(ltop)
          istop=1
          iprnt=.true.
          exit iter
        endif
!     
!...convert theta to t...
!     
        do nk=1,ltop
          exn(nk)=(p00/p0(nk))**(0.2854*(1.-0.28*qg(nk)))
          tg(nk)=thtag(nk)/exn(nk)
          tvg(nk)=tg(nk)*(1.+0.608*qg(nk))
        enddo
        if(ishall.eq.1)then
          exit iter
        endif
!     
!*******************************************************************
!                                                                  *
!     compute new cloud and change in available buoyant energy.    *
!                                                                  *
!*******************************************************************
!     
!...the following computations are similar to that for updraft
!     
!        thmix=0.
          tmix=0.
          qmix=0.
!
!...find the thermodynamic characteristics of the layer by
!...mass-weighting the characteristics of the individual model
!...layers...
!
          do nk=lc,kpbl
            tmix=tmix+dp(nk)*tg(nk)
            qmix=qmix+dp(nk)*qg(nk)  
          enddo
          tmix=tmix/dpthmx
          qmix=qmix/dpthmx
          es=aliq*exp((tmix*bliq-cliq)/(tmix-dliq))
          qss=0.622*es/(pmix-es)
!     
!...remove supersaturation for diagnostic purposes, if necessary...
!     
          if(qmix.gt.qss)then
            rl=xlv0-xlv1*tmix
            cpm=cp*(1.+0.887*qmix)
            dssdt=qss*(cliq-bliq*dliq)/((tmix-dliq)*(tmix-dliq))
            dq=(qmix-qss)/(1.+rl*dssdt/cpm)
            tmix=tmix+rl/cp*dq
            qmix=qmix-dq
            tlcl=tmix
          else
            qmix=amax1(qmix,0.)
            emix=qmix*pmix/(0.622+qmix)
            astrt=1.e-3
            binc=0.075
            a1=emix/aliq
            tp=(a1-astrt)/binc
            indlu=int(tp)+1
            value=(indlu-1)*binc+astrt
            aintrp=(a1-value)/binc
            tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
            tdpt=(cliq-dliq*tlog)/(bliq-tlog)
            tlcl=tdpt-(.212+1.571e-3*(tdpt-t00)-4.36e-4*(tmix-t00))
     1             *(tmix-tdpt)
            tlcl=amin1(tlcl,tmix)
          endif
          tvlcl=tlcl*(1.+0.608*qmix)
          zlcl = zmix+(tlcl-tmix)/gdry
          do nk = lc,kl
            klcl=nk
            if(zlcl.le.z0(nk))then
              exit 
            endif
          enddo
          k=klcl-1
          dlp=(zlcl-z0(k))/(z0(klcl)-z0(k))
!     
!...estimate environmental temperature and mixing ratio at the lcl...
!     
          tenv=tg(k)+(tg(klcl)-tg(k))*dlp
          qenv=qg(k)+(qg(klcl)-qg(k))*dlp
          tven=tenv*(1.+0.608*qenv)
          plcl=p0(k)+(p0(klcl)-p0(k))*dlp
          theteu(k)=tmix*(1.e5/pmix)**(0.2854*(1.-0.28*qmix))* 
     1            exp((3374.6525/tlcl-2.5403)*qmix*(1.+0.81*qmix))
!     
!...compute adjusted abe(abeg).
!     
          abeg=0.
          do nk=k,ltopm1
            nk1=nk+1
            theteu(nk1) = theteu(nk)
!
            call tpmix2dd(p0(nk1),theteu(nk1),tgu(nk1),qgu(nk1))
!
            tvqu(nk1)=tgu(nk1)*(1.+0.608*qgu(nk1)-qliq(nk1)-qice(nk1))
            if(nk.eq.k)then
              dzz=z0(klcl)-zlcl
              dilbe=((tvlcl+tvqu(nk1))/(tven+tvg(nk1))-1.)*dzz
            else
              dzz=dza(nk)
              dilbe=((tvqu(nk)+tvqu(nk1))/(tvg(nk)+tvg(nk1))-1.)*dzz
            endif
            if(dilbe.gt.0.)abeg=abeg+dilbe*g
!
!...dilute by entrainment by the rate as original updraft...
!
            call envirtht(p0(nk1),tg(nk1),qg(nk1),thteeg(nk1),
     1                aliq,bliq,cliq,dliq)
            theteu(nk1)=theteu(nk1)*ddilfrc(nk1)+thteeg(nk1)
     1                *(1.-ddilfrc(nk1))
          enddo
!     
!...assume at least 90% of cape (abe) is removed by convection during
!...the period timec...
!     
          if(noitr.eq.1)then
!         write(98,*)' '
!         write(98,*)'tau, i, j, =',ntsd,i,j
!         write(98,1060)fabe
!          goto 265
          exit iter
          endif
          dabe=amax1(abe-abeg,0.1*abe)
          fabe=abeg/abe
          if(fabe.gt.1. .and. ishall.eq.0)then
!          write(98,*)'updraft/downdraft couplet increases cape at this
!     *grid point; no convection allowed!'
            return  
          endif
          if(ncount.ne.1)then
            if(abs(ainc-aincold).lt.0.0001)then
              noitr=1
              ainc=aincold
              cycle iter
            endif
            dfda=(fabe-fabeold)/(ainc-aincold)
            if(dfda.gt.0.)then
              noitr=1
              ainc=aincold
              cycle iter
            endif
          endif
          aincold=ainc
          fabeold=fabe
          if(ainc/aincmx.gt.0.999.and.fabe.gt.1.05-stab)then
!           write(98,*)' '
!           write(98,*)'tau, i, j, =',ntsd,i,j
!           write(98,1055)fabe
!            goto 265
            exit
          endif
          if((fabe.le.1.05-stab.and.fabe.ge.0.95-stab) .or.
     1           ncount.eq.10)then
            exit iter
          else
            if(ncount.gt.10)then
!             write(98,*)' '
!             write(98,*)'tau, i, j, =',ntsd,i,j
!             write(98,1060)fabe
!             goto 265
              exit
            endif
!     
!...if more than 10% of the original cape remains, increase the convecti
!...mass flux by the factor ainc:
!     
            if(fabe.eq.0.)then
              ainc=ainc*0.5
            else
              if(dabe.lt.1.e-4)then
                noitr=1
                ainc=aincold
                cycle iter
              else
                ainc=ainc*stab*abe/dabe
              endif
            endif
!           ainc=amin1(aincmx,ainc)
            ainc=amin1(aincmx,ainc)
!...if ainc becomes very small, effects of convection ! jsk mods
!...will be minimal so just ignore it...              ! jsk mods
            if(ainc.lt.0.05)then
              return                          ! jsk mods
            endif
!            ainc=amax1(ainc,0.05)                        ! jsk mods
            tder=tder2*ainc
            pptflx=pptfl2*ainc
!           if (xtime.lt.10.)then
!           write(98,1080)lfs,ldb,ldt,timec,tadvec,nstep,ncount,
!          *              fabeold,aincold 
!           endif
            do nk=1,ltop
              umf(nk)=umf2(nk)*ainc
              dmf(nk)=dmf2(nk)*ainc
              detlq(nk)=detlq2(nk)*ainc
              detic(nk)=detic2(nk)*ainc
              udr(nk)=udr2(nk)*ainc
              uer(nk)=uer2(nk)*ainc
              der(nk)=der2(nk)*ainc
              ddr(nk)=ddr2(nk)*ainc
            enddo
!     
!...go back up for another iteration...
!     
          endif
        enddo iter
!     
!...compute hydrometeor tendencies as is done for t, qv...
!     
!...frc2 is the fraction of total condensate      !  ppt fb mods
!...generated that goes into precipitiation       !  ppt fb mods
!
!  redistribute hydormeteors according to the final mass-flux values:
!
        if(cpr.gt.0.)then 
          frc2=pptflx/(cpr*ainc)                    !  ppt fb mods
        else
           frc2=0.
        endif
        do nk=1,ltop
          qlpa(nk)=ql0(nk)
          qipa(nk)=qi0(nk)
          qrpa(nk)=qr0(nk)
          qspa(nk)=qs0(nk)
          rainfb(nk)=pptliq(nk)*ainc*fbfrc*frc2   !  ppt fb mods
          snowfb(nk)=pptice(nk)*ainc*fbfrc*frc2   !  ppt fb mods
        enddo
        do ntc=1,nstep
!     
!...assign hydrometeors concentrations at the top and bottom of each lay
!...based on the sign of omega...
!     
          do nk=1,ltop
            qlfxin(nk)=0.
            qlfxout(nk)=0.
            qifxin(nk)=0.
            qifxout(nk)=0.
            qrfxin(nk)=0.
            qrfxout(nk)=0.
            qsfxin(nk)=0.
            qsfxout(nk)=0.
          enddo   
          do nk=2,ltop
            if(omg(nk).le.0.)then
              qlfxin(nk)=-fxm(nk)*qlpa(nk-1)
              qifxin(nk)=-fxm(nk)*qipa(nk-1)
              qrfxin(nk)=-fxm(nk)*qrpa(nk-1)
              qsfxin(nk)=-fxm(nk)*qspa(nk-1)
              qlfxout(nk-1)=qlfxout(nk-1)+qlfxin(nk)
              qifxout(nk-1)=qifxout(nk-1)+qifxin(nk)
              qrfxout(nk-1)=qrfxout(nk-1)+qrfxin(nk)
              qsfxout(nk-1)=qsfxout(nk-1)+qsfxin(nk)
            else
              qlfxout(nk)=fxm(nk)*qlpa(nk)
              qifxout(nk)=fxm(nk)*qipa(nk)
              qrfxout(nk)=fxm(nk)*qrpa(nk)
              qsfxout(nk)=fxm(nk)*qspa(nk)
              qlfxin(nk-1)=qlfxin(nk-1)+qlfxout(nk)
              qifxin(nk-1)=qifxin(nk-1)+qifxout(nk)
              qrfxin(nk-1)=qrfxin(nk-1)+qrfxout(nk)
              qsfxin(nk-1)=qsfxin(nk-1)+qsfxout(nk)
            endif
          enddo   
!     
!...update the hydrometeor concentration values at each level...
!     
          do nk=1,ltop
            qlpa(nk)=qlpa(nk)+(qlfxin(nk)+detlq(nk)-qlfxout(nk))
     1              *dtime*emsd(nk)
            qipa(nk)=qipa(nk)+(qifxin(nk)+detic(nk)-qifxout(nk))
     1              *dtime*emsd(nk)
            qrpa(nk)=qrpa(nk)+(qrfxin(nk)-qrfxout(nk)+rainfb(nk))
     1              *dtime*emsd(nk)         !  ppt fb mods
            qspa(nk)=qspa(nk)+(qsfxin(nk)-qsfxout(nk)+snowfb(nk))
     1              *dtime*emsd(nk)         !  ppt fb mods
          enddo     
        enddo
        do nk=1,ltop
          qlg(nk)=qlpa(nk)
          qig(nk)=qipa(nk)
          qrg(nk)=qrpa(nk)
          qsg(nk)=qspa(nk)
        enddo   
!
!...clean things up, calculate convective feedback tendencies for this
!...grid point...
!     
!     if (xtime.lt.10.)then
!     write(98,1080)lfs,ldb,ldt,timec,tadvec,nstep,ncount,fabe,ainc 
!     endif
       if(iprnt)then  
         write(98,1080)lfs,ldb,ldt,timec,tadvec,nstep,ncount,fabe,ainc
!        call flush(98)   
       endif  
!     
!...send final parameterized values to output files...
!     
!297   if(iprnt)then 
       if(iprnt)then 
!    if(i.eq.16 .and. j.eq.41)then
!      if(istop.eq.1)then
         write(98,*)
!        write(98,*)'at t(h), i, j =',float(ntsd)*72./3600.,i,j
         write(98,*)'p(lc), dtp, wkl, wklcl =',p0(lc)/100.,       
     1               tlcl+dtlcl+dtrh-tenv,wkl,wklcl
         write(98,*)'tlcl, dtlcl, dtrh, tenv =',tlcl,dtlcl,       
     1                dtrh,tenv   
         write(98,1025)klcl,zlcl,dtlcl,ltop,p0(ltop),iflag,       
     1   tmix-t00,pmix,qmix,abe
         write(98,1030)p0(let)/100.,p0(ltop)/100.,vmflcl,plcl/100., 
     1   wlcl,cldhgt(lc)
         write(98,1035)pef,pefcbh,lc,let,wkl,vws 
         write(98,*)'precip efficiency =',peff 
      write(98,1080)lfs,ldb,ldt,timec,tadvec,nstep,ncount,fabe,ainc
!      endif
!!!!! here !!!!!!!
        write(98,1070)'  p  ','   dp ',' dt k/d ',' dr k/d ','   omg  ',
     1 ' domgdp ','   umf  ','   uer  ','   udr  ','   dmf  ','   der  '
     2 ,'   ddr  ','   ems  ','    w0  ','  detlq ',' detic '
           write(98,*)'just before do 300...'
!          call flush(98)
           do nk=1,ltop
             k=ltop-nk+1
             dtt=(tg(k)-t0(k))*86400./timec
             rl=xlv0-xlv1*tg(k)
             dr=-(qg(k)-q0(k))*rl*86400./(timec*cp)
             udfrc=udr(k)*timec*emsd(k)
             uefrc=uer(k)*timec*emsd(k)
             ddfrc=ddr(k)*timec*emsd(k)
             defrc=-der(k)*timec*emsd(k)
       write(98,1075)p0(k)/100.,dp(k)/100.,dtt,dr,omg(k),domgdp(k)*1.e4,
     1    umf(k)/1.e6,uefrc,udfrc,dmf(k)/1.e6,defrc,ddfrc,ems(k)/1.e11, 
     2    w0avg1d(k)*1.e2,detlq(k)*timec*emsd(k)*1.e3,detic(k)*     
     3    timec*emsd(k)*1.e3
           enddo
           write(98,1085)'k','p','z','t0','tg','dt','tu','td','q0','qg',
     1            'dq','qu','qd','qlg','qig','qrg','qsg','rh0','rhg'
           do nk=1,kl
             k=kx-nk+1
             dtt=tg(k)-t0(k)
             tuc=tu(k)-t00
             if(k.lt.lc.or.k.gt.ltop)tuc=0.
             tdc=tz(k)-t00
             if((k.lt.ldb.or.k.gt.ldt).and.k.ne.lfs)tdc=0.
             if(t0(k).lt.t00)then
               es=aliq*exp((bliq*tg(k)-cliq)/(tg(k)-dliq))
             else
               es=aliq*exp((bliq*tg(k)-cliq)/(tg(k)-dliq))
             endif  
             qgs=es*0.622/(p0(k)-es)
             rh0=q0(k)/qes(k)
             rhg=qg(k)/qgs
           write(98,1090)k,p0(k)/100.,z0(k),t0(k)-t00,tg(k)-t00,dtt,tuc,
     1       tdc,q0(k)*1000.,qg(k)*1000.,(qg(k)-q0(k))*1000.,qu(k)*     
     2       1000.,qd(k)*1000.,qlg(k)*1000.,qig(k)*1000.,qrg(k)*1000.,  
     3       qsg(k)*1000.,rh0,rhg
           enddo
!     
!...if calculations above show an error in the mass budget, print out a
!...to be used later for diagnostic purposes, then abort run...
!     
!         if(istop.eq.1 .or. ishall.eq.1)then

!         if(ishall.ne.1)then
!            write(98,4421)i,j,iyr,imo,idy,ihr,imn
!           write(98)i,j,iyr,imo,idy,ihr,imn,kl
! 4421       format(7i4)
!            write(98,4422)kl
! 4422       format(i6) 
            do 310 nk = 1,kl
              k = kl - nk + 1
              write(98,4455) p0(k)/100.,t0(k)-273.16,q0(k)*1000.,   
     1                 u0(k),v0(k),w0avg1d(k),dp(k),tke(k)
!             write(98) p0,t0,q0,u0,v0,w0,dp,tke
!           write(98,1115)z0(k),p0(k)/100.,t0(k)-273.16,q0(k)*1000.,
!    *               u0(k),v0(k),dp(k)/100.,w0avg(i,j,k)
310        continue
            if(istop.eq.1)then
            endif
!         endif
4455    format(8f11.3) 
       endif
        cndtnf=(1.-eqfrc(lfs))*(qliq(lfs)+qice(lfs))*dmf(lfs)
        raincv(i,j)=dt*pptflx*(1.-fbfrc)/dxsq     !  ppt fb mods
!        raincv(i,j)=.1*.5*dt*pptflx/dxsq               !  ppt fb mods
!         rnc=0.1*timec*pptflx/dxsq
        rnc=raincv(i,j)*nic
       if(ishall.eq.0.and.iprnt)write (98,909)i,j,rnc

!     write(98,1095)cpr*ainc,tder+pptflx+cndtnf
!     
!  evaluate moisture budget...
!     

        qinit=0.
        qfnl=0.
        dpt=0.
        do 315 nk=1,ltop
          dpt=dpt+dp(nk)
          qinit=qinit+q0(nk)*ems(nk)
          qfnl=qfnl+qg(nk)*ems(nk)
          qfnl=qfnl+(qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)
315     continue
        qfnl=qfnl+pptflx*timec*(1.-fbfrc)       !  ppt fb mods
!        qfnl=qfnl+pptflx*timec                 !  ppt fb mods
        err2=(qfnl-qinit)*100./qinit
       if(iprnt)write(98,1110)qinit,qfnl,err2
      if(abs(err2).gt.0.05 .and. istop.eq.0)then 
!       write(99,*)'!!!!!!!! moisture budget error in kfpara !!!'
!       write(99,1110)qinit,qfnl,err2
        iprnt=.true.
        istop=1
            write(98,4422)kl
4422    format(i6)
            do 311 nk = 1,kl
              k = kl - nk + 1
!             write(99,4455) p0(k)/100.,t0(k)-273.16,q0(k)*1000.,   
!     1                u0(k),v0(k),w0avg1d(k),dp(k)
!             write(98) p0,t0,q0,u0,v0,w0,dp,tke
!           write(98,1115)p0(k)/100.,t0(k)-273.16,q0(k)*1000.,      
!     1              u0(k),v0(k),w0avg1d(k),dp(k)/100.,tke(k)
            write(98,4456)p0(k)/100.,t0(k)-273.16,q0(k)*1000.,     
     1              u0(k),v0(k),w0avg1d(k),dp(k)/100.,tke(k)
311        continue
!           call flush(98)

!        goto 297
!         stop 'qverr'
      endif
1115    format (2x,f7.2,2x,f5.1,2x,f6.3,2(2x,f5.1),2x,f7.2,2x,f7.4)
4456    format(8f12.3)
        if(pptflx.gt.0.)then
          relerr=err2*qinit/(pptflx*timec)
        else
          relerr=0.
        endif
      if(iprnt)then
        write(98,1120)relerr
        write(98,*)'tder, cpr, trppt =',tder,cpr*ainc,trppt*ainc
      endif
!     
!...feedback to resolvable scale tendencies.
!     
!...if the advective time period (tadvec) is less than specified minimum
!...timec, allow feedback to occur only during tadvec...
!     
        if(tadvec.lt.timec)nic=nint(tadvec/dt)
        nca(i,j)=float(nic)
        if(ishall.eq.1)then
          timec = 2400.
          nca(i,j) = float(ntst)
          nshall = nshall+1
        endif 
        do k=1,kx
#ifdef CLD1
!         if(imoist(inest).ne.2)then
!
!...if hydrometeors are not allowed, they must be evaporated or sublimat
!...and fed back as vapor, along with associated changes in temperature.
!...note:  this will introduce changes in the convective temperature and
!...water vapor feedback tendencies and may lead to supersaturated value
!...of qg...
!
           rlc=xlv0-xlv1*tg(k)
           rls=xls0-xls1*tg(k)
           cpm=cp*(1.+0.887*qg(k))
           tg(k)=tg(k)-(rlc*(qlg(k)+qrg(k))+rls*(qig(k)+qsg(k)))/cpm
           qg(k)=qg(k)+(qlg(k)+qrg(k)+qig(k)+qsg(k))
           dqcdt(k)=0.
           dqidt(k)=0.
           dqrdt(k)=0.
           dqsdt(k)=0.
!         else
#else 
!
!...if ice phase is not allowed, melt all frozen hydrometeors...
!
          if(p_qi .lt. p_first_scalar .and. warm_rain)then

            cpm=cp*(1.+0.887*qg(k))
            tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
            dqcdt(k)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
            dqidt(k)=0.
            dqrdt(k)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
            dqsdt(k)=0.
          elseif(p_qi .lt. p_first_scalar .and. .not. warm_rain)then
!
!...if ice phase is allowed, but mixed phase is not, melt frozen hydrome
!...below the melting level, freeze liquid water above the melting level
!
            cpm=cp*(1.+0.887*qg(k))
            if(k.le.ml)then
              tg(k)=tg(k)-(qig(k)+qsg(k))*rlf/cpm
            elseif(k.gt.ml)then
              tg(k)=tg(k)+(qlg(k)+qrg(k))*rlf/cpm
            endif
            dqcdt(k)=(qlg(k)+qig(k)-ql0(k)-qi0(k))/timec
            dqidt(k)=0.
            dqrdt(k)=(qrg(k)+qsg(k)-qr0(k)-qs0(k))/timec
            dqsdt(k)=0.
          elseif(p_qi .ge. p_first_scalar) then
!
!...if mixed phase hydrometeors are allowed, feed back convective tenden
!...of hydrometeors directly...
!
            dqcdt(k)=(qlg(k)-ql0(k))/timec
            dqidt(k)=(qig(k)-qi0(k))/timec
            dqrdt(k)=(qrg(k)-qr0(k))/timec
            if (p_qs .ge. p_first_scalar) then
               dqsdt(k)=(qsg(k)-qs0(k))/timec
            else
               dqidt(k)=dqidt(k)+(qsg(k)-qs0(k))/timec
            endif
          else
!        print *,'this combination of imoist, iexice, iice not allowed!'
          endif
#endif
          dtdt(k)=(tg(k)-t0(k))/timec
          dqdt(k)=(qg(k)-q0(k))/timec
        enddo
        raincv(i,j)=dt*pptflx*(1.-fbfrc)/dxsq     !  ppt fb mods
!        raincv(i,j)=.1*.5*dt*pptflx/dxsq               !  ppt fb mods
!         rnc=0.1*timec*pptflx/dxsq
        rnc=raincv(i,j)*nic
909     format('at i, j =',i3,1x,i3,' convective rainfall =',
     1           f8.4,' mm')
!      write (98,909)i,j,rnc
!      write (6,909)i,j,rnc
!      write(98,*)'at ntsd =',ntsd,',no. of kf points activated =',
!     *            nccnt
!      call flush(98)
1000  format(' ',10a8)
1005  format(' ',f6.0,2x,f6.4,2x,f7.3,1x,f6.4,2x,4(f6.3,2x),2(f7.3,1x))
1010  format(' ',' vertical velocity is negative at ',f4.0,' mb')
1015  format(' ','all remaining mass detrains below ',f4.0,' mb')
1025  format(5x,' klcl=',i2,' zlcl=',f7.1,'m',                         
     1  ' dtlcl=',f5.2,' ltop=',i2,' p0(ltop)=',-2pf5.1,'mb frz lv=',  
     2  i2,' tmix=',0pf4.1,1x,'pmix=',-2pf6.1,' qmix=',3pf5.1,         
     3  ' cape=',0pf7.1)
1030   format(' ',' p0(let) = ',f6.1,' p0(ltop) = ',f6.1,' vmflcl =',  
     1 e12.3,' plcl =',f6.1,' wlcl =',f6.3,' cldhgt =',f8.1) 
1035   format(1x,'pef(ws)=',f4.2,'(cb)=',f4.2,'lc,let=',2i3,'wkl='   
     1 ,f6.3,'vws=',f5.2)
!1055  format('*** degree of stabilization =',f5.3,   
!      ', no more mass flux is allowed!')
!1060    format(' iteration does not converge to give the specified    
!     1 degree of stabilization!  fabe= ',f6.4) 
 1070 format (16a8) 
 1075 format (f8.2,3(f8.2),2(f8.3),f8.2,2f8.3,f8.2,6f8.3) 
 1080 format(2x,'lfs,ldb,ldt =',3i3,' timec, tadvec, nstep=',
     1        2(1x,f5.0),i3,'ncount, fabe, ainc=',i2,1x,f5.3,f6.2) 
 1085 format (a3,16a7,2a8) 
 1090 format (i3,f7.2,f7.0,10f7.2,4f7.3,2f8.3) 
 1095 format(' ','  ppt production rate= ',f10.0,
     1        ' total evap+ppt= ',f10.0)
1105   format(' ','net latent heat release =',e12.5,' actual heating =',
     1 e12.5,' j/kg-s, difference = ',f9.3,'%')
1110   format(' ','initial water =',e12.5,' final water =',e12.5,    
     1 ' total water change =',f8.2,'%')
!1115 format(2x,f6.0,2x,f7.2,2x,f5.1,2x,f6.3,2(2x,f5.1),2x,f7.2,2x,f7.4)
1120   format(' ','moisture error as function of total ppt =',f9.3,'%')
!
      end subroutine  kf_eta_para
