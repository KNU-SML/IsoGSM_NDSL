      subroutine rasv2(k,   dt, pdd
     *,                ncrnd, krmin, krmax, kfmax, frac, rasal, botop
     *,                revap, max_neg_bouy, alfint, cumfrc
!    *,                revap, max_neg_bouy, alfint, alfinq
     *,                cp,  alhl, alhf, grav, rkap, rannum
     *,                mct, kctop
!    *,                mct, kctop, sgh, sig
!
     *,                poi, qoi, uvi, kbl
!    *,                qli, qii
     *,                qli, qii, clw
     *,                rain, cln,    q1,  q2, pcu, flx, rcu
     *,                cdrag, dsfc
     *,                prs,   psj, wrkfun, calkbl, crtfun, updret)
!
!
!
!*********************************************************************
!*********************************************************************
!******************** relaxed arakawa-schubert ***********************
!************************ parameterization ***************************
!********************* plug compatible driver  **********************
!************************* 27 february 1998 ***************************
!*********************************************************************
!************************** developed by *****************************
!**************************              *****************************
!************************ shrinivas moorthi **************************
!************************                *****************************
!************************      emc/ncep  *****************************
!*********************************************************************
!*********************************************************************
!
!
!
      parameter (icm=100, one=1.0, daylen=86400.0, pfac=1.0/450.0)
!
      integer k, mct, ncrnd
      real rasal(k-1), uvi(k,2), poi(k), qoi(k)
!    *,    qli(k),     qii(k)
     *,    qli(k),     qii(k), clw(k), cln(k)
     *,    cdrag,      wfnc,   rcu(k,2), dsfc
      real dt, pdd,    frac,   max_neg_bouy
     *,    cp, alhl,   alhf,   grav, rkap
      integer kbl, kctop(mct+1), krmin, krmax, kfmax
!
      real q1(k), q2(k), flx(k), pcu(k), rain
      real prs(k+1), psj(k+1),  alfint(k), alfinq(k)
      real rannum(ncrnd)
!
      integer  ic(icm), idiag(100)
      logical  botop, revap, lprnt, wrkfun, calkbl, crtfun
     *,        dndrft, updret, cumfrc
!
      real cfac, savt, savq, savw(2), st, sq, sw(2), sprs, sprj, wr1
     *,    rhfacs, rhfacl, tem
!
      integer kcr, kfx, ncmx, nc, ktem, i, l, lct, lcb, ntrc
!
cfpp$ noconcur r
!
!
      rhfacs = 0.70
      rhfacl = 0.70
!
!     lend = max(1, len / 200)
!     if (lend .le. 100) then
!        do i=1,lend
!           idiag(i) = (i-1)*100 + 1
!        enddo
!     endif
!
      kcr   = min(mct,krmax)

      ktem  = min(mct,kfmax)
      kfx   = ktem - kcr
      ncmx  = kfx + ncrnd
!
      if (kfx .gt. 0) then
         if (botop) then 
            do nc=1,kfx
              ic(nc) = ktem + 1 - nc
            enddo
         else   
            do nc=kfx,1,-1
             ic(nc) = ktem + 1 - nc
            enddo
         endif  
      endif
!
      if (ncrnd .gt. 0) then
         do i=1,ncrnd
           irnd = (rannum(i)-0.0005)*(kcr-krmin+1)
           ic(kfx+i) = irnd + krmin
         enddo
      endif
!
      ntrc = 0
      if (cumfrc) ntrc = 2
      do 100 nc=1,ncmx
!
         ib = ic(nc)
!
         lct  = kctop(ib)
         lcb  = kctop(ib+1) - 1
!
!        if (me .eq. 0) then
!        write(6,*) ' lct=',lct,' lcb=',lcb,' ib=',ib
!    *,' kctop=',kctop(ib),kctop(ib+1),' prs=',prs(1,lcb),prs(1,lcb+1)
!        write(6,*) ' ib=',ib,' sgh=',sgh(ib),sgh(ib+1)
!    *,' kctop=',kctop(ib),kctop(ib+1),' sig=',sig(lcb),sig(lcb+1)
!        endif
!
         if (lct .lt. lcb) then
!          tem1    = 1.0 / (sgh(ib+1) - sgh(ib))
!          tem     = (sig(lcb+1) - sig(lcb)) * tem1
!          tem2    = tem
           wr1     = 1.0 / (prs(lcb+1) - prs(lct))
           tem     = (prs(lcb+1) - prs(lcb)) * wr1
!
           savt    = poi(lcb)
           savq    = qoi(lcb)
           savw(1) = uvi(lcb,1)
           savw(2) = uvi(lcb,2)
           sprs    = prs(lcb)
           sprj    = psj(lcb)
!
           poi(lcb)   = poi(lcb)   * tem
           qoi(lcb)   = qoi(lcb)   * tem
           uvi(lcb,1) = uvi(lcb,1) * tem
           uvi(lcb,2) = uvi(lcb,2) * tem
!          prs(lcb)   = prs(lct)
!          psj(lcb)   = psj(lct)
           do l=lct, lcb-1
!            tem  = (sig(l+1) - sig(l)) * tem1 
!            tem2 = tem2 + tem
             tem        = (prs(l+1) - prs(l))   * wr1
             poi(lcb)   = poi(lcb)   + poi(l)   * tem
             qoi(lcb)   = qoi(lcb)   + qoi(l)   * tem
             uvi(lcb,1) = uvi(lcb,1) + uvi(l,1) * tem
             uvi(lcb,2) = uvi(lcb,2) + uvi(l,2) * tem
!        if (me .eq. 0) then
!           print *,' tem=',tem,' poi=',poi(i,lcb),poi(i,l)
!    *,' l=',l,' lcb=',lcb,' lct=',lct
!        endif
           enddo
!          write(0,*) ' tem2=',tem2
           st       = poi(lcb)
           sq       = qoi(lcb)
           sw(1)    = uvi(lcb,1)
           sw(2)    = uvi(lcb,2)
           prs(lcb) = prs(lct)
           psj(lcb) = psj(lct)
         endif
!
         dndrft = pdd .gt. 0.0
!
!     write(0,*)' ib=',ib,' lct=',lct,' lcb=',lcb,' kdd=',kdd,' dndrft='
!    *,       dndrft
!
!     if (me .eq. 0) then
!     print *,' calling cloud type ib=',ib,' kbl=',kbl
!    *,' kpbl=',kpbl,' alfint=',alfint,' frac=',frac
!    *,' ntrc=',ntrc
!
!     if (lprnt) then
!     ix = 1
!     print *,' poi=',(poi(ix,l),l=1,k)
!     print *,' qoi=',(qoi(ix,l),l=1,k),' kbl=',kbl
!     endif
!
      do i=1,k-lcb+1
!        q1(ib)   = 0.0
!        q2(ib)   = 0.0
!        pcu(ib)  = 0.0
         flx(lcb) = 0.0
      enddo
      wfnc = 0.0
!
      do i=lcb,k-1
         alfinq(i) = alfint(i)
      enddo
      alfinq(k) = 0.0
!     alfinq = 1.0     ! 96/09/20
!
!     if (me .eq. 0) then
!     print *,' calling cloud type ib= ', ib,' dt=',dt
!     print *,' poi=',(poi(1,l),l=1,k)
!     print *,' qoi=',(qoi(1,l),l=1,k)
!     endif
!     print *,' alft=',alfint
!
      tla = -10.0
!
          call cloudv2(k, lcb, ntrc
     *,              rasal(lcb), frac,  max_neg_bouy
     *,              alfint, alfinq, rhfacl, rhfacs
     *,              cp,  alhl, alhf, rkap, grav
!
     *,              poi, qoi, uvi, prs, psj
     *,              qli, qii, kbl, dsfc
     *,              cdrag
     *,              idiag
     *,              q1, q2, rcu, pcu, flx
     *,              rain, revap, dt
     *,              wfnc, wrkfun, calkbl, crtfun, tla, dndrft, pdd)
!    *,              wfnc, wrkfun, calkbl, crtfun, tla, dndrft, updret)
!     if (me .eq. 0) then
!     print *,' after calling cloud type ib= ', ib
!     print *,' poi=',(poi(l),l=1,k)
!     print *,' qoi=',(qoi(l),l=1,k)
!     endif
! 
!     compute cloud amounts for the goddard radiation
!
       if (flx(kbl) .gt. 0.0) then
          pl = 0.5 * (prs(lct) + prs(lcb+1))
          cfac = min(1.0, max(0.0, (850.0-pl)*pfac))
       else
          cfac = 0.0
       endif
!
      if (lct .lt. lcb) then
        st        = poi(lcb)   - st
        sq        = qoi(lcb)   - sq
        sw(1)     = uvi(lcb,1) - sw(1)
        sw(2)     = uvi(lcb,2) - sw(2)
!
        poi(lcb)   = savt
        qoi(lcb)   = savq
        uvi(lcb,1) = savw(1)
        uvi(lcb,2) = savw(2)
        prs(lcb)   = sprs
        psj(lcb)   = sprj
        
        do l=lct,lcb
          poi(l)   = poi(l)     + st
          qoi(l)   = qoi(l)     + sq
          uvi(l,1) = uvi(l,1)   + sw(1)
          uvi(l,2) = uvi(l,2)   + sw(2)
          cln(l)   = min(cln(l) + cfac, 1.0)
        enddo
        do l=lct,lcb-1
          qli(l) = qli(lcb)
          qii(l) = qii(lcb)
        enddo
      else
cdir$ ivdep
        cln(lcb) = min(cln(lcb) + cfac, 1.0)
      endif
!
#ifdef DBG
!     print *,' rain=',rain(1),' ib=',ib
!     tem1 = 0.0
!     tem2 = 0.0
!     tem3 = 0.0
!     tem4 = 0.0
!     do l=1,k
!        cf1a = (100.0*cp) / (alhl*grav)
!        cf2a = 100.0 / grav
!        print *,' cp=',cp,' alhl=',alhl,' l=',l
!    *,' cf1a=',cf1a,' cf2a=',cf2a
!        print *,' tcu=',(tcu(i,l),i=1,len)
!        print *,' qcu=',(qcu(i,l),i=1,len)
!        do i=1,len
!           tem = prs(i,l+1)-prs(i,l)
!           cf1 = cf1a * tem
!           cf2 = cf2a * tem
!           tem1 = tem1 + q1(i,l)*cf1
!           tem2 = tem2 + (q2(i,l) - (alhf/alhl)*qii(i,l))*cf2
!           tem4 = tem4 + (qli(i,l)+(1.0+alhf/alhl)*qii(i,l))*cf2

!           tem1 = tem1 + (q1(i,l)-(alhf/cp)*qii(i,l))*cf1
!           tem2 = tem2 + q2(i,l)*cf2
!           tem4 = tem4 + (qli(i,l)+qii(i,l))*cf2
!        enddo
!     enddo
!     print *,' tem1=',tem1,' tem2=',tem2,' tem3=',tem3
!     do i=1,len
!        tem3 = tem3 + rain(i)
!     enddo
!     tem4 = tem4 + tem3
!     print *,' after cloud: tem1=',tem1,' tem2=',tem2,' tem3=',tem3
!    *,                    ' tem4=',tem4
!     print *,' rain=',(rain(i),i=1,len)
!
!     if (lprnt) then
!     print *,' pcu=',pcu(ix)
!     print *,' poi=',(poi(ix,l),l=1,k)
!     print *,' qoi=',(qoi(ix,l),l=1,k)
!     endif
!
!     print *,' rain=',rain(ix)
#endif

!   warining!!!!
!   by doing the following, cloud does not contain environmental
!   condensate!
!
        do l=1,k
          clw(l ) = clw(l) + qli(l) + qii(l)
          qli(l)  = 0.0
          qii(l)  = 0.0
!
!         q1(l) = 0.0
!         q2(l) = 0.0
        enddo
!       rain = 0.0
!
  100 continue
!
!        cln(ib) = min(cln(ib), 1.0)
!
      return
      end
