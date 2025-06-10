#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NTRAC ntrac
#define PI pi
#define RD rd
#define RV rv
#include <define.h>
      subroutine cldjms(imx2,imx22,kmx,nbin,mcld,
     1           ps,q,t,vvel,cv,cvt,cvb,si,sl,
     2           slmsk,cld,mtop,mbot,cldary,ivva,invr,rhmax,
#ifdef DBG
     3           xlatrd,rhcld,istrat,ntotal,lat)
#else
     3           xlatrd,rhcld,istrat,ntotal)
#endif
#ifdef CRAY_THREAD
cfpp$ expand(fpvs,fpvs0)
ctune
cfpp$ noconcur r
#endif
c   from yh.rad.mdl93(cldnew28).......
c     later updated from yh.rad.mdl94(cldmul28)...22jan94
c     later updated from yh.rad.mdl94(cldml28a)... 1feb94
c     later updated from yh.rad.mdl94(cldml28b)... 5feb94
c               subr cldprp replaced
c               added vertical interp of cld-rh relations(istrat gt 1)
c     later updated from yh.rad.mdl94(cldml28e)... 11mar94
c               subr cldprp replaced,gcl adjusted
c     later updated from cloud6................... 24mar94
c              subr cldprp , low enhanced to old value..0.14..
c              subr gclnew , llyr calculation adj to old valu(kl-1)
c                             llyrl was ok.. ive removed it and
c                             replaced it by its equivalent, klowb
c     later updated from cloud6................... 30mar94
c               subr cldprp , low and middle (not cv) enhanced=0.10
c---------------------------------------------------------------------
c     nov., 1992 - y.h., k.a.c., and a.k.
c        cloud parameterization patterned after slingo and slingos
c        work (jgr, 1991).
c     stratiform clouds are allowed in any layer except the surface
c        and upper stratosphere.  the relative humidity criterion may
c        vary in different model layers.
c
c     output cloud amounts are in cldary(i,k), k=1 is the lowest
c        model layer, stratiform (str) and convective (cnv) types of
c        cloud are compressed into one word: camt = str + 1.0e4*cnv
c        low marine stratus amts are flaged by adding 2.
c
c   for istrat = 0, there is rh-cld relation for each layer..
c                      crit rh computed within..
c   for istrat = 1, rh-cld relation from tables created using
c                     mitchell-hahn tuning technique (a.f. rtneph obs)
c                  ...stratus computed similar to old opnl cldjms.....
c                      except no cloud below layer=klowb..approx 955mb
c
c     convective clouds are from model convective scheme and are
c        no longer broken into .75,.25,.25..rather cc itself is used..
c        convective still takes precedence over stratiform in radfs
c         but here cv+st merge exits in cldary...(in radiation use of
c        cc gives improvement to tropical middle cld (as did st+cv))
c
c     clouds are also divided into 3 atmospheric domains (l,m,h) for
c        diagnostic purposes.  they are computed from random overlap
c        assumption for separated cloud layers and maximum overlap
c        for adjacent cloud layers.  a total cloud fraction is also
c        computed.
c
c     h,m,l domain pressure tops 'ptop1(k)' vary linearly from
c        'ptopc(k,1)' at 45deg to 'ptopc(k,2)' at the pole
c
c     input variables:
c        ps (cb)       - surface pressure
c        q  (kg/kg)    - specific humidity
c        t  (deg k)    - absolute temperature
c        vvel(cb/sec)  - vertical velocity
c        cv,cvt,cvb    - conv cld fraction, top, bottom layer from
c                        kuo scheme
c        si,sl         - mdl sigma interface and layer mean
c        slmsk         - sea/land mask array(sea:0.,land:1.,snow:2.)
c        ivva          - flag to control vertical velocity adj.
c                        =1: with, =0: without
c        invr          - flag to control lapse rate inversion cld
c                        =1: with, =0: without
c        rhmax         - upper limit of relative humidity to
c                        form overcast cloud (cld fractn = 1.)
c
c  modify to as an array (h.-m. h. juang)
c       xlatrd        - current latitude in radians (1st data pt)
c                         for models with diff lat at each pt, need to
c                         use the lat of all points....careful.....
c        rhcld         - cloud-rh relations from mitchell+hahn,
c                        using a.f. rtneph analyses
c        istrat        - 0 or 1:for default or 'rhcld' tables
c                        in the stratiform cloud calculation
c
c    output variables:
c       cldary         - vertical column array of cloud fraction
c                        profile
c       cld            - cld fraction in 3 types of domains (l,m,h)
c                          and total in 4th layer
c       mtop,mbot      - top, bottom layers of clouds (l,m,h)
c
c--------------------------------------------------------------------
c
#include <paramodel.h>
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
#ifndef ICE
        parameter(nxx=7501)
        dimension tbpvs0(nxx)
        common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
#endif
                      p a r a m e t e r
     1 ( rd=rd_, rv=rv_, eps=rd/rv, epsm1=rd/rv-1.0, pi=pi_ )
                      d i m e n s i o n
     1  ps(imx22), cv(imx22),  cvt(imx22),  cvb(imx22)
     2, slmsk(imx22),si(imx22,kmx+1),  sl (imx22,kmx)
     3, t(imx22,kmx), vvel(imx22,kmx),   q(imx22,kmx*ntotal)
     4, cld(imx22,4),    mtop(imx22,3),      mbot(imx22,3)
     5, cldary(imx22,kmx),xlatrd(imx22)
c
c   rh-cld relationships for each point
c
      dimension rhcld(imx22,nbin,mcld)
c
c   ptopc(k,l): top presure of each cld domain (k=1-4 are sfc,l,m,h;
c       l=1,2 are low-lat (<45 degree) and pole regions)
c
#include <comcd1.h>
c
c  workspace
c
      logical bitx(ILOTS), 
     *        bity(ILOTS), 
     *        bitz(ILOTS), 
     *        bit1, bit2,
     *        bitm(ILOTS)
      dimension rhrh (ILOTS,levs_),
     *          prsly(ILOTS,levs_),
     *          dthdp(ILOTS,levs_),
     *          theta(ILOTS,levs_),
     *          kcut (ILOTS),kbase(ILOTS),
     *          kbt1 (ILOTS),kth1 (ILOTS),
     *          cl1  (ILOTS),
     *          kbt2 (ILOTS),kth2 (ILOTS),
     *          cl2  (ILOTS),
     *          kcvb (ILOTS),kcvt (ILOTS),
     *          omeg (ILOTS),
     *          ksave(ILOTS),
     *          ptop1(ILOTS,4)
c
c   begin here 
c
      kdim=kmx
      kdimp=kmx+1
      levm1=kmx-1
      levm2=kmx-2
      icloud = ntrac_ + 1
c
c  find top pressure for each cloud domain
c
      do 4 k=1,4
      do 4 i=1,imx2
       fac = max (0.0e0, 4.0e0*abs(xlatrd(i))/pi-1.0e0)
       ptop1(i,k) = ptopc(k,1) + (ptopc(k,2)-ptopc(k,1)) * fac
    4 continue
c
c  low cloud top sigma level, computed for each lat cause
c       domain definition changes with latitude...
      klow=kdim
      do 10 k=kdim,1,-1
      do 10 i=1,imx2
        if (si(i,k) .lt. ptop1(i,2) * 1.0e-3) klow = min(klow,k)
  10  continue
c
c  potential temp and layer relative humidity
c
      do 40 k=1,kdim
      do 40 i=1,imx2
        cldary(i,k) = 0.0e0
        prsly(i,k) = ps(i) * sl(i,k) * 10.0e0
        exnr = (prsly(i,k)*0.001e0) ** (-rocp)
        theta(i,k) = exnr * t(i,k)
#ifdef ICE
        es = fpvs(t(i,k))
#else
        xj1=min(max(c1xpvs0+c2xpvs0*t(i,k),1.),float(nxx))
        jx1=min(xj1,nxx-1.)
        fpvs01=tbpvs0(jx1)+(xj1-jx1)*(tbpvs0(jx1+1)-tbpvs0(jx1))
        es = fpvs01
c       es = fpvs0(t(i,k))
#endif
        qs = eps * es / (sl(i,k)*ps(i) + epsm1*es)
        rhrh(i,k) = max (0.0e0, min (1.0e0, q(i,k)/qs))
  40  continue
c
c   potential temp lapse rate
c
      do 50 k=1,levm1
      do 50 i=1,imx2
        dthdp(i,k) = (theta(i,k+1) - theta(i,k)) /
     1               (prsly(i,k+1) - prsly(i,k))
  50  continue
c 
c     find the stratosphere cut off layer for high cloud. it
c      is assumed to be above the layer with dthdp less than
c      -0.25 in the high cloud domain (from looking at 1 case).
c
      do 60 i=1,imx2
        kcut(i) = levm2
  60  continue
      do 80 k=klow+1,levm2
        bit1 = .false.
        do 70 i=1,imx2
          if (kcut(i).eq.levm2 .and. prsly(i,k).le.ptop1(i,3) .and.
     1        dthdp(i,k).lt.-0.25e0) then
            kcut(i) = k
          end if
          bit1    = bit1 .or. kcut(i).eq.levm2
  70    continue
        if (.not. bit1) go to 85
  80  continue
  85  continue
c
c ------------------------------------------------------------------
c
      if (istrat.le.0) then
c
c ------------------------------------------------------------------
c     default scheme ....tuned for 28 lyrs by y-t hou.
c     calculate stratiform cloud and put into array 'cldary'
c     the relative humidity criteria are preset for each model
c     sigma level, (1) for ocean points, and (2) for land points.
c ------------------------------------------------------------------
c
      do 130 k=klowb,levm2
        bit1 = .false.
        do 90 i=1,imx2
          bitx(i) = k.le.kcut(i)
          bit1 = bit1 .or. bitx(i)
  90    continue
        if (.not. bit1) go to 130
        do 100 i=1,imx2
          spnt = max (0.6e0, min (0.85e0, 0.96e0-0.6e0*sl(i,k)))
          cr1se1 = (0.41e0*sl(i,k) - 0.71e0)**2 + 0.52e0
          cr1sl1 = 0.8e0 - 0.167e0*sl(i,k)
          if (slmsk(i).eq.1.0e0) then
            cr1 = cr1se1
          else
            cr1 = cr1sl1
          end if
          xcrh1 = spnt * (rhmax - cr1)
          cr2 = cr1 + xcrh1
          xcrh2 = rhmax - cr2
          cl1(i) = max (0.e0, (rhrh(i,k)-cr1)/xcrh1) ** 3
          if (cl1(i).gt.1.0e0)
     1      cl1(i) = 1.e0 + sqrt((rhrh(i,k)-cr2)/xcrh2)
 100    continue
        do 120 i=1,imx2
          if (bitx(i)) then
            cldary(i,k) = min (1.0e0, 0.5e0*cl1(i))
          end if
 120    continue
 130  continue
c
c  special treatment on low clouds
c
      dvvcld = vvcld(1) - vvcld(2)
      rclap = 1.0e0 / (0.8e0 - crhrh)
      do 180 i=1,imx2
        kbase(i) = 0
 180  continue
c
      do 350 k=klowb,klowt
c
        do 190 i=1,imx2
          omeg(i) = 10.0e0 * vvel(i,k)
          cl1 (i) = 0.0e0
 190    continue
        if (ivva .le. 0) go to 250
c
c  vertical velocity adjustment on low clouds
c
        bit1 = .false.
        do 210 i=1,imx2
          bitx(i) = prsly(i,k).ge.ptop1(i,2) .and. cldary(i,k).gt.0.0e0
          bit1 = bit1 .or. bitx(i)
 210    continue
        if (.not. bit1) go to 250
        do 220 i=1,imx2
          if (bitx(i)) then
            if (omeg(i).ge.vvcld(1)) then
              cldary(i,k) = 0.0e0
            else if(omeg(i).gt.vvcld(2)) then
              cr1 = (vvcld(1) - omeg(i)) / dvvcld
c             cldary(i,k) = cldary(i,k) * cr1
              cldary(i,k) = cldary(i,k) * sqrt(cr1)
            endif
          endif
 220    continue
c
c  t inversion related stratus clouds
c
 250    if (invr .lt. 1) go to 350
        bit1 = .false.
        do 260 i=1,imx2
          bitx(i) = prsly(i,k).ge.pstrt .and. slmsk(i).le.0.0
     1              .and. dthdp(i,k).le.clapse
          bit1 = bit1 .or. bitx(i)
 260    continue
        if (.not. bit1) go to 350
        do 270 i=1,imx2
          if ( kbase(i).eq.0  .and. rhrh(i,k).gt.crhrh .and. bitx(i) )
     1      kbase(i) = k
 270    continue
        do 300 i=1,imx2
          if (kbase(i).gt.0 .and. bitx(i) .and. cldary(i,k+1).le.0.1e-1
     1        .and. cldary(i,k+2).le.0.1e-1) then
            cr1 = min (1.0e0,
     1            max (0.0e0,  16.67e0*(clapse-dthdp(i,k)) ))
            if(rhrh(i,kbase(i)).lt.0.8e0) then
              cr1 = cr1 * (rhrh(i,kbase(i))-crhrh) * rclap
            endif
c
c  for t inversion type cloud, add flag value of 2.0
c
            cldary(i,k) = max (cldary(i,k), cr1) + 2.0e0
          endif
 300    continue
 350  continue
c
      end if
c
      if (istrat.gt.0) then
c
c ------------------------------------------------------------------
c     calculate stratiform cloud and put into array 'cldary' using
c       the cloud-rel.humidity relationship from table look-up..where
c       tables obtained using k.mitchell frequency distribution tuning
c        (observations are daily means from us af rtneph).....k.a.c.
c       tables created without lowest 10 percent of atmos.....k.a.c.
c ------------------------------------------------------------------
c  this loop to retrieve cloud from rh rewritten 950113 -mi
c
      do klev=klowb,levm2
        do i=1,imx2
          kbase(i)=0
          bitx(i)=.false.
        enddo
        do kc=mcld,1,-1
          do i=1,imx2
            if(prsly(i,klev).ge.ptop1(i,kc+1)) kbase(i)=kc
          enddo
        enddo
        nx=0
        nhalf=(nbin+1)/2
        do i=1,imx2
          if(kbase(i).le.0.or.klev.gt.kcut(i)) then
            cldary(i,klev)=0.
          elseif(rhrh(i,klev).le.rhcld(i,1,kbase(i))) then
            cldary(i,klev)=0.
          elseif(rhrh(i,klev).ge.rhcld(i,nbin,kbase(i))) then
            cldary(i,klev)=1.
          else
            bitx(i)=.true.
            ksave(i)=nhalf
            nx=nx+1
          endif
        enddo
        dowhile(nx.gt.0)
          nhalf=(nhalf+1)/2
          do i=1,imx2
            if(bitx(i)) then
              crk=rhrh(i,klev)
              cr1=rhcld(i,ksave(i),kbase(i))
              cr2=rhcld(i,ksave(i)+1,kbase(i))
              if(crk.le.cr1) then
                ksave(i)=max(ksave(i)-nhalf,1)
              elseif(crk.gt.cr2) then
                ksave(i)=min(ksave(i)+nhalf,nbin-1)
              else
#ifndef INTERACTIVE_STRATUS
                  cldary(i,klev)=0.01*(ksave(i)+(crk-cr1)/(cr2-cr1))
#else
                  qsum = 0.
                  do ic = icloud,ntotal
                    kk = klev + levs_ * (ic-1)
                    qsum = q(i,kk) + qsum
                  enddo
                  onemrh = max(1.e-10,1.-rhrh(i,klev))
                  value = max(min(1000.*qsum/onemrh,50.),0.)
                  cldary(i,klev) = max(rhrh(i,klev)*(1.-exp(-value)),0.)
#endif
                bitx(i)=.false.
                nx=nx-1
              endif
            endif
          enddo
        enddo
      enddo
c
c     clean out not-suspected marine stratus regions...
c      cause tuning procedure not carried out down to lyr3 and we
c      get too much lo cloud if we don t clean it out..
c
      do 831 i=1,imx2
       bitm(i) = .true.
  831 continue
      do 833 k=klowb,llyr
       do 832 i=1,imx2
        if(bitm(i)) then
          bitm(i) = prsly(i,k).lt.pstrt
     1         .or. slmsk(i).gt.0.0 .or. dthdp(i,k).gt.clapkc.or.
     2        rhrh (i,k+1).gt.0.60e0.or.rhrh (i,k+2).gt.0.60e0
          kbase(i) = k
        endif
  832  continue
  833 continue
      do 835 k=1,llyr
       do 834 i=1,imx2
        if(bitm(i)) cldary(i,k) = 0.0e0
  834  continue
  835 continue
c
c  special treatment on low clouds
c
      dvvcld = vvcld(1) - vvcld(2)
c
      do 950 k=klowb,klow
c
        do 904 i=1,imx2
          omeg(i) = 10.0e0 * vvel(i,k)
          cl1 (i) = 0.0e0
 904    continue
        if (ivva .le. 0) go to 920
c
c  vertical velocity adjustment on low clouds
c
        bit1 = .false.
        do 906 i=1,imx2
          bitx(i) = prsly(i,k).ge.ptop1(i,2) .and. cldary(i,k).gt.0.0e0
          bit1 = bit1 .or. bitx(i)
 906    continue
        if (.not. bit1) go to 920
        if(k.gt.llyr) then
          do 910 i=1,imx2
            if (bitx(i)) then
              if(omeg(i).ge.vvcld(1)) then
                cldary(i,k) = 0.0e0
              else if(omeg(i).gt.vvcld(2)) then
                cr1 = (vvcld(1) - omeg(i)) / dvvcld
c               cldary(i,k) = cldary(i,k) * cr1
                cldary(i,k) = cldary(i,k) * sqrt(cr1)
              endif
            endif
 910      continue
        else
          do 915 i=1,imx2
c
c  no vvel filter for marine stratus region
c
            if (bitm(i)) then
              if (bitx(i)) then
                if(omeg(i).ge.vvcld(1)) then
                  cldary(i,k) = 0.0e0
                else if(omeg(i).gt.vvcld(2)) then
                  cr1 = (vvcld(1) - omeg(i)) / dvvcld
                  cldary(i,k) = cldary(i,k) * sqrt(cr1)
                endif
              endif
            endif
 915      continue
        endif
c
c t inversion related stratus clouds
c
 920    if (invr .lt. 1) go to 950
        if (k.gt.llyr) go to 950
        bit1 = .true.
        do 930 i=1,imx2
          bit1 = bit1 .and. bitm(i)
 930    continue
        if (bit1) go to 950
        do 940 i=1,imx2
         if (.not.bitm(i)) then
          if (dthdp(i,kbase(i)).gt.clpse) then
c
c   smooth transition for cloud when dthdp between
c       clapse and clapse+dclps  (-0.05 and -0.06)
c
           cfiltr = 1.0e0 - ((clpse - dthdp(i,kbase(i))) / dclps)
           cldary(i,k) = cldary(i,k)*cfiltr
          end if
c
c  for t inversion type cloud, add flag value of 2.0
c
          cldary(i,k) = cldary(i,k)+2.0e0
         end if
 940    continue
 950  continue
c 
      end if
c
c ------------------------------------------------------------------
c     add convective cloud into 'cldary', no merge at this point..
c     two types of clouds are separated by a factor of 1.0e+4
c ------------------------------------------------------------------
c
 360  bit1 = .false.
      do 380 i=1,imx2
        bitx(i) = cv(i).gt.0.0e0 .and. cvt(i).ge.cvb(i)
        bit1 = bit1 .or. bitx(i)
 380  continue
      if (.not. bit1) go to 550
      do 390 i=1,imx2
        if (bitx(i)) then
          kcvb(i) = nint(cvb(i))
          kcvt(i) = min(levm2, nint(cvt(i)))
        else
          kcvb(i) = 1
          kcvt(i) = 1
        end if
 390  continue
      do 450 k=klowb,levm2
        bit2 = .false.
        do 400 i=1,imx2
          bity(i) = bitx(i) .and. kcvb(i).le.k .and. kcvt(i).ge.k
          bit2 = bit2 .or. bity(i)
 400    continue
        if (.not. bit2) go to 450
        do 420 i=1,imx2
          if (bity(i)) cldary(i,k) = cldary(i,k)
     1               + 10.0e0 * aint(1.0e3 * cv(i))
 420    continue
 450  continue
c
c     if mean cvt layer higher than 400mb add anvil cirrus
c
      bit2 = .false.
      do 460 i=1,imx2
        bitz(i) = bitx(i) .and. prsly(i,kcvt(i)).le.cvtop
        bit2 = bit2 .or. bitz(i)
 460  continue
      if (.not. bit2) go to 500
      do 480 i=1,imx2
        if (bitz(i)) then
          kk = kcvt(i)
          cr1 = max (0.0e0, min (1.0e0, 2.0e0*(cv(i)-0.3e0)))
c
c    get stratus back before doing anvil calculation
c
          cr2 =  mod(cldary(i,kk),10.e0)
          cldary(i,kk) = cr2 + 10.0e0*aint(1.0e3*cr1)
        end if
 480  continue
 500  continue
 550  continue
c
c  introduce adjustment of clouds (mk)
c
#ifdef CLDADJ
#ifdef DBG
      if(lat.eq.10) then
        print *,'before cldadj lat=',lat
        print *,(cldary(20,k),k=1,levs_)
      endif
      if(lat.eq.30) then
        print *,'before cldadj lat=',lat
        print *,(cldary(20,k),k=1,levs_)
      endif
#endif
      call cldadj(cldary,imx22,imx2,kmx,ptop1,ILOTS,si,xlatrd)
#ifdef DBG
      if(lat.eq.10) then
        print *,'after cldadj lat=',lat
        print *,(cldary(20,k),k=1,levs_)
      endif
      if(lat.eq.30) then
        print *,'after cldadj lat=',lat
        print *,(cldary(20,k),k=1,levs_)
      endif
#endif
#endif
c
c -------------------------------------------------------------------
c     separate clouds into 3 pressure domains (l,m,h).  within each
c     of the domains, assume separated cloud layers are randomly
c     overlapped and adjacent cloud layers are maximum overlapped.
c     vertical location of each type of cloud is determined by
c     the thickest continuing cloud layers in the domain.
c -------------------------------------------------------------------
c  loop over 3 cloud domains (l,m,h)
c
      do 750 l=1,3
c
        do 580 i=1,imx2
          cld (i,l) = 0.0e0
          mtop(i,l) = 1
          mbot(i,l) = 1
          cl1 (i) = 0.0e0
          cl2 (i) = 0.0e0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
 580    continue
c
        do 700 k=2,levm2
          bit1 = .false.
          do 600 i=1,imx2
            bitx(i) = (prsly(i,k).ge.ptop1(i,l+1)) .and.
     1        (prsly(i,k).lt.ptop1(i,l)) .and. (cldary(i,k).gt.0.0e0)
            bit1 = bit1 .or. bitx(i)
 600      continue
          if (.not. bit1) go to 700
          do 630 i=1,imx2
            cr1  =  mod(cldary(i,k), 2.0e0)
            cr2  = float(int(cldary(i,k)) / 10) * 1.0e-3
            if (bitx(i)) then
              if(kth2(i).le.0) then
c
c kth2 le 0 : 1st cld layer.
c
                kbt2(i) = k
                kth2(i) = 1
              else
c
c  kth2 gt 0 : consecutive cld layer.
c
                kth2(i) = kth2(i) + 1
              endif
c
c   physical cloud as seen by radiation..conv takes precedence
c   except anvil cirrus not random overlapped with cv tower as
c   in radiation code(so hi may be slight underestimate)....
c
              if (cr2.gt.0.0e0) then
               cl2 (i) = max (cl2(i), cr2)
              else
               cl2 (i) = max (cl2(i), cr1)
              end if
            endif
 630      continue
          bit2 = .false.
c
c  bity=true if next lyr=clear or we change cloud domains..
c
          do 640 i=1,imx2
            bity(i) = bitx(i) .and. (cldary(i,k+1).le.0.0e0
     1                         .or.  prsly(i,k+1).lt.ptop1(i,l+1) )
            bit2 = bit2 .or. bity(i)
 640      continue
          if (.not. bit2) go to 700
c
c     at the domain boundary or separated cld lyrs, random overlap.
c     choose the thickest or the largest fraction amt as the cld
c     layer in that domain
c
          do 650 i=1,imx2
            if (bity(i)) then
              if (cl1(i).gt.0.0e0) then
                kbt1(i) = int( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))
     1                       / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))
     1                        / (cl1(i) + cl2(i)) ) + 1
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)
              else
                kbt1(i) = kbt2(i)
                kth1(i) = kth2(i)
                cl1 (i) = cl2 (i)
              endif
              kbt2(i) = 1
              kth2(i) = 0
              cl2 (i) = 0.0e0
            endif
 650      continue
 700    continue
c
c  finish one domain, save effective clouds
c
        do 720 i=1,imx2
          cld(i,l) =  cl1(i)
          mtop(i,l) = max(kbt1(i), kbt1(i)+kth1(i)-1)
          mbot(i,l) = kbt1(i)
 720    continue
 750  continue
c
c      calculate total cloud from the multi-lyr cloud array
c      in a manner as seen by the radiation code.
c      where, max overlap is used for vertically adjacent cloud layers
c      a clear layer separates two contiguously layered cloud types.
c      where, for convection any anvil is considered a separate
c         randomly overlapped cloud..
c      ilow=0,1 if no,yes preceeding model layer was cloudy..
c      clow contains the cloudiness of preceeding separate layered cld
c
      do 780 i=1,imx2
       cld(i,4) = 0.e0
       icvec = 0
       ilow = 0
       clow = 0.e0
       do 780 k=1,kdim
        ccldy =  mod(cldary(i,k), 2.0e0)
        ccvec = float(int(cldary(i,k)) / 10) * 1.0e-3
        if (ccvec.gt.0.e0) then
         ccldy = ccvec
         icvec = 1
        end if
        if (ccldy.gt.0.e0) then
         if (ilow.eq.0) then
          clow = ccldy
          ilow = 1
         else
          if (icvec.gt.0) then
           if (ccldy.ne.clow) then
c
c  if convective and an adjacent lyr=stratiform (ie ccldy changes),
c  then random overlap the preceeding cloud tower...
c
            cld(i,4) = cld(i,4) + (1.-cld(i,4))*clow
            clow = ccldy
           end if
          else
c
c max overlap for non convective adjacent cld layers...
c
           clow = max(ccldy,clow)
          end if
         end if
        else
         if (ilow.eq.1) then
c
c  if this is first clear layer in a gap betwixt cldlyrs, then
c   random overlap the preceeding clouds with the ones below..
c
          cld(i,4) = cld(i,4) + (1.-cld(i,4))*clow
          ilow = 0
         end if
        end if
 780  continue
      return
      end
