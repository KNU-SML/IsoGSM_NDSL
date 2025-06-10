#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine cldprp(idimt,prslv,t,cldary,ibeg,ipts,xlatrd,
     1                  ktop,kbtm,nclds,cldlw,taucl,cfac,cldsw,
     2                  taulw,cld0)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c  mar., 2000 - fanglin yang
c  to use chou_uiuc's lw scheme, save
c   taulw(i,l)   - cloud optical depth for each layer, k=1-->top layer
c   cld0 (i,l) - unweighted cloud fraction, k=1-->top layer
c---------------------------------------------------------------------
c     feb., 1993 - y.h.
c        cloud radiative properties calculations after davis (1982)
c        and harshvardhan et al. (1987).
c     nov., 1995 - y.h.
c        modified to provide mixed cloud overlapping scheme for
c        chous sw radiation scheme (1992).
c--------------------------------------------------------------------
c     input variables:
c        prslv(i,k)    - level pressure (mb) k=1 is toa
c        t (i,k)       - absolute temperature, k=1 is top layer (k)
c        cldary(i,k)   - cloud array contains compressed cloud
c                        fractions of 3 types (stratiform, conv
c                        and stratus), k=1 is the mdl sfc layer
c        ibeg,ipts     - indices for the beginnig no. and the
c                        total no. of array elements to be processed
c ---   modify xlatrd to general for regional and global (h.-m.h. juangk
c ***    xlatrd        - current latitude in radians (1st data pt)
c                         for models with diff lat at each pt, need to
c                         use the lat of each point....careful.....
c    output variables:
c --- clouds for lw rad. k=1 is the sfc, k=2 is the
c      lowest cloud layer, and so on
c        ktop,kbtm(i,k)- cloud top and bottom indeces, ktop and
c                        kbtm values from 1 to l model layers,
c                        with value of 1 being the top mdl layer
c        nclds(i)      - no. of separated cloud layers in a column
c        cldlw(i,k)    - cloud fractions for lw, emissivity adjusted
c        emis(i,k)     - cloud emissivity
c  ***   ityp(i,k)     - type of clouds, ityp=1, and 2 are for
c                        the rh, and conv types
c --- clouds for sw rad. k=1 is the top layer or level
c        taucl(i,k)    - cloud optical depth in every model layer
c                        k=1 is the top layer
c        cfac(i,k)     - feaction of clear sky view at the layer interface
c        cldsw(i,k)    - layer cloud fractions for sw
c
c--------------------------------------------------------------------
c
#include <paramodel.h>
#include <rdparm.h>
                      d i m e n s i o n
     1  prslv(imbx,lp1),   cldary(idimt,l),   xlatrd(idimt)
     2, taucl(imbx,l),     t   (imbx,lp1),    cldlw(imbx,lp1)
     3, ktop(imbx,lp1),    kbtm(imbx,lp1),    cldsw(imbx,l)
     4, ityp(imbx,lp1),    nclds(imax),       cfac(imbx,lp1)
     5, taulw(imbx,l),     cld0(imbx,l)
c
#include <comcd1.h>
c
c ---  workspace ---
                      d i m e n s i o n
     1  mtyp (imax),    xamt (imax),    tauc (imax)
     2, kcld (imax),    mbtm (imax),    cl1  (imax)
c97  3, cl2  (imax),    ptopd(imax),    alfa (imax)
     3, cl2  (imax),                    alfa (imax)
                        l o g i c a l
     1  bitx(imax),   bitw(imax),  bit1,   bit2
c===>    begin here ................................................
      do 20 i=1,imax
        kcld(i) = 2
        mbtm(i) = 1
        mtyp(i) = 0
        xamt(i) = 0.0e0
        ityp(i,1) = 0
        ktop(i,1) = lp1
        kbtm(i,1) = lp1
        cldlw(i,1)= 1.0e0
        cfac(i,1) = 1.0e0
  20  continue
      do 30 k=2,lp1
      do 30 i=1,imax
        ityp(i,k) = 0
        ktop(i,k) = 1
        kbtm(i,k) = 1
        cldlw(i,k)   = 0.0e0
c       cld0 (i,k)   = 0.0e0
        cldsw(i,k-1) = 0.0e0
        taucl(i,k-1) = 0.0e0
        taulw(i,k-1) = 0.0e0
        cld0 (i,k-1)   = 0.0e0
        cfac(i,k) = 1.0e0
  30  continue
c
c --- loop over mdl layers (bottom up)
c
      do 200 k=2,l
c
        bit1 = .false.
        do 60 i=1,ipts
          ir = i + ibeg - 1
          bitx(i) = cldary(ir,k).gt.0.0e0
          bit1 = bit1 .or. bitx(i)
  60    continue
        if (.not. bit1) go to 200
c --- decompress cloud array
        do 70 i=1,ipts
          cl1(i) = 0.0e0
          cl2(i) = 0.0e0
          bitw(i) = bitx(i)
  70    continue
        do 80 i=1,ipts
          if (bitx(i)) then
            ir = i + ibeg - 1
            cl1(i) = amod(cldary(ir,k), 2.0e0)
            cltemp = amod(cldary(ir,k), 10.0e0)
            cl2(i) = 1.0e-4 * (cldary(ir,k) - cltemp)
c --- mtyp=1,2 for rh+stratus, and conv cloud types
            if (cl2(i) .gt. 0.0e0) then
              mtyp(i) = 2
            else
              mtyp(i) = 1
            end if
          end if
  80    continue
        if(k.lt.l) then
          do 100 i=1,ipts
            ir = i + ibeg - 1
            if(bitw(i)) then
              bitw(i) = cldary(ir,k+1).le.0.0e0
            endif
 100      continue
        endif
        bit2 = .false.
        do 110 i=1,ipts
          bit2 = bit2 .or. bitw(i)
          if (bitx(i)) then
            if(ityp(i,kcld(i)).eq.0) then
              ityp(i,kcld(i)) = mtyp(i)
              xamt(i) = cl1(i)
              if (mtyp(i) .eq. 2) xamt(i) = cl2(i)
              mbtm(i) = k
            else if(ityp(i,kcld(i)).ne.mtyp(i) .or.
     1             (mtyp(i).eq.2 .and. xamt(i).ne.cl2(i)) ) then
              cldlw(i,kcld(i)) = xamt(i)
              ktop(i,kcld(i)) = lp1 - (k - 1)
              kbtm(i,kcld(i)) = lp1 - mbtm(i)
              ityp(i,kcld(i)+1) = mtyp(i)
              mbtm(i) = k
              xamt(i) = cl1(i)
              if (mtyp(i).eq.2) xamt(i) = cl2(i)
              kcld(i) = kcld(i) + 1
            else if(mtyp(i).eq.1) then
              xamt(i) = amax1(xamt(i), cl1(i))
            endif
          end if
 110    continue
        if (.not. bit2) go to 200
        do 160 i=1,ipts
          if (bitw(i)) then
            cldlw(i,kcld(i)) = xamt(i)
            ktop(i,kcld(i)) = lp1 - k
            kbtm(i,kcld(i)) = lp1 - mbtm(i)
            kcld(i) = kcld(i) + 1
            mtyp(i) = 0
            mbtm(i) = 1
            xamt(i) = 0.0e0
          end if
 160    continue
c
 200  continue
c --- record num of cld lyrs and find max num of cld lyrs
      mclds = 0
      do 220 i=1,ipts
        nclds(i) = kcld(i) - 2
        mclds = max(mclds, nclds(i))
 220  continue
c
c --- estimate cloud optical properties from t and q
c     (top down)
c
      do 400 nncld=1,mclds
        nc = mclds - nncld + 2
c
        do 230 i=1,ipts
          tauc(i) = 0.0e0
          bitx(i) = cldlw(i,nc) .gt. 0.0e0
          bitw(i) = bitx(i)
cconv - reduce conv cloud amount for sw rad
          if (ityp(i,nc) .eq. 2) then
            alfa(i) = amax1(0.25e0,
     1                1.0e0-0.125e0*(kbtm(i,nc)-ktop(i,nc)))
          else
            alfa(i) = 1.0e0
          end if
 230    continue
c --- find top pressure for mid cloud (3) domain=function of latitude
        minktp=lp1
        maxkbt=1
        do 240 i=1,ipts
          if (bitx(i)) then
            minktp = min(minktp,ktop(i,nc))
            maxkbt = max(maxkbt,kbtm(i,nc))
          end if
 240    continue
       if (nncld .eq. 1) kstrt = minktp
c --- calc cld thickness delp and mean temp (celsius)
        do 260 kk=minktp,maxkbt
        do 260 i=1,ipts
          if (kk.ge.ktop(i,nc) .and. kk.le.kbtm(i,nc) .and.bitx(i)) then
            delp = prslv(i,kk+1) - prslv(i,kk)
            tcld = t(i,kk) - 273.16e0
c --- convective cloud
            if (ityp(i,nc) .eq. 2) then
              tau0 = delp * 0.06e0
c --- rh clouds
            else
              if (tcld .le. -10.0e0) then
                tau0 = delp
     1               * amax1(0.1e-3, 2.00e-6*(tcld+82.5e0)**2)
              else
                tau0 = delp*amin1(0.08e0, 6.949e-3*tcld+0.08e0)
              end if
            end if
            tauc(i) = tauc(i) + tau0
            cldsw(i,kk) = cldlw(i,nc)
            taucl(i,kk) = tau0 * alfa(i) * cldlw(i,nc)
            cld0 (i,kk) = cldlw(i,nc)
            taulw(i,kk) = tau0
            if (bitw(i)) then
              cfac(i,kk+1) = cfac(i,kk) * (1.0e0 - cldsw(i,kk))
              bitw(i) = .false.
            else
              cfac(i,kk+1) = cfac(i,kk)
            end if
          elseif (kk.gt.kbtm(i,nc) .and. bitx(i)) then
            cfac(i,kk+1) = cfac(i,kk)
          end if
 260    continue
        mkbtp1 = maxkbt + 1
        do 280 k=mkbtp1,l
        do 280 i=1,ipts
          if (bitx(i)) cfac(i,k+1) = cfac(i,mkbtp1)
 280    continue
c --- calc cld emis
        do 320 i=1,ipts
          if (bitx(i))
     1      cldlw(i,nc) = cldlw(i,nc)*(1.0e0-exp(-0.75e0*tauc(i)))
 320    continue
c
 400  continue
c --- cloud scaled for sw
      do 420 kk=1,l
      do 420 i=1,ipts
        if (cfac(i,lp1) .lt. 1.0e0)
     1    taucl(i,kk) = taucl(i,kk) / (1.0e0 - cfac(i,lp1))
 420  continue
      if (ipts .eq. imax) go to 565
      ipts1 = ipts + 1
      do 520 i=ipts1,imax
        nclds(i) = nclds(ipts)
 520  continue
      do 540 k=1,lp1
      do 540 i=ipts1,imax
        cldlw(i,k) = cldlw(ipts,k)
        ktop(i,k) = ktop(ipts,k)
        kbtm(i,k) = kbtm(ipts,k)
        cfac(i,k) = cfac(ipts,k)
 540  continue
      do 560 k=1,l
      do 560 i=ipts1,imax
        taulw(i,k) = taulw(ipts,k)
        cld0 (i,k) = cld0 (ipts,k)
        taucl(i,k) = taucl(ipts,k)
        cldsw(i,k) = cldsw(ipts,k)
 560  continue
 565  continue
c
      return
      end
