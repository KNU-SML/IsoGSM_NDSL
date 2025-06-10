#define G g
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2 latg2
#define LATG latg
#define LEVS levs
#define LONF2 lonf2
#define LONF lonf
#define LONFP lonfp
#define PI pi
#define RD rd
#include <define.h>
      subroutine cldprpi
     *        (idimt,sigl,prssl,prslv,t,cldtot,cldcnv,ibeg,ipts,lat
     1,                 icwp,clwp,slmsk,xlat,ktop,kbtm,nclds
     2,                 cldlw,taucl,cfac,cldsw
     3,                 cwp,cip,rew,rei,fice
     4,                 taulw,cld0)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c---------------------------------------------------------------------
c     aug., 2005 - a.s.
c        implement effective radius for ice clouds, Wyser(1998,JC,p1793-1802)
c        and for water clouds, Bower et al.(1994,JAS,p2722-2732) 
c        modify some parts like cldprp
c        consideration about in-cloud (cloud water/ice devide by cloudiness)
c     feb., 1993 - y.h.
c        cloud radiative properties calculations after davis (1982)
c        and harshvardhan et al. (1987).
c     nov., 1995 - y.h.
c        modified to provide mixed cloud overlapping scheme for
c        chou's sw radiation scheme (1992)_.
c     aug., 1998 - y.h.
c        modified to estimate effective radius of water/ice cloud
c        drop, fraction of ice water content, and cloud emissivity
c        for lw radiation calculation (kiehl et al. 1998,j.clim).
c     mar., 2000 - fanglin yang
c        to use chou_uiuc's lw scheme, save
c        taulw(i,l)   - cloud optical depth for each layer, k=1-->top layer
c        cld0 (i,l) - unweighted cloud fraction, k=1-->top layer
c
c--------------------------------------------------------------------
c     input variables:
c        sigl(i,k)     - mdl sigma layer mean value, k=1 is sfc layer
c        prslv(i,k)    - interface level pressure (mb) k=1 is toa
c        prssl(i,k)    - layer mean pressure (mb)      k=1 is mdl top layer
c                                                        l is mdl sfc layer
c        t (i,k)       - absolute temperature(K)  k=1 is mdl top layer
c                                                   l is mdl sfc layer
c                                                   lp1 is sfc
c        cldtot(i,k)   - stratiform cloud      k=1 is mdl sfc layer
c        cldcnv(i,k)   - convective cloud      k=1 is mdl sfc layer
c        ibeg,ipts     - indices for the beginnig no. and the
c                        total no. of array elements to be processed
c ---   modify xlatrd to general for regional and global (h.-m.h. juangk
c ***    xlat        - current latitude in radians (1st data pt)
c                         for models with diff lat at each pt, need to
c                         use the lat of each point....careful.....
c        icwp          - flag indicates the method used for cloud
c                        properties, =0 use t-p; =1 use clwp.
c        clwp          - layer cloud water+ice path (g/m**2), k=1 is toa
c        slmsk         - land/sea/ice mask (0:sea.1:land,2:ice)
c    output variables:
c --- clouds for lw rad. k=1 is the sfc, k=2 is the
c      lowest cloud layer, and so on
c        ktop,kbtm(i,k)- cloud top and bottom indeces, ktop and
c                        kbtm values from 1 to l model layers,
c                        with value of 1 being the top mdl layer
c        nclds(i)      - no. of separated cloud layers in a column
c        cldlw(i,k)    - cloud fractions for lw, emissivity adjusted
c  ***   ityp(i,k)     - type of clouds, ityp=1, and 2 are for
c                        the rh, and conv types
c --- clouds for sw rad. k=1 is the top layer or level
c        taucl(i,k)    - cloud optical depth in every model layer
c                        k=1 is the top layer (used for icwp=0 only)
c        cfac(i,k)     - fraction of clear sky view at the layer interfa
c        cldsw(i,k)    - layer cloud fractions for sw
c  - - following are output variables for icwp=1
c        cwp           - layer cloud water path (g/m**2)
c        cip           - layer cloud ice path (g/m**2)
c        rew (i,k)     - effective layer cloud water drop size (micron)
c        rei (i,k)     - effective layer cloud ice drop size (micron)
c        fice(i,k)     - fraction of cloud ice content
c
c--------------------------------------------------------------------
c
#include <paramodel.h>
#include <rdparm.h>
                      d i m e n s i o n
     1  cldtot(idimt,l),   cldcnv(idimt,l)
     2, sigl(idimt,l),           slmsk(idimt),  xlat(idimt)
     3, prslv(imax,lp1),   t (imax,lp1)
     4, prssl(imax,l)   

c - output array
                      d i m e n s i o n
     1  taucl(imax,l),     cldlw(imax,lp1),   cldsw(imax,l)
     2, ktop(imax,lp1),    kbtm(imax,lp1),    cfac(imax,lp1)
     3, ityp(imax,lp1),    nclds(imax)
     4, taulw(imbx,l),     cld0(imbx,l)
     5, rew (imax,l),      rei (imax,l),      fice(imax,l)
     6, clwp(imax,l),      cwp (imax,l),      cip (imax,l)
c
#include <comcd1.h>
c
        parameter (gr  = g_ / rd_, pi = pi_)

#define REI_W98    ! effective radius in ice clouds. Wyser(1998),JC,p1793-1802
#define REW_B94    ! effective radius in water clouds. Bower et al.(1994),JAS,p2722-2732

#ifdef REI_W98
        parameter (aiwc0=50.0,riwc0=1.0/aiwc0)     ! Wyser(1998)
#endif
#ifdef REW_B94
        parameter (c1 = 1.0/3.0)                   ! Bower et al.(1994)
        parameter (aocean = 150.0, aconti = 600.0) ! Bower et al.(1994) 
        parameter (cocean = 3.0/4.0/pi/aocean)     ! Bower et al.(1994)
        parameter (cconti = 3.0/4.0/pi/aconti)     ! Bower et al.(1994)
#endif

c ---  workspace ---
                      d i m e n s i o n
     1  mtyp (imax),    xamt (imax),    tauc (imax)
     2, kcld (imax),    mbtm (imax),    cl1  (imax)
     3, cl2  (imax),    alfa (imax)
                        l o g i c a l
     1  bitx(imax),  bitw(imax),  bit1,   bit2

c===>    begin here ................................................
      do i=1,imax
        kcld(i)    = 2
        mbtm(i)    = 1
        xamt(i)    = 0.0e0
        ityp(i,1)  = 0
        ktop(i,1)  = lp1
        kbtm(i,1)  = lp1
        cldlw(i,1) = 1.0e0
        cfac(i,1)  = 1.0e0
      end do

      do k=2,lp1
        do i=1,imax
          ityp(i,k)  = 0
          ktop(i,k)  = 1
          kbtm(i,k)  = 1
          cldlw(i,k) = 0.0e0
          cfac(i,k)  = 1.0e0
        end do
      end do

      do k=1,l
        do i=1,imax
          cldsw(i,k) = 0.0e0
          taucl(i,k) = 0.0e0
          taulw(i,k) = 0.0e0
          cld0(i,k)  = 0.0e0
          cwp (i,k)  = 0.0e0
          cip (i,k)  = 0.0e0
          rew (i,k)  = 10.0e0
          rei (i,k)  = 10.0e0
          fice(i,k)  = 0.0e0
        end do
      end do

c0499 if (icwp .eq. 1) then ! 0898 - use cloud water/ice content
        do k=1,l
          do i=1,ipts
            ir   = i + ibeg -1
            delt = 263.16 - t(i,k)
!--- fraction of ice
            fice(i,k) = amin1(1.0, amax1(0.0e0, delt/20.0e0 ))
                                        ! - fraction of ice in cloud
            cip (i,k) = clwp(i,k) * fice(i,k) 
            cwp (i,k) = clwp(i,k) - cip(i,k)

            rdz = gr*prssl(i,k)/(prslv(i,k+1)-prslv(i,k))/t(i,k)  ! 1./dz
!--- effective radius of water cloud
#ifndef REW_B94
            if (nint(slmsk(ir)) .eq. 1) then
              rew(i,k) = 5.0e0 + 5.0e0*amin1(1.0e0, amax1(0.0e0,
     1                   delt*0.05e0 )) ! - effective radius for water
            end if
#else
            if (fice(i,k).lt.1.0) then  ! -> fwater=1-fice>0
              aa = cwp(i,k) * rdz                  ! lwc[g/m**3]
              if (aa.gt.1.e-2) then
                if (nint(slmsk(ir)) .eq. 1) then    ! land
                  aa = cconti*aa
                else                                ! ocean
                  aa = cocean*aa
                endif
                aa = 100.0 * (aa)**c1
                rew(i,k) = max(3.0,min(aa,20.0))
              else
                rew(i,k) = 3.0
              endif
            endif
#endif

!--- effective radius of ice cloud
#ifndef REI_W98
            wgt = amin1(1.0e0, amax1(0.0e0,
     1           (sigl(i,lp1-k) ) / 0.4e0 ))
!    1           (sigl(i,lp1-k) - 0.1e0) / 0.7e0 ))
!     1           (sigl(i,lp1-k) - 0.3e0) / 0.5e0 ))
!    1           (sigl(i,lp1-k) - 0.4e0) / 0.4e0 ))
           rei(i,k) = 30.0 - 20.0*wgt  ! - effective radius for ice (default)
!            rei(i,k) = 50.0 - 40.0*wgt  ! - effective radius for ice
!           rei(i,k) = 90.0 - 70.0*wgt  ! - effective radius for ice
#else   ! Wyser(1998)
        ! dz=(R/g)*(dp/p)*t
            if (fice(i,k).gt.0.0) then
              aa = cip(i,k) * rdz                  ! iwc[g/m**3]
              if (aa.gt.1.e-6) then
                aa = aa*riwc0                      ! iwc/iwc0
                delt = 273.16 - t(i,k)
                b = -2.0 + 1.e-3*delt*sqrt(delt)*log10(aa)
                bb = b*b
                aa = 377.4 + 203.3*b + 37.91*bb + 2.3696*b*bb
                rei(i,k) = max(10.0,min(aa,100.0))
              endif
            endif
#endif
          enddo
        enddo
c0499 end if
c
c --- loop over mdl layers (bottom up)
c
      do k=2,l
!
        do i=1,ipts
          ir      = i + ibeg - 1
          cl1(i)  = cldtot(ir,k)
          cl2(i)  = cldcnv(ir,k)
        end do
!
        bit1 = .false.
        do i=1,ipts
          bitx(i) = cl1(i).gt.0.001e0 .or. cl2(i).gt.0.001e0
          bit1    = bit1 .or. bitx(i)
        end do
        if (bit1) then
          do i=1,ipts
c --- mtyp=1,2 for rh+stratus, and conv cloud types
            if (cl2(i) .gt. 0.001e0) then
              mtyp(i) = 2
            else if (cl1(i) .gt. 0.001e0) then
              mtyp(i) = 1
            else
              mtyp(i) = 0
            end if
          end do
!
          if (k .lt. l) then
            do i=1,ipts
              ir = i + ibeg - 1
c --- set bitw for clear gap above cloud
              bitw(i) = bitx(i) .and. cldtot(ir,k+1).le.0.001e0
     1                          .and. cldcnv(ir,k+1).le.0.001e0
            end do
          end if
          bit2 = .false.
          do i=1,ipts
            bit2 = bit2 .or. bitw(i)
            if (bitx(i)) then
              kkcl = kcld(i)
              if(ityp(i,kkcl) .eq. 0) then
                ityp(i,kkcl) = mtyp(i)
                xamt(i) = cl1(i)
                if (mtyp(i) .eq. 2) xamt(i) = cl2(i)
                mbtm(i) = k
              else if(ityp(i,kkcl).ne.mtyp(i) .or.
     1               (mtyp(i).eq.2 .and. xamt(i).ne.cl2(i)) ) then
                cldlw(i,kkcl)  = xamt(i)
                ktop(i,kkcl)   = lp1 - (k - 1)
                kbtm(i,kkcl)   = lp1 - mbtm(i)
                ityp(i,kkcl+1) = mtyp(i)
                mbtm(i)        = k
                xamt(i)        = cl1(i)
                if (mtyp(i).eq.2) xamt(i) = cl2(i)
                kcld(i) = kkcl + 1
              else if(mtyp(i).eq.1) then
                xamt(i) = amax1(xamt(i), cl1(i))
              endif
            end if
          end do
          if (bit2) then
            do i=1,ipts
              if (bitw(i)) then
                kkcl = kcld(i)
                cldlw(i,kkcl) = xamt(i)
                ktop(i,kkcl)  = lp1 - k
                kbtm(i,kkcl)  = lp1 - mbtm(i)
                kcld(i)       = kkcl + 1
                mtyp(i)       = 0
                mbtm(i)       = 1
                xamt(i)       = 0.0e0
              end if
            end do
          endif                             ! bit2 endif
        endif                               ! bit1 endif
      enddo                                 ! the k loop ends here!

c --- record num of cld lyrs and find max num of cld lyrs
      mclds = 0
      do i=1,ipts
        nclds(i) = kcld(i) - 2
        mclds    = max(mclds, nclds(i))
      end do

c     write(6,221) mclds, ibeg
c221  format(' in cldprp: maxclds =',i4,' ibeg=',i4)
c     if (mclds .eq. 0) return
c
c --- estimate cloud optical properties from t and q  -- (top down)
c
      do nncld=1,mclds
        nc = mclds - nncld + 2
c
        do i=1,ipts
          bitx(i) = cldlw(i,nc) .ge. 0.001e0
c         bity(i) = bitx(i) .and. ityp(i,nc).eq.2
          bitw(i) = bitx(i)
        end do

c --- find top pressure for mid cloud (3) domain=function of latitude
        minktp=lp1
        maxkbt=1
        do i=1,ipts
          if (bitx(i)) then
            minktp = min(minktp,ktop(i,nc))
            maxkbt = max(maxkbt,kbtm(i,nc))
          end if
        end do
c       write(6,241) nc,minktp, maxkbt
c241    format(3x,'nc, minktp, maxkbt =',3i6/3x,'bitx, bitw :')

c --- find clear sky view at each levels

        do kk=minktp,maxkbt
          do i=1,ipts
            if (kk.ge.ktop(i,nc) .and.
     1          kk.le.kbtm(i,nc) .and. bitx(i)) then
              cldsw(i,kk) = cldlw(i,nc)
              if (bitw(i)) then
                cfac(i,kk+1) = cfac(i,kk) * (1.0e0 - cldsw(i,kk))
                bitw(i) = .false.
              else
                cfac(i,kk+1) = cfac(i,kk)
              end if
            elseif (kk.gt.kbtm(i,nc) .and. bitx(i)) then
              cfac(i,kk+1) = cfac(i,kk)
            end if
          end do
        end do
c
        mkbtp1 = maxkbt + 1
        do k=mkbtp1,l
          do i=1,ipts
            if (bitx(i)) cfac(i,k+1) = cfac(i,mkbtp1)
          end do
        end do
c
        do i=1,ipts
          tauc(i) = 0.0e0
        end do

        if (icwp .ne. 1) then

cconv - reduce conv cloud amount for sw rad
          do i=1,ipts
            if (ityp(i,nc) .eq. 2) then
c0799         alfa(i) = cldlw(i,nc)
              alfa(i) = 
     1            amax1(0.25e0, 1.0e0-0.125e0*(kbtm(i,nc)-ktop(i,nc)))
            else
!             alfa(i) = sqrt(cldlw(i,nc))
              alfa(i) = 1.0
            end if
          end do
!
c --- calc cld thickness delp and mean temp (celsius)
          do kk=minktp,maxkbt
            do i=1,ipts
              if (kk.ge.ktop(i,nc) .and.
     1            kk.le.kbtm(i,nc) .and. bitx(i)) then
                delp = prslv(i,kk+1) - prslv(i,kk)
                tcld = t(i,kk) - 273.16e0
c --- convective cloud
                if (ityp(i,nc) .eq. 2) then
                  tau0 = delp * 0.05e0
c    comment out like cldprp
c0499 - if conv cld, set to water cloud only
c                  fice(i,kk) = 0.0e0
c --- rh clouds
                else
                  if (tcld .le. -10.0e0) then
                    tau0 = delp
     1                   * amax1(0.1e-3, 2.00e-6*(tcld+82.5e0)**2)
                  else
                    tau0 = delp
     1                   * amin1(0.08e0, 6.949e-3*tcld+0.08e0)
                  end if
                end if
                taulw(i,kk) = tau0
                tauc(i)     = tauc(i) + tau0
                taucl(i,kk) = tau0 * alfa(i) * cldlw(i,nc)
                cld0(i,kk)  = cldlw(i,nc)
              end if
            end do
          end do

c --- calc cld emis and effective cloud cover for lw
           do i=1,ipts
            if (bitx(i))
     1        cldlw(i,nc) = cldlw(i,nc)*(1.0e0-exp(-0.75e0*tauc(i)))
          enddo
!
        else                         ! prognostic liquid water
!
c --- calc tauc and emis
          do kk=minktp,maxkbt
            do i=1,ipts
              if (kk .ge. ktop(i,nc) .and.
     1            kk .le. kbtm(i,nc) .and. bitx(i)) then
c    comment out like cldprp
c0499 - if conv cld, set to liquid content only
c                if (ityp(i,nc) .eq. 2) then
c                  fice(i,kk) = 0.0e0
c                  cwp (i,kk) = cwp(i,kk) + cip(i,kk)
c                  cip (i,kk) = 0.0e0
c                end if
                awlw    = 0.090361e0 * cwp(i,kk)
                ailw    = (0.005e0 + 1.0e0/rei(i,kk)) * cip(i,kk)
                tau0 = awlw + ailw
                taulw(i,kk) = tau0
                tauc(i) = tauc(i) + tau0
                cld0(i,kk)  = cldlw(i,nc)
!
              endif
            end do
          end do
          do i=1,ipts
            if (bitx(i))
     1        cldlw(i,nc) = cldlw(i,nc) * amin1(1.0e0,
     2                      1.0e0 - exp(-1.66e0*tauc(i)) )
          enddo
        end if
!
      end do               ! end of nncld loop!
!
      if (icwp .eq. 1) then
        do k=1,l
        do i=1,ipts
          if (cldsw(i,k) .lt. 0.001e0) then
            taucl(i,k) = 0.0e0
            cldsw(i,k) = 0.0e0
            cwp  (i,k) = 0.0e0
            cip  (i,k) = 0.0e0
          end if
        end do
        end do
      end if
!
      if (ipts .lt. imax) then     ! --- fill up arrays
        ipts1 = ipts + 1
        do i=ipts1,imax
          nclds(i) = nclds(ipts)
        end do
        do k=1,lp1
          do i=ipts1,imax
            cldlw(i,k) = cldlw(ipts,k)
            ktop(i,k) = ktop(ipts,k)
            kbtm(i,k) = kbtm(ipts,k)
            cfac(i,k) = cfac(ipts,k)
          end do
        end do
        do k=1,l
          do i=ipts1,imax
            taulw(i,k) = taulw(ipts,k)
            cld0 (i,k) = cld0 (ipts,k)
            taucl(i,k) = taucl(ipts,k)
            cldsw(i,k) = cldsw(ipts,k)
            cwp  (i,k) = cwp  (ipts,k)
            cip  (i,k) = cip  (ipts,k)
            rew  (i,k) = rew  (ipts,k)
            rei  (i,k) = rei  (ipts,k)
            fice (i,k) = fice (ipts,k)
          end do
        end do
      endif

#ifdef DBG
#ifndef MP
#ifdef RSM
      latp = jgrd12_
      ir = igrd12_/4   ! center point
#else
      latp = latg2_
      ir = lonf2_/4   ! (180,0) point
#endif
#ifdef SCM
      i = 1
      ir = 1
#else
      i   = ir - ibeg +1
      if(i.ge.1.and.i.le.ipts.and.lat.eq.latp) then
#endif
       write(6,501) ir,lat
       write(6,502)
501   format(//,'=== cloud properties in cldprpi at (i,lat) ',
     1    2i10,'   ===',/)
502   format('    k','      clwp','      fice','       cwp',
     1'       cip','       rew', '       rei','      ityp')
503   format(i5,6f10.5,i10)
       do k=1,l
         write(6,503) k,clwp(i,k),fice(i,k),cwp(i,k),
     1         cip(i,k),rew(i,k),rei(i,k),ityp(i,k)
       enddo
c
       write(6,505)
505   format(//,'    k','   tem (C)','     cldlw',' ktop',' kbtm',
     1'     cldsw','      cfac','    cldtot','    cldcnv')
506   format(i5,2f10.5,2i5,4f10.5)
       do k=1,l
         write(6,506) k,t(i,k)-273.16,cldlw(i,k),ktop(i,k),kbtm(i,k),
     1     cldsw(i,k),cfac(i,k),cldtot(ir,l-k+1),cldcnv(ir,l-k+1)
       enddo
#ifndef SCM
      endif
#endif
#endif
#endif
c
      return
      end
