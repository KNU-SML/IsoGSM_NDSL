#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define PI pi
#define T0C t0c
#include <define.h>
      subroutine albaer(imx2,jmx2,
     1                  mon,slmsk,snowf,zorlf,coszf,tseaf,hprif,jsno,
#ifndef SW99UPDATE
     2                  alvsf,alnsf,alvwf,alnwf,facsf,facwf,paerf,
     3                  xlat,alvbr,alnbr,alvdr,alndr,paerr,latdef)
#else
     2                  alvsf,alnsf,alvwf,alnwf,facsf,facwf,
     3                  xlat,alvbr,alnbr,alvdr,alndr,latdef
     4,                 iaer,nfaer,kprfg,idxcg,cmixg,denng)
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c*******************************************************************
c   this program computes four components of surface albedos (i.e.
c     vis-nir, direct-diffused) based on brieglebs scheme. and
c     bilinearly interpolates albedo and aerosol distribution to
c     radiation grid.
c
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
c
c
#ifdef SW99UPDATE
c   updates to opac aerosol algorithm (1998)  .....y.h...sep99
c
#endif
c   input variables:
c     mon     - month of the year
c     slmsk   - sea(0),land(1),ice(2) mask on fcst model grid
c     snowf   - snow depth water equivalent in mm
c     zorlf   - surface roughness in cm
c     coszf   - cosin of solar zenith angle
c     tseaf   - sea surface temperature in k
c1197 hprif   - topographic sdv in m
c1197 jsno    - lat at 70 deg - indicating extent of perm snow cover
c     alvsf   - mean vis albedo with strong cosz dependency
c     alnsf   - mean nir albedo with strong cosz dependency
c     alvwf   - mean vis albedo with weak cosz dependency
c     alnwf   - mean nir albedo with weak cosz dependency
#ifndef SW99UPDATE
c     paerf   - aerosol distribution factor on fcst grid
#else
c     nfaer   - file number for input aerosol data
#endif
c
c   output variables:   (all on radiation grid)
c     alvbr   - vis beam surface albedo
c     alnbr   - nir beam surface albedo
c     alvdr   - vis diff surface albedo
c     alndr   - nir diff surface albedo
#ifndef SW99UPDATE
c     paerr   - aerosol distribution factor
#else
c     kprfg     - aerosol profile type index
c     idxcg     - aerosol components types indices
c     cmixg     - aerosol components mixing ratioes
c     denng     - first two layers aerosol number densities
#endif
c******************************************************************
#include <paramodel.h>
#ifdef SW99UPDATE
#include <aerparm.h>
#endif
c --- input
                       d i m e n s i o n
     1  slmsk(imx2,jmx2),snowf(imx2,jmx2),zorlf(imx2,jmx2)
     2, tseaf(imx2,jmx2),coszf(imx2,jmx2),hprif(imx2,jmx2)
#ifdef ALB4
     3, alvsf(imx2,jmx2), alnsf(imx2,jmx2)
     4, alvwf(imx2,jmx2), alnwf(imx2,jmx2)
#else
     3, alvsf(imx2,jmx2,4), alnsf(imx2,jmx2,4)
     4, alvwf(imx2,jmx2,4), alnwf(imx2,jmx2,4)
#endif
     5, facsf(imx2,jmx2),   facwf(imx2,jmx2)
#ifndef SW99UPDATE
     6, paerf(imx2,jmx2,5)
#endif
     7, xlat(imx2,jmx2)
c --- output
                       d i m e n s i o n
     1  alvbr(imx2,jmx2),alnbr(imx2,jmx2),alvdr(imx2,jmx2)
#ifndef SW99UPDATE
#ifdef RSMVECTORIZE
     2, alndr(imx2,jmx2),paerr(imx2,jmx2,5),latdef(jmx2)
#else
     2, alndr(imx2,jmx2),paerr(imx2,5,jmx2),latdef(jmx2)
#endif
#else
     2, alndr(imx2,jmx2),latdef(jmx2)
     3, idxcg(nxc,imxae,jmxae),  cmixg(nxc,imxae,jmxae)
     4, denng(ndn,imxae,jmxae),  kprfg(imxae,jmxae)
#endif
c --- internal variables
#ifdef MP
#ifdef RSM
#define ILOTS igrd12p_
#define JLOTS jgrd12p_
#else
#define ILOTS lonf2p_
#define JLOTS latg2p_
#endif
#else
#ifdef RSM
#define ILOTS igrd12_
#define JLOTS jgrd12_
#else
#define ILOTS lonf2_
#define JLOTS latg2_
#endif
#endif
      dimension alvbf(ILOTS,JLOTS),
     *          alnbf(ILOTS,JLOTS),
     *          alvdf(ILOTS,JLOTS),
     *          alndf(ILOTS,JLOTS),
     *          asnvb(ILOTS),asnnb(ILOTS),
     *          asnvd(ILOTS),asnnd(ILOTS),
     *          asevb(ILOTS),asenb(ILOTS),
     *          asevd(ILOTS),asend(ILOTS),
     *          fsno (ILOTS),fsea (ILOTS),
     *          rfcs (ILOTS),rfcw (ILOTS),
#ifndef SW99UPDATE
     *          flnd (ILOTS), mm(4)
#else
     *          flnd (ILOTS), mm(4),
     *          idxc (nxc),   cmix (nxc)
      character cline*80, ctyp*3, aerosol_file*40
#endif
c
#ifndef ALB4
      data mm /4, 7, 10, 13/
#endif
c
      logical snochk
      parameter(snodeg=70.* pi_ /180.)
c
      ifcs=imx2
      jfcs=jmx2
      irad=imx2
      jrad=jmx2
c
#ifndef ALB4
      if (mon.lt.1 .or. mon.gt.12) then
        write(6,4) mon
   4    format(3x,'error in month specification - mon=',i4)
        stop 4
      end if
      k=1
  10  if (mon .lt. mm(k)) go to 12
      k = k + 1
      go to 10
  12  k1 = k
      k2 = k + 1
      if (k2 .gt. 4) k2 = 1
      f12 = float(mon - mm(k1) + 3) / 3.0e0
#endif
#ifdef SW99UPDATE
c ==============================================
c ===  first section defines surface albedo  ===
c ==============================================
#endif
      do 100 j=1,jfcs
#ifdef SCM
        nxlat = nint(xlat(1,j)*180./pi_)
        jj = 90 - nxlat
        if (jj.gt.90) jj = 180 - jj
#else
        jj=latdef(j)
#endif
c
        do 20 i=1,ifcs
          snochk = .false.
          if(jj.le.jsno.or.(jsno.eq.0.and.abs(xlat(i,j)).gt.snodeg))
     1       snochk = .true.
c --- modified snow albedo scheme - units convert to m
c     (originally snowf in mm; zorlf in cm)
          asnow = 0.02*snowf(i,j)
          argh  = min(1.0, max(.025, 0.01*zorlf(i,j)))
          fsno0 = asnow / (argh + asnow)
          if (slmsk(i,j).eq.0.0 .and. tseaf(i,j).gt.271.2)
     1      fsno0 = 0.0
          fsno1 = 1.0 - fsno0
          flnd0 = facsf(i,j) + facwf(i,j)
          fsea0 = max(0.0, 1.0 - flnd0)
          fsno (i) = fsno0
          fsea (i) = fsea0 * fsno1
          flnd (i) = flnd0 * fsno1
c --- diffused sea surface albedo
          if (slmsk(i,j).eq.2.0 .or. tseaf(i,j) .lt. 271.2) then
            asevd(i) = 0.70
            asend(i) = 0.65
          else
            asevd(i) = 0.06
            asend(i) = 0.06
          end if
c --- diffused snow albedo
          if (slmsk(i,j).eq.1.0 .and.
     1       ((.not.snochk).or.tseaf(i,j).gt.271.2)) then
              hfac = max(0.10, min(1.0, 1.225-1.125e-3*hprif(i,j)))
          else
              hfac = max(0.75, min(1.0, 1.0625-0.3125e-3*hprif(i,j)))
          end if
          asnvd(i) = 0.90 * hfac
          asnnd(i) = 0.75 * hfac
  20    continue
c
        do 40 i=1,ifcs
c --- direct snow albedo
          if (coszf(i,j) .lt. 0.5) then
            csnow = 0.5 * (3.0 / (1.0+4.0*coszf(i,j)) - 1.0)
            asnvb(i) = min( 0.98, asnvd(i)+(1.0-asnvd(i))*csnow )
            asnnb(i) = min( 0.98, asnnd(i)+(1.0-asnnd(i))*csnow )
          else
            asnvb(i) = asnvd(i)
            asnnb(i) = asnnd(i)
          end if
c --- direct sea surface albedo
          if (coszf(i,j) .gt.0.0) then
            rfcs(i) = 1.4 / (1.0 + 0.4*coszf(i,j))
            rfcw(i) = 1.1 / (1.0 + 0.2*coszf(i,j))
            if (tseaf(i,j) .ge. t0c_) then
              asevb(i) = max(0.055, 0.026/(coszf(i,j)**1.7+0.065)
     1                 + 0.15 * (coszf(i,j)-0.1) * (coszf(i,j)-0.5)
     2                 * (coszf(i,j)-1.0))
              asenb(i) = asevb(i)
            else
              asevb(i) = asevd(i)
              asenb(i) = asend(i)
            end if
          else
            rfcs(i) = 1.0
            rfcw(i) = 1.0
            asevb(i) = asevd(i)
            asenb(i) = asend(i)
          end if
  40    continue
        do 60 i=1,ifcs
#ifdef ALB4
          a1   = alvsf(i,j) * facsf(i,j)
          b1   = alvwf(i,j) * facwf(i,j)
          a2   = alnsf(i,j) * facsf(i,j)
          b2   = alnwf(i,j) * facwf(i,j)
#else
          a1   = (alvsf(i,j,k1) + (alvsf(i,j,k2)-alvsf(i,j,k1))*f12)
     1         * facsf(i,j)
          b1   = (alvwf(i,j,k1) + (alvwf(i,j,k2)-alvwf(i,j,k1))*f12)
     1         * facwf(i,j)
          a2   = (alnsf(i,j,k1) + (alnsf(i,j,k2)-alnsf(i,j,k1))*f12)
     1         * facsf(i,j)
          b2   = (alnwf(i,j,k1) + (alnwf(i,j,k2)-alnwf(i,j,k1))*f12)
     1         * facwf(i,j)
#endif
          alvbf(i,j) = (a1*rfcs(i) + b1*rfcw(i))*flnd(i)
     1               + asevb(i)*fsea(i) + asnvb(i)*fsno(i)
          alvdf(i,j) = (a1         + b1        )*flnd(i)
     1               + asevd(i)*fsea(i) + asnvd(i)*fsno(i)
          alnbf(i,j) = (a2*rfcs(i) + b2*rfcw(i))*flnd(i)
     1               + asenb(i)*fsea(i) + asnnb(i)*fsno(i)
          alndf(i,j) = (a2         + b2        )*flnd(i)
     1               + asend(i)*fsea(i) + asnnd(i)*fsno(i)
  60    continue
c
 100  continue
c
c ... for one grid, no interpolation ...
c
      do j=1,jfcs
        do i=1,ifcs
          alvbr(i,j)=alvbf(i,j)
          alvdr(i,j)=alvdf(i,j)
          alnbr(i,j)=alnbf(i,j)
          alndr(i,j)=alndf(i,j)
        enddo
      enddo
#ifndef SW99UPDATE
      do k = 1,5
        do j=1,jfcs
          do i=1,ifcs
#ifndef RSMVECTORIZE
            paerr(i,k,j)=min(1.0,paerf(i,j,k))
            if (paerr(i,k,j).lt. 0.01) paerr(i,k,j) = 0.0
#else
            paerr(i,j,k)=min(1.0,paerf(i,j,k))
            if (paerr(i,j,k).lt. 0.01) paerr(i,j,k) = 0.0
#endif
          enddo
        enddo
      enddo
c
c....  final check to make total is one
      do 200 j=1,jrad
      do 200 i=1,irad
        if (slmsk(i,j).eq.0.0 .or. slmsk(i,j).eq.2.0) then
#ifndef RSMVECTORIZE
          psea = paerr(i,3,j) + paerr(i,5,j)
          paerr(i,1,j) = 0.0
          paerr(i,2,j) = 0.0
          paerr(i,4,j) = 0.0
          if (psea .gt. 1.0) then
            paerr(i,3,j) = min(1.0, paerr(i,3,j))
            paerr(i,5,j) = 1.0 - paerr(i,3,j)
          else if (psea .lt. 1.0) then
            paerr(i,3,j) = 1.0 - paerr(i,5,j)
          end if
        else
          plnd = paerr(i,1,j) + paerr(i,2,j) + paerr(i,4,j)
          paerr(i,5,j) = 0.0
          if (plnd .gt. 1.0) then
            paerr(i,3,j) = 0.0
            paerr(i,2,j) = min(1.0, paerr(i,2,j))
            if (paerr(i,1,j) .gt. 0.0) then
              paerr(i,4,j) = 0.0
              paerr(i,1,j) = 1.0 - paerr(i,2,j)
            else
              paerr(i,1,j) = 0.0
              paerr(i,4,j) = 1.0 - paerr(i,2,j)
            end if
          else if (plnd .lt. 1.0) then
c --- use mar-i as the background fill
            paerr(i,3,j) = 1.0 - plnd
#else
          psea = paerr(i,j,3) + paerr(i,j,5)
          paerr(i,j,1) = 0.0
          paerr(i,j,2) = 0.0
          paerr(i,j,4) = 0.0
          if (psea .gt. 1.0) then
            paerr(i,j,3) = min(1.0, paerr(i,j,3))
            paerr(i,j,5) = 1.0 - paerr(i,j,3)
          else if (psea .lt. 1.0) then
            paerr(i,j,3) = 1.0 - paerr(i,j,5)
          end if
        else
          plnd = paerr(i,j,1) + paerr(i,j,2) + paerr(i,j,4)
          paerr(i,j,5) = 0.0
          if (plnd .gt. 1.0) then
            paerr(i,j,3) = 0.0
            paerr(i,j,2) = min(1.0, paerr(i,j,2))
            if (paerr(i,j,1) .gt. 0.0) then
              paerr(i,j,4) = 0.0
              paerr(i,j,1) = 1.0 - paerr(i,j,2)
            else
              paerr(i,j,1) = 0.0
              paerr(i,j,4) = 1.0 - paerr(i,j,2)
            end if
          else if (plnd .lt. 1.0) then
c --- use mar-i as the background fill
            paerr(i,j,3) = 1.0 - plnd
#endif
          end if
        end if
 200  continue
#else
 
c ============================================================
c ===  second section defines global aerosol distribution  ===
c ============================================================

      if (iaer .eq. 1) then

      write(aerosol_file,101) mon
  101 format('aeropac3a.m',i2.2)
      print *,' aerosol_file=',aerosol_file,' imon=',imon
     *,' imxae=',imxae,' jmxae=',jmxae,' nxc=',nxc,' ndn=',ndn
c
      open(unit=nfaer,file=aerosol_file,status='OLD',form='FORMATTED')
c
      read(nfaer,102) cline
 102  format(a80)
c
      do j=1,jmxae
        do i=1,imxae
          read(nfaer,110) (idxc(k),cmix(k),k=1,nxc),kprf,denn,nc,ctyp
 110      format(5(i2,e11.4),i2,f8.2,i3,1x,a3)
c
          kprfg(i,j)     = kprf
          denng(1,i,j)   = denn       ! num density of 1st layer
          if (kprf .ge. 6) then
            denng(2,i,j) = cmix(nxc)  ! num density of 2dn layer
          else
            denng(2,i,j) = 0.0
          end if
c
          temp = 1.0
          do k=1,nxc-1
            idxcg(k,i,j) = idxc(k)    ! component index
            cmixg(k,i,j) = cmix(k)    ! component mixing ratio
            temp         = temp - cmix(k)
          end do
          idxcg(nxc,i,j) = idxc(nxc)
          cmixg(nxc,i,j) = temp       ! to make sure all add to 1.
        end do
      end do
c
      close (nfaer)

      endif
#endif
c
      return
      end
