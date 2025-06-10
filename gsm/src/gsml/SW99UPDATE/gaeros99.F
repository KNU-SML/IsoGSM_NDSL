#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONFP lonfp
      subroutine gaeros99(si,sl)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c********************************************************************
c  setup common block 'swaer' for aerosols and rayleigh scattering
c  optical properties in eight uv+vis bands and four nir bands.
c   band: 1. 0.175-0.225 (uv-c)     2. 0.225-0.245;0.260-0.280 (uv-c)
c         3. 0.245-0.260 (uv-c)     4. 0.280-0.295 (uv-b)
c         5. 0.295-0.310 (uv-b)     6. 0.310-0.320 (uv-b)
c         7. 0.320-0.400 (uv-a)     8. 0.400-0.700 (par)
c         9. 2.27 - 4.0  (nir)     10. 1.22 - 2.27 (nir)
c        11. 0.70 - 1.22 (nir)     12. 0.70 - 4.0  (nir)
c  ref: wmo report wcp-112 (1986)
c
c  there are seven typical vertical structures:
c     1.antarctic, 2.arctic, 3.continent, 4.maritime, 5.desert,
c     6.maritime with mineral overlay, 7.continent with mineral overlay
c
c     sigref  - ref. sigma level                   n/d   ndm*nae
c  arrays in the common block:
c     haer    - scale height of aerosols           km    ndm*nae
c     hh      - atmospheric scale height           km    l
c     hz      - level height                       km    l+1
c     dz      - layer thickness                    km    l
c     idm     - domain index                       n/d   l*nae
c     taur    - rayleigh scattering optical depth  n/d   l*nbd
c********************************************************************
c
#include <paramodel.h>
#include <rdparm99.h>
#include <comswaer99.h>
      dimension si(lp1),sl(l),sigln(lp1)
c
c===> ... compute layer distributions of rayleigh scattering
c
      do 10 n=1,nbd
      do 10 k=1,l
        taur(k,n) = tauray(n) * (si(lp1-k)-si(lp1-k+1))
  10  continue
c
c===> ... setup log sigma array (set toa sigma=0.0001)
c         rem: si,sl k=1 is sfc; but in radiation k=1 is toa
      sigln(1) = alog(1.0e-4)
      do 20 k=1,l
        sigln(k+1) = alog(si(lp1-k))
  20  continue
      do 30 k=1,l
        hh(k) = 6.05e0 + 2.5e0 * sl(lp1-k)
        dz(k) = hh(k) * (sigln(k+1)-sigln(k))
  30  continue
      hz(lp1) = 0.0
      do k=l,1,-1
        hz(k) = hz(k+1) + dz(k)
      end do
c
      do 60 iaer=1,nae
        idom = 1
        do 50 k=l,1,-1
          if (si(lp1-k+1) .lt. sigref(idom,iaer)) then
            idom = idom + 1
            if (idom.eq.2 .and.
     1          sigref(2,iaer).eq.sigref(3,iaer)) then
              idom = 3
            end if
          end if
          idm (k,iaer) = idom
  50    continue
  60  continue
c
      return
      end
