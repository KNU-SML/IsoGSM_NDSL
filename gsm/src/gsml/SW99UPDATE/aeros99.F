#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONFP lonfp
      subroutine aeros99(ipts,ib,kaer,kprf,idxc,cmix,denn,rh
     1,                 tau,ssa,asy,taurs)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c********************************************************************
c  compute aerosols optical properties in eight uv+vis bands and
c  four nir bands. there are seven different vertical prifile
c  structures. in the troposphere, aerosol distribution at each
c  grid point is composed from up to six components out of a total
c  of ten different substances.
c  1999-10-13  y.h.  updated to opac data (1998)
c
c   band: 1. 0.175-0.225 (uv-c)     2. 0.225-0.245;0.260-0.280 (uv-c)
c         3. 0.245-0.260 (uv-c)     4. 0.280-0.295 (uv-b)
c         5. 0.295-0.310 (uv-b)     6. 0.310-0.320 (uv-b)
c         7. 0.320-0.400 (uv-a)     8. 0.400-0.700 (par)
c         9. 2.27 - 4.0  (nir)     10. 1.22 - 2.27 (nir)
c        11. 0.70 - 1.22 (nir)     12. 0.70 - 4.0  (nir)
c
c  input parameters:
c     ib   - spectral band index                   -    1
c     kaer - =0 do not compute aerosols            -    1
c            =1 compute aerosol profiles
c     kprf - indecies of aerosol prof structures   -    imax
c     idxc - indecies of aerosol components        -    nxc*imax
c     cmix - mixing ratioes of aerosol components  -    ncx*imax
c     denn - aerosol number densities              -    ndn*imax
c     rh   - relative humidity in fraction         -    imax*l
c
c  output parameters:
c     tau  - optical depth                         n/d
c     ssa  - single scattering albedo              n/d
c     asy  - asymmetry parameter                   n/d
c     taurs- rayleigh scattering optical depth     n/d
c
c  variables in common block:
c     haer - scale height of aerosols              km   ndm*nae
c     idm  - aerosol domain index                  -    l*nae
c     dz   - layer thickness                       km   l
c     hz   - level height                          km   l+1
c     taur - rayleigh scattering optical depth     -    l*nbd
c********************************************************************
c
#include <paramodel.h>
#include <rdparm99.h>
#include <aerparm.h>
#include <comswaer99.h>
c --- input
      dimension  idxc(nxc,imax), cmix(nxc,imax), denn(ndn,imax)
     1,          rh(imbx,l),     kprf(imax)
c --- output
                          d i m e n s i o n
     1  tau(imbx,l),  ssa(imbx,l),  asy(imbx,l),  taurs(l)
      data  crt1,crt2 / 30.0, 0.03333 /
      save crt1, crt2
c
c===> ... layer optical depth due to rayleigh scattering
c
      do k=1,l
        taurs(k) = taur(k,ib)
      end do
c
      do k=1,l
      do i=1,ipts
        ssa(i,k) = 0.0e0
        asy(i,k) = 0.0e0
        tau(i,k) = 0.0e0
      end do
      end do
c
      if (kaer .lt. 1) return
c
      do i=1,ipts
c
        kpf = kprf(i)
        do k=1,l
          idom = idm(k,kpf)
          drh = rh(i,k) - 0.5
c
          if (idom .eq. 1) then
c --- 1st domain - mixing layer
            ext1 = 0.0
            sca1 = 0.0
            ssa1 = 0.0
            asf1 = 0.0
            do icmp=1,nxc
              ic = idxc(icmp,i)
              if (ic .gt. ncm1) then
                ic1 = ic - ncm1
                drh1 = exp(abpw(ic1,ib)*drh)
                drh2 = drh * drh
                ex00 = aext(1,ic1,ib) + aext(2,ic1,ib)*drh
     1               + aext(3,ic1,ib)*drh1
                sc00 = bsca(1,ic1,ib) + bsca(2,ic1,ib)*drh
     1               + bsca(3,ic1,ib)*drh1
                ss00 = cssa(1,ic1,ib) + cssa(2,ic1,ib)*drh
     1               + cssa(3,ic1,ib)*drh2
                as00 = dasf(1,ic1,ib) + dasf(2,ic1,ib)*drh
     1               + dasf(3,ic1,ib)*drh2
              else if (ic .gt. 0) then
                ex00 = ext0(ic,ib)
                sc00 = sca0(ic,ib)
                ss00 = ssa0(ic,ib)
                as00 = asf0(ic,ib)
              else
                ex00 = 0.0
                sc00 = 0.0
                ss00 = 0.0
                as00 = 0.0
              end if
              ext1 = ext1 + cmix(icmp,i) * ex00
              sca1 = sca1 + cmix(icmp,i) * sc00
              ssa1 = ssa1 + cmix(icmp,i) * ss00 * ex00
              asf1 = asf1 + cmix(icmp,i) * as00 * sc00
            end do
            ext2 = ext1 * denn(1,i)
            ssa2 = ssa1 / ext1
            asf2 = asf1 / sca1
          else if (idom .eq. 2) then
c --- 2nd domain - mineral transport layers
            ext2 = ext0(6,ib) * denn(2,i)
            ssa2 = ssa0(6,ib)
            asf2 = asf0(6,ib)
          else if (idom .eq. 3) then
c --- 3rd domain - free tropospheric layers
c   1:inso 0.17E-3; 2:soot 0.4; 7:wasO 0.59983; n:730
            drh1 = exp(abpw(1,ib)*drh)
            drh2 = drh * drh
            ex01 = ext0(1,ib)
            sc01 = sca0(1,ib)
            ss01 = ssa0(1,ib)
            as01 = asf0(1,ib)
            ex02 = ext0(2,ib)
            sc02 = sca0(2,ib)
            ss02 = ssa0(2,ib)
            as02 = asf0(2,ib)
            ex03 = aext(1,1,ib) + aext(2,1,ib)*drh + aext(3,1,ib)*drh1
            sc03 = bsca(1,1,ib) + bsca(2,1,ib)*drh + bsca(3,1,ib)*drh1
            ss03 = cssa(1,1,ib) + cssa(2,1,ib)*drh + cssa(3,1,ib)*drh2
            as03 = dasf(1,1,ib) + dasf(2,1,ib)*drh + dasf(3,1,ib)*drh2
            ext1 = 0.17E-3*ex01 + 0.4*ex02 + 0.59983*ex03
            sca1 = 0.17E-3*sc01 + 0.4*sc02 + 0.59983*sc03
            ssa1 = 0.17E-3*ss01*ex01 + 0.4*ss02*ex02 + 0.59983*ss03*ex03
            asf1 = 0.17E-3*as01*sc01 + 0.4*as02*sc02 + 0.59983*as03*sc03
            ext2 = ext1 * 730.0
            ssa2 = ssa1 / ext1
            asf2 = asf1 / sca1
          else if (idom .eq. 4) then
c --- 4th domain - stratospheric layers
            ext2 = estr(ib)
            ssa2 = 0.9
            asf2 = 0.6
          else
c --- upper stratosphere assume no aerosols
            ext2 = 0.0
            ssa2 = 1.0
            asf2 = 0.0
          end if
c
          hd = haer(idom,kpf)
          if (hd .gt. 0.0e0) then
            hd1 = 1.0 / hd
            sig0u = exp(-hz(k)  *hd1)
            sig0l = exp(-hz(k+1)*hd1)
            tau(i,k) = ext2 * hd*(sig0l - sig0u)
          else
            tau(i,k) = ext2 * dz(k)
          end if
          ssa(i,k) = ssa2
          asy(i,k) = asf2
        end do
      end do
c===> ... smooth profile at domain boundaries
      do k=2,l
      do i=1,ipts
        ratio = 1.0e0
        if (tau(i,k) .gt. 0.0e0) ratio = tau(i,k-1) / tau(i,k)
        tt = tau(i,k) + tau(i,k-1)
        if (ratio .gt. crt1) then
          tau(i,k) = 0.2e0 * tt
          tau(i,k-1) = tt - tau(i,k)
        else if (ratio .lt. crt2) then
          tau(i,k) = 0.8e0 * tt
          tau(i,k-1) = tt - tau(i,k)
        end if
      end do
      end do
c
      return
      end
