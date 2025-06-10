#include <define.h>
      subroutine rasv (len, lenc,   k,   dt
     *,                ncrnd, krmin, krmax, frac, rasal, botop
     *,                cp,  alhl, grav, crmsfx, rannum
     *,                sgo, onbdsg, prh, prj, fpk, hpk, sge
     *,                pii, pik,    pki,  gam
     *,                uoi, voi,    poi, qoi
     *,                pcu, cln,    q1,  q2
     *,                cd
#ifdef ISOTOPE
     *,                qoi1,qoi2,pcui1,pcui2,cndi1,cndi2,cnd
     *,                eqfo,eqfd
#endif
     *,                cwm,cmf,detr)
c
c*********************************************************************
c*********************************************************************
c******************** relaxed arakawa-schubert ***********************
c************************ parameterization ***************************
c********************* plug compatible version  **********************
c*********************  for sigma=p/ps models  ***********************
c************************* 04 august 1992 ****************************
c*********************************************************************
c************************** developed by *****************************
c**************************              *****************************
c************************ shrinivas moorthi **************************
c************************       and         **************************
c************************  max j. suarez *****************************
c************************                *****************************
c******************** laboratory for atmospheres *********************
c****************** nasa/gsfc, greenbelt, md 20771 *******************
c*********************************************************************
c*********************************************************************
c
      parameter (icm=100, one=1.0, daylen=86400.0, rknob=10.0)
      parameter (ibmax=1)
c
      dimension pii(len),   pik(len),   pki(len), rasal(k-1)
     *,         uoi(len,k), voi(len,k), poi(len,k), qoi(len,k)
     *,         cd(len)
#ifdef ISOTOPE
      dimension qoi1(len,k),qoi2(len,k),
     $          pcui1(len),pcui2(len),
     $          cndi1(len),cndi2(len),cnd(len)
      dimension eqfo(len),eqfd(len)
      integer kk
      real    ff
#endif
c
      dimension q1(len,k), q2(len,k),  cln(len,k), pcu(len)
c
      real  gam(len,k)
      real  rannum(ncrnd*3)
c
      dimension sgo(k),    prj(k),    prh(k),    fpk(k),    hpk(k)
     *,         sge(k),    eps(icm),  del(icm),  sdu(icm),  sdd(icm)
     *,         epss(icm)
      dimension onbdsg(k), ic(icm)
      dimension rhfacl(ibmax), rhfacs(ibmax)
      logical   botop
      real      cwm(len,k)

      real      cmf(len,k)
      real      detr(len,k)
c
      km1   = k  - 1
      albycp = alhl / cp
c
c-clm rhfacs(1) = 0.85
c
      rhfacs(1) = 0.70
c     rhfacs(2) = 0.70
c     rhfacs(3) = 0.70
c
c-clm rhfacl(1) = 0.85
c
      rhfacl(1) = 0.70
c     rhfacl(2) = 0.70
c     rhfacl(3) = 0.70
c
      do l=1,km1
         eps(l)   = (prj(l+1)-prh(l  ))/prh(l  )
         del(l+1) = (prh(l+1)-prj(l+1))/prh(l+1)
c
         tem = cp / (prh(l+1)-prh(l))
         sdu(l+1) = prh(l+1) * (prj(l+1)-prh(l  )) * tem
         sdd(l  ) = prh(l  ) * (prh(l+1)-prj(l+1)) * tem
      enddo
c
      do ibl=1,ibmax
         kpbl1 = k - ibl
         kpbl  = kpbl1 + 1
         kcr   = min(kpbl1,krmax)
         kfx   = kpbl1 - kcr
         ncmx  = kfx + ncrnd
c
         if (kfx .gt. 0) then
            if (botop) then
               do 20 nc=1,kfx
               ic(nc) = kpbl - nc
   20         continue
            else
               do 25 nc=kfx,1,-1
               ic(nc) = kpbl - nc
   25          continue
            endif
         endif
c
c   warning!! if ibl> 3, this wont work!!!
c   ----------------------------------------
c   rannum dimension not sufficient!
c
         if (ncrnd .gt. 0) then
            ii = (ibl-1)*ncrnd
            do 30 i=1,ncrnd
            irnd = (rannum(i+ii)-0.0005)*(kcr-krmin+1)
            ic(kfx+i) = irnd + krmin
   30       continue
         endif
c
         fracs  = (sge(kpbl) - sge(kpbl1)) * frac
c
         do 100 nc=1,ncmx
c
         ib = ic(nc)
c
          do i=ib,k
             epss(i) = eps(i)
          enddo
          epss(ib) = eps(ib) * 0.5
c
          call cloud(len, lenc, kpbl-ib+1, k-ib+1
     *,      dt, rasal(ib), fracs, rhfacl(ibl), rhfacs(ibl)
     *,      cp,  alhl, grav, crmsfx, prj(ib)
     *,      cd
     *,      sgo(ib),   onbdsg(ib), prh(ib),   fpk(ib), hpk(ib), sge(ib)
     *,      epss(ib),  del(ib),    sdu(ib),   sdd(ib)
     *,      pii,       pik,     pki, gam(1,ib)
     *,      uoi(1,ib), voi(1,ib),  poi(1,ib), qoi(1,ib)
     *,      pcu,       cln(1,ib),  q1(1,ib),  q2(1,ib)
#ifdef ISOTOPE
     *,      qoi1(1,ib),qoi2(1,ib),pcui1,pcui2,cndi1,cndi2,cnd
     *,      eqfo,eqfd
#endif
     *,      cwm(1,ib) ,cmf(1,ib), detr(1,ib))
c
  100    continue
      enddo
c
      return
      end
