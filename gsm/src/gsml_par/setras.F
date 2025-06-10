#include <paramodel.h>
      subroutine setras(lev, si, sl, del, cp, rgas, dt, nsphys, fh
     *,                 sig, sgb, prh, prj, hpk, fpk, ods, prns
     *,                 rasal, lm, krmin, krmax, nstrp
     *,                 ncrnd, rannum, afac, ufac)
c
#ifdef INTEL
      USEIFPORT
#endif
c
      parameter (actp=1.7,   facm=1.50)
c
      dimension si(lev+1), sl(lev), del(lev)
      dimension sig(lev+1), prj(lev+1), prh(lev),fpk(lev), hpk(lev)
     *,         sgb(lev),   ods(lev),   rasal(lev), prns(lev/2)
     *,         rannum(200,nsphys)
      dimension ph(15), a(15), ac(15), ad(15), au(15)
#ifdef MAC
      integer, allocatable::nrnd(:)
#endif
c
#ifdef SX6
      integer seed 
      parameter (seed=1.0)
#endif
#ifdef NEC
      integer seed 
      parameter (seed=1.0)
#endif
#ifdef ES
      integer seed
      parameter (seed=1.0)
#endif

c
      common /rasacr/ a, ac, ad, actop
c
      dimension adata(15)
c
      data ph/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0
     *,       550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/
c
      data adata/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677
     *,       0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664
     *,       0.0553, 0.0445, 0.0633/
c
      logical first
      data iseed/0/, first/.true./, fh0/0./
      save iseed, first
c
#ifdef NEC
      real*16 rseed
#endif
#ifdef ES
      real*16 rseed
#endif
#ifdef SX6
      real*16 rseed
#endif
c
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c     set critical workfunction arrays
c
      if (first) then
         do i=1,15
           a(i)=adata(i)
         enddo
         first   = .false.
         actop   = actp*facm
         do 5 l=1,15
         a(l) = a(l)*facm
    5    continue
         do 10 l=2,15
         tem   = ph(l) - ph(l-1)
         au(l) = a(l-1) / tem
         ad(l) = a(l)   / tem
         ac(l) = ph(l)*au(l) - ph(l-1)*ad(l)
         ad(l) = ad(l) - au(l)
   10    continue
         call setes
      endif
c
c      set other parameters
c
         if (fh .eq. fh0) then
            iseed = iseed + 1
         else
            iseed = fh*3600.0 + 0.0001
            iseed = max(iseed, 1)
            fh0   = fh
         endif
c
         rkap = rgas / cp
c
         nstrp = 1
         do l=2, lev
            if (sl(l) .gt. 0.95) nstrp = nstrp + 1
         enddo
c        nstrp = 1
         do l=1,nstrp+1
           prj(l) = (si(l)/1000.0) ** rkap
         enddo
c
         afac = 0.0
         do l=1,nstrp
            prns(l) =  (si(l)*prj(l)-si(l+1)*prj(l+1)) /
     *                ((si(l) - si(l+1)) * (rkap+1.0) )
            afac    = afac + prns(l) * del(l)
         enddo
         afac = afac / (si(1) - si(nstrp+1))
c
         lm   = lev - nstrp + 1
         lmm1 = lm - 1
         do l=1,lm
            sig(l) = si(lev-l+2)
            prj(l) = (sig(l)/1000.0) ** rkap
         enddo
         sig(lm+1) = si(1)
         prj(lm+1) = (sig(lm+1)/1000.0) ** rkap
c
         krmin = 1
         krmax = lmm1
         do l=1,lmm1
            sgb(l) = sl(lev-l+1)
            if (sgb(l) .le. 0.067) krmin = l
            if (sgb(l) .le. 0.76) krmax = l
         enddo
         sgb(lm) = 0.5 * (sig(lm) + sig(lm+1))
c
         ncrnd   = 48 * (dt/3600) + 0.50001
         ncrnd   = max(ncrnd, 1)
         rasalf  = 0.30
c
         if (ncrnd*3 .gt. 100) then
            print *,' dimension of rannum too small -- job terminated'
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
         endif
#define DEFAULT
#ifdef RANF
#undef DEFAULT
         call ranset(iseed)
#endif
#ifdef DEFAULT
#ifdef NEC
c         call xsrandom(iseed)
#else
#ifdef ES
c         call xsrandom(iseed)
#else
#ifdef SX6
c         call xsrandom(iseed)
#else
#ifdef sun
          r=rand(i)
#else
#ifndef MAC
         call srand(iseed)
#endif
#endif
#endif
#endif
#endif
#endif

#ifndef MAC
         do nn=1,nsphys
           do i=1,ncrnd*3
#define DEFAULT
#ifdef RANF
#undef DEFAULT
             rannum(i,nn) = ranf()
#endif
#ifdef DEFAULT
#ifdef NEC
             rseed=iseed        
             rannum(i,nn) = frand(rseed)
#else
#ifdef ES
             rseed=iseed
             rannum(i,nn) = frand(rseed)
#else
#ifdef SX6
             rseed=iseed
             rannum(i,nn) = frand(rseed)
#else
             rannum(i,nn) = rand()
#endif
#endif
#endif
#endif
           enddo
         enddo
#endif

#ifdef MAC
         call random_seed(size=krsize)
         allocate ( nrnd(krsize) )
         nrnd=iseed
         call random_seed(put=nrnd)
         call random_number(rannum)
         deallocate ( nrnd )
#endif
c
         critsig=0.95
c
         do l=1,lm
            prh(l  ) = (sig(l+1)*prj(l+1)-sig(l)*prj(l)) /
     *               ( (sig(l+1)-sig(l))*(rkap+1.0) )
            fpk(l  ) = cp * (prj(l+1)-prj(l))
            hpk(l  ) = cp * (prj(l+1)-prh(l))
            ods(l) = sig(l+1) - sig(l)
            rasal(l) = rasalf
c
            if(sig(l).gt.critsig) rasal(l)=0.
c
         enddo
c
         do l=1,lm
           ods(l) = 1.0 / ods(l)
         enddo
         afac = 1.0 / (afac - prh(lmm1))
         ufac = 0.0
c
      return
      end
