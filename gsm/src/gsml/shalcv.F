#define CP cp
#define G g
#define HVAP hvap
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
#define RD rd
      subroutine shalcv(ims2,imx2,kmx,kmx2,ndiff,dt,
     $                  kuo,ps,q,t,
     *                  del,sl,slk,si)
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpkap)
#endif
c-----------------------------------------------------------------------
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:   shalcv       computes shallow convective heating and moi
c   prgmmr: peter caplan     org: w/nmc23    date: 91-03-19
c
c abstract: sub-grid-scale shallow convective cloud parameterization.
c   this routine computes the effects of shallow convection
c   based on tiedtke (1984), ecmwf workshop on convection in
c   large-scale numerical models.
c   tapered k profile in cloud developed by caplan and long.
c   originally coded by r. kistler and p. caplan, converted to standard
c   fortran for cray by h.juang. h. pan modified it and kuo91 to
c   do mstadb after compression and the rest of the computation
c   in compressed arrays. tidy up of moist processes by m. iredell.
c
c program history log:
c   91-03-19  hua-lu pan
c   91-05-07  iredell             arguments changed, tridi2 split off
c   02-09-01  yifeng cui          rsm mpi added
c   06-01-25  shimpo              add tracer like moninp.F
c
c usage:    call shalcv(im,km,dt,del,si,sl,slk,kuo,ps,q,t)
c
c   input argument list:
c     im       - integer number of points
c     km       - integer number of levels
c     ndiff    - number of passive variables including moisture
c     dt       - real time step in seconds
c     del      - real (km) sigma layer thickness
c     sl       - real (km) sigma values
c     slk      - real (km) sigma values to the kappa
c     ps       - real (im) surface pressure in kilopascals (cb)
c     q        - real (im,km) current specific humidity in kg/kg
c     t        - real (im,km) current temperature in kelvin
c
c   output argument list:
c     q        - real (im,km) adjusted specific humidity in kg/kg
c     t        - real (im,km) adjusted temperature in kelvin
c
c subprograms called:
c   mstadb   - computes moist adiabat and returns cloud values
c   tridin   - solves tridiagonal matrix problem for t, q, and tracers
cx  tridi2   - solves tridiagonal matrix problem
c
c remarks: nonstandard automatic arrays are used.
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
c
c$$$
#include <paramodel.h>
      dimension del(imx2,kmx),si(imx2,kmx+1),sl(imx2,kmx),
     $          slk(imx2,kmx),kuo(imx2),ps(imx2),
     &          q(imx2,kmx*ndiff),t(imx2,kmx)
c  physical parameters
      parameter(g=g_,rd=rd_,
     &          cp=cp_,hvap=hvap_)
      parameter(gocp=g/cp)
c  bounds of parcel origin
      parameter(kliftl=2,kliftu=2)
c  local variables and arrays
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
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      logical lshc
      dimension lshc(ILOTS)
      dimension index2(ILOTS)
      dimension klcl(ILOTS),kbot(ILOTS),
     *          ktop(ILOTS)
      dimension ps2(ILOTS),
     &          q2(ILOTS*levs_*ndiff),
     *          t2(ILOTS*levs_),
     &          al(ILOTS*(levs_-1)),
     *          ad(ILOTS*levs_),
     &          au(ILOTS*(levs_-1))
      dimension sl2(ILOTS*levs_)
      dimension slk2(ILOTS*levs_)
      dimension del2(ILOTS*levs_)
      dimension si2(ILOTS*levs_+1)
c hmhj
      im=ims2
      km=kmx2
c-----------------------------------------------------------------------
c  compress fields to points with no deep convection
c  and moist static instability.
      do i=1,im
        lshc(i)=.false.
      enddo
      do k=1,km-1
        do i=1,im
          if(kuo(i).eq.0) then
            eldq=hvap*(q(i,k)-q(i,k+1))
            cpdt=cp*(t(i,k)-t(i,k+1))
            rtdls=(sl(i,k)-sl(i,k+1))/si(i,k+1)*rd*0.5*(t(i,k)+t(i,k+1))
            dmse=eldq+cpdt-rtdls
            lshc(i)=lshc(i).or.dmse.gt.0.
          endif
        enddo
      enddo
      n2=0
      do i=1,im
        if(lshc(i)) then
          n2=n2+1
          index2(n2)=i
        endif
      enddo
      if(n2.eq.0) return
      do i=1,n2
        ps2(i)=ps(index2(i))
      enddo
      do k=1,km
#ifdef CRAY_THREAD
cfpp$ select(vector)
#endif
        do i=1,n2
          ik=(k-1)*n2+i
          q2(ik)=q(index2(i),k)
          t2(ik)=t(index2(i),k)
          sl2(ik)=sl(index2(i),k)
          slk2(ik)=slk(index2(i),k)
          si2(ik)=si(index2(i),k)
          del2(ik)=del(index2(i),k)
        enddo
      enddo
      do i=1,n2
        k=km+1
        ik=(k-1)*n2+i
        si2(ik)=si(index2(i),k)
      enddo
         
      if(ndiff.ge.2) then
        do kk = 2, ndiff
          is = (kk-1) * km
          is2= (kk-1) * kmx
          do k = 1, km
            do i = 1, n2
              ik=(k+is-1)*n2+i
              q2(ik)=q(index2(i),k+is2)  ! tracer(s)
            enddo
          enddo
        enddo
      endif
c-----------------------------------------------------------------------
c  compute moist adiabat and determine limits of shallow convection.
c  check for moist static instability again within cloud.
      call mstadb(n2,km-1,kliftl,kliftu,sl2,slk2,ps2,t2,q2,
     &            klcl,kbot,ktop,al,au)
      do i=1,n2
        kbot(i)=klcl(i)-1
        ktop(i)=ktop(i)+1
        lshc(i)=.false.
      enddo
      do k=1,km-1
        do i=1,n2
          if(k.ge.kbot(i).and.k.lt.ktop(i)) then
            ik=(k-1)*n2+i
            iku=k*n2+i
            eldq=hvap*(q2(ik)-q2(iku))
            cpdt=cp*(t2(ik)-t2(iku))
            rtdls=(sl2(ik)-sl2(iku))/si2(iku)
     $           *rd*0.5*(t2(ik)+t2(iku))
            dmse=eldq+cpdt-rtdls
            lshc(i)=lshc(i).or.dmse.gt.0.
            au(ik)=g/rtdls
          endif
        enddo
      enddo
      k1=km+1
      k2=0
      do i=1,n2
        if(.not.lshc(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        k1=min(k1,kbot(i))
        k2=max(k2,ktop(i))
      enddo
      kt=k2-k1+1
      if(kt.lt.2) return
c     use all levels
      k1=1
      k2=km
      kt=k2-k1+1
c-----------------------------------------------------------------------
c  set eddy viscosity coefficient cku at sigma interfaces.
c  compute diagonals and rhs for tridiagonal matrix solver.
c  expand final fields.
      do i=1,n2
        ik=(k1-1)*n2+i
        ad(ik)=1.
      enddo
      do k=k1,k2-1
        do i=1,n2
          ik=(k-1)*n2+i
          iku=k*n2+i
          dtodsl=2.*dt/del2(ik)
          dtodsu=2.*dt/del2(iku)
          dsig=sl2(ik)-sl2(iku)
          if(k.eq.kbot(i)) then
            ck=1.5
          elseif(k.eq.ktop(i)-1) then
            ck=1.
          elseif(k.eq.ktop(i)-2) then
            ck=3.
          elseif(k.gt.kbot(i).and.k.lt.ktop(i)-2) then
            ck=5.
          else
            ck=0.
          endif
          dsdz1=ck*dsig*au(ik)*gocp
          dsdz2=ck*dsig*au(ik)*au(ik)
          au(ik)=-dtodsl*dsdz2
          al(ik)=-dtodsu*dsdz2
          ad(ik)=ad(ik)-au(ik)
          ad(iku)=1.-al(ik)
          t2(ik)=t2(ik)+dtodsl*dsdz1
          t2(iku)=t2(iku)-dtodsu*dsdz1
        enddo
      enddo
      ik1=(k1-1)*n2+1
      call tridin(n2,n2,kt,ndiff,al(ik1),ad(ik1),
     &            au(ik1),t2(ik1),q2(ik1),
     &            au(ik1),t2(ik1),q2(ik1))
      do k=k1,k2
#ifdef CRAY_THREAD
cfpp$ select(vector)
#endif
        do i=1,n2
          ik=(k-1)*n2+i
          q(index2(i),k)=q2(ik)
          t(index2(i),k)=t2(ik)
        enddo
      enddo
      if(ndiff.ge.2) then
        do kk = 2, ndiff
          is = (kk-1) * km
          is2= (kk-1) * kmx
          do k = k1, k2
            do i = 1, n2
              ik=(k+is-1)*n2+i
              q(index2(i),k+is2)=q2(ik)  ! tracer(s)
            enddo
          enddo
        enddo
      endif
c-----------------------------------------------------------------------
      return
      end
