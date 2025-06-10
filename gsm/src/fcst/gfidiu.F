#define CP cp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NTOTAL ntotal
#define NTRAC ntrac
#define OMEGA omega
#define RD rd
#define RERTH rerth

#include <define.h>
      subroutine gfidiu(lons2,
#ifdef SCM
     * dg,tg,rqg,qg,
     1 del,sl,rdel2,
#ifdef DBG
     2 dtdt,drdt,vadt,adbc)
#else
     2 dtdt,drdt)
#endif
#else			/* not SCM */
     * dg,tg,zg,ug,vg,rqg,dphi,dlam,
     1 rcl,del,rdel2,ci,p1,p2,h1,h2,tov,spdmax,
     2 dtdf,dtdl,drdf,drdl,dudl,dvdl,dudf,dvdf,
#ifndef NISLQ
     2 dqdt,dtdt,drdt,dudt,dvdt)
#else
     2 dqdt,dtdt,drdt,dudt,dvdt,dot)
#endif
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c input variables
c
#include <paramodel.h>
#ifdef MP
#define LONF22S lonf22p_
#else
#define LONF22S lonf22_
#endif
      dimension
     1  dg(LONF22S,levs_),tg(LONF22S,levs_), zg(LONF22S,levs_),
     2  ug(LONF22S,levs_),vg(LONF22S,levs_),
     3 rqg(LONF22S,levs_,ntotal_)
#ifdef SCM
#ifdef DBG
     + ,vadt(LONF22S,levs_),adbc(LONF22S,levs_)
#endif
     +  ,qg(LONF22S), sl(levs_)
#else
     3  ,dphi(LONF22S),dlam(LONF22S)
#endif
c
#ifndef SCM
      dimension
     1  dtdf(LONF22S,levs_),dtdl(LONF22S,levs_),
     1  drdf(LONF22S,levs_,ntotal_),
     1  drdl(LONF22S,levs_,ntotal_),
     1  dudl(LONF22S,levs_),dvdl(LONF22S,levs_),
     1  dudf(LONF22S,levs_),dvdf(LONF22S,levs_)
#endif
c output variables
c
      dimension spdmax(levs_),
#ifndef SCM
     1  dudt(LONF22S,levs_),dvdt(LONF22S,levs_),
#endif
     1  dtdt(LONF22S,levs_),drdt(LONF22S,levs_,ntotal_)
#ifndef SCM
     1  ,dqdt(LONF22S)
#endif
c
c constant arrays
c
       dimension
     1 del(levs_),rdel2(levs_),
     2 ci(levp1_),tov(levs_)
#ifndef SCM
     3 ,p1(levs_),p2(levs_),h1(levs_),h2(levs_)
#endif
c
c local variables
c
      dimension
     1  cg (LONF22S,levs_), db(LONF22S,levs_),cb(LONF22S,levs_),
     2  dot(LONF22S,levp1_),dup(LONF22S,levs_),dvp(LONF22S,levs_),
     3  dum(LONF22S,levs_ ),dvm(LONF22S,levs_), ek(LONF22S,levs_),
     4  rmu(levs_ ),rnu(levs_),rho(levs_),si(levp1_)
c
      rk= rd_ /cp_
#ifndef SCM
      sinra=sqrt(1.-1./rcl)
      fnor=2.*omega_*sinra
      fsou=-fnor
      sinra=sinra/rerth_
#endif
c
      si(1)=1.0
      do 4 k=1,levs_
      si(k+1)=si(k)-del(k)
4     continue
c
      do 1 k=1,levm1_
      rho(k)=alog(si(k)/si(k+1))
1     continue
      rho(levs_)=0.
c
      do 2 k=1,levs_
      rmu(k)=1.-si(k+1)*rho(k)/del(k)
2     continue
c
      do 3 k=1,levm1_
      rnu(k+1)=-1.+si(k)*rho(k)/del(k)
3     continue
      rnu(1)=0.
c
#ifndef SCM
      do 1234 k=1,levs_
      spdmax(k)=0.
1234  continue
      rcl2=.5e0*rcl
c
      do 140 k=1,levs_
      do 140 j=1,lons2
      ek(j,k)=(ug(j,k)*ug(j,k)+vg(j,k)*vg(j,k))*rcl
  140 continue
c
      do 10 k=1,levs_
      do 10 j=1,lons2
      if (ek(j,k) .gt. spdmax(k))  spdmax(k)=ek(j,k)
   10 continue
c
c     compute c=v(true)*del(ln(ps)).divide by cos for del, cos for v
c
      do 150 j=1,lons2
      dphi(j)=dphi(j)*rcl
      dlam(j)=dlam(j)*rcl
  150 continue
      do 180 k=1,levs_
      do 180 j=1,lons2
      cg(j,k)=ug(j,k)*dlam(j)+vg(j,k)*dphi(j)
  180 continue
c
      do 190 j=1,lons2
      db(j,1)=del(1)*dg(j,1)
      cb(j,1)=del(1)*cg(j,1)
  190 continue
      do 210 k=1,levm1_
      do 210 j=1,lons2
      db(j,k+1)=db(j,k)+del(k+1)*dg(j,k+1)
      cb(j,k+1)=cb(j,k)+del(k+1)*cg(j,k+1)
  210 continue
c
c   store integral of cg in dlax
c
      do 220 j=1,lons2
      dqdt(j)= -cb(j,levs_)
  220 continue
c
c   sigma dot computed only at interior interfaces.
c
      do 230 j=1,lons2
      dot(j,1)=0.e0
      dvm(j,1)=0.e0
      dum(j,1)=0.e0
      dot(j,levp1_)=0.e0
      dvp(j,levs_ )=0.e0
      dup(j,levs_ )=0.e0
  230 continue
c
      do 240 k=1,levm1_
      do 240 j=1,lons2
      dot(j,k+1)=dot(j,k)+
     1                 del(k)*(db(j,levs_)+cb(j,levs_)-
     2                 dg(j,k)-cg(j,k))
  240 continue
c
c
c
      do 260 k=1,levm1_
      do 260 j=1,lons2
      dvp(j,k  )=vg(j,k+1)-vg(j,k)
      dup(j,k  )=ug(j,k+1)-ug(j,k)
      dvm(j,k+1)=vg(j,k+1)-vg(j,k)
      dum(j,k+1)=ug(j,k+1)-ug(j,k)
  260 continue
      do j=1,lons2
       dphi(j)=dphi(j)/rcl
       dlam(j)=dlam(j)/rcl
      enddo
      do k=1,levs_
       do j=1,lons2
        dudt(j,k)=dudt(j,k)
     &            -ug(j,k)*dudl(j,k)-vg(j,k)*dudf(j,k)
     1 -rdel2(k)*(dot(j,k+1)*dup(j,k)+dot(j,k)*dum(j,k))
     2 -rd_*tg(j,k)*dlam(j)
c
        dvdt(j,k)=dvdt(j,k)
     &            -ug(j,k)*dvdl(j,k)-vg(j,k)*dvdf(j,k)
     1 -rdel2(k)*(dot(j,k+1)*dvp(j,k)+dot(j,k)*dvm(j,k))
     2 -rd_*tg(j,k)*dphi(j)
       enddo
      enddo
c
      lons=lons2/2
      do k=1,levs_
       do j=1,lons
        dudt(j,k)=dudt(j,k)+vg(j,k)*fnor
        dudt(j+lons,k)=dudt(j+lons,k)+vg(j+lons,k)*fsou

        dvdt(j,k)=dvdt(j,k)-ug(j,k)*fnor
     1 -sinra*ek(j,k)
        dvdt(j+lons,k)=dvdt(j+lons,k)-ug(j+lons,k)*fsou
     1 +sinra*ek(j+lons,k)
       enddo
      enddo
      do k=1,levs_
       do j=1,lons2
        dudt(j,k)=dudt(j,k)*rcl
        dvdt(j,k)=dvdt(j,k)*rcl
       enddo
      enddo
c
#endif
c
#ifdef SCM
      do j=1,lons2
c
       dot(j,1)=0.e0
       dot(j,levp1_)=0.e0
       do k=1,levs_-1
         psfc = exp(qg(j))  ! cb
         dot(j,k+1) = (dg(j,k)+dg(j,k+1))/(2.*psfc) ! cb/sec / cb
cyh         if (j.eq.1)
cyh     +      print 1963,j,k,dg(j,k),dg(j,k+1),psfc,dot(j,k+1)
cyh1963     format('DOT',2i5,4e13.5)
       enddo
c
c ... dT/dt ...
c
       k = 1
       difm = 0.0
       sdotm = 0.0
       difp=tg(j,k+1) - tg(j,k)
       sdotp = dot(j,k+1)
       delk = si(k+1) - 1.0
       vadvT = -1.*(sdotp*difp + sdotm*difm)/delk
       pres_k = sl(k) * exp(qg(j))
       omega_k = dg(j,k)
       adiabatic = rk*tg(j,k)*omega_k/pres_k
       dtdt(j,k) = vadvT + adiabatic
#ifdef DBG
       vadt(j,k) = vadvT
       adbc(j,k) = adiabatic
       if (j.eq.1) print 1967, j,k,delk,sdotp,difp,sdotm,difm,
     +             vadvT,pres_k,omega_k,adiabatic,dtdt(j,k)
#endif
c
       do k=2,levs_-1
         difm = tg(j,k) - tg(j,k-1)
         sdotm = dot(j,k-1)
         difp = tg(j,k+1) - tg(j,k)
         sdotp = dot(j,k+1)
         delk = si(k+1) - si(k-1)
         vadvT = -1.*(sdotp*difp + sdotm*difm)/delk
         pres_k = sl(k) * exp(qg(j))
         omega_k = dg(j,k)
         adiabatic = rk*tg(j,k)*omega_k/pres_k
         dtdt(j,k) = vadvT + adiabatic
#ifdef DBG
         vadt(j,k) = vadvT
         adbc(j,k) = adiabatic
         if (j.eq.1) print 1967, j,k,delk,sdotp,difp,sdotm,difm,
     +               vadvT,pres_k,omega_k,adiabatic,dtdt(j,k)
1967     format(2i5,10e13.5)
#endif
       enddo
c
       k = levs_
       difm = tg(j,k) - tg(j,k-1)
       sdotm = dot(j,k-1)
       difp =0.0
       sdotp = 0.0
       delk = 0.0 - si(k-1)
       vadvT = -1.*(sdotp*difp + sdotm*difm)/delk
       pres_k = sl(k) * exp(qg(j))
       omega_k = dg(j,k)
       adiabatic = rk*tg(j,k)*omega_k/pres_k
#ifdef DBG
       vadt(j,k) = vadvT
       adbc(j,k) = adiabatic
       dtdt(j,k) = vadvT + adiabatic
       if (j.eq.1) print 1967, j,k,delk,sdotp,difp,sdotm,difm,
     +             vadvT,pres_k,omega_k,adiabatic,dtdt(j,k)
#endif
c
c ... dr/dt ...
c
       do n = 1, ntotal_
c
        k = 1
        difm = 0.0
        sdotm = 0.0
        difp = rqg(j,k+1,n) - rqg(j,k,n)
        sdotp = dot(j,k+1)
        delk = si(k+1) - 1.0
        drdt(j,k,n) = -1.*(sdotp*difp + sdotm*difm)/delk
c
        do k = 2, levs_-1
         difm = rqg(j,k,n) - rqg(j,k-1,n)
         sdotm = dot(j,k-1)
         difp = rqg(j,k+1,n) - rqg(j,k,n)
         sdotp = dot(j,k+1)
         delk = si(k+1) - si(k-1)
         drdt(j,k,n) = -1.*(sdotp*difp + sdotm*difm)/delk
         if (j.eq.1) print 1966, j,k,delk,sdotp,difp,sdotm,difm,
     +               drdt(j,k,n)
1966     format(2i5,6e13.5)
        enddo
c
        k = levs_
        difm = rqg(j,k,n) - rqg(j,k-1,n)
        sdotm = dot(j,k-1)
        difp = 0.0
        sdotp = 0.0
        delk = 0.0 - si(k-1)
        drdt(j,k,n) = -1.*(sdotp*difp + sdotm*difm)/delk
c
       enddo
c
      enddo
#else			/* SCM */
      do 280 k=1,levm1_
      do 280 j=1,lons2
cecmwf:
      dup(j,k  )=tg(j,k+1)+tov(k+1)-tg(j,k)-tov(k)+2.*rk*rnu(k+1)*
     .(tg(j,k)+tov(k))
cecmwf:
      dum(j,k+1)=tg(j,k+1)+tov(k+1)-tg(j,k)-tov(k)+2.*rk*rmu(k+1)*
     .(tg(j,k+1)+tov(k+1))
  280 continue
c
c
      do k=1,levs_
       do j=1,lons2
        dtdt(j,k)=dtdt(j,k)
     &  -ug(j,k)*dtdl(j,k)-vg(j,k)*dtdf(j,k)
     1 -rdel2(k)*(dot(j,k+1)*dup(j,k)+dot(j,k)*dum(j,k))
       enddo
      enddo
c
      do k=1,levs_
       do j=1,lons2
        dtdt(j,k)=dtdt(j,k)
     1  +rk*(tov(k)+tg(j,k))*(cg(j,k)-cb(j,levs_)-db(j,levs_))
       enddo
      enddo
c
#ifndef NISLQ
      do 330 n=1,ntotal_
      do 300 k=1,levm1_
      do 300 j=1,lons2
      dup(j,k  )=rqg(j,k+1,n)-rqg(j,k,n)
      dum(j,k+1)=rqg(j,k+1,n)-rqg(j,k,n)
  300 continue
      do 310 j=1,lons2
      dup(j,levs_)=0.e0
  310 continue
      do 320 k=1,levs_
      do 320 j=1,lons2
      drdt(j,k,n)=drdt(j,k,n)
     & -ug(j,k)*drdl(j,k,n)-vg(j,k)*drdf(j,k,n)
     1 -rdel2(k)*(dot(j,k+1)*dup(j,k)+dot(j,k)*dum(j,k))
  320 continue
  330 continue
#endif
#endif
c
      return
      end
