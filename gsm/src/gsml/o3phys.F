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
#define PI pi
      subroutine o3phys (im,ix,kx,deltim,pis,ozi,ozo,xlat,sl)
#include <paramodel.h>
c
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
c
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
c
      parameter (jo3=18,ko3=46,blat=-85.0,dphi=10.0)
c
      dimension prdin(jo3,ko3),disin(jo3,ko3)
      common /o3pcom/prdin,disin
      dimension po3(ko3), slog(ILOTS,levs_)
     *,         prdout(ILOTS,ko3),
     *          disout(ILOTS,ko3)
      dimension ozi(ix,kx), ozo(ix,kx),xlat(ix),sl(ix,kx),pis(ix)
      dimension wk1(ILOTS), wk2(ILOTS), 
     *          wk3(ILOTS), 
     *          wk4(ILOTS)
     *,         wk5(ILOTS), wkp(ILOTS,2), 
     *          wkd(ILOTS,2)
     *,         prod(ILOTS),dist(ILOTS)
     *,         ddy(ILOTS),jindx1(ILOTS),
     *          jindx2(ILOTS)
c
cfpp$ noconcur r
c
      dt=deltim * 2
c
        dp = 0.2844
        do n = 1,ko3
           po3(n) = 101.30 * exp(-dp*(float(n)-0.5))
           po3(n) = log(po3(n))
        enddo
        do i=1,im
        do k = 1,levs_
           slog(i,k) = log(sl(i,k))
        enddo
        enddo
        dphii = 1.0 / dphi
        elat  = blat + (jo3-1)*dphi
c
        do i = 1,im
          xdeg = xlat(i)*180./pi_
          ddy(i) = (xdeg - blat) * dphii + 1.0
          jindx1(i)  = ddy(i)
          jindx2(i)  = jindx1(i) + 1
          ddy(i)     = ddy(i) - jindx1(i)
c         if(xdeg .le. blat) then
c           jindx1(i) = 1
c           jindx2(i) = 1
c         endif
c         if (xdeg .ge. elat) then
c           jindx1(i) = jo3
c           jindx2(i) = jo3
c         endif
          jindx1(i) = min(max (jindx1(i),1),jo3)
          jindx2(i) = min(max (jindx2(i),1),jo3)
c         print*,' o3phys i xdeg ddy j1 j2 ',i,xdeg,ddy(i),
c    1             jindx1(i),jindx2(i)
        enddo
c
      do l = 1,ko3
        do i = 1,im
          j1 = jindx1(i)
          j2 = jindx2(i)
          tem= 1.0 - ddy(i)
          prdout(i,l) = (tem * prdin(j1,l) + ddy(i) * prdin(j2,l))
     *                * 1.655
          disout(i,l) = tem * disin(j1,l) + ddy(i) * disin(j2,l)
        enddo
      enddo
c
      do i = 1,im
        wkp(i,1)     = prdout(i,1)
        wkp(i,2)     = prdout(i,ko3)
        wkd(i,1)     = disout(i,1)
        wkd(i,2)     = disout(i,ko3)
      enddo
c
      psmin =  1.0e10
      psmax = -1.0e10
      do i = 1,im
         psmin = min(psmin,pis(i))
         psmax = max(psmax,pis(i))
      enddo
c
      do l = 1,levs_
c
         pmin = psmin + slog(1,l)
         pmax = psmax + slog(1,l)
         do i=1,im
           pmin = min(pmin,psmin + slog(i,l))
           pmax = max(pmax,psmax + slog(i,l))
         enddo
         kmax = 1
         kmin = 1
         do k = 1,ko3-1
            if (pmin .lt. po3(k)) kmax = k
            if (pmax .lt. po3(k)) kmin = k
         enddo
c        kmax = min(kmax+1,ko3-1)
c
c        print *,' pmin=',pmin,' pmax=',pmax,' kmin=',kmin
c    *, ' kmax=',kmax,' po3=',po3(kmin),po3(kmax),' l=',l
c
         do i = 1,im
            wk1(i) = pis(i) + slog(i,l)
            prod(i) = 0.0
            dist(i) = 0.0
         enddo
c
c        tem = 10.0 * exp(wk1(1))
c        print *,' l=',l,' wk1=',wk1(1),' tem=',tem
c

         do k = kmin,kmax
            do i = 1,im
              wk2(i)     = prdout(i,k)
              wk3(i)     = prdout(i,k+1)
              wk4(i)     = disout(i,k)
              wk5(i)     = disout(i,k+1)
            enddo
            temp = 1.0 / (po3(k) - po3(k+1))
            do i = 1,im
               if (wk1(i) .lt. po3(k) .and. wk1(i) .ge. po3(k+1)) then
                  tem      = (wk1(i) - po3(k+1)) * temp
                  prod(i)  = tem * wk2(i) + (1.0-tem) * wk3(i)
                  dist(i)  = tem * wk4(i) + (1.0-tem) * wk5(i)
               endif
            enddo
         enddo
c
         do i = 1,im
            if (wk1(i) .lt. po3(ko3)) then
               prod(i) = wkp(i,2)
               dist(i) = wkd(i,2)
            endif
            if (wk1(i) .ge. po3(1)) then
               prod(i) = wkp(i,1)
               dist(i) = wkd(i,1)
            endif
            ozib = ozi(i,l)  ! no filling
            ozo(i,l) = (ozib + prod(i)*dt) / (1.0 + dist(i)*dt)
         enddo
c
      enddo
c
      return
      end
