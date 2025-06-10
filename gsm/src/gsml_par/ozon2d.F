#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define PI pi

      subroutine ozon2d(ims2,imx2,kmx,
     1                  sigl,qo3,sfcp,xlat,rsin1,rcos1,rcos2)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c ...  compute model lyr o3 profile from the original gfdl data
#include <paramodel.h>
      parameter (nl=81,nlp1=nl+1,lngth=37*nl)
      dimension sigl(imx2,kmx)
      dimension qo3(imx2,kmx),sfcp(imx2),xlat(imx2)
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
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      dimension jjrow(ILOTS),tthan(ILOTS)
      dimension qo3o3(ILOTS,nl),psm(ILOTS)
c...    **************************************************************
c--   seasonal climatologies of o3 from o3intn
      common /seaso3/
c-       ...winter....  ...spring....  ...summer....  ....fall.....
     1   dduo3n(37,nl), ddo3n2(37,nl), ddo3n3(37,nl), ddo3n4(37,nl)
     2  , prgfdl(nl)
c...    **************************************************************
c...       begin here .....
      rndg = 180./pi_
      do 10 i=1,ims2
       th2=0.2e0*xlat(i)*rndg
       jjrow(i)=19.001e0-th2
       tthan(i)=(19-jjrow(i))-th2
   10 continue
c....   seasonal and spatial interpolation done below.
      do 20 k=1,nl
      do 20 i=1,ims2
        do3v = dduo3n(jjrow(i),k) + rsin1*ddo3n2(jjrow(i),k)
     1             +rcos1*ddo3n3(jjrow(i),k)
     2             +rcos2*ddo3n4(jjrow(i),k)
        do3vp = dduo3n(jjrow(i)+1,k) + rsin1*ddo3n2(jjrow(i)+1,k)
     1             +rcos1*ddo3n3(jjrow(i)+1,k)
     2             +rcos2*ddo3n4(jjrow(i)+1,k)
c...   now latitudinal interpolation, and
c          convert o3 into mass mixing ratio(original data mpy by 1.e4)
c      flip vertical coordinate too...
        qo3o3(i,nl+1-k) = 1.e-4 * (do3v+tthan(i)*(do3vp-do3v))
   20 continue
c...    vertical (linear in ln p) interpolate for each gridpoint
      numitr = 0
      ilog = nl
   21 continue
      ilog = (ilog+1)/2
        if(ilog.eq.1) go to 22
        numitr = numitr + 1
        go to 21
   22 continue
      do 60 k=1,kmx
        nhalf=(nl+1)/2
        do 30 i=1,ims2
          jjrow(i) = nhalf
          psm(i) = sfcp(i) * sigl(i,k)
   30   continue
        do 40 it=1,numitr
          nhalf=(nhalf+1)/2
          do 40 i=1,ims2
            if(psm(i).lt.prgfdl(jjrow(i))) then
              jjrow(i) = jjrow(i) + nhalf
            else if(psm(i).ge.prgfdl(jjrow(i)-1)) then
              jjrow(i) = jjrow(i) - nhalf
            endif
            jjrow(i) = min(jjrow(i),nl)
            jjrow(i) = max(jjrow(i),2)
   40   continue
        do 50 i=1,ims2
          if(psm(i).gt.prgfdl(1)) then
            qo3(i,k) = qo3o3(i,1)
          else if(psm(i).lt.prgfdl(nl)) then
            qo3(i,k) = qo3o3(i,nl)
          else
            aplo = alog(prgfdl(jjrow(i)-1))
            apup = alog(prgfdl(jjrow(i)))
            qo3(i,k) = qo3o3(i,jjrow(i)) + (alog(psm(i))-apup) /
     1                 (aplo-apup) *
     2                 (qo3o3(i,jjrow(i)-1)-qo3o3(i,jjrow(i)))
          endif
   50   continue
   60 continue
      return
      end
