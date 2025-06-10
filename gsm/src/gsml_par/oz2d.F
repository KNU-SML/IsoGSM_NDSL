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
      subroutine oz2d(ims2,imx2,kmx,sigl,qo3,sfcp,xlat)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c ...  compute model lyr o3 profile from the original nasa data
#include <paramodel.h>
      parameter (loz=17,j5=37)
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
      dimension wgts(ILOTS),jlats(ILOTS),
     *          psm(ILOTS)
      dimension jjrow(ILOTS)
      dimension qo3o3(ILOTS,loz)
c...    **************************************************************
c--   monthly climatology of o3 from o3sbuv
      common /sbuv/ psnasa(loz),o3nasa(j5,loz)
c...    **************************************************************
c...       begin here .....
      rndg = 180./pi_
      do 10 i=1,ims2
       rlat=xlat(i)*rndg
       do 5 j=1,j5
        lato3=90-(j-1)*5
        xlato3=lato3
        if (xlato3.lt.rlat) then
         jlats(i)=j
         wgts(i) = (rlat-xlato3)/5.
cc         go to 1013
         go to 10
        end if
    5  continue
   10 continue
c....   latitudinal (linear) interpolation done below.
c      flip vertical coordinate too...
      do 20 k=1,loz
       do 20 i=1,ims2
        qo3o3(i,loz+1-k)=wgts(i)
     1                   * (o3nasa(jlats(i)-1,k)-o3nasa(jlats(i),k))
     2                   +  o3nasa(jlats(i),k)
   20 continue
c...    vertical (linear in ln p) interpolate for each gridpoint
      numitr = 0
      ilog = loz
   21 continue
      ilog = (ilog+1)/2
        if(ilog.eq.1) go to 22
        numitr = numitr + 1
        go to 21
   22 continue
      do 60 k=1,kmx
        nhalf=(loz+1)/2
        do 30 i=1,ims2
          jjrow(i) = nhalf
          psm(i) = sfcp(i) * sigl(i,k)
   30   continue
        do 40 it=1,numitr
          nhalf=(nhalf+1)/2
          do 40 i=1,ims2
            if(psm(i).lt.psnasa(jjrow(i))) then
              jjrow(i) = jjrow(i) + nhalf
            else if(psm(i).ge.psnasa(jjrow(i)-1)) then
              jjrow(i) = jjrow(i) - nhalf
            endif
cc            jjrow(i) = min(jjrow(i),nl)
            jjrow(i) = min(jjrow(i),loz)
            jjrow(i) = max(jjrow(i),2)
   40   continue
        do 50 i=1,ims2
          if(psm(i).gt.psnasa(1)) then
            qo3(i,k) = qo3o3(i,1)
          else if(psm(i).lt.psnasa(loz)) then
            qo3(i,k) = qo3o3(i,loz)
          else
            aplo = alog(psnasa(jjrow(i)-1))
            apup = alog(psnasa(jjrow(i)))
            qo3(i,k) = qo3o3(i,jjrow(i)) + (alog(psm(i))-apup) /
     1                 (aplo-apup) *
     2                 (qo3o3(i,jjrow(i)-1)-qo3o3(i,jjrow(i)))
          endif
   50   continue
   60 continue
      return
      end
