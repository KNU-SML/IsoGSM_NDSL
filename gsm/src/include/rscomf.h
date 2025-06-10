c
c  begin two rloop (comf)
c     in-code nesting regional to global
c
c fillowing common has to be in order with the call of sums in
c     rloopa, rloopb, rloopr, rloopz
c
       common /rcomfi/ idate(4)
       common /rcomf/ 
     1    z(lnwav_)
     2 ,  y(lnwav_,levs_)
     3 , rt(lnwav_,levs_)
     4 , di(lnwav_,levs_)
     5 ,  q(lnwav_)
     6 , te(lnwav_,levs_)
     7 , rq(lnwav_,levs_)
     8 ,tem(lnwav_,levs_)
     9 , qm(lnwav_)
     4 , rm(lnwav_,levs_)
     a , gz(lnwav_)
     b ,  x(lnwav_,levs_)
     & ,dpdlam(lnwav_)
     c , uu(lnwav_,levs_)
     d ,uum(lnwav_,levs_)
     e ,  w(lnwav_,levs_)
     & ,dpdphi(lnwav_)
     f , vv(lnwav_,levs_)
     g ,vvm(lnwav_,levs_)
c
       common /rcomf1/
     4 flat(lngrd_), flon(lngrd_), 
     5  fm2(lngrd_), fm2x(lngrd_), fm2y(lngrd_)
c
      common /rcomf2i/ lastep
      logical lastep
      common /rcomf2/
     *  slmsk(igrd12_,jgrd12_),hprime(igrd12_,jgrd12_),
#ifdef RSMVECTORIZE
     *         swh(igrd12p_,jgrd12p_,levs_),
     *         hlw(igrd12p_,jgrd12p_,levs_),
#else
     * swh(igrd12_,levs_,jgrd12_),hlw(igrd12_,levs_,jgrd12_),
#endif
     * sfcnsw(igrd12_,jgrd12_),sfcdlw(igrd12_,jgrd12_),
#ifdef VIC
     * sfcusw(igrd12_,jgrd12_),sfcdsw(igrd12_,jgrd12_)   
#endif
     * sdec,cdec,slag,solhr,clstp,
     * sinlar(igrd12_,jgrd12_),coslar(igrd12_,jgrd12_),
     * albed (igrd12_,jgrd12_),coszer(igrd12_,jgrd12_),
     *     cv(igrd12_,jgrd12_),   cvt(igrd12_,jgrd12_),
     *    cvb(igrd12_,jgrd12_), tsflw(igrd12_,jgrd12_)
c
      common /rcomf3/ 
     *       dusfc(igrd12_,jgrd12_), dvsfc(igrd12_,jgrd12_),
     *       dtsfc(igrd12_,jgrd12_), dqsfc(igrd12_,jgrd12_),
     *      dlwsfc(igrd12_,jgrd12_),ulwsfc(igrd12_,jgrd12_),
     *      geshem(igrd12_,jgrd12_),  tsea(igrd12_,jgrd12_),
#ifdef NOAHYDRO
     *      geshem2(igrd12_,jgrd12_),
#endif
     *        f10m(igrd12_,jgrd12_),
     *       u10m(igrd12_,jgrd12_),v10m(igrd12_,jgrd12_),
     *       t2m(igrd12_,jgrd12_),q2m(igrd12_,jgrd12_),
     *       dugwd(igrd12_,jgrd12_), dvgwd(igrd12_,jgrd12_),
     *       psurf(igrd12_,jgrd12_),psmean(igrd12_,jgrd12_),
     *      dtflux
c
      common /rcomf4/
     *        tg3(igrd12_,jgrd12_),zorl(igrd12_,jgrd12_),
     *        plantr(igrd12_,jgrd12_),
     *        sheleg(igrd12_,jgrd12_),bengsh(igrd12_,jgrd12_),
     *        gflux(igrd12_,jgrd12_),  slrad(igrd12_,jgrd12_),
     *        smc(igrd12_,jgrd12_,lsoil_),
	 *        stc(igrd12_,jgrd12_,lsoil_),
     *        canopy(igrd12_,jgrd12_),runoff(igrd12_,jgrd12_),
     *        tmpmax(igrd12_,jgrd12_),tmpmin(igrd12_,jgrd12_),
     *        ep(igrd12_,jgrd12_),cldwrk(igrd12_,jgrd12_),
     *        hpbl(igrd12_,jgrd12_),pwat(igrd12_,jgrd12_)
#include <rscomver.h>
#include <rscomio.h>
#include <rscompln.h>
#include <rscomrad.h>
#include <rscomloc.h>
