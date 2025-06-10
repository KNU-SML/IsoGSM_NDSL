#ifdef MP
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#define LNWAVS lnwavp_
#define LNGRDS lngrdp_
#else
#define IGRD12S igrd12_ 
#define JGRD12S jgrd12_
#define LNWAVS lnwav_
#define LNGRDS lngrd_
#endif
c
c
c  following common has to be in order with the call of sums in
c    rloopa, rloopb, rloopr, rloopz
c
c
       integer idate
       common /rcomfi/ idate(4)
       real z,y,rt,di,q,te,rq,tem,qm,rm,gz,x,dpdlam,uu,uum,w,dpdphi,
     1      vv,vvm
       common /rcomf/ 
     1    z(LNWAVS)
     2 ,  y(LNWAVS,levs_)
     3 , rt(LNWAVS,levh_)
     4 , di(LNWAVS,levs_)
     5 ,  q(LNWAVS)
     6 , te(LNWAVS,levs_)
     7 , rq(LNWAVS,levh_)
     8 ,tem(LNWAVS,levs_)
     9 , qm(LNWAVS)
     4 , rm(LNWAVS,levh_)
     a , gz(LNWAVS)
     b ,  x(LNWAVS,levs_)
     & ,dpdlam(LNWAVS)
     c , uu(LNWAVS,levs_)
     d ,uum(LNWAVS,levs_)
     e ,  w(LNWAVS,levs_)
     & ,dpdphi(LNWAVS)
     f , vv(LNWAVS,levs_)
     g ,vvm(LNWAVS,levs_)
c
#ifdef MP
       real za,ya,rta,dia,qa,tea,rqa,tema,qma,rma,gza,xa,dpdlama,uua
     1      uuma,wa,dpdphia,vva,vvma
       common /rcomfp/ 
     1       za(llwavp_)
     2 ,     ya(llwavp_,levsp_)
     3 ,    rta(llwavp_,levhp_)
     4 ,    dia(llwavp_,levsp_)
     5 ,     qa(llwavp_)
     6 ,    tea(llwavp_,levsp_)
     7 ,    rqa(llwavp_,levhp_)
     8 ,   tema(llwavp_,levsp_)
     9 ,    qma(llwavp_)
     4 ,    rma(llwavp_,levhp_)
     a ,    gza(llwavp_)
     b ,     xa(llwavp_,levsp_)
     & ,dpdlama(llwavp_)
     c ,    uua(llwavp_,levsp_)
     d ,   uuma(llwavp_,levsp_)
     e ,     wa(llwavp_,levsp_)
     & ,dpdphia(llwavp_)
     f ,    vva(llwavp_,levsp_)
     g ,   vvma(llwavp_,levsp_)
#endif
c
       real flat,flon,fm2,fm2x,fm2y
       common /rcomf1/
     4 flat(LNGRDS), flon(LNGRDS), 
     5  fm2(LNGRDS), fm2x(LNGRDS), fm2y(LNGRDS)
c
#include <comsfc.h>
c
      logical lastep
      common /rcomf2i/ lastep
      real swh,hlw,ozon,solin,sfcnsw,sfcdlw,sfcusw,sfcdsw,coszer,
     1     coszdg,sdec,cdec,slag,solhr,clstp,sinlar,coslar,albed,
     2     tsflw
      common /rcomf2/
#ifdef RSMVECTORIZE
     *         swh(igrd12p_,jgrd12p_,levs_),
     *         hlw(igrd12p_,jgrd12p_,levs_),
     *        ozon(igrd12p_,jgrd12p_,levs_),
     *        solin(igrd12p_,jgrd12p_,levs_),
#else
     *    swh(IGRD12S,levs_,JGRD12S),
     *    hlw(IGRD12S,levs_,JGRD12S),
     *   ozon(IGRD12S,levs_,JGRD12S),
     *   solin(IGRD12S,levs_,JGRD12S),
#endif
     * sfcnsw(IGRD12S,JGRD12S),sfcdlw(IGRD12S,JGRD12S),
#ifdef VIC
     * sfcusw(IGRD12S,JGRD12S),sfcdsw(IGRD12S,JGRD12S),
#endif
     * coszer(IGRD12S,JGRD12S),coszdg(IGRD12S,JGRD12S),
     * sdec,cdec,slag,solhr,clstp,
     * sinlar(IGRD12S,JGRD12S),coslar(IGRD12S,JGRD12S),
     * albed (IGRD12S,JGRD12S),
     *  tsflw(IGRD12S,JGRD12S)
c
#include <comrflx.h>
#ifdef ISOTOPE
#include <comriso.h>
#endif
#include <rscomver.h>
#include <rscomio.h>
#include <rscompln.h>
#include <rscomrad.h>
#include <rscomloc.h>
