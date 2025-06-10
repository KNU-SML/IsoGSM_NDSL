      subroutine rascnv(len, lenc, km, lm, nstrp, dt, si, sl
     *,           krmin, krmax, ncrnd, afac, rannum, ufac
     *,           rgas, cp, grav, alhl
     *,           tin, qin, uin, vin, pin
     *,           evapc, kbot, ktop, kuo, lat, cd
     *,           sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns 
#ifdef ISOTOPE
     *,           ntrac,evapci,condi,cond
     *,           eqfo,eqfd
#endif
     *,           cwm ,cmf,detr)
c
#ifdef ISOTOPE
      dimension tin(lenc,km), qin(lenc,km*ntrac), uin(lenc,km)
#else
      dimension tin(lenc,km), qin(lenc,km), uin(lenc,km)
#endif
     *,         vin(lenc,km), pin(len), cd(len)
     *,         si(km+1), sl(km)
     *,         evapc(len), kbot(len), ktop(len), kuo(len)
     *,         rannum(ncrnd*3), cwm(lenc,km), cmf(lenc,km)
     *,         detr(lenc,km)
#ifdef ISOTOPE
     *,         evapci(lenc,2),condi(lenc,2),cond(len)
      dimension eqfo(len),eqfd(len)
      integer   ntrac
#endif

c
      dimension sig(lm+1),  prj(lm+1), prh(lm),  fpk(lm), hpk(lm)
     *,         sgb(lm),  ods(lm), rasal(lm), prns(nstrp)
c
c     common /rascmr/ sig, sgb, prh, prj, hpk, fpk, ods, prns
c    *,               rasal, rannum, afac, ufac
c     common / rascmi/krmin, krmax, nstrp, ncrnd
c
      call rascnp(len, lenc, km, lm, nstrp, dt, si, sl
     *,           krmin, krmax, ncrnd, afac, rannum, ufac
     *,           rgas, cp, grav, alhl
     *,           tin, qin, uin, vin, pin
     *,           evapc, kbot, ktop, kuo, lat, cd
     *,           sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns
#ifdef ISOTOPE
     *,           ntrac,evapci,condi,cond
     *,           eqfo,eqfd
#endif
     *,           cwm, cmf,detr )
c
      return
      end
