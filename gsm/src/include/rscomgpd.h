#ifdef RKN
c
      parameter(nvrken=rslvark_+rmlvark_*levs_,nptken=rlpnt_)
      parameter(nstken=rltstp_)
       common/rcomgpd/ svdata(nvrken,nptken,nstken)
       common/rcomgpdi/ 
     &               igpd (nptken),jgpd (nptken),
     &               igpdr(nptken),jgpdr(nptken),
     &               itnum,npoint,isave,isshrt,ilshrt,ikfreq,
     &               imodk,irstkn
#endif
