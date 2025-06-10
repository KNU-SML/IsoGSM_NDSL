c
      integer nvrken,nptken,nstken,igrd,jgrd
     1       ,itnum,npoint,isave,isshrt,ilshrt,ikfreq
      real svdata
      parameter(nvrken=slvark_+mlvark_*levs_,nptken=lpnt_)
      parameter(nstken=ltstp_)
#ifdef KEN
      common/comgpd/ svdata(nvrken,nptken,nstken),
     1               igrd(nptken),jgrd(nptken),
     2               itnum,npoint,isave,isshrt,ilshrt,ikfreq
#endif
