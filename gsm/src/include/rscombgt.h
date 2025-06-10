#ifdef A
c  include budget common
#endif
#ifdef T
      common/rbgtt/ wt(lngrd_,levs_,nt_+1)
#endif
#ifdef Q
      common/rbgtq/ wq(lngrd_,levs_,nq_+1)
#endif
#ifdef U
      common/rbgtu/ wu(lngrd_,levs_,nu_+1)
#endif
#ifdef V
      common/rbgtv/ wv(lngrd_,levs_,nv_+1)
#endif
#ifdef P
      common/rbgtp/ wp(lngrd_,1,np_+1)
#endif
#ifdef A
      common/rbgta/ tmpbgt(lngrd_,levs_)
#endif
