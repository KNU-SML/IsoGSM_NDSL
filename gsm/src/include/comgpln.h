#ifdef MP
#define LNT2S  lln2p_
#define TWOJ1S lcapp_
#else
#define LNT2S  lnt2_
#define TWOJ1S twoj1_
#endif
c
      common/comgpln/ qtt(LNT2S,latg2_),qvv(TWOJ1S,latg2_),
     1                qdd(LNT2S,latg2_),qww(LNT2S,latg2_)
#undef LNT2S
#undef TWOJ1S
