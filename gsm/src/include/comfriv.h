#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
c
c  begin comfriv (kei)
c
      real trunof,imap,gdriv,rflow,roff
      common /comfriv/
     *          trunof(LONF2S,LATG2S),
     *          imap(io2_,jo2_),
     *          gdriv(io2_,jo2_),
     *          rflow(io2_,jo2_),
     *          roff(io2_,jo2_)
c
#ifdef ISOTOPE
      real trunof1,gdriv1,rflow1,roff1
      real trunof2,gdriv2,rflow2,roff2
      common /comfriviso/
     *          trunof1(LONF2S,LATG2S),
     *          trunof2(LONF2S,LATG2S),
     *          gdriv1(io2_,jo2_),
     *          gdriv2(io2_,jo2_),
     *          rflow1(io2_,jo2_),
     *          rflow2(io2_,jo2_),
     *          roff1(io2_,jo2_),
     *          roff2(io2_,jo2_)
#endif
c
