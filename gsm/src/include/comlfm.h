#ifdef MP
#define LNT22S lnt22p_
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LNT22S lnt22_
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
c
      common/comlfm/ weight(1000)
      common/comlfm/ fq (LNT22S)
      common/comlfm/ fte(LNT22S,levs_)
      common/comlfm/ fdi(LNT22S,levs_)
      common/comlfm/ fze(LNT22S,levs_)
      common/comlfm/ frq(LNT22S,levh_)
      common/comlfm/
     *         ftsea(LONF2S,LATG2S),   fsmc(LONF2S,LATG2S,lsoil_),
     1       fsheleg(LONF2S,LATG2S),   fstc(LONF2S,LATG2S,lsoil_),
     2          ftg3(LONF2S,LATG2S),  fzorl(LONF2S,LATG2S),
     3       fplantr(LONF2S,LATG2S),    fcv(LONF2S,LATG2S),
     4          fcvb(LONF2S,LATG2S),   fcvt(LONF2S,LATG2S),
     5       falbedo(LONF2S,LATG2S,lalbd_), fslmsk(LONF2S,LATG2S),
     6         ff10m(LONF2S,LATG2S),fcanopy(LONF2S,LATG2S)
      common/comlfm/ 
     *        islmsk(LONF2S,LATG2S,3),
     1          wcvb(LONF2S,LATG2S),   wcvt(LONF2S,LATG2S)
      common/comlfm/ filtwin,klenp,nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo
