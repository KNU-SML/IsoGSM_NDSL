#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAPP lcapp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define MLVARK mlvark
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTRAC ntrac
#define SLVARK slvark
#define TWOJ1 twoj1
#define NTOTAL ntotal
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#define G g
#define PI pi





#ifdef NISLQ
c
c  common block for radioactive transport in semi-lagrangian mode
c
       real slq_q1, slq_q2, slq_q3, slq_psfc2, slq_u2, 
     1      slq_v2, slq_w2, racos, racos2, glat, gglat, 
     1      gglati, gglon, ggloni, slq_qtemp
#ifdef NISLQ_MASS
     1,     slq_psfc1, slq_psfc3
#endif /* end NISLQ_MASS */
c
       common /nislq/
     1    slq_q1   (LONF2S,levh_ ,LATG2S)
     2 ,  slq_q2   (LONF2S,levh_ ,LATG2S)
     3 ,  slq_q3   (LONF2S,levh_ ,LATG2S)
     3 ,  slq_qtemp (LONF2S,levh_ ,LATG2S)
     4 ,  slq_psfc2(LONF2S       ,LATG2S)
#ifdef NISLQ_MASS
     4 ,  slq_psfc1(LONF2S       ,LATG2S)
     4 ,  slq_psfc3(LONF2S       ,LATG2S)
#endif /* end NISLQ_MASS */
     5 ,  slq_u2   (LONF2S,levs_ ,LATG2S)
     6 ,  slq_v2   (LONF2S,levs_ ,LATG2S)
     7 ,  slq_w2   (LONF2S,levp1_,LATG2S)
     8 ,  racos(LATG2S)
     9 ,  racos2(LATG2S)
     1 ,  glat(latg_), gglat(latg_*2)
     1 ,  gglati(latg_*2+1), gglon(lonf_), ggloni(lonf_+1)
c
       integer lonfull, lonpart, latfull, latpart
     1        ,mysllonlen, mysllatlen, jlist1
     1        ,sllonstr, sllonlen ,sllatstr ,sllatlen 
     1        ,nx, my, my_max, lev, nlevs, nlevsp
     1        ,lonhalf, lathalf, truej, shflj
       common /nislqdim/ lonfull,lonpart,latfull,latpart
     1                  ,mysllonlen, mysllatlen, nx, my
     1                  , my_max, lev, nlevs
     1                  , nlevsp, lonhalf, lathalf, truej(latg_)
     1		        , shflj(latg_), jlist1(LATG2S)
#ifdef MP
     1                  ,sllonstr(nrow_), sllonlen(nrow_)
     1                  ,sllatstr(nrow_), sllatlen(nrow_)
#else
     1                  ,sllonstr(1), sllonlen(1)
     1                  ,sllatstr(1), sllatlen(1)
#endif /* end MP */
c
       real wgtb, wgtm
       common /nislqbdy/ wgtb(LONF2S, LATG2S)
     1                  ,wgtm(LONF2S, LATG2S)
c
#endif /* end NISLQ */
