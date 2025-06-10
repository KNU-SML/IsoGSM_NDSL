#define IGRD12P igrd12p
#define JGRD12P jgrd12p
#define IGRD12 igrd12
#define JGRD12 jgrd12
#ifdef MP
#define IGRD12S igrd12p_
#define JGRD12S jgrd12p_
#else
#define IGRD12S igrd12_ 
#define JGRD12S jgrd12_
#endif
c   begin commap
      common/mapcom/  rlat(IGRD12S,JGRD12S),  rlon(IGRD12S,JGRD12S)
     1              ,  xm2(IGRD12S,JGRD12S),  xm2p(IGRD12S,JGRD12S)
     2              ,xm2px(IGRD12S,JGRD12S), xm2py(IGRD12S,JGRD12S)
     2              , gzdx(IGRD12S,JGRD12S),  gzdy(IGRD12S,JGRD12S)
     3              , corf(IGRD12S,JGRD12S),    xm(IGRD12S,JGRD12S)
     4              , xm2m
