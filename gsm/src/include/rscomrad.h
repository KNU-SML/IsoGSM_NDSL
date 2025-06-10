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
c
      common /rradia1/ fluxr(IGRD12S,JGRD12S,26)
      common /rradia2/ cvavg(IGRD12S,JGRD12S)

