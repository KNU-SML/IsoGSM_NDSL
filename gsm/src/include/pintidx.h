#ifdef MP
#define IGRD12S igrd12p_
#else
#define IGRD12S igrd12_
#endif
      parameter(ncntmax=border_*2)
      common/pintidx/ncount,irra(IGRD12S,ncntmax),
     1                      ibba(IGRD12S,ncntmax),
     2                      iwa (IGRD12S,ncntmax),
     3                      ina (IGRD12S,ncntmax)
