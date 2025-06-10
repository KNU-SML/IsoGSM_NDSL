      integer idxffcb,idxffcc
      real idxffca1,idxffca2
#ifdef MP
#define IWAV1S iwav1p_
#define LEVHS levhp_
#define LEVSS levsp_
#else
#define IWAV1S iwav1_
#define LEVHS levh_
#define LEVSS levs_
#endif
      common/sumfidx/
     1 idxffca1(IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),jwav1_,7),
     2 idxffca2(IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),jwav1_,7),
     3 idxffcb (IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),3,7),
     4 idxffcc (IWAV1S*2*jgrd12_*(LEVHS*2+LEVSS*2+2),jwav1_,7)
