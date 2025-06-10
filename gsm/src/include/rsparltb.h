c
c   begin parltb 
c
#ifdef MP
#define IGRD1S  igrd1p_
#define JGRD12S jgrd12p_
#else
#define IGRD1S  igrd1_
#define JGRD12S jgrd12_
#endif
      parameter(ibgd1=(IGRD1S-1)/bgf_+2*border_+1,     ! for safety, we add 1
     1          jbgd1=((JGRD12S-1)/bgf_+2*border_+1)*2,  ! for safety, we add 2
     2          lngrdb=ibgd1*jbgd1)
      common/comparltb/ ib1,ib2,jb1,jb2,jb3,jb4,jbx,istr,ilen
