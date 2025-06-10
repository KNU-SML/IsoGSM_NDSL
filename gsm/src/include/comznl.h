#ifdef MP
#define LATG2S latg2p_
#else
#define LATG2S latg2_
#endif
c
      parameter(nrm=23,nrs=32)
      parameter(nlb=6,nst=6)
      parameter(nmu=1,nmv=2,nmtv=3,nmq=4,nmvot2=5,
     &          nmdiv2=6,nmomega=7,nmt=8,nmrh=9,nmke=10,
     &          nmtconv=11,nmtlarg=12,nmtshal=13,nmtvrdf=14,nmqconv=15,
     &          nmqshal=16,nmqvrdf=17,nmuvrdf=18,nmvvrdf=19,nmthsw=20,
     &          nmthlw=21,nmtcld=22,nmtccv=23)
      parameter(nsrain=1,nsrainc=2,nstsfc=3,nsqsfc=4,nsusfc=5,
     &          nsvsfc=6,nsrcov=7,nsrcovc=8,nsps=9,nstskin=10,
     &          nswet=11,nssnow=12,nstg1=13,nstg2=14,nstg3=15,
     &          nssfcsw=16,nssfclw=17,nsrhs=18,nstvs=19,nsts=20,
     &          nsqs=21,nszorl=22,nsslmsk=23,nsugwd=24,nsvgwd=25,
     &          nsuasfc=26,nsuagwd=27,nsuamtn=28,nsua=29,nsuap=30,
     &          nsep=31,nscldwrk=32)
      common /comznl/ zdm(2,levs_,nrm,LATG2S),zwm(2,LATG2S),zhm(nrm)
      common /comznl/ zds(2,nst,nrs,LATG2S),zws(2,nst,LATG2S)
