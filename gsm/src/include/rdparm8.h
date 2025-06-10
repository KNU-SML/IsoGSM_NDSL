#ifdef REDUCE_GRID
#define NVECT 2
#else
#ifdef MP
#ifdef RSM
#ifdef NEC
#define NVECT igrd12p_
#else
#define NVECT min(igrd12p_,64)
#endif
#else
#define NVECT lonfp_*2
#endif
#else
#ifdef RSM
#ifdef NEC
#define NVECT igrd12_
#else
#define NVECT 64
#endif
#else
#ifdef SCM
#define NVECT 2
#else
#define NVECT 64
#endif
#endif
#endif
#endif
c
c     parameter settings for the longwave and shortwave radiation code:
c          imax   =  no. points along the lat. circle used in calcs.
c          l      =  no. vertical levels (also layers) in model
c   note: the user normally will modify only the imax and l parameters
c          nblw   =  no. freq. bands for approx computations. see
c                      bandta for definition
c          nblx   =  no. freq bands for approx cts computations
c          nbly   =  no. freq. bands for exact cts computations. see
c                      bdcomb for definition
c          inlte  =  no. levels used for nlte calcs.
c          nnlte  =  index no. of freq. band in nlte calcs.
c          nb,ko2 are shortwave parameters; other quantities are derived
c                    from the above parameters.
      parameter (l=levs_)
      parameter (imax=NVECT)
      parameter (imbx=imax)
      parameter (nblw=163,nblx=47,nbly=15)
      parameter (nblm=nbly-1)
      parameter (lp1=l+1,lp2=l+2,lp3=l+3)
      parameter (lm1=l-1,lm2=l-2,lm3=l-3)
      parameter (ll=2*l,llp1=ll+1,llp2=ll+2,llp3=ll+3)
      parameter (llm1=ll-1,llm2=ll-2,llm3=ll-3)
      parameter (lp1m=lp1*lp1,lp1m1=lp1m-1)
      parameter (lp1v=lp1*(1+2*l/2))
      parameter (lp121=lp1*nbly)
      parameter (ll3p=3*l+2)
      parameter (nb=12)
      parameter (inlte=3,inltep=inlte+1,nnlte=56)
      parameter (nb1=nb-1)
      parameter (ko2=12)
      parameter (ko21=ko2+1,ko2m=ko2-1)
      parameter (nbd=12,nvb=8,nrb=4,nk0=10,nae=6,nsrc=5)
