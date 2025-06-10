#define JCAP1 jcap1
#define JCAP jcap
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT2 lnt2
#define RERTH rerth
      subroutine dellnp(q,dpdphs,dpdtop,dpdla,llstr,llens,lwvdef)
c
#include <paramodel.h>
      parameter (jcap=jcap_)
      parameter (lnep=(jcap+2)*(jcap+3)/2)
c
c input q is in ibm triang. order
c output  is in ibm triang. order
c
#ifdef RSM
#undef MP
#endif
#ifdef MP
#define LLN2S lln2p_
#define LLN22S lln22p_
#else
#define LLN2S lnt2_
#define LLN22S lnt22_
#endif
      dimension         q(LLN22S)
      dimension    dpdphs(LLN22S)
      dimension    dpdtop(2,jcap1_)
      dimension     dpdla(LLN22S)
      dimension     lwvdef(jcap1_)
      dimension         e(lnep)
c
      save e, ifirst
c
      je(n,l) =((jcap+2)*(jcap+3)-(jcap+2-l)*(jcap+3-l))/2+n-l
c
      jc(n,l) = (jcap+1)*(jcap+2)-(jcap+1-l)*(jcap+2-l)+2*(n-l)
c
      data ifirst/1/
      if(ifirst.ne.1)go to 280
c --------------------- initial setting---------------
      do 220 l=0,jcap
      n=l
      ie=je(n,l)+1
      e(ie)=0.
220   continue
      do 260 l=  0,jcap
      do 240 n=l+1,jcap+1
      rn=n
      rl=l
      a=(rn*rn-rl*rl)/(4.*rn*rn-1.)
      ie=je(n,l)+1
      e(ie)=sqrt(a)
240   continue
260   continue
      ifirst=0
280   continue
c ------------------ end of initial setting -----------
c
      aa=1./rerth_
c
      lntx=0
      do 10000 ll=1,llens
        l=lwvdef(llstr+ll)
        lnt0=jc(l,l)-lntx
        lntx=lntx+jc(l+1,l+1)-jc(l,l)
c
c ----- l = 0, jcap
c
        rl=l
        do 320 n=l,jcap
          icr=jc(n,l)+1 - lnt0
          ici=jc(n,l)+2 - lnt0
          dpdla(ici)= rl*q(icr)
          dpdla(icr)=-rl*q(ici)
320     continue
c
c ------ l = 0, jcap-1
c
        if( l.le.jcap-1 ) then
          do 420 n=l,jcap-1
             ie=je(n+1,l)+1
            icr=jc(n,l)+1 - lnt0
            ici=jc(n,l)+2 - lnt0
             rn=n
            dpdphs(icr)=(rn+2.)*e(ie)*q(icr+2)
            dpdphs(ici)=(rn+2.)*e(ie)*q(ici+2)
420       continue
        endif
c
c ------ l = 0, jcap
c
        n=  jcap
        icr=jc(n,l)+1 - lnt0
        ici=jc(n,l)+2 - lnt0
        dpdphs(icr)=0.0
        dpdphs(ici)=0.0
c
c ------ l = 0, jcap-1
c
        if( l.le.jcap-1 ) then
          do 620 n=l+1,jcap
             ie=je(n,l)+1
            icr=jc(n,l)+1 - lnt0
            ici=jc(n,l)+2 - lnt0
             rn=n
            dpdphs(icr)=dpdphs(icr)+(1.-rn)*e(ie)*q(icr-2)
            dpdphs(ici)=dpdphs(ici)+(1.-rn)*e(ie)*q(ici-2)
620       continue
        endif
c
c ------ l = 0, jcap
c
        n=jcap+1
        rn=n
         ie=je(n,l)+1
        icr=jc(n,l)+1 - lnt0
        ici=jc(n,l)+2 - lnt0
        dpdtop(1,ll)=(1.-rn)*e(ie)*q(icr-2)
        dpdtop(2,ll)=(1.-rn)*e(ie)*q(ici-2)
        dpdtop(1,ll)=dpdtop(1,ll)*aa
        dpdtop(2,ll)=dpdtop(2,ll)*aa
c
10000 continue
c
      do 820 j=1,LLN2S
       dpdla(j)= dpdla(j)*aa
      dpdphs(j)=dpdphs(j)*aa
820   continue
c
      return
      end
