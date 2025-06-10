#define G g
#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAP jcap
#define JCAPP jcapp
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAP22P lcap22p
#define LCAPP lcapp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NPES npes
#define RERTH rerth
#define TWOJ1 twoj1

      subroutine mtntrq(dt,snnp1,colrad,zs,z00,psmean,dumtn)                    
#include <paramodel.h>
#include <comgpln.h>
#ifdef MP
#define JCAP1S jcap1p_
#define LNT22S lnt22p_
#define LATG2S latg2p_
#define LONF2S lonf2p_
#define LCAP22S lcap22p_
#else
#define JCAP1S jcap1_
#define LNT22S lnt22_
#define LATG2S latg2_
#define LONF2S lonf2_
#define LCAP22S lonf22_
#endif
#ifdef MP
#include <commpi.h>
      dimension zsa(lln22p_),syn(LCAP22S,2,latg2_)
      dimension syf(lonf22_,2,latg2p_),grs(lonf22p_,2,latg2p_)
#else
      real syn(lonf22_,2)                                                       
#endif
      real snnp1(lnt2_),colrad(latg2_)                                          
      real zs(LNT22S),psmean(LONF2S,LATG2S),dumtn(LONF2S,LATG2S) 
      real spc(LNT22S,3),spctop(2,JCAP1S)
      asq=rerth_*rerth_/g_                                                      
c
#ifdef MP
      llstr=lwvstr(mype)
      llens=lwvlen(mype)
      lats2=latlen(mype)
      lons2=lonlen(mype)*2
      lnps2=lnplen(mype)*2
#else
      llstr=0
      llens=jcap1_
#endif
c
#ifdef MP
      call mpnn2n (zs,lnt22p_,zsa,lln22p_,1)
      do j=1,lnps2
        jj=lnpstr(mype)*2+j
        if(snnp1(jj).ne.0.0) spc(j,1)=zsa(j)*asq/snnp1(jj)
      enddo                                                                     
      if( mype.eq.master ) then
       spc(1,1)=z00
       spc(2,1)=0.
      endif
#else
      do j=3,lnt2_
        spc(j,1)=zs(j)*asq/snnp1(j)
      enddo
      spc(1,1)=z00
      spc(2,1)=0.
#endif
      call dellnp(spc(1,1),spc(1,3),spctop,spc(1,2),
     1            llstr,llens,lwvdef)
c.................................................................
#ifdef MP
#ifdef ORIGIN_THREAD
c$doacross share(colrad,spc,syn,llstr,llens,lwvdef),
c$&        local(lat,qtt)                                 
#endif
#ifdef CRAY_THREAD
cmic$ do all                                                                    
cmic$& shared(colrad,spc,syn,llstr,llens,lwvdef)
cmic$& private(lat,qtt)                                   
#endif
#ifdef OPENMP
!$omp parallel do private(lat,qtt)
#endif
      do lat=1,latg2_                                                           
        call sums2i(spc,syn(1,1,lat),qtt(1,lat),
     1               llstr,llens,lwvdef,2)
      enddo
      call mpnl2ny(syn,lcap22p_,latg2_,
     1             syf,lonf22_,latg2p_,2,1,2)
#ifdef ORIGIN_THREAD
c$doacross share(syf,lats2),local(lat)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syf,lats2)
cmic$& private(lat)
#endif
#ifdef OPENMP
!$omp parallel do private(lat)
#endif
      do lat=1,lats2
        call ftilonf(syf(1,1,lat),syf(1,1,lat),2*2,1)
      enddo
      call mpx2nx (syf,lonf22_,2,grs,lonf22p_,2,latg2p_,1,1,2)
#ifdef ORIGIN_THREAD
c$doacross share(grs,dt,psmean,dumtn,lats2,lons2),local(lat,i,psref)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(grs,dt,psmean,dumtn,lats2,lons2)
cmic$& private(lat,i,psref)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,psref)
#endif
      do lat=1,lats2
        do i=1,lons2
          psref=dt*101.325*(1.-grs(i,1,lat)/8.e3)
          dumtn(i,lat)=1.e3*rerth_*(psmean(i,lat)-psref)*grs(i,2,lat)
        enddo                                                                   
      enddo                                                                     
c .......................................................................
#else

#ifdef ORIGIN_THREAD
c$doacross share(colrad,spc,dt,psmean,dumtn),
c$&        local(lat,qtt,syn,i,psref)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& private(lat,qtt,syn,i,psref)
cmic$& shared(colrad,spc,dt,psmean,dumtn)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,qtt,syn,i,psref)
#endif
      do lat=1,latg2_
        call sums2i(spc,syn,qtt(1,lat),2)
        call ftilonf(syn,syn,2*2,1)
        do i=1,lonf2_
          psref=dt*101.325*(1.-syn(i,1)/8.e3)
          dumtn(i,lat)=1.e3*rerth_*(psmean(i,lat)-psref)*syn(i,2)
        enddo
      enddo
c ......................................................................
#endif
      return                                                                    
      end                                                                       
