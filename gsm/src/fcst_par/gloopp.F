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
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NPES npes
#define TWOJ1 twoj1

      subroutine gloopp(q,avprs0)                                               
c                                                                               
c  surface pressure correction for long integrations                            
c                                                                               
#include <paramodel.h>
#include <comfgrid.h>
#include <comgpln.h>
#ifdef MP
#include <commpi.h>
      dimension qa(lln22p_),syf(lonf22_,latg2p_),avprsjp(latg2p_)
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
c                                                                               
#ifdef MP
#define JCAP1S jcap1p_
#define LNT22S lnt22p_
#define LLN22S lln22p_
#define LCAP22S lcap22p_
#else
#define JCAP1S jcap1_
#define LNT22S lnt22_
#define LLN22S lnt22_
#define LCAP22S lonf22_
#endif
      dimension q(LNT22S)
      dimension prs(LCAP22S,latg2_)
      dimension qnew(LLN22S,latg2_)
      dimension flp(2,JCAP1S,latg2_),flm(2,JCAP1S,latg2_)
      dimension avprsj(latg2_)                                                  
c=============================================================
#ifdef MP
      llstr=lwvstr(mype)
      llens=lwvlen(mype)
c                                                                               
c  loop to find average pressure                                                
c                                                                               
      avprs=0.                                                                  
      sumwgt=0.                                                                 
c
      call mpnn2n (q,lnt22p_,qa,lln22p_,1)
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(qtt,colrad,qa,prs,llstr,llens,lwvdef,lcapd,lcapdp,mype)
c$& local(lat,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,llensd)
#endif
      do lat=1,latg2_
#ifdef REDUCE_GRID
        llensd=lcapdp(lat,mype)
#else
        llensd=llens
#endif
        call sums2i(qa,prs(1,lat),qtt(1,lat),
     &              llstr,llensd,lwvdef,1)                                  
      enddo
c
      call mpnl2ny(prs,lcap22p_,latg2_,
     1                 syf,lonf22_,latg2p_,1,1,1)
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(syf,latlen,latstr,mype,avprsjp,wgb,lcapd,lonfd,latdef)
c$&        local(lan,lat,i,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lan,lat,lat,i,lcapf,lonff)
#endif
      do lan=1,latlen(mype)
        lat=latstr(mype)+lan-1
#ifdef REDUCE_GRID
        lcapf=lcapd(latdef(lat))
        lonff=lonfd(latdef(lat))
#else
        lcapf=jcap1_
        lonff=lonf_
#endif
        call ftilonf(syf(1,lan),syf(1,lan),2,lcapf,lonff,latdef(lat),1)
        avprsjp(lan)=0.
        do i=1,lonff*2
          avprsjp(lan)=avprsjp(lan)+syf(i,lan)
        enddo                                                                   
        avprsjp(lan)=avprsjp(lan)/float(lonff*2)
      enddo
c
      call mplatall(avprsjp,latg2p_,avprsj,latg2_,1)
c
c ******** no multithread or multi-node for next loop
c          otherwise non-reproduciable *****************
c
      do lat=1,latg2_
        sumwgt=sumwgt+wgt(lat)
        avprs=avprs+avprsj(lat)*wgt(lat)
      enddo                                                                     
      avprs=avprs/sumwgt                                                        
#ifndef NOPRINT
      if(mype.eq.master) then
        print *,'sumwgt=',sumwgt,' avprs=',avprs,' avprs0=',avprs0 
      endif
#endif
c                                                                               
c  pressure correction loop                                                     
c                                                                               
      if(avprs0.gt.0) then                                                      
c
#ifdef ORIGIN_THREAD
c$doacross share(qa) local(i)
#endif
#ifdef OPENMP
!$omp parallel do private(i)
#endif
        do i=1,lln22p_
          qa(i)=0.0
        enddo                                                                   
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(syf,avprs,avprs0,latlen,mype,lcapd,lonfd,latdef)                         
c$& local(lan,lat,i,lcapf,lonff) 
#endif
#ifdef OPENMP
!$omp parallel do private(lan,lat,i,lcapf,lonff)
#endif
c                                                                               
        do lan=1,latlen(mype)
          lat=latstr(mype)+lan-1
#ifdef REDUCE_GRID
          lcapf=lcapd(latdef(lat))
          lonff=lonfd(latdef(lat))
#else
          lcapf=jcap1_
          lonff=lonf_
#endif
          do i=1,lonff*2
            syf(i,lan)=syf(i,lan)-avprs+avprs0
          enddo                                                                 
          call ftilonf(syf(1,lan),syf(1,lan),2,lcapf,lonff,
     1                                         latdef(lat),-1)
        enddo
        call mpny2nl(syf,lonf22_,latg2p_,
     1               prs,lcap22p_,latg2_,1,1,1)
c
#ifdef ORIGIN_THREAD
c$doacross share(prs,flp,flm,qa,llens,
c$&              llstr,llens,lwvdef,lcapdp,mype)
c$& local(lat,i,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,llensd)
#endif
        do lat=1,latg2_
#ifdef REDUCE_GRID
          llensd=lcapdp(lat,mype)
#else
          llensd=llens
#endif
          call flpflm(flp(1,1,lat),flm(1,1,lat),prs(1,lat),llensd,1)
        enddo
c
c ******** no multithread for the next loop
c          otherwise non-reproduciable *****************
c
        do lat=1,latg2_
#ifdef REDUCE_GRID
        llensd=lcapdp(lat,mype)
#else
        llensd=llens
#endif
          call fl2ip(flp(1,1,lat),flm(1,1,lat),qa,qww(1,lat),
     &               llstr,llensd,lwvdef,1)
        enddo                                                                   
c                                                                               
        call mpn2nn(qa,lln22p_,q,lnt22p_,1)
c                                                                               
      else                                                                      
c                                                                               
        avprs0=avprs                                                            
c                                                                               
      endif                                                                     
c                                                                               
#else
c ================================================================
c ============= following is non MP ==============================
c ================================================================
      llstr=0
      llens=jcap1_
c                                                                               
c  loop to find average pressure                                                
c                                                                               
      avprs=0.                                                                  
      sumwgt=0.                                                                 
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(qtt,colrad,q,prs,avprsj,wgt,
c$&              llstr,llens,lwvdef,lcapd,lonfd),
c$&        local(lat,i,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,lcapf,lonff)
#endif
      do lat=1,latg2_                                                           
#ifdef REDUCE_GRID
        lcapf=lcapd(lat)
        lonff=lonfd(lat)
#else
        lcapf=jcap1_
        lonff=lonf_
#endif
        call sums2i(q,prs(1,lat),qtt(1,lat),
     &              llstr,lcapf,lwvdef,1)                                  
        call ftilonf(prs(1,lat),prs(1,lat),2,lcapf,lonff,lat,1)
        avprsj(lat)=0.                                                          
        do i=1,lonff*2
          avprsj(lat)=avprsj(lat)+prs(i,lat)
        enddo                                                                   
        avprsj(lat)=avprsj(lat)/float(lonff*2)
      enddo
c
c ******** no multithread or multi-node computational
c          otherwise non-reproduciable *****************
c
      do lat=1,latg2_
        sumwgt=sumwgt+wgt(lat)
        avprs=avprs+avprsj(lat)*wgt(lat)
      enddo                                                                     
      avprs=avprs/sumwgt                                                        
#ifndef NOPRINT
      print *,'sumwgt=',sumwgt,' avprs=',avprs,' avprs0=',avprs0                
#endif
c                                                                               
c  pressure correction loop                                                     
c                                                                               
      if(avprs0.gt.0) then                                                      
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(prs,avprs,avprs0,flp,flm,qnew,qww,
c$&              llstr,llens,lwvdef,lcapd,lonfd),
c$&        local(lat,i,lcapf,lonff)  
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,lcapf,lonff)
#endif
c                                                                               
        do lat=1,latg2_                                                         
#ifdef REDUCE_GRID
         lcapf=lcapd(lat)
         lonff=lonfd(lat)
#else
         lcapf=jcap1_
         lonff=lonf_
#endif
          do i=1,lnt2_                                                          
            qnew(i,lat)=0.e0                                                    
          enddo                                                                 
          do i=1,lonff*2
            prs(i,lat)=prs(i,lat)-avprs+avprs0
          enddo                                                                 
          call ftilonf(prs(1,lat),prs(1,lat),2,lcapf,lonff,lat,-1)                              
          call flpflm(flp(1,1,lat),flm(1,1,lat),prs(1,lat),lcapf,1)
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
          call fl2i(flp(1,1,lat),flm(1,1,lat),qnew(1,lat),qww(1,lat),
     &              llstr,lcapf,lwvdef,1)         
#endif
#ifdef DEFAULT
          call fl2ip(flp(1,1,lat),flm(1,1,lat),qnew(1,lat),qww(1,lat),
     &               llstr,lcapf,lwvdef,1)        
#endif
        enddo                                                                   
c                                                                               
#ifdef ORIGIN_THREAD
c$doacross share(q,qnew) local(i,lat)                                           
#endif
#ifdef OPENMP
!$omp parallel do private(i,lat)
#endif
        do i=1,lnt2_                                                            
          q(i)=0.0                                                              
          do lat=1,latg2_                                                       
            q(i)=q(i)+qnew(i,lat)                                               
          enddo                                                                 
        enddo                                                                   
c                                                                               
      else                                                                      
c                                                                               
        avprs0=avprs                                                            
c                                                                               
      endif                                                                     
c===================================================================
#endif
      return                                                                    
      end                                                                       
