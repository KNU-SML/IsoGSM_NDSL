#define CP cp
#define LEVS levs
#define RD rd
#define RERTH rerth
      subroutine ambmsv_hyb(kmx,ak5,bk5,am,bm,sv,gv,cm)
#include <paramodel.h>

      real
     &     pk5ref(levs_+1),beta,dpkref(levs_),
     &     tref(levs_),psref,kappa,factor,
     &     alfaref(levs_),
     &     vecm(levs_),   yecm(levs_,levs_),tecm(levs_,levs_)
      integer k,kk,j,irow,icol,icolbeg,icolend
      real am(levs_,levs_),bm(levs_,levs_),sv(levs_),
     $     gv(levs_),cm(levs_,levs_)
      real ak5(levs_+1),bk5(levs_+1)
      integer kmx
 
!     print *,' enter get_am_bm_hyb_fd ' 		! hmhj

      do k=1,levs_
       tref(k)=300.
      enddo
      psref=80.
      beta=1.
      kappa=rd_/cp_

cprint print *,' begin  am_bm '
 
      do k=1,levs_+1
       pk5ref(k)=ak5(k)+bk5(k)*psref
cprint print 100,k,ak5(k),bk5(k),pk5ref(k)
      enddo
 
100   format('k=',i2,2x,'ak5=',e10.3,2x,'bk5=',e10.3,2x,'pk5ref=',e10.3)
cprintprint*,'-------------------------------------------------------'
 
      do k=1,levs_
       dpkref(k)=pk5ref(k+1)-pk5ref(k)
       gv(k)=beta*rd_*tref(k)/(rerth_*rerth_)
cprint print 110,k,dpkref(k),gv(k)
      enddo
 
110   format('k=',i2,2x,' in am_bm dpkref=',e11.4,2x,'gv=',e11.4)
cprintprint*,'-------------------------------------------------------'
cprintprint*,'-  calculate alfaref  watch alfaref(1)     '
 
      alfaref(1)=log(2.) ! could also be=1.  but watch for layer values
 
cprintprint*,'-------------------------------------------------------'
      do k=2,levs_
       alfaref(k)=1.-(pk5ref(k)/dpkref(k))*log(pk5ref(k+1)/pk5ref(k))
cprintprint 210,k,k,k,k,k
210   format('alfa(',i2,')=1.-(pk5(',i2,')/dpk(',i2,'))*log(pk5(',i2,
     & '+1)/pk5(',i2,'))')
      enddo
 
!sela print 125,alfaref(1)
125   format('worry --- alfaref(1)=',e10.3)
      do k=1,levs_
cprint print 130,k,alfaref(k)
      enddo
130   format('k=',i2,2x,'alfaref',e16.8)
cprintprint*,'---- begin matrices computation -----------'
 
c     print*,'333333333333333333333333333333333333333333'
cprintprint 144
144   format('begin yecm computation')
       yecm=0.
       do irow=1,levs_
          yecm(irow,irow)=alfaref(irow)*rd_
          icolbeg=irow+1
          if(icolbeg.le.levs_)then
           do icol=icolbeg,levs_
            yecm(irow,icol)=rd_*log( pk5ref(icol+1)/pk5ref(icol) )
           enddo
          endif
       enddo
150    format('yecm(',i2,',',i2,')=rd_*log( pk5ref(',i2,
     &        '+1)/pk5ref(',i2,'))')
c     print*,'-----------------1234567------------------'
160    format('yecm=',4(1x,e10.3))
 
      tecm=0.
 
      do irow=1,levs_
c       print*,' doing row ...............................',irow
         tecm(irow,irow)=kappa*tref(irow)*alfaref(irow)
         icolend=irow-1
 
      do icol=1,icolend
      factor=(kappa*tref(irow)/
     &                   dpkref(irow))*log(pk5ref(irow+1)/pk5ref(irow))
      tecm(irow,icol)=factor*dpkref(icol)
      enddo
      enddo
165    format('irow=',i2,2x,'factor=',e16.8,2x,'icolend=',i2)
166    format('factor=(kappa*tref/dpkref(',i2,'))*log(pk5ref(',i2,
     & '+1)/pk5ref(',i2,'))')
167    format('innerlup irow=',i2,2x,'icol=',i2,2x,'tecm(ir,ic)=',e12.4)
c     print*,'4444444  print yecm      44444444444444444'
 
       do irow=1,levs_
c       print*,'yecm row=',irow,'levs_=',levs_
c       print 1700,(yecm(irow,j),j=1,levs_/2)
c       print 1701,(yecm(irow,j),j=levs_/2+1,levs_)
       enddo
1700   format('  a  ',10(1x,e10.3))
1701   format('  b  ',10(1x,e10.3))
 
c     print*,'5555555  print tecm      55555555555555555'
 
       do irow=1,levs_
c       print*,'tecm row=',irow,'levs_=',levs_
c       print 1700,(tecm(irow,j),j=1,levs_/2)
c       print 1701,(tecm(irow,j),j=levs_/2+1,levs_)
       enddo
 
c     print*,'666666666666666666666666666666666666666666'
cprintprint 171
171   format('begin vvec dcomputation')
       do icol=1,levs_
        vecm(icol)=dpkref(icol)/psref
       enddo
       do icol=1,levs_
cprint  print 175,icol,vecm(icol)
       enddo
175    format('icol=',i2,2x,'vecm=',e16.8)
 
 
 
      do j=1,levs_
       sv(j)=vecm(levs_+1-j)
      do k=1,levs_
        am(k,j)=yecm(levs_+1-k,levs_+1-j)
        bm(k,j)=tecm(levs_+1-k,levs_+1-j)
      enddo
      enddo
 
      do j=1,levs_
      do k=1,levs_
        am(k,j)=am(k,j)*beta/(rerth_*rerth_)
      enddo
      enddo

!     print *,' leave get_am_bm_hyb_fd ' 		! hmhj
 
      do j=1,levs_
        do k=1,levs_
          cm(k,j)=gv(k)*sv(j)
        enddo
        do k=1,levs_
          do i=1,levs_
            cm(k,j)=cm(k,j)+am(k,i)*bm(i,j)
          enddo
        enddo
      enddo
c
      return
      end
