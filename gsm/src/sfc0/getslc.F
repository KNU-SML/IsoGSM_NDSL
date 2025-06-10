#define G g
#define HFUS hfus
#define LSOIL lsoil
#define T0C t0c
      subroutine  getslc(slmsk,stype,stc,smc,ijdim,lsoil,slc)
c
      implicit none
#include <paramodel.h>
c
c  initialie slc from smc and stc
c
      integer ijdim,lsoil
      real slmsk(ijdim),stype(ijdim),stc(ijdim,lsoil),smc(ijdim,lsoil)
      real slc(ijdim,lsoil)
c
      real blim,bx,fk,frh2o
      real psis(9),beta(9),smcmax(9)
      integer isltpk
      data blim/5.5/
      data psis/0.04,0.62,0.47,0.14,0.10,0.26,0.14,0.36,0.04/
      data beta/4.26,8.72,11.55,4.74,10.73,8.17,6.77,5.25,4.26/
      data smcmax/0.421,0.464,0.468,0.434,0.406,
     &            0.465,0.404,0.439,0.421/
c
      integer ij,k
c
      do ij=1,ijdim
        if(slmsk(ij) .ne. 0.) then
          isltpk=int(stype(ij)+.5)
          do k=1,lsoil_
            if (stc(ij,k).lt.273.149) then
c
c sh2o <= smc for t < 273.149k (-0.001c)
c
c first guess following explicit solution for flerchinger eqn from koren
c et al, jgr, 1999, eqn 17 (kcount=0 in function frh2o).
c
              bx = beta(isltpk)
              if ( beta(isltpk) .gt. blim ) bx = blim
                fk=(((hfus_/(g_*(-psis(isltpk))))*
     1           ((stc(ij,k)-t0c_)/stc(ij,k)))**
     1           (-1/bx))*smcmax(isltpk)
              if (fk .lt. 0.02) fk = 0.02
              slc(ij,k) = min ( fk, smc(ij,k) )
c
c now use iterative solution for liquid soil water content using
c function frh2o with the initial guess for sh2o from above explicit
c first guess.
c
              slc(ij,k)=frh2o(stc(ij,k),
     +             smc(ij,k), slc(ij,k),
     .             smcmax(isltpk),beta(isltpk),
     .             psis(isltpk))
            else
c
c sh2o = smc for t => 273.149k (-0.001c)
c
              slc(ij,k)=smc(ij,k)
            endif
          enddo
        endif
      enddo
c
      return
      end
