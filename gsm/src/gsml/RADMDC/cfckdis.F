      subroutine cfckdis(ipts,m,np,k,cfcexp,tcfc,tran)
c**********************************************************************
c   compute cfc-(11,12,22) transmittances between levels k1 and k2 for m x n
c   soundings using the k-distribution method with linear pressure scaling.
c
c---- input parameters
c   number of grid intervals in zonal direction (m)
c   number of grid intervals in meridional direction (n)
c   number of levels (np)
c   current level (k)
c   exponentials for cfc absorption (cfcexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to cfc absorption
c     for the various values of the absorption coefficient (tcfc)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer m,n,np,k,i,ipts

c---- input parameters -----

      real cfcexp(m,np)

c---- updated parameters -----

      real tcfc(m),tran(m)

c-----tcfc is the exp factors between levels k1 and k2. 

          do i=1,ipts   

            tcfc(i)=tcfc(i)*cfcexp(i,k)
            tran(i)=tran(i)*tcfc(i)

          enddo

      return
      end
