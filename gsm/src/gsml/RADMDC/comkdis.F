      subroutine comkdis(ipts,ib,m,np,k,comexp,tcom,tran)
c**********************************************************************
c   compute co2-minor transmittances between levels k1 and k2 for m x n
c   soundings using the k-distribution method with linear pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals in zonal direction (m)
c   number of grid intervals in meridional direction (n)
c   number of levels (np)
c   current level (k)
c   exponentials for co2-minor absorption (comexp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to co2-minor absorption
c     for the various values of the absorption coefficient (tcom)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np,k,i,ipts

c---- input parameters -----

      real comexp(m,np,2)

c---- updated parameters -----

      real tcom(m,2),tran(m)

c---- temporary arrays -----

      real xc

c-----tcom is the 2 exp factors between levels k1 and k2. 
c     xc is the total co2-minor transmittance

          do i=1,ipts   

c-----band 4

           if (ib.eq.4) then

            tcom(i,1)=tcom(i,1)*comexp(i,k,1)
            xc=   0.972025*tcom(i,1)
            tcom(i,2)=tcom(i,2)*comexp(i,k,2)
            xc=xc+0.027975*tcom(i,2)

c-----band 5

           else

            tcom(i,1)=tcom(i,1)*comexp(i,k,1)
            xc=   0.961324*tcom(i,1)
            tcom(i,2)=tcom(i,2)*comexp(i,k,2)
            xc=xc+0.038676*tcom(i,2)

           endif

            tran(i)=tran(i)*xc

          enddo

      return
      end
