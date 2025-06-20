      subroutine ch4kdis(ipts,ib,m,np,k,ch4exp,tch4,tran)
c**********************************************************************
c   compute ch4 transmittances between levels k1 and k2 for m x n soundings
c   using the k-distribution method with linear pressure scaling.
c
c---- input parameters
c   spectral band (ib)
c   number of grid intervals in zonal direction (m)
c   number of grid intervals in meridional direction (n)
c   number of levels (np)
c   current level (k)
c   exponentials for ch4 absorption (ch4exp)
c
c---- updated parameters
c   transmittance between levels k1 and k2 due to ch4 absorption
c     for the various values of the absorption coefficient (tch4)
c   total transmittance (tran)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np,k,i,ipts

c---- input parameters -----

      real ch4exp(m,np,4)

c---- updated parameters -----

      real tch4(m,4),tran(m)

c---- temporary arrays -----

      real xc

c-----tch4 is the 2 exp factors between levels k1 and k2. 
c     xc is the total ch4 transmittance

         do i=1,ipts   

c-----band 6

          if (ib.eq.6) then

           tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
           xc= tch4(i,1)

c-----band 7

          else

           tch4(i,1)=tch4(i,1)*ch4exp(i,k,1)
           xc=   0.610650*tch4(i,1)

           tch4(i,2)=tch4(i,2)*ch4exp(i,k,2)
           xc=xc+0.280212*tch4(i,2)

           tch4(i,3)=tch4(i,3)*ch4exp(i,k,3)
           xc=xc+0.107349*tch4(i,3)

           tch4(i,4)=tch4(i,4)*ch4exp(i,k,4)
           xc=xc+0.001789*tch4(i,4)

          endif

           tran(i)=tran(i)*xc

         enddo

      return
      end
