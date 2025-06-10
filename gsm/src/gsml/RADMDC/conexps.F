      subroutine conexps(ipts,ib,m,np,dcont,xke,ne,conexp)
c**********************************************************************
c   compute exponentials for continuum absorption in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer scaled water vapor amount for continuum absorption (dcont)
c  absorption coefficients for the first k-distribution function
c     due to water vapor continuum absorption (xke)
c  number of terms used in each band to compute h2o continuum
c     transmittance (ne)
c
c---- output parameters
c  1 or 3 exponentials for each layer (conexp)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np,i,k,iq,ipts

c---- input parameters ------

      real dcont(m,np)

c---- updated parameters -----

      real conexp(m,np,3)

c---- static data -----

      integer ne(9)
      real xke(9)

c**********************************************************************

        do k=1,np
          do i=1,ipts
           conexp(i,k,1) = exp(-dcont(i,k)*xke(ib))
          enddo
        enddo

       if (ib .eq. 3) then

c-----the absorption coefficients for sub-bands 3b (iq=2) and 3a (iq=3)
c     are, respectively, double and quadruple that for sub-band 3c (iq=1)
c     (table 6).

        do iq=2,3
         do k=1,np
           do i=1,ipts
            conexp(i,k,iq) = conexp(i,k,iq-1) *conexp(i,k,iq-1)
           enddo
         enddo
        enddo

       endif

      return
      end
