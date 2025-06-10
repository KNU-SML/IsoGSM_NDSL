      subroutine comexps(ipts,ib,m,np,dcom,dt,comexp)
c**********************************************************************
c   compute co2-minor exponentials for individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer co2 amount (dcom)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c  2 exponentials for each layer (comexp)
c**********************************************************************
      implicit none
      integer ib,m,n,np,i,k,ipts

c---- input parameters -----

      real dcom(m,np),dt(m,np)

c---- output parameters -----

      real comexp(m,np,2)

c---- temporary arrays -----

      real xc,xc1,xc2

c**********************************************************************

       do k=1,np
         do i=1,ipts   

c-----two exponentials by powers of 60 for band 4

          if (ib.eq.4) then

           xc=dcom(i,k)*(1.+(3.5775e-2+4.0447e-4*dt(i,k))*dt(i,k))
           comexp(i,k,1)=exp(-xc*1.95404e-5)

           xc=comexp(i,k,1)*comexp(i,k,1)*comexp(i,k,1)
           xc=xc*xc
           xc1=xc*xc
           xc=xc1*xc1
           xc=xc*xc
           comexp(i,k,2)=xc*xc1

c-----two exponentials by powers of 44 for band 5

          else

           xc=dcom(i,k)*(1.+(3.4268e-2+3.7401e-4*dt(i,k))*dt(i,k))
           comexp(i,k,1)=exp(-xc*4.25830e-5)

           xc=comexp(i,k,1)*comexp(i,k,1)
           xc1=xc*xc
           xc2=xc1*xc1
           xc=xc2*xc2
           xc=xc*xc
           comexp(i,k,2)=xc1*xc2*xc

          endif

         enddo
       enddo

      return
      end
