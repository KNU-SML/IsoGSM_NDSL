      subroutine b10exps(ipts,m,np,dh2o,dcont,dco2,dn2o,pa,dt
     *          ,h2oexp,conexp,co2exp,n2oexp)
c**********************************************************************
c   compute band3a exponentials for individual layers.
c
c---- input parameters
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer h2o amount for line absorption (dh2o)
c  layer h2o amount for continuum absorption (dcont)
c  layer co2 amount (dco2)
c  layer n2o amount (dn2o)
c  layer pressure (pa)
c  layer temperature minus 250k (dt)
c
c---- output parameters
c
c  exponentials for each layer (h2oexp,conexp,co2exp,n2oexp)
c
c**********************************************************************
      implicit none
      integer m,n,np,i,k,ipts

c---- input parameters -----

      real dh2o(m,np),dcont(m,np),dn2o(m,np)
      real dco2(m,np),pa(m,np),dt(m,np)

c---- output parameters -----

      real h2oexp(m,np,6),conexp(m,np,3),co2exp(m,np,6,2)
     *      ,n2oexp(m,np,4)

c---- temporary arrays -----

      real xx,xx1,xx2,xx3

c**********************************************************************

        do k=1,np
          do i=1,ipts   

c-----compute the scaled h2o-line amount for the sub-band 3a
c     table 1, chou et al. (jas, 1993)

           xx=dh2o(i,k)*(pa(i,k)/500.0)
     1           *(1.+(0.0149+6.20e-5*dt(i,k))*dt(i,k))

c-----six exponentials by powers of 8

           h2oexp(i,k,1)=exp(-xx*0.10624)

           xx=h2oexp(i,k,1)*h2oexp(i,k,1)
           xx=xx*xx
           h2oexp(i,k,2)=xx*xx

           xx=h2oexp(i,k,2)*h2oexp(i,k,2)
           xx=xx*xx
           h2oexp(i,k,3)=xx*xx

           xx=h2oexp(i,k,3)*h2oexp(i,k,3)
           xx=xx*xx
           h2oexp(i,k,4)=xx*xx

           xx=h2oexp(i,k,4)*h2oexp(i,k,4)
           xx=xx*xx
           h2oexp(i,k,5)=xx*xx

           xx=h2oexp(i,k,5)*h2oexp(i,k,5)
           xx=xx*xx
           h2oexp(i,k,6)=xx*xx

c-----compute the scaled co2 amount for the sub-band 3a
c     table 1, chou et al. (jas, 1993)

           xx=dco2(i,k)*(pa(i,k)/300.0)**0.5
     1           *(1.+(0.0179+1.02e-4*dt(i,k))*dt(i,k))

c-----six exponentials by powers of 8

           co2exp(i,k,1,1)=exp(-xx*2.656e-5)

           xx=co2exp(i,k,1,1)*co2exp(i,k,1,1)
           xx=xx*xx
           co2exp(i,k,2,1)=xx*xx

           xx=co2exp(i,k,2,1)*co2exp(i,k,2,1)
           xx=xx*xx
           co2exp(i,k,3,1)=xx*xx

           xx=co2exp(i,k,3,1)*co2exp(i,k,3,1)
           xx=xx*xx
           co2exp(i,k,4,1)=xx*xx

           xx=co2exp(i,k,4,1)*co2exp(i,k,4,1)
           xx=xx*xx
           co2exp(i,k,5,1)=xx*xx

           xx=co2exp(i,k,5,1)*co2exp(i,k,5,1)
           xx=xx*xx
           co2exp(i,k,6,1)=xx*xx

c-----one exponential of h2o continuum for sub-band 3a

           conexp(i,k,1)=exp(-dcont(i,k)*1.04995e+2)

c-----compute the scaled n2o amount for sub-band 3a

           xx=dn2o(i,k)*(1.+(1.4476e-3+3.6656e-6*dt(i,k))*dt(i,k))

c-----two exponential2 by powers of 58

           n2oexp(i,k,1)=exp(-xx*0.22646)

           xx=n2oexp(i,k,1)*n2oexp(i,k,1)
           xx1=xx*xx
           xx1=xx1*xx1
           xx2=xx1*xx1
           xx3=xx2*xx2
           n2oexp(i,k,2)=xx*xx1*xx2*xx3

          enddo
        enddo

      return
      end
