      subroutine h2oexps(ipts,ib,m,np,dh2o,pa,dt,xkw,aw,bw,pm,mw,h2oexp)
c**********************************************************************
c   compute exponentials for water vapor line absorption
c   in individual layers.
c
c---- input parameters
c  spectral band (ib)
c  number of grid intervals in zonal direction (m)
c  number of grid intervals in meridional direction (n)
c  number of layers (np)
c  layer water vapor amount for line absorption (dh2o) 
c  layer pressure (pa)
c  layer temperature minus 250k (dt)
c  absorption coefficients for the first k-distribution
c     function due to h2o line absorption (xkw)
c  coefficients for the temperature and pressure scaling (aw,bw,pm)
c  ratios between neighboring absorption coefficients for
c     h2o line absorption (mw)
c
c---- output parameters
c  6 exponentials for each layer  (h2oexp)
c
c**********************************************************************
      implicit none
      integer ib,m,n,np,i,k,ik,ipts

c---- input parameters ------

      real dh2o(m,np),pa(m,np),dt(m,np)

c---- output parameters -----

      real h2oexp(m,np,6)

c---- static data -----

      integer mw(9)
      real xkw(9),aw(9),bw(9),pm(9)

c---- temporary arrays -----

      real xh,xh1

c**********************************************************************
c    note that the 3 sub-bands in band 3 use the same set of xkw, aw,
c    and bw.  therefore, h2oexp for these sub-bands are identical.
c**********************************************************************

        do k=1,np
          do i=1,ipts   

c-----xh is   the scaled water vapor amount for line absorption
c     computed from (27).

           xh = dh2o(i,k)*(pa(i,k)/500.)**pm(ib)
     1        * ( 1.+(aw(ib)+bw(ib)* dt(i,k))*dt(i,k) )

c-----h2oexp is the water vapor transmittance of the layer (k2-1)
c     due to line absorption

           h2oexp(i,k,1) = exp(-xh*xkw(ib))

          enddo
        enddo

        do ik=2,6

         if (mw(ib).eq.6) then

          do k=1,np
            do i=1,ipts   
             xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             h2oexp(i,k,ik) = xh*xh*xh
            enddo
          enddo

        elseif (mw(ib).eq.8) then

          do k=1,np
            do i=1,ipts   
             xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             xh = xh*xh
             h2oexp(i,k,ik) = xh*xh
            enddo
          enddo

        elseif (mw(ib).eq.9) then

          do k=1,np
            do i=1,ipts   
             xh=h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             xh1 = xh*xh
             h2oexp(i,k,ik) = xh*xh1
            enddo
          enddo

        else

          do k=1,np
            do i=1,ipts   
             xh = h2oexp(i,k,ik-1)*h2oexp(i,k,ik-1)
             xh = xh*xh
             xh = xh*xh
             h2oexp(i,k,ik) = xh*xh
            enddo
          enddo

        endif
       enddo

      return
      end
