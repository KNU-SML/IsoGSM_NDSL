      subroutine column (ipts,m,np,pa,dt,sabs0,sabs,spre,stem)
c
c**************************************************************************
c compute column-integrated (from top of the model atmosphere)
c     absorber amount (sabs), absorber-weighted pressure (spre) and
c     temperature (stem).
c     computations of spre and stem follows eqs. (37) and (38).
c
c input parameters
c   number of soundings in zonal direction (m)
c   number of soundings in meridional direction (n)
c   number of atmospheric layers (np)
c   layer pressure (pa)
c   layer temperature minus 250k (dt)
c   layer absorber amount (sabs0)
c
c output parameters
c   column-integrated absorber amount (sabs)
c   column absorber-weighted pressure (spre)
c   column absorber-weighted temperature (stem)
c
c units of pa and dt are mb and k, respectively.
c    units of sabs are g/cm**2 for water vapor and (cm-atm)stp for co2 and o3
c
      implicit none
      integer m,n,np,i,k,ipts
c
c input parameters
c
      real pa(m,np),dt(m,np),sabs0(m,np)
c
c  output parameters
c
      real sabs(m,np+1),spre(m,np+1),stem(m,np+1)
c
         do i=1,ipts
          sabs(i,1)=0.0
          spre(i,1)=0.0
          stem(i,1)=0.0
         enddo

        do k=1,np
          do i=1,ipts
           sabs(i,k+1)=sabs(i,k)+sabs0(i,k)
           spre(i,k+1)=spre(i,k)+pa(i,k)*sabs0(i,k)
           stem(i,k+1)=stem(i,k)+dt(i,k)*sabs0(i,k)
          enddo
        enddo
       return
       end
