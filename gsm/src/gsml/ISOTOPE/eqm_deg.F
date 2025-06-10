!
!----------------------------------------------
      subroutine eqm_deg(id,temp,prcp,p,dz,fe)

!  calculate the degree of equilibration
!----------------------------------------------------
      implicit none
      integer nr                ! number of bin
      parameter(nr=20)

!------------------------------Arguments--------------------------------
      real temp                 ! tempeature in K
      real p                    ! pressur in Pa
      real pi
      real prcp                 ! precipitation rate in kg/m2/s
      real dz                   ! thickness of the layer in m
      integer id
!     id=2:H2O18, id=3:HDO

      real fe                   ! degree of equilibration
                                ! 1 means 100% equilibrium
!--------------------working variables--------------------------------
      integer i
      real r2(nr)               ! diameter-(bin size) of the raindrop
      real r(nr)                ! radius of the raindrop
      real rt
      real mr(nr)               ! mass of the raindrop
      real atime(nr)            ! equilibration time of the raindrop
      real vt(nr)               ! velocity of raindrop
      real prate                ! precipitation rate prcp*3600
      real prsum                ! sum of raindrop frin M-P relationship
                                ! should be identical to prate, but
                                ! small difference exists 
      real xl,xn,prdt           ! needed parameter
      real eqm(nr)              ! 
      real xnd(nr)              ! # of raindrop/sec
      real xndsum               ! # of total raindrop
      real fesum                ! 
      real time
      real dd
      real hr
      real rhol
      real xn0
      data hr/3600./
      data rhol/1.e3/
      data dd/0.0002/
      data pi/3.141592/
      data xn0/8.e3/
!-----------------------------------------------------------------------
!
      prate=prcp*hr
      do i=1,nr
! for Marshall and Palmer distribution N(D)d(D) is the number of drops
! per unit volume with diamters between D and D+dD where D is the
! diameter of the raindrop
         r2(i)=0.0004*(i-1)+0.0002    ! m
         r(i)=(r2(i)+0.0002)/2.
         rt=r(i)*1000.
         vt(i)=-0.1927+4.96255*rt-0.90441*rt**2+0.05658*rt**3
!        call velocity(r(i),temp,p,vt(i))
         mr(i)=4./3.*pi*r(i)**3*rhol
         call equil_time(id,r(i),temp,p,vt(i),mr(i),atime(i))
      enddo
!------------------------------------------------------------------
!
!    Marshall-Palmer distribution
!
! N(D)=N0*exp(-L*D)           ! N0=0.08 cm-4
! L(R)=41*R**(-0.21)          ! cm-1, R=mm/h
! conversion to mks unit
! N(D)=8e6*exp(-4.1e3*(R**(-0.21))*D)
!
! for Marshall and Palmer distribution N(D)d(D) is the number of drops
! per unit volume with diamters between D and D+dD where D is the
! diameter of the raindrop
!
!
!------------------------------------------------------------------
!    xmass_r: calculte raindrop mass
      prsum=0.
      fesum=0.
      xndsum=0.
      do i=1,nr
         xl=4.1e3*prate**(-1.*0.21)               ! m-1
         xn=xn0*exp(-1.*xl*r2(i))*dd             ! m-3: N(D)*dD

! number of raindrop/m2 = N(D)*dD*v      ! #/m2/h

         xnd(i)=xn*vt(i)                         ! #/m2/s
         time=dz/vt(i)
         eqm(i)=exp(-1.*time/atime(i))*xnd(i)
         fesum=fesum+eqm(i)*mr(i)
         xndsum=xndsum+xnd(i)*mr(i)
!        write(*,*) 'EQM',r(i)*1000.,atime(i),eqm(i),xnd(i),prate,time
      enddo

      fe=max(fesum/xndsum,0.)
      fe=1.-fe
!        write(*,*) 'FE',fe
      fe=min(fe,1.)

      return
      end subroutine eqm_deg
!--------------------------------------------
      subroutine equil_time(id,r,temp,p,vt,mr,etime)
!  calculate the time to reach isotopic equilibrium
!-----------------------------------------------------------------
      implicit none
!------------------------------Arguments--------------------------------
      real temp                 ! tempeature in K
      real r                    ! radius of the raindrop
      real p                    ! pressure in Pa
      real mr                   ! mass of the raindrop
      real vt                   ! velocity of the raindrop
      real etime                ! equilibration time
      integer id                !id, 1:H2O, 2:H2O18, 3:HDO
!------------------------------working variable-----------------------

      real es                    ! saturation vapor pressure
      real tc                    ! t-273.15
      real ws
      real di
      real f
      real alpha
      real vent
      real vis
      real d0
      real p0                   ! pressure at STP
      real t0                   ! 273.15oC
      real t1                   ! 
      real a(5),b(4)
! viscosity
      real a1,a05,tr,c,vis0,dvis
      real h
      real rho_air              ! density of the air
      real rhor                 ! density ratio 
      real rho0                 ! 
      real r2
      real re
      real sc
      real x

      integer i
      data a/-1.,-0.709661,0.662534,-0.197846,0.00770147/
      real pi
      data pi/3.141592/
      real freq                 ! function
!-----------------------------------------------------------------------
      tc=temp-273.15
      es=611.2*exp(17.67*tc/(tc+243.5))
      ws=max(es/461./temp,1.e-12)

! calculate diffusivity
      d0=0.211*1.e-4
      t0=273.15
      p0=1.e5
      di=d0*(temp/273.15)**1.94*(p0/p) ! diffusivity of H2O16

      if (id .eq. 2) then
! Merlivat (1970)
         di=0.9723*di
! based on molecular weight
!     di=0.9843*di
!     di=di/1.0319
      elseif (id .eq. 3) then
! Merlivat (1970)
         di=0.9755*di
! based on molecular weight
!     di=0.9691*di
!     di=di/1.016
      else
         di=di
      endif
      alpha=freq(id-1,temp)
c      call frac_coeff_l(id,temp,alpha)

!calculate viscosity of the air
      h=6.1609*1.e-6
      
      a1=0.128517
      a05=2.60661
      t1=132.5
      tr=temp/t1

      c=0.
      do i=1,5
         c=c+a(i)/tr**float(i-1)
      enddo

      vis0=a1*tr+a05*tr**.5+c

      rho0=314.3
      rho_air=p/287./temp
      rhor=rho_air/rho0

      dvis=0.
      do i=1,5
         dvis=dvis+b(i)*rhor**float(i)
      enddo

      vis=h*(vis0+dvis)

      r2=r*2.
      sc=vis/rho_air/di
      re=rho_air*r2*vt/vis
      x=sc**(1./3.)*re**0.5

      if (x .lt. 1.4) then
         vent=1.+0.108*x**2
      else
         vent=0.78+0.308*x
      endif


      etime=alpha*mr/di/ws/4./pi/r/vent

      return
      end subroutine equil_time
