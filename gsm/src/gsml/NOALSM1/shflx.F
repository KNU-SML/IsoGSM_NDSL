      subroutine shflx (ssoil,stc,smc,smcmax,nsoil,t1,dt,yy,zz1,zsoil,
     &                  tbot,zbot,smcwlt,psisat,sh2o,bexp,f1,df1,ice,
     &                  quartz,csoil)
c      
      implicit none
c      
c ----------------------------------------------------------------------
c subroutine shflx
c ----------------------------------------------------------------------
c update the temperature state of the soil column based on the thermal
c diffusion equation and update the frozen soil moisture content based
c on the temperature.
c ----------------------------------------------------------------------
c ** clu_rev:  fix 'variable nsold has undefined type' error here
      integer nsold
      parameter(nsold = 20)
c
      integer i
      integer ice
      integer ifrz
      integer nsoil
!clu  integer nsold
      real ai(nsold)
      real bi(nsold)
      real ci(nsold)
c
      real bexp
      real csoil
      real df1
      real dt
      real f1
      real psisat
      real quartz
      real rhsts(nsold)
      real ssoil
      real sh2o(nsoil)
      real smc(nsoil)
      real smcmax
      real smcwlt
      real stc(nsoil)
      real stcf(nsold)
      real t0
      real t1
      real tbot
      real yy
      real zbot
      real zsoil(nsoil)
      real zz1
c
      parameter(t0 = 273.15)
c ----------------------------------------------------------------------
c hrt routine calcs the right hand side of the soil temp dif eqn
c ----------------------------------------------------------------------
      if (ice.eq.1) then
c ----------------------------------------------------------------------
c sea-ice case
c ----------------------------------------------------------------------
         call hrtice (rhsts,stc,nsoil,zsoil,yy,zz1,df1,ai,bi,ci)
         call hstep (stcf,stc,rhsts,dt,nsoil,ai,bi,ci)
      else
c ----------------------------------------------------------------------
c land-mass case
c ----------------------------------------------------------------------
         call hrt (rhsts,stc,smc,smcmax,nsoil,zsoil,yy,zz1,tbot,
     &             zbot,psisat,sh2o,dt,
     &             bexp,f1,df1,quartz,csoil,ai,bi,ci)
         call hstep (stcf,stc,rhsts,dt,nsoil,ai,bi,ci)
      endif
      do i = 1,nsoil
         stc(i) = stcf(i)
      end do
c ----------------------------------------------------------------------
c in the no snowpack case (via routine nopac branch,) update the grnd
c (skin) temperature here in response to the updated soil temperature 
c profile above.  (note: inspection of routine snopac shows that t1
c below is a dummy variable only, as skin temperature is updated
c differently in routine snopac) 
c ----------------------------------------------------------------------
      t1 = (yy + (zz1 - 1.0) * stc(1)) / zz1
c ----------------------------------------------------------------------
c calculate surface soil heat flux
c ----------------------------------------------------------------------
      ssoil = df1 * (stc(1) - t1) / (0.5 * zsoil(1))
c ----------------------------------------------------------------------
c end subroutine shflx
c ----------------------------------------------------------------------
      return
      end
