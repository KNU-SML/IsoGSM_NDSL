      subroutine hrt (rhsts,stc,smc,smcmax,nsoil,zsoil,yy,zz1,
     &                tbot,zbot,psisat,sh2o,dt,bexp,
     &                f1,df1,quartz,csoil,ai,bi,ci)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine hrt
c ----------------------------------------------------------------------
c calculate the right hand side of the time tendency term of the soil
c thermal diffusion equation.  also to compute ( prepare ) the matrix
c coefficients for the tri-diagonal matrix of the implicit time scheme.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      logical itavg
c
      integer i
      integer k
      integer nsoil
c ----------------------------------------------------------------------
c declare work arrays needed in tri-diagonal implicit solver
c ----------------------------------------------------------------------
      real ai(nsold)
      real bi(nsold)
      real ci(nsold)
c ----------------------------------------------------------------------
c declarations
c ----------------------------------------------------------------------
      real bexp
      real cair
      real ch2o
      real cice
      real csoil
      real ddz
      real ddz2
      real denom
      real df1
      real df1n
      real df1k
      real dt
      real dtsdz
      real dtsdz2
      real f1
      real hcpct
      real psisat
      real quartz
      real qtot
      real rhsts(nsoil)
      real ssoil
      real sice
      real smc(nsoil)
      real sh2o(nsoil)
      real smcmax
      real snksrc
      real stc(nsoil)
      real t0
      real tavg
      real tbk
      real tbk1
      real tbot
      real zbot
      real tsnsr
      real tsurf
      real yy
      real zsoil(nsoil)
      real zz1
c
      parameter(t0 = 273.15)
c ----------------------------------------------------------------------
c set specific heat capacities of air, water, ice, soil mineral       
c ----------------------------------------------------------------------
      parameter(cair = 1004.0)
      parameter(ch2o = 4.2e6)
      parameter(cice = 2.106e6)
c note: csoil now set in routine redprm and passed in
c      parameter(csoil = 1.26e6)
c ----------------------------------------------------------------------
c initialize logical for soil layer temperature averaging.
c ----------------------------------------------------------------------
      itavg = .true.
c      itavg = .false.
c ----------------------------------------------------------------------
c begin section for top soil layer
c ----------------------------------------------------------------------
c calc the heat capacity of the top soil layer
c ----------------------------------------------------------------------
      hcpct = sh2o(1)*ch2o + (1.0-smcmax)*csoil + (smcmax-smc(1))*cair
     &        + ( smc(1) - sh2o(1) )*cice
c ----------------------------------------------------------------------
c calc the matrix coefficients ai, bi, and ci for the top layer
c ----------------------------------------------------------------------
      ddz = 1.0 / ( -0.5 * zsoil(2) )
      ai(1) = 0.0
      ci(1) = (df1 * ddz) / (zsoil(1) * hcpct)
      bi(1) = -ci(1) + df1 / (0.5 * zsoil(1) * zsoil(1)*hcpct*zz1)
c ----------------------------------------------------------------------
c calculate the vertical soil temp gradient btwn the 1st and 2nd soil
c layers.  then calculate the subsurface heat flux. use the temp
c gradient and subsfc heat flux to calc "right-hand side tendency
c terms", or "rhsts", for top soil layer.
c ----------------------------------------------------------------------
      dtsdz = (stc(1) - stc(2)) / (-0.5 * zsoil(2))
      ssoil = df1 * (stc(1) - yy) / (0.5 * zsoil(1) * zz1)
      rhsts(1) = (df1 * dtsdz - ssoil) / (zsoil(1) * hcpct)
c ----------------------------------------------------------------------
c next capture the vertical difference of the heat flux at top and
c bottom of first soil layer for use in heat flux constraint applied to
c potential soil freezing/thawing in routine snksrc.
c ----------------------------------------------------------------------
      qtot = ssoil - df1*dtsdz
c ----------------------------------------------------------------------
c if temperature averaging invoked (itavg=true; else skip):
c set temp "tsurf" at top of soil column (for use in freezing soil
c physics later in function subroutine snksrc).  if snowpack content is
c zero, then tsurf expression below gives tsurf = skin temp.  if
c snowpack is nonzero (hence argument zz1=1), then tsurf expression
c below yields soil column top temperature under snowpack.  then
c calculate temperature at bottom interface of 1st soil layer for use
c later in function subroutine snksrc
c ----------------------------------------------------------------------
      if (itavg) then 
        tsurf = (yy + (zz1-1) * stc(1)) / zz1
        call tbnd (stc(1),stc(2),zsoil,zbot,1,nsoil,tbk)
      endif
c ----------------------------------------------------------------------
c calculate frozen water content in 1st soil layer. 
c ----------------------------------------------------------------------
      sice = smc(1) - sh2o(1)
c ----------------------------------------------------------------------
c if frozen water present or any of layer-1 mid-point or bounding
c interface temperatures below freezing, then call snksrc to
c compute heat source/sink (and change in frozen water content)
c due to possible soil water phase change
c ----------------------------------------------------------------------
      if ( (sice   .gt. 0.) .or. (tsurf .lt. t0) .or.
     &     (stc(1) .lt. t0) .or. (tbk   .lt. t0) ) then

        if (itavg) then 
          call tmpavg(tavg,tsurf,stc(1),tbk,zsoil,nsoil,1)
        else
          tavg = stc(1)
        endif
        tsnsr = snksrc (tavg,smc(1),sh2o(1), 
     &    zsoil,nsoil,smcmax,psisat,bexp,dt,1,qtot)
c
        rhsts(1) = rhsts(1) - tsnsr / ( zsoil(1) * hcpct )
      endif
c ----------------------------------------------------------------------
c this ends section for top soil layer.
c ----------------------------------------------------------------------
c initialize ddz2
c ----------------------------------------------------------------------
      ddz2 = 0.0
c ----------------------------------------------------------------------
c loop thru the remaining soil layers, repeating the above process
c (except subsfc or "ground" heat flux not repeated in lower layers)
c ----------------------------------------------------------------------
      df1k = df1
      do k = 2,nsoil
c ----------------------------------------------------------------------
c calculate heat capacity for this soil layer.
c ----------------------------------------------------------------------
        hcpct = sh2o(k)*ch2o +(1.0-smcmax)*csoil +(smcmax-smc(k))*cair
     &        + ( smc(k) - sh2o(k) )*cice
        if (k .ne. nsoil) then
c ----------------------------------------------------------------------
c this section for layer 2 or greater, but not last layer.
c ----------------------------------------------------------------------
c calculate thermal diffusivity for this layer.
c ----------------------------------------------------------------------
          call tdfcnd (df1n,smc(k),quartz,smcmax,sh2o(k))
c ----------------------------------------------------------------------
c calc the vertical soil temp gradient thru this layer
c ----------------------------------------------------------------------
          denom = 0.5 * ( zsoil(k-1) - zsoil(k+1) )
          dtsdz2 = ( stc(k) - stc(k+1) ) / denom
c ----------------------------------------------------------------------
c calc the matrix coef, ci, after calc'ng its partial product
c ----------------------------------------------------------------------
          ddz2 = 2. / (zsoil(k-1) - zsoil(k+1))
          ci(k) = -df1n * ddz2 / ((zsoil(k-1) - zsoil(k)) * hcpct)
c ----------------------------------------------------------------------
c if temperature averaging invoked (itavg=true; else skip):  calculate
c temp at bottom of layer.
c ----------------------------------------------------------------------
          if (itavg) then 
            call tbnd (stc(k),stc(k+1),zsoil,zbot,k,nsoil,tbk1)
          endif
        else
c ----------------------------------------------------------------------
c special case of bottom soil layer:  calculate thermal diffusivity for
c bottom layer.
c ----------------------------------------------------------------------
          call tdfcnd (df1n,smc(k),quartz,smcmax,sh2o(k))
c ----------------------------------------------------------------------
c calc the vertical soil temp gradient thru bottom layer.
c ----------------------------------------------------------------------
          denom = .5 * (zsoil(k-1) + zsoil(k)) - zbot
          dtsdz2 = (stc(k)-tbot) / denom
c ----------------------------------------------------------------------
c set matrix coef, ci to zero if bottom layer.
c ----------------------------------------------------------------------
          ci(k) = 0.
c ----------------------------------------------------------------------
c if temperature averaging invoked (itavg=true; else skip):  calculate
c temp at bottom of last layer.
c ----------------------------------------------------------------------
          if (itavg) then 
            call tbnd (stc(k),tbot,zsoil,zbot,k,nsoil,tbk1)
          endif 
        endif
c ----------------------------------------------------------------------
c this ends special loop for bottom layer.
c ----------------------------------------------------------------------
c calculate rhsts for this layer after calc'ng a partial product.
c ----------------------------------------------------------------------
        denom = ( zsoil(k) - zsoil(k-1) ) * hcpct
        rhsts(k) = ( df1n * dtsdz2 - df1k * dtsdz ) / denom
        qtot = -1.0*denom*rhsts(k)
        sice = smc(k) - sh2o(k)
        if ( (sice .gt. 0.) .or. (tbk .lt. t0) .or.
     &     (stc(k) .lt. t0) .or. (tbk1 .lt. t0) ) then
          if (itavg) then 
            call tmpavg(tavg,tbk,stc(k),tbk1,zsoil,nsoil,k)
          else
            tavg = stc(k)
          endif
          tsnsr = snksrc(tavg,smc(k),sh2o(k),zsoil,nsoil,
     &                   smcmax,psisat,bexp,dt,k,qtot)
          rhsts(k) = rhsts(k) - tsnsr / denom
        endif 
c ----------------------------------------------------------------------
c calc matrix coefs, ai, and bi for this layer.
c ----------------------------------------------------------------------
        ai(k) = - df1 * ddz / ((zsoil(k-1) - zsoil(k)) * hcpct)
        bi(k) = -(ai(k) + ci(k))
c ----------------------------------------------------------------------
c reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.
c ----------------------------------------------------------------------
        tbk   = tbk1
        df1k  = df1n
        dtsdz = dtsdz2
        ddz   = ddz2
      end do
c ----------------------------------------------------------------------
c end subroutine hrt
c ----------------------------------------------------------------------
      return
      end
