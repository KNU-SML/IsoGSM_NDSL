      subroutine hrtice (rhsts,stc,nsoil,zsoil,yy,zz1,df1,ai,bi,ci)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine hrtice
c ----------------------------------------------------------------------
c calculate the right hand side of the time tendency term of the soil
c thermal diffusion equation in the case of sea-ice pack.  also to
c compute (prepare) the matrix coefficients for the tri-diagonal matrix
c of the implicit time scheme.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      integer k
      integer nsoil
c
      real ai(nsold)
      real bi(nsold)
      real ci(nsold)
c
      real ddz
      real ddz2
      real denom
      real df1
      real dtsdz
      real dtsdz2
      real hcpct
      real rhsts(nsoil)
      real ssoil
      real stc(nsoil)
      real tbot
      real yy
      real zbot
      real zsoil(nsoil)
      real zz1
c
      data tbot /271.16/
c ----------------------------------------------------------------------
c set a nominal universal value of the sea-ice specific heat capacity,
c hcpct = 1880.0*917.0.
c ----------------------------------------------------------------------
      parameter(hcpct = 1.72396e+6)
c ----------------------------------------------------------------------
c the input argument df1 is a universally constant value of sea-ice
c thermal diffusivity, set in routine snopac as df1 = 2.2.
c ----------------------------------------------------------------------
c set ice pack depth.  use tbot as ice pack lower boundary temperature
c (that of unfrozen sea water at bottom of sea ice pack).  assume ice
c pack is of n=nsoil layers spanning a uniform constant ice pack
c thickness as defined by zsoil(nsoil) in routine sflx.
c ----------------------------------------------------------------------
      zbot = zsoil(nsoil)
c ----------------------------------------------------------------------
c calc the matrix coefficients ai, bi, and ci for the top layer
c ----------------------------------------------------------------------
      ddz = 1.0 / ( -0.5 * zsoil(2) )
      ai(1) = 0.0
      ci(1) = (df1 * ddz) / (zsoil(1) * hcpct)
      bi(1) = -ci(1) + df1/(0.5 * zsoil(1) * zsoil(1) * hcpct * zz1)
c ----------------------------------------------------------------------
c calc the vertical soil temp gradient btwn the top and 2nd soil layers.
c recalc/adjust the soil heat flux.  use the gradient and flux to calc
c rhsts for the top soil layer.
c ----------------------------------------------------------------------
      dtsdz = ( stc(1) - stc(2) ) / ( -0.5 * zsoil(2) )
      ssoil = df1 * ( stc(1) - yy ) / ( 0.5 * zsoil(1) * zz1 )
      rhsts(1) = ( df1 * dtsdz - ssoil ) / ( zsoil(1) * hcpct )
c ----------------------------------------------------------------------
c initialize ddz2
c ----------------------------------------------------------------------
      ddz2 = 0.0
c ----------------------------------------------------------------------
c loop thru the remaining soil layers, repeating the above process
c ----------------------------------------------------------------------
      do k = 2,nsoil
        if (k .ne. nsoil) then
c ----------------------------------------------------------------------
c calc the vertical soil temp gradient thru this layer.
c ----------------------------------------------------------------------
          denom = 0.5 * ( zsoil(k-1) - zsoil(k+1) )
          dtsdz2 = ( stc(k) - stc(k+1) ) / denom
c ----------------------------------------------------------------------
c calc the matrix coef, ci, after calc'ng its partial product.
c ----------------------------------------------------------------------
          ddz2 = 2. / (zsoil(k-1) - zsoil(k+1))
          ci(k) = -df1 * ddz2 / ((zsoil(k-1) - zsoil(k)) * hcpct)
        else
c ----------------------------------------------------------------------
c calc the vertical soil temp gradient thru the lowest layer.
c ----------------------------------------------------------------------
          dtsdz2 = (stc(k)-tbot)/(.5 * (zsoil(k-1) + zsoil(k))-zbot)
c ----------------------------------------------------------------------
c set matrix coef, ci to zero.
c ----------------------------------------------------------------------
          ci(k) = 0.
        endif
c ----------------------------------------------------------------------
c calc rhsts for this layer after calc'ng a partial product.
c ----------------------------------------------------------------------
        denom = ( zsoil(k) - zsoil(k-1) ) * hcpct
        rhsts(k) = ( df1 * dtsdz2 - df1 * dtsdz ) / denom
c ----------------------------------------------------------------------
c calc matrix coefs, ai, and bi for this layer.
c ----------------------------------------------------------------------
        ai(k) = - df1 * ddz / ((zsoil(k-1) - zsoil(k)) * hcpct)
        bi(k) = -(ai(k) + ci(k))
c ----------------------------------------------------------------------
c reset values of dtsdz and ddz for loop to next soil lyr.
c ----------------------------------------------------------------------
        dtsdz = dtsdz2
        ddz   = ddz2
      end do
c ----------------------------------------------------------------------
c end subroutine hrtice
c ----------------------------------------------------------------------
      return
      end
