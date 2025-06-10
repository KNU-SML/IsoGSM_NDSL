      subroutine srt (rhstt,edir,et,sh2o,sh2oa,nsoil,pcpdrp,
     &                zsoil,dwsat,dksat,smcmax,bexp,runoff1, 
     &                runoff2,dt,smcwlt,slope,kdt,frzx,sice,ai,bi,ci)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine srt
c ----------------------------------------------------------------------
c calculate the right hand side of the time tendency term of the soil
c water diffusion equation.  also to compute ( prepare ) the matrix
c coefficients for the tri-diagonal matrix of the implicit time scheme.
c ----------------------------------------------------------------------
      integer nsold
      parameter(nsold = 20)
c
      integer cvfrz      
      integer ialp1
      integer iohinf
      integer j
      integer jj      
      integer k
      integer ks
      integer nsoil
c
      real acrt
      real ai(nsold)
      real bexp
      real bi(nsold)
      real ci(nsold)
      real dd
      real ddt
      real ddz
      real ddz2
      real denom
      real denom2
      real dice
      real dksat
      real dmax(nsold)
      real dsmdz
      real dsmdz2
      real dt
      real dt1
      real dwsat
      real edir
      real et(nsoil)
      real fcr
      real frzx
      real infmax
      real kdt
      real mxsmc
      real mxsmc2
      real numer
      real pcpdrp
      real pddum
      real px
      real rhstt(nsoil)
      real runoff1
      real runoff2
      real sh2o(nsoil)
      real sh2oa(nsoil)
      real sice(nsoil)
      real sicemax
      real slope
      real slopx
      real smcav
      real smcmax
      real smcwlt
      real sstt
      real sum
      real val
      real wcnd
      real wcnd2
      real wdf
      real wdf2
      real zsoil(nsoil)
c ----------------------------------------------------------------------
c frozen ground version:
c reference frozen ground parameter, cvfrz, is a shape parameter of
c areal distribution function of soil ice content which equals 1/cv.
c cv is a coefficient of spatial variation of soil ice content.  based
c on field data cv depends on areal mean of frozen depth, and it close
c to constant = 0.6 if areal mean frozen depth is above 20 cm.  that is
c why parameter cvfrz = 3 (int{1/0.6*0.6}).
c current logic doesn't allow cvfrz be bigger than 3
c ----------------------------------------------------------------------
        parameter(cvfrz = 3)
c ----------------------------------------------------------------------
c determine rainfall infiltration rate and runoff.  include the
c infiltration formule from schaake and koren model.
c modified by q duan
c ----------------------------------------------------------------------
      iohinf=1
c ----------------------------------------------------------------------
c let sicemax be the greatest, if any, frozen water content within soil
c layers.
c ----------------------------------------------------------------------
      sicemax = 0.0
      do ks=1,nsoil
       if (sice(ks) .gt. sicemax) sicemax = sice(ks)
      end do
c ----------------------------------------------------------------------
c determine rainfall infiltration rate and runoff
c ----------------------------------------------------------------------
      pddum = pcpdrp
      runoff1 = 0.0
      if (pcpdrp .ne. 0.0) then
c ----------------------------------------------------------------------
c modified by q. duan, 5/16/94
c ----------------------------------------------------------------------
c        if (iohinf .eq. 1) then
        dt1 = dt/86400.
        smcav = smcmax - smcwlt
        dmax(1)=-zsoil(1)*smcav
c ----------------------------------------------------------------------
c frozen ground version:
c ----------------------------------------------------------------------
        dice = -zsoil(1) * sice(1)
c          
        dmax(1)=dmax(1)*(1.0 - (sh2oa(1)+sice(1)-smcwlt)/smcav)
        dd=dmax(1)
c
        do ks=2,nsoil
c ----------------------------------------------------------------------
c frozen ground version:
c ----------------------------------------------------------------------
          dice = dice + ( zsoil(ks-1) - zsoil(ks) ) * sice(ks)
          dmax(ks) = (zsoil(ks-1)-zsoil(ks))*smcav
          dmax(ks) = dmax(ks)*(1.0 - (sh2oa(ks)+sice(ks)-smcwlt)/smcav)
          dd = dd+dmax(ks)
        end do
c ----------------------------------------------------------------------
c val = (1.-exp(-kdt*sqrt(dt1)))
c in below, remove the sqrt in above
c ----------------------------------------------------------------------
        val = (1.-exp(-kdt*dt1))
        ddt = dd*val
        px = pcpdrp*dt  
        if (px .lt. 0.0) px = 0.0
        infmax = (px*(ddt/(px+ddt)))/dt
c ----------------------------------------------------------------------
c frozen ground version:
c reduction of infiltration based on frozen ground parameters
c ----------------------------------------------------------------------
        fcr = 1. 
        if (dice .gt. 1.e-2) then 
          acrt = cvfrz * frzx / dice 
          sum = 1.
          ialp1 = cvfrz - 1 
          do j = 1,ialp1
            k = 1
            do jj = j+1,ialp1
              k = k * jj
            end do
            sum = sum + (acrt ** ( cvfrz-j)) / float (k) 
          end do
          fcr = 1. - exp(-acrt) * sum 
        endif 
        infmax = infmax * fcr
c ----------------------------------------------------------------------
c correction of infiltration limitation:
c if infmax .le. hydrolic conductivity assign infmax the value of
c hydrolic conductivity
c ----------------------------------------------------------------------
c         mxsmc = max ( sh2oa(1), sh2oa(2) ) 
        mxsmc = sh2oa(1)
c
        call wdfcnd (wdf,wcnd,mxsmc,smcmax,bexp,dksat,dwsat,
     &               sicemax)
c
        infmax = max(infmax,wcnd)
        infmax = min(infmax,px)
c
        if (pcpdrp .gt. infmax) then
          runoff1 = pcpdrp - infmax
          pddum = infmax
        endif
      endif
c ----------------------------------------------------------------------
c to avoid spurious drainage behavior, 'upstream differencing' in line
c below replaced with new approach in 2nd line:
c 'mxsmc = max(sh2oa(1), sh2oa(2))'
c ----------------------------------------------------------------------
      mxsmc = sh2oa(1)
      call wdfcnd (wdf,wcnd,mxsmc,smcmax,bexp,dksat,dwsat,
     &             sicemax)
c ----------------------------------------------------------------------
c calc the matrix coefficients ai, bi, and ci for the top layer
c ----------------------------------------------------------------------
      ddz = 1. / ( -.5 * zsoil(2) )
      ai(1) = 0.0
      bi(1) = wdf * ddz / ( -zsoil(1) )
      ci(1) = -bi(1)
c ----------------------------------------------------------------------
c calc rhstt for the top layer after calc'ng the vertical soil moisture
c gradient btwn the top and next to top layers.
c ----------------------------------------------------------------------
      dsmdz = ( sh2o(1) - sh2o(2) ) / ( -.5 * zsoil(2) )
      rhstt(1) = (wdf * dsmdz + wcnd - pddum + edir + et(1))/zsoil(1)
      sstt = wdf * dsmdz + wcnd + edir + et(1)
c ----------------------------------------------------------------------
c initialize ddz2
c ----------------------------------------------------------------------
      ddz2 = 0.0
c ----------------------------------------------------------------------
c loop thru the remaining soil layers, repeating the abv process
c ----------------------------------------------------------------------
      do k = 2,nsoil
        denom2 = (zsoil(k-1) - zsoil(k))
        if (k .ne. nsoil) then
          slopx = 1.
c ----------------------------------------------------------------------
c again, to avoid spurious drainage behavior, 'upstream differencing' in
c line below replaced with new approach in 2nd line:
c 'mxsmc2 = max (sh2oa(k), sh2oa(k+1))'
c ----------------------------------------------------------------------
          mxsmc2 = sh2oa(k)
          call wdfcnd (wdf2,wcnd2,mxsmc2,smcmax,bexp,dksat,dwsat,
     &                 sicemax)
c ----------------------------------------------------------------------
c calc some partial products for later use in calc'ng rhstt
c ----------------------------------------------------------------------
          denom = (zsoil(k-1) - zsoil(k+1))
          dsmdz2 = (sh2o(k) - sh2o(k+1)) / (denom * 0.5)
c ----------------------------------------------------------------------
c calc the matrix coef, ci, after calc'ng its partial product
c ----------------------------------------------------------------------
          ddz2 = 2.0 / denom
          ci(k) = -wdf2 * ddz2 / denom2
        else
c ----------------------------------------------------------------------
c slope of bottom layer is introduced
c ----------------------------------------------------------------------
          slopx = slope
c ----------------------------------------------------------------------
c retrieve the soil water diffusivity and hydraulic conductivity for
c this layer
c ----------------------------------------------------------------------
          call wdfcnd (wdf2,wcnd2,sh2oa(nsoil),smcmax,bexp,dksat,dwsat,
     &                 sicemax)
c ----------------------------------------------------------------------
c calc a partial product for later use in calc'ng rhstt
c ----------------------------------------------------------------------
          dsmdz2 = 0.0
c ----------------------------------------------------------------------
c set matrix coef ci to zero
c ----------------------------------------------------------------------
          ci(k) = 0.0
        endif
c ----------------------------------------------------------------------
c calc rhstt for this layer after calc'ng its numerator
c ----------------------------------------------------------------------
        numer = (wdf2 * dsmdz2) + slopx * wcnd2 - (wdf * dsmdz) 
     &    - wcnd + et(k)
        rhstt(k) = numer / (-denom2)
c ----------------------------------------------------------------------
c calc matrix coefs, ai, and bi for this layer
c ----------------------------------------------------------------------
        ai(k) = -wdf * ddz / denom2
        bi(k) = -( ai(k) + ci(k) )
c ----------------------------------------------------------------------
c reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr
c runoff2:  sub-surface or baseflow runoff
c ----------------------------------------------------------------------
        if (k .eq. nsoil) then
          runoff2 = slopx * wcnd2
        endif
        if (k .ne. nsoil) then
          wdf = wdf2
          wcnd = wcnd2
          dsmdz = dsmdz2
          ddz = ddz2
        endif
      end do
c ----------------------------------------------------------------------
c end subroutine srt
c ----------------------------------------------------------------------
      return
      end
