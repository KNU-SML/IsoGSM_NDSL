      subroutine runoff(msl,  dtime,    ppt,  binf,
     &                   ws,     ds,    dsm,   cef,
     &                  smx,    smr,    kst,  expt,
     &                  dph,    bub,
     &                 fsmc,   fsic,  sflow, bflow)
#include <paramodel.h>

#include <vartyp.h>

************************************************************************
* source file:       runoff.f
* prgmmr:            ji chen   
* org:               ecpc/crd/sio/ucsd 
* date:              june, july, august & september 2003
* prgm history:      modified from 4.0.3 uw vic (runoff.c)
*
* subroutine documentation block
*
* abstract: this subroutine calculates infiltration and runoff from 
*            the surface, gravity driven drainage between all soil 
*            layers, and generates bflow from the bottom layer.
*
* sublayer indecies are always is the current vic model moisture layer
*    1 = thawed sublayer, 2 = frozen sublayer, and 3 = unfrozen sublayer.
************************************************************************

* ------------------- input variables ----------------------------------
      integer msl          ! number of soil layer
      real dtime           ! time step (second)
      real ppt             ! incoming precip and snow melt (m/time step)
      real binf            ! vic parameter
      real ws              ! fraction of max soilm nonlinear baseflow
      real ds              ! fraction of dsm begin nonlinar bflow
      real dsm             ! maximum velocity of baseflow (mm/day)
      real cef             ! exponent used in infiltration curve (n/a)
      real fsic(msl)       ! soil layer ice (mm)
      real kst(msl)        ! saturated hydraulic conductivity (mm/day)
      real expt(msl)       ! parameter describing kst with soil moisture
      real dph(msl)        ! soil layer thickness (m)
      real bub(msl)        ! bubbling pressure of soil layer (cm)
      real smx(msl)        ! soil layer maximum soil moisture (mm)
      real smr(msl)        ! soil layer residue (mm)
* ----------------------------------------------------------------------

* ------------------- modified variables -------------------------------
      real fsmc(msl)       ! soil layer moisture (mm)
* ----------------------------------------------------------------------

* ------------------- output variables ---------------------------------
      real sflow           ! surface runoff (mm/time step)
      real bflow           ! base flow (mm/time step)
* ----------------------------------------------------------------------

* ---------------------- local variables -------------------------------
      integer k, tmp

      real frcday          ! fractional day of the model time step
      real inflow, max_infil, ex, a, i_0, basis
      real moist(msl), q12(msl-1)
      real top_moist, tmp_moist, tmp_smcn
      real top_smx, tmp_inflow, tmpws, frac

      real avgmatri, avgsmr, avgsmx, avgbub, avgb, avgkst, tmpdph
      real matri(msl), dif(msl), b(msl)
      real smcon(msl), smrcon(msl), smxcon(msl)

      real d1, d2, d3, d4

      real secpday
      parameter (secpday=86400.0)
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* ppt = amount of liquid water coming to the surface
* ----------------------------------------------------------------------
      inflow = ppt*1000.0

      frcday = dtime / secpday
* ----------------------------------------------------------------------
* initialize variables
* ----------------------------------------------------------------------

      do k = 1, msl
* ----------------------------------------------------------------------
* set layer unfrozen moisture content
* ----------------------------------------------------------------------
         moist(k) = fsmc(k) - fsic(k)
         if(moist(k).lt.0) moist(k) = 0.0
      end do
      
* ----------------------------------------------------------------------
* surf_runoff based on soil moisture level of upper layers
* ----------------------------------------------------------------------

      top_moist = 0.
      top_smx   = 0.

      do k = 1, msl-1
         top_moist = top_moist + fsmc(k)
         top_smx   = top_smx + smx(k)
      end do

      if(top_moist.gt.top_smx) top_moist=top_smx
      
* ----------------------------------------------------------------------
* calculate surf_runoff from surface--surf_runoff calculations 
*   for top layer only
*   a and i_0 as in wood et al. in jgr 97, d3, 1992 equation (1)
* ----------------------------------------------------------------------
      if(binf.le.0) then
         print *,'error runoff vic binf must be a positive value '
         print *,' binf = ',binf
#ifdef MP
#ifdef RSM
         call rmpabort
#else
         call mpabort
#endif
#else
         call abort
#endif
      endif

      max_infil = (1.0 + binf) * top_smx
      ex        = binf / (1.0 + binf)

      if(abs(top_moist-top_smx).lt.1.e-6) then
         a   = 1.0
         i_0 = max_infil
      else
         a   = 1.0 - (1.0 - top_moist / top_smx)**ex
         i_0 = max_infil*(1.0-(1.0 - a)**(1.0/binf))
      endif

* ----------------------------------------------------------------------
* maximum inflow ---- equation (3a) wood et al.
* ----------------------------------------------------------------------

      if(inflow.le.0.0) sflow = 0.0

      if(inflow.gt.0) then
         if (max_infil.eq.0.0) then
            sflow = inflow
         else if((i_0 + inflow) .ge. max_infil) then
            sflow = inflow - top_smx + top_moist
         else 
            basis = 1.0 - (i_0 + inflow) / max_infil
            sflow = inflow - top_smx + top_moist + top_smx
     &                   *(basis**(1.0*(1.0 + binf)))
         endif
      endif

      if(sflow.lt.0.) sflow = 0.0
      if(sflow.gt.inflow) sflow = inflow

      inflow = inflow - sflow

#ifdef DBGVIC
      print *,'in runoff before computing bflow, sflow=',sflow
#endif

* ----------------------------------------------------------------------
* compute flow between soil layers 
* ----------------------------------------------------------------------
#ifdef RUNSMDIF
      do k = 1, msl
         if(expt(k).gt.3.0) then
            b(k) = (expt(k)-3.0)/2.0
         else
            b(k) = 1.e-6
         endif
         if(moist(k).gt.smr(k))then
            matri(k)=10*bub(k)*((moist(k)-smr(k))/(smx(k)-smr(k)))**
     &                  (-b(k))
         else
            matri(k)=1.e6
         endif
         smcon(k) =max(1.e-6, moist(k)/1000.0/dph(k))   ! mm -> unitless
         smrcon(k)=max(1.e-6, smr(k)/1000.0/dph(k))
         smxcon(k)=smx(k)/1000.0/dph(k)
      enddo
#endif

* ----------------------------------------------------------------------
* compute drainage between sublayers 
* ----------------------------------------------------------------------

      do k = 1, msl-1

* ----------------------------------------------------------------------
* brooks & corey relation for hydraulic conductivity
* ----------------------------------------------------------------------

#ifdef RUNSMDIF
         if(moist(k).gt.smr(k))then
            tmpdph = dph(k)+dph(k+1)
            avgmatri=10**((dph(k+1)*log10(matri(k))+dph(k)*
     &                    log10(matri(k+1)))/tmpdph)
            avgsmr=10**((dph(k+1)*log10(smrcon(k))+dph(k)*
     &                    log10(smrcon(k+1)))/tmpdph)
            avgsmx=10**((dph(k+1)*log10(smxcon(k))+dph(k)*
     &                    log10(smxcon(k+1)))/tmpdph)
            avgbub=10**((dph(k+1)*log10(bub(k))+dph(k)*
     &                    log10(bub(k+1)))/tmpdph)
            avgb  =10**((dph(k+1)*log10(b(k))+dph(k)*
     &                    log10(b(k+1)))/tmpdph)
            avgkst=10**((dph(k+1)*log10(kst(k))+dph(k)*
     &                    log10(kst(k+1)))/tmpdph)

            tmp_smcn = avgsmr + (avgsmx-avgsmr)*
     &                  ((avgmatri/avgbub)**(-1.0/avgb))

            q12(k)= kst(k) * ((moist(k)-smr(k))/  
     &             (smx(k)-smr(k)))**expt(k)

            if(abs(smcon(k)-smcon(k+1)).lt.1.e-6) then
               dif(k) = 0.0
            else
               dif(k) = (avgb*avgkst*avgbub/avgsmx)*(((tmp_smcn-
     &                   avgsmr)/(avgsmx-avgsmr))**(avgb+2))
               dif(k) = dif(k)*(smcon(k)-smcon(k+1))/(tmpdph/2.0)
            endif
         else
            q12(k) = 0.0
            dif(k) = 0.0
         end if
#else
         dif(k) = 0.0
         if(moist(k).gt.smr(k))then
            q12(k) = kst(k)* ((moist(k)-smr(k))/  
     &           (smx(k)-smr(k)))**expt(k)
         else
            q12(k) = 0.0
         endif

c         if(fsic(k).lt.1.e-6 .and. q12(k).lt.1.e-6)then
c         if(expt(k).gt.13)then
c            print 101,k,q12(k),kst(k),moist(k)/dph(k)/1000.0,
c     &               fsic(k)/dph(k)/1000.0,smr(k)/dph(k)/1000.0,
c     &               smx(k)/dph(k)/1000.0,expt(k),dph(k)
c 101        format('k=',i2,8(1x,e10.4))
c         endif
#endif
      end do

* -- q12(k) (mm/day) -> (mm/time step)
      
      do k = 1, msl-1
         q12(k) = q12(k) * frcday
         dif(k) = dif(k) * frcday
      end do

* ----------------------------------------------------------------------
* solve for current soil layer moisture, and check versus maximum 
* and minimum moisture contents.  
* ----------------------------------------------------------------------

      do k = 1, msl-1
         tmp_inflow = 0.

* ----------------------------------------------------------------------
* update soil layer moisture content
* ----------------------------------------------------------------------

         moist(k) = moist(k) + inflow - q12(k) - dif(k)

* ----------------------------------------------------------------------
* verify that soil layer moisture is less than maximum
* ----------------------------------------------------------------------

         if((moist(k)+fsic(k)).gt.smx(k))then
            tmp_inflow=moist(k)+fsic(k)-smx(k)
            moist(k) = smx(k)-fsic(k)
            
            if(k.eq.1) then
               q12(k) = q12(k) + tmp_inflow
               tmp_inflow = 0
            else
               tmp = k
               do while (tmp_inflow.gt.0)
                  tmp=tmp-1
                  if(tmp.eq.0) then
* ----------------------------------------------------------------------
* if top layer saturated, add to surf_runoff
* ----------------------------------------------------------------------
                     sflow = sflow + tmp_inflow
                     tmp_inflow = 0.0
                  else
* ----------------------------------------------------------------------
* else add excess soil moisture to next higher layer
* ----------------------------------------------------------------------
                     moist(tmp) = moist(tmp)+tmp_inflow
                     if((moist(tmp)+fsic(tmp)).gt.smx(tmp))then
                        tmp_inflow=(moist(tmp)+fsic(tmp))-smx(tmp)
                        moist(tmp)=smx(tmp)-fsic(tmp)
                     else
                        tmp_inflow= 0
                     endif
                  endif
               enddo    ! end check if excess moisture in top layer
            endif
         endif

* ----------------------------------------------------------------------
* verify that current layer moisture is greater than minimum
* ----------------------------------------------------------------------

         if ((moist(k)+fsic(k)).lt.smr(k))then

* ----------------------------------------------------------------------
* moisture cannot fall below residual moisture content
* ----------------------------------------------------------------------
            q12(k)  = q12(k)+moist(k)+fsic(k)-smr(k)
            moist(k)=smr(k)-fsic(k)
            if(moist(k).lt.0) moist(k) = 0.0
         end if
         inflow = q12(k)+dif(k)
      end do                    ! end loop through (msl-1) soil layers 

* ----------------------------------------------------------------------
* compute bflow
*    arno model for the bottom soil layer (based on bottom
*    soil layer moisture from previous time step)
* ----------------------------------------------------------------------
      
      k = msl

#ifdef VICGLB
*
* when using vic global soil data, please check the following paper for
* the detail of computing baseflow algorithm.
* nijssen b., g.m. odonnell, d.p. lettenmaier, d. lohmann, and e.f. wood,
* predicting the discharge of global rivers. j. climate, 3307-3322,2001.
*
* d1=ds, d2=dsm, d3=ws, d4=cef
*
      if(moist(k).ge.smr(k))then
         d1 = ds
         d2 = dsm
         d3 = ws
         d4 = cef

         bflow = d1*moist(k)
         if(moist(k).gt.d3) then
            bflow = bflow + d2*(moist(k)-d3)**d4
         endif
         bflow = bflow * frcday
      else
         bflow = 0.0
      endif
#else
      if(ws.gt.1.0.or.ws.lt.0.0) then
         print *,'error in runoff_diffu wrong ws: ',ws
#ifdef MP
#ifdef RSM
         call rmpabort
#else
         call mpabort
#endif
#else
         call abort
#endif
      endif
      if(moist(k).ge.smr(k))then
         tmpws = ws*(smx(k)-smr(k))
         frac  = frcday * ds * dsm /tmpws
         bflow = frac * (moist(k)-smr(k))             ! mm/time step

         if (moist(k).gt.tmpws) then
            frac = (moist(k)-tmpws)/(smx(k)-tmpws)
            bflow= bflow+frcday*dsm*(1.0-ds/ws)*(frac**cef) ! mm/time step
         end if
      else
         bflow = 0.0
      endif
#endif

* turn off bflow
*      bflow = 0
      
* ----------------------------------------------------------------------
* extract bflow from the bottom soil layer 
* ----------------------------------------------------------------------
      moist(k) = moist(k) + inflow - bflow
c      print 100,inflow,bflow,moist(k)
c 100  format('runof ',3f12.7)
* ----------------------------------------------------------------------
* check lower sub-layer moistures 
* ----------------------------------------------------------------------
      tmp_moist = 0
	
      if((moist(k)+fsic(k)).lt.smr(k)) then

* ----------------------------------------------------------------------
* soil moisture is below minimum 
* ----------------------------------------------------------------------
         bflow = bflow+moist(k)+fsic(k)-smr(k)
         moist(k) = smr(k)-fsic(k)
      end if
      if(bflow .lt. 0) bflow = 0

      if((moist(k)+fsic(k)) .gt. smx(k)) then

* ----------------------------------------------------------------------
* soil moisture above maximum 
* ----------------------------------------------------------------------
         tmp_moist= moist(k)+fsic(k)-smx(k)
         moist(k) = smx(k) - fsic(k)
         tmp = k
         do while (tmp_moist .gt. 0)
            tmp = tmp - 1
            if(tmp.eq.0) then
* ----------------------------------------------------------------------
* if top layer saturated, add to surf_runoff
* ----------------------------------------------------------------------
               sflow = sflow + tmp_moist
               tmp_moist = 0
            else 
* ----------------------------------------------------------------------
* else if sublayer exists, add excess soil moisture
* ----------------------------------------------------------------------
               moist(tmp) = moist(tmp) + tmp_moist 
               if((moist(tmp)+fsic(tmp)).gt.smx(tmp))then
                  tmp_moist = moist(tmp) + fsic(tmp) - smx(tmp)
                  moist(tmp)= smx(tmp)-fsic(tmp)
               else
                  tmp_moist=0
               end if
            endif
         enddo
      end if

      do k = 1, msl
         fsmc(k) = moist(k)+fsic(k)
c         fsmc(k) = fsmc(k)
      end do

      return
      end 
