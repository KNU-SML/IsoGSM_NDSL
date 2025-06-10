#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define NSOIL nsoil
c
      subroutine vic1ini(ijdim,im)
c
      implicit none
c
      integer ijdim,im
c
*-----------------------------------------------------------------------
*  SUBROUTINE DOCUMENTATION BLOCK
*  
*  SUBPROGRAM vic1ini.F
*  PRGMMR: Ji Chen   
*  ORG: ECPC/CRD/SIO/UCSD DATE: SEPTEMBER 2003
*
*  ABSTRACT: initialize the rest variables and parameters after merging
*            grib data, analysis data, and other land surface initial data
*
*  The following parameters and variables should be initialized after 
*  calling sfc
*    
*    jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,
*    jsli,jveg,jcpy,jf10,jvet,jrot,jalf,just,jffm,jffh,
*    jds ,jdsm,jws ,jcef,jexp,jkst,jbub,jqrt,jbkd,jsld,
*    jwcr,jwpw,jsmr,
*
*  We need initialize the following parameters and variables here:
*
*           Number        Comments
*  1. jsmx:  lsoil        maximum soil moisture (mm)
*  2. jdpn:  nsoil        thickness of soil node (m)
*  3. jsxn:  nsoil        maximum sm at soil node (m3/m3)
*  4. jepn:  nsoil        Para the vari of Ksat at soil node (N/A)
*  5. jbbn:  nsoil        bubbling pressure at soil node (cm)
*  6. japn:  nsoil        para alpha at soil node
*  7. jbtn:  nsoil        para beta at soil node
*  8. jgmn:  nsoil        para gamma at soil node
*  9. jlai:      1        leaf area index (monthly)
* 10. jsic:  lsoil        soil ice content (mm)
* 11. jcsn:      1        canopy snow  (m h2o)
* 12. jrsn:      1        snow density (kg/m^3)
* 13. jtsn:      1        snow surface temperature (K)
* 14. jtpk:      1        snow pack temperature (K)
* 15. jsfw:      1        surface snow water equivalent (m h2o)
* 16. jpkw:      1        snow pack snow water equivalent (m h2o)
* 17. jlst:      1        time step since last snow fall
*-----------------------------------------------------------------------
c
#include <paramodel.h>
c
#ifdef VICLSM1
c
#include <varsfc.h>
#include <comsfcmrg.h>
#include <vic_veglib.h>
c
      integer lsoil, nsoil
      parameter (lsoil=lsoil_, nsoil=nsoil_)

      integer ij, m, n

      real vdep(lsoil), vdzn(nsoil), sumd, sumdp, totdep
      real smx(lsoil), bub(lsoil), expt(lsoil)
      real smmaxnd(nsoil),expnd(nsoil),bubnd(nsoil)
      real alpha(nsoil),beta(nsoil),gamma(nsoil)
c
      real tmp_T, tmp_smc, smmax, buble, exple, max_unfwat
c
      real t0c
      parameter (t0c = 273.15)
c
*------------------- end of variable definition ------------------------
c
**********************************************************************
* in VICLSM1, some parameters are constant over the study regions,   *
* so we just give the fixed values of them here. In the next version *
* VICLSM, we will read their grib files                              *
**********************************************************************
c
      do ij = 1, ijdim
         if(sfcfcs(ij,jsli).eq.1.) then

            sfcfcs(ij,jlai) = veg_lai(im,sfcfcs(ij,jvet))
            sfcfcs(ij,jcsn) =   0.0
            sfcfcs(ij,jrsn) =  50.0
            if(sfcfcs(ij,jsno).gt.0.) then
               sfcfcs(ij,jlst) =  0.0
               sfcfcs(ij,jstc) = sfcfcs(ij,jstc+1)  ! 6/2004
            else
               sfcfcs(ij,jlst) = -1.0
            endif
            sfcfcs(ij,jtsn) = min(sfcfcs(ij,jtsf),t0c)
            sfcfcs(ij,jtpk) = min(sfcfcs(ij,jtsf),t0c)
            sfcfcs(ij,jsfw) = 0.0
            sfcfcs(ij,jpkw) = 0.0
            
*--- the Wcr and Wpwp are the fraction of maximum moisture 05/25/2004

            do m = 1, lsoil
               if(sfcfcs(ij,jsld+m-1).gt.0) then
                  sfcfcs(ij,jsmx+m-1) = (sfcfcs(ij,jsld+m-1)-
     &                 sfcfcs(ij,jbkd+m-1))/sfcfcs(ij,jsld+m-1) ! fraction
               else
                  sfcfcs(ij,jsmx+m-1) = 0.0
               end if
               if(sfcfcs(ij,jsmc+m-1).gt.sfcfcs(ij,jsmx+m-1)) then
                  sfcfcs(ij,jsmc+m-1) = sfcfcs(ij,jsmx+m-1)
               endif
               sfcfcs(ij,jwcr+m-1)=sfcfcs(ij,jwcr+m-1)*
     &                             sfcfcs(ij,jsmx+m-1)
               sfcfcs(ij,jwpw+m-1)=sfcfcs(ij,jwpw+m-1)*
     &                             sfcfcs(ij,jsmx+m-1)
               sfcfcs(ij,jsmx+m-1)=sfcfcs(ij,jsmx+m-1)*
     &                 sfcfcs(ij,jdph+m-1)*1000.0 ! maximum sm (mm)
            end do

*---check soil layer thickness effects on the VIC performancs 5/22/2004
c            do m = 1, lsoil_
c               if(sfcfcs(ij,jdph+m-1).gt.0)then
c                  sfcfcs(ij,jwcr+m-1) = sfcfcs(ij,jwcr+m-1) / 
c     &                                  sfcfcs(ij,jdph+m-1)
c                  sfcfcs(ij,jwpw+m-1) = sfcfcs(ij,jwpw+m-1) / 
c     &                                  sfcfcs(ij,jdph+m-1)
c                  sfcfcs(ij,jsmr+m-1) = sfcfcs(ij,jsmr+m-1) / 
c     &                                  sfcfcs(ij,jdph+m-1)
c                  sfcfcs(ij,jsmx+m-1) = sfcfcs(ij,jsmx+m-1) / 
c     &                                  sfcfcs(ij,jdph+m-1)
c               endif
c            enddo
c            sfcfcs(ij,jdph  ) = 0.1
c            sfcfcs(ij,jdph+1) = 0.3
c            sfcfcs(ij,jdph+2) = 1.0
c            do m = 1, lsoil_
c               sfcfcs(ij,jwcr+m-1) = sfcfcs(ij,jwcr+m-1) *
c     &                                  sfcfcs(ij,jdph+m-1)
c               sfcfcs(ij,jwpw+m-1) = sfcfcs(ij,jwpw+m-1) *
c     &                                  sfcfcs(ij,jdph+m-1)
c               sfcfcs(ij,jsmr+m-1) = sfcfcs(ij,jsmr+m-1) *
c     &                                  sfcfcs(ij,jdph+m-1)
c               sfcfcs(ij,jsmx+m-1) = sfcfcs(ij,jsmx+m-1) *
c     &                                  sfcfcs(ij,jdph+m-1)
c            enddo
*---end check

            do m = 1, lsoil
               tmp_T= sfcfcs(ij,jstc+m)
               smmax= sfcfcs(ij,jsmx+m-1)
               buble= sfcfcs(ij,jbub+m-1)
               exple = sfcfcs(ij,jexp+m-1)
               if(tmp_T.lt.t0c)then
                  sfcfcs(ij,jsic+m-1)=sfcfcs(ij,jsmc+m-1)*
     &                           sfcfcs(ij,jdph+m-1)*1000.0-
     &                           max_unfwat(tmp_T,smmax,buble,exple)
                  if(sfcfcs(ij,jsic+m-1).lt.0.0)sfcfcs(ij,jsic+m-1)=0.0
               else
                  sfcfcs(ij,jsic+m-1) = 0.0
               endif
               tmp_smc=sfcfcs(ij,jsmc+m-1)*sfcfcs(ij,jdph+m-1)*1000.0
               if(sfcfcs(ij,jsic+m-1).gt.tmp_smc) then
                  sfcfcs(ij,jsic+m-1) = tmp_smc
               endif
            enddo
         endif
      end do
      
c  define the depth of soil nodes, the number of soil nodes >= 3
c  currently, we are using a very simple way to find the soil temperture
c  profile for each soil node (Ji 2003)

      print *,'in vic1ini for computing node parameters'

      do ij = 1, ijdim
         if(sfcfcs(ij,jsli).eq.1.) then
            do m = 1, lsoil
               vdep(m) = sfcfcs(ij,jdph+m-1)
            end do
            if(vdep(1) .gt.0.0) then
               vdzn(1) = 0.0
               vdzn(2) = vdep(1)
               sumd = 0
               do n = 3, nsoil
                  vdzn(n) = (n-1)*0.1
                  sumd = sumd + vdzn(n)
               end do
               
               sumdp = 0
               do m = 2, lsoil
                  sumdp = sumdp + vdep(m)
               end do
               
               do n = 3, nsoil
                  if(sumd.gt.0) then
                     vdzn(n) = vdzn(n)/sumd*sumdp
                  else
                     vdzn(n) = 0
                  end if
               end do
            else
               do n = 1, nsoil
                  vdzn(n) = 0.0
               end do
            end if
            do n = 1, nsoil
               sfcfcs(ij,jdpn+n-1) = vdzn(n)
            end do
         endif
      end do

*-- compute thermal node soil parameters ---
c
      do ij = 1, ijdim
         if(sfcfcs(ij,jsli).eq.1.) then
            totdep = 0
            do m = 1, lsoil
               vdep(m) =  sfcfcs(ij,jdph+m-1)
               totdep  = totdep + vdep(m)
            end do
            do n = 1, nsoil
               vdzn(n) = sfcfcs(ij,jdpn+n-1)
            end do
c     
            if(totdep.gt.0.0) then
               do m = 1, lsoil
                  smx(m)  = sfcfcs(ij,jsmx+m-1)
                  bub(m)  = sfcfcs(ij,jbub+m-1)
                  expt(m) = sfcfcs(ij,jexp+m-1)
               end do
               
               call nodepara(lsoil, nsoil,  vdep,    vdzn,
     1                         smx,   bub,  expt, smmaxnd,
     2                       expnd, bubnd, alpha,    beta,
     3                       gamma)

               do n = 1, nsoil
                  sfcfcs(ij,jsxn+n-1) = smmaxnd(n)
                  sfcfcs(ij,jepn+n-1) = expnd(n)
                  sfcfcs(ij,jbbn+n-1) = bubnd(n)
                  sfcfcs(ij,japn+n-1) = alpha(n)
                  sfcfcs(ij,jbtn+n-1) = beta(n)
                  sfcfcs(ij,jgmn+n-1) = gamma(n)
               end do
            else
               do n = 1, nsoil
                  sfcfcs(ij,jsxn+n-1) = 0.0
                  sfcfcs(ij,jepn+n-1) = 0.0
                  sfcfcs(ij,jbbn+n-1) = 0.0
                  sfcfcs(ij,japn+n-1) = 0.0
                  sfcfcs(ij,jbtn+n-1) = 0.0
                  sfcfcs(ij,jgmn+n-1) = 0.0
               end do
            end if
c
         endif
      end do
#endif
      return
      end
