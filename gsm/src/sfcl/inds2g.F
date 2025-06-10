#define LALBD lalbd
#define LSOIL lsoil
#define MSUB msub
#define NSOIL nsoil
      subroutine inds2g(is2g)
c
      implicit none
#include <paramodel.h>
#include <varsfc.h>
#include <vargrb.h>
c
      integer is2g(numsfcs),n
c
      do n=1,numsfcs
        is2g(n)=0
      enddo
c
      is2g(jtsf)=itsf  ! sst and land surface temperature
      is2g(jsno)=isno  ! snow depth
      is2g(jtg3)=itg3  ! deep soil temperature
      is2g(jzor)=izor  ! surface roughness
      is2g(jcv )=9999  ! convective cloud cover
      is2g(jcvb)=9999  ! convective cloud base
      is2g(jcvt)=9999  ! convective cloud top
      is2g(jsli)=iais  ! sea ice mask (to be converted to slimsk)
      is2g(jcpy)=9999  ! canopy water content
      is2g(jf10)=9999  ! 10m and lowest sigma level coversion factor
#ifndef OSULSM1
      is2g(jalf)=ialf  ! albedo fraction
      is2g(just)=9999  ! fricational velocity
      is2g(jffm)=9999  ! coefficient for momentun exchange
      is2g(jffh)=9999  ! coefficient for heat exchange
#endif

#ifdef OSULSM1
      is2g(jsmc)=iso2  ! soil moisture
      is2g(jstc)=ito2  ! soil temperature
      is2g(jalb)=iab1  ! albedo
      is2g(jplr)=iplr  ! stomata resistance
#endif

#ifdef OSULSM2
      is2g(jsmc)=iso2  ! soil moisture
      is2g(jstc)=ito2  ! soil temperature
      is2g(jalb)=iab4  ! albedo
      is2g(jveg)=iveg  ! vegetation cover
      is2g(jvet)=ivet  ! vegetation type
      is2g(jsot)=isot  ! soil type
#endif

#ifdef NOALSM1
      is2g(jsmc)=isn1  ! soil moisture
      is2g(jstc)=itn1  ! soil temperature
      is2g(jalb)=iab4  ! albedo
      is2g(jveg)=iveg  ! vegetation cover
      is2g(jvet)=ivet  ! vegetation type
      is2g(jsot)=isot  ! soil type
      is2g(jprc)=9999  ! precip for noah
      is2g(jsrf)=9999  ! sr flag for noah
      is2g(jsnd)=9999  ! snow depth for noah
      is2g(jslc)=9999  ! slc for noah
      is2g(jsmn)=ismn  ! min vegetation cover for noah
      is2g(jsmx)=ismx  ! max vegetation cover for noah
      is2g(jslo)=islo  ! slope type for noah
      is2g(jsna)=isna  ! snow albedo for noah
#endif

#ifdef VICLSM1
      is2g(jsmc)=isv1  ! soil moisture
      is2g(jstc)=itv1  ! soil temperature
      is2g(jalb)=iab4  ! albedo
      is2g(jveg)=ivgv  ! vegetation cover
      is2g(jvet)=ivtv  ! vegetation type
      is2g(jrot)=ivrt  ! vegetation root
      is2g(jprc)=9999  ! precip for noah
      is2g(jsrf)=9999  ! sr flag for noah
      is2g(jbif)=ibif  ! Variable infil curve parameter (N/A)
      is2g(jds )=ids   ! Fract of Dsm nonlinear baseflow begins
      is2g(jdsm)=idsm  ! Maximum velocity of baseflow (mm/day)
      is2g(jws )=iws   ! Fract maxi sm nonlinear baseflow occurs
      is2g(jcef)=icef  ! c
      is2g(jexp)=iexp  ! Para the vari of Ksat with sm
      is2g(jkst)=ikst  ! Saturated hydrologic conductivity (mm/day)
      is2g(jdph)=idph  ! soil layer thickness (m)
      is2g(jbub)=ibub  ! Bubbling pressure of soil layer (cm)
      is2g(jqrt)=iqrt  ! Quartz content of soil layer (fraction)
      is2g(jbkd)=ibkd  ! Bulk density of soil layer (kg/m3)
      is2g(jsld)=isld  ! Soil density of soil layer (kg/m3)
      is2g(jwcr)=iwcr  ! sm content at the critical point (mm)
      is2g(jwpw)=iwpw  ! sm content wilting point (mm)
      is2g(jsmr)=ismr  ! Soil moisture residual moisture (mm)
      is2g(jsmx)=9999  ! maximum soil moisture (mm)
      is2g(jdpn)=9999  ! maximum sm at soil node (m3/m3)
      is2g(jsxn)=9999  ! maximum sm at soil node (m3/m3)
      is2g(jepn)=9999  ! Para the vari of Ksat at soil node (N/A)
      is2g(jbbn)=9999  ! bubbling pressure at soil node (cm)
      is2g(japn)=9999  ! para alpha at soil node
      is2g(jbtn)=9999  ! para beta at soil node
      is2g(jgmn)=9999  ! para gamma at soil node
      is2g(jlai)=9999  ! lai
      is2g(jslz)=islz  ! soil roughness (m)
      is2g(jsnz)=isnz  ! snow roughness (m)
      is2g(jsic)=9999  ! soil ice content (mm)
      is2g(jcsn)=9999  ! canopy snow  (mm h2o)
      is2g(jrsn)=9999  ! snow density (kg/m^3)
      is2g(jtsn)=9999  ! snow surface temperature (K)
      is2g(jtpk)=9999  ! snow pack temperature (K)
      is2g(jsfw)=9999  ! surface snow water equivalent (mm h2o)
      is2g(jpkw)=9999  ! snow pack snow water equivalent (mm h2o)
      is2g(jlst)=9999  ! time step since last snow fall
#endif

#ifdef VICLSM2
      is2g(jsmc)=isv2  ! soil moisture
      is2g(jstc)=itv2  ! soil temperature
      is2g(jalb)=iab4  ! albedo
      is2g(jveg)=ivgv  ! vegetation cover
      is2g(jvet)=ivtv  ! vegetation type
      is2g(jrot)=ivrt  ! vegetation root
      is2g(jprc)=9999  ! precip for noah
      is2g(jsrf)=9999  ! sr flag for noah
      is2g(jbif)=ibif  ! Variable infil curve parameter (N/A)
      is2g(jds )=ids   ! Fract of Dsm nonlinear baseflow begins
      is2g(jdsm)=idsm  ! Maximum velocity of baseflow (mm/day)
      is2g(jws )=iws   ! Fract maxi sm nonlinear baseflow occurs
      is2g(jcef)=icef  ! c
      is2g(jexp)=iexp  ! Para the vari of Ksat with sm
      is2g(jkst)=ikst  ! Saturated hydrologic conductivity (mm/day)
      is2g(jdph)=idph  ! soil layer thickness (m)
      is2g(jbub)=ibub  ! Bubbling pressure of soil layer (cm)
      is2g(jqrt)=iqrt  ! Quartz content of soil layer (fraction)
      is2g(jbkd)=ibkd  ! Bulk density of soil layer (kg/m3)
      is2g(jsld)=isld  ! Soil density of soil layer (kg/m3)
      is2g(jwcr)=iwcr  ! sm content at the critical point (mm)
      is2g(jwpw)=iwpw  ! sm content wilting point (mm)
      is2g(jsmr)=ismr  ! Soil moisture residual moisture (mm)
      is2g(jsmx)=9999  ! maximum soil moisture (mm)
      is2g(jdpn)=9999  ! maximum sm at soil node (m3/m3)
      is2g(jsxn)=9999  ! maximum sm at soil node (m3/m3)
      is2g(jepn)=9999  ! Para the vari of Ksat at soil node (N/A)
      is2g(jbbn)=9999  ! bubbling pressure at soil node (cm)
      is2g(japn)=9999  ! para alpha at soil node
      is2g(jbtn)=9999  ! para beta at soil node
      is2g(jgmn)=9999  ! para gamma at soil node
      is2g(jlai)=9999  ! lai
      is2g(jslz)=islz  ! soil roughness (m)
      is2g(jsnz)=isnz  ! snow roughness (m)
      is2g(jsic)=9999  ! soil ice content (mm)
      is2g(jcsn)=9999  ! canopy snow  (mm h2o)
      is2g(jrsn)=9999  ! snow density (kg/m^3)
      is2g(jtsn)=9999  ! snow surface temperature (K)
      is2g(jtpk)=9999  ! snow pack temperature (K)
      is2g(jsfw)=9999  ! surface snow water equivalent (mm h2o)
      is2g(jpkw)=9999  ! snow pack snow water equivalent (mm h2o)
      is2g(jlst)=9999  ! time step since last snow fall
      is2g(jlai)=ilai  ! leaf area index
#endif
c
#ifdef ADJSOILW
      is2g(jobp)=iobp  ! observed precipitation
#endif
c
      return
      end
