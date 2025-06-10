#define LSOIL lsoil
#define NSOIL nsoil
#define MSUB  msub
#define LALBD lalbd
c
c  varsfc.h
c
c  define sfc file variables and their properties.
c
c  this subroutine should list sfc file variables for
c  each land surface program
c  
c  if you are introducing new land surface model, you need to 
c  introduce new numsfcs and numsfcv and add subsequent 3 parameters 
c  for each new field.
c
c  special note:
c   it is not necessary to match the length of the comsfc array
c   to numsfcs.  In this case, make sure that the array pointed by
c   the index specified here is not referenced anywhere in the program.
c   you may still reference it in the grbfld array in the sfcl programs.
c   this mechanism is used for implementation of soil moisture
c   adjustment using observed precipitation
c
      integer numsfcs,numsfcv
c
c  note that numsfcv is one less than actual number of variables
c  for all cases except OSULSM1.  this is because frictional velocity
c  array is treated as just a single field.
c
c
#ifdef OSULSM1
c  osulsm1
c
#ifdef ADJSOILW
      parameter(numsfcs=11+lsoil_*2+lalbd_+1)
      parameter(numsfcv=15)
#else
      parameter(numsfcs=11+lsoil_*2+lalbd_)
      parameter(numsfcv=14)
#endif
      integer jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,jsli,
     1        jplr,jcpy,jf10
#ifdef ADJSOILW
      integer jobp
#endif

      parameter(jtsf=     1     )  ! sst and land surface temperature
      parameter(jsmc=jtsf+1     )  ! soil moisture
      parameter(jsno=jsmc+lsoil_)  ! snow depth
      parameter(jstc=jsno+1     )  ! soil temperature
      parameter(jtg3=jstc+lsoil_)  ! deep soil temperature
      parameter(jzor=jtg3+1     )  ! surface roughness
      parameter(jcv =jzor+1     )  ! convective cloud cover
      parameter(jcvb=jcv +1      ) ! convective cloud base
      parameter(jcvt=jcvb+1     )  ! convective cloud top
      parameter(jalb=jcvt+1     )  ! albedo
      parameter(jsli=jalb+lalbd_)  ! land/sea ice/snow mask
      parameter(jplr=jsli+1     )  ! stomato resistance
      parameter(jcpy=jplr+1     )  ! canopy
      parameter(jf10=jcpy+1     )  ! 10m and sig1 conv factor
#ifdef ADJSOILW
      parameter(jobp=jf10+1     )  ! observed precip
#endif
#endif
c
#ifdef OSULSM2
c  osulsm2
c
#ifdef ADJSOILW
      parameter(numsfcs=16+lsoil_*2+lalbd_+2+1)
      parameter(numsfcv=21)
#else
      parameter(numsfcs=16+lsoil_*2+lalbd_+2)
      parameter(numsfcv=20)
#endif
      integer jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,jsli,
     1        jveg,jcpy,jf10,jvet,jsot,jalf,just,jffm,jffh
#ifdef ADJSOILW
      integer jobp
#endif
      parameter(jtsf=     1     )  ! sst and land surface temperature
      parameter(jsmc=jtsf+1     )  ! soil moisture
      parameter(jsno=jsmc+lsoil_)  ! snow depth
      parameter(jstc=jsno+1     )  ! soil temperature
      parameter(jtg3=jstc+lsoil_)  ! deep soil temperature
      parameter(jzor=jtg3+1     )  ! surface roughness
      parameter(jcv =jzor+1     )  ! convective cloud cover
      parameter(jcvb=jcv +1     )  ! convective cloud base
      parameter(jcvt=jcvb+1     )  ! convective cloud top
      parameter(jalb=jcvt+1     )  ! albedo
      parameter(jsli=jalb+lalbd_)  ! land/sea ice/snow mask
      parameter(jveg=jsli+1     )  ! vegetation cover
      parameter(jcpy=jveg+1     )  ! canopy water content
      parameter(jf10=jcpy+1     )  ! 10m factor
      parameter(jvet=jf10+1     )  ! vegetation type
      parameter(jsot=jvet+1     )  ! soil type
      parameter(jalf=jsot+1     )  ! albedo fraction
      parameter(just=jalf+2     )  ! frictional velocity
      parameter(jffm=just+1     )  ! coefficient for momentun exchange
      parameter(jffh=jffm+1     )  ! coefficient for heat exchange
#ifdef ADJSOILW
      parameter(jobp=jffh+1     )  ! observed precipitation
#endif
#endif
c
#ifdef NOALSM1
c  noalsm1
c
#ifdef ADJSOILW
      parameter(numsfcs=23+lsoil_*3+lalbd_+2+1)
      parameter(numsfcv=29)
#else
      parameter(numsfcs=23+lsoil_*3+lalbd_+2)
      parameter(numsfcv=28)
#endif
      integer jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,jsli,
     1        jveg,jcpy,jf10,jvet,jsot,jalf,just,jffm,jffh,jprc,jsrf,
     2        jsnd,jslc,jsmn,jsmx,jslo,jsna
#ifdef ADJSOILW
      integer jobp
#endif
      parameter(jtsf=     1     )  ! sst and land surface temperature
      parameter(jsmc=jtsf+1     )  ! soil moisture
      parameter(jsno=jsmc+lsoil_)  ! snow depth
      parameter(jstc=jsno+1     )  ! soil temperature
      parameter(jtg3=jstc+lsoil_)  ! deep soil temperature
      parameter(jzor=jtg3+1     )  ! surface roughness
      parameter(jcv =jzor+1     )  ! convective cloud cover
      parameter(jcvb=jcv +1     )  ! convective cloud base
      parameter(jcvt=jcvb+1     )  ! convective cloud top
      parameter(jalb=jcvt+1     )  ! albedo
      parameter(jsli=jalb+lalbd_)  ! land/sea ice/snow mask
      parameter(jveg=jsli+1     )  ! vegetation cover
      parameter(jcpy=jveg+1     )  ! canopy water content
      parameter(jf10=jcpy+1     )  ! 10m and sig1 conv factor
      parameter(jvet=jf10+1     )  ! vegetation type
      parameter(jsot=jvet+1     )  ! soil type
      parameter(jalf=jsot+1     )  ! albedo fraction
      parameter(just=jalf+2     )  ! frictional velocity
      parameter(jffm=just+1     )  ! coefficient for momentun exchange
      parameter(jffh=jffm+1     )  ! coefficient for heat exchange
      parameter(jprc=jffh+1     )  ! precip for noah
      parameter(jsrf=jprc+1     )  ! sr flag for noah
      parameter(jsnd=jsrf+1     )  ! snow depth for noah
      parameter(jslc=jsnd+1     )  ! slc for noah
      parameter(jsmn=jslc+lsoil_)  ! min vegetation cover for noah
      parameter(jsmx=jsmn+1     )  ! max vegetation cover for noah
      parameter(jslo=jsmx+1     )  ! slope type for noah
      parameter(jsna=jslo+1     )  ! snow albedo for noah
#ifdef ADJSOILW
      parameter(jobp=jsna+1     )  ! observed precipitation
#endif
#endif
c
#ifdef VICLSM1
c  viclsm1
c
#ifdef ADJSOILW
      parameter(numsfcs=32+lsoil_*14+nsoil_*8+lalbd_+2+1)
      parameter(numsfcv=57)
#else
      parameter(numsfcs=32+lsoil_*14+nsoil_*8+lalbd_+2)
      parameter(numsfcv=56)
#endif
      integer jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,
     2        jsli,jveg,jcpy,jf10,jvet,jrot,jalf,just,jffm,jffh,
     3        jprc,jsrf,jbif,jds ,jdsm,jws ,jcef,jexp,jkst,jdph,
     4        jbub,jqrt,jbkd,jsld,jwcr,jwpw,jsmr,jsmx,jdpn,jsxn,
     5        jepn,jbbn,japn,jbtn,jgmn,jlai,jslz,jsnz,jsic,jcsn,
     6        jrsn,jtsn,jtpk,jsfw,jpkw,jlst
#ifdef ADJSOILW
      integer jobp
#endif
      parameter(jtsf=     1     )  ! sst and land surface temperature
      parameter(jsmc=jtsf+1     )  ! soil moisture
      parameter(jsno=jsmc+lsoil_)  ! snow depth
      parameter(jstc=jsno+1     )  ! soil temperature
      parameter(jtg3=jstc+nsoil_)  ! deep soil temperature
      parameter(jzor=jtg3+1     )  ! surface roughness
      parameter(jcv =jzor+1     )  ! convective cloud cover
      parameter(jcvb=jcv +1     )  ! convective cloud base
      parameter(jcvt=jcvb+1     )  ! convective cloud top
      parameter(jalb=jcvt+1     )  ! albedo
      parameter(jsli=jalb+lalbd_)  ! land/sea ice/snow mask
      parameter(jveg=jsli+1     )  ! vegetation cover
      parameter(jcpy=jveg+1     )  ! canopy water content
      parameter(jf10=jcpy+1     )  ! 10m and sig1 conv factor
      parameter(jvet=jf10+1     )  ! vegetation type
      parameter(jrot=jvet+1     )  ! root content
      parameter(jalf=jrot+lsoil_)  ! albedo fraction
      parameter(just=jalf+2     )  ! frictional velocity
      parameter(jffm=just+1     )  ! coefficient for momentun exchange
      parameter(jffh=jffm+1     )  ! coefficient for heat exchange
      parameter(jprc=jffh+1     )  ! precip for vic
      parameter(jsrf=jprc+1     )  ! sr flag for vic
      parameter(jbif=jsrf+1     )  ! Variable infil curve parameter (N/A)
      parameter(jds =jbif+1     )  ! Fract of Dsm nonlinear baseflow begins
      parameter(jdsm=jds +1     )  ! Maximum velocity of baseflow (mm/day)
      parameter(jws =jdsm+1     )  ! Fract maxi sm nonlinear baseflow occurs
      parameter(jcef=jws +1     )  ! c
      parameter(jexp=jcef+1     )  ! Para the vari of Ksat with sm
      parameter(jkst=jexp+lsoil_)  ! Saturated hydrologic conductivity (mm/day)
      parameter(jdph=jkst+lsoil_)  ! thickness of soil layer (m)
      parameter(jbub=jdph+lsoil_)  ! Bubbling pressure of soil layer (cm)
      parameter(jqrt=jbub+lsoil_)  ! Quartz content of soil layer (fraction)
      parameter(jbkd=jqrt+lsoil_)  ! Bulk density of soil layer (kg/m3)
      parameter(jsld=jbkd+lsoil_)  ! Soil density of soil layer (kg/m3)
      parameter(jwcr=jsld+lsoil_)  ! sm content at the critical point (mm)
      parameter(jwpw=jwcr+lsoil_)  ! sm content wilting point (mm)
      parameter(jsmr=jwpw+lsoil_)  ! Soil moisture residual moisture (mm)
      parameter(jsmx=jsmr+lsoil_)  ! maximum soil moisture (mm)
      parameter(jdpn=jsmx+lsoil_)  ! thickness of soil node (m)
      parameter(jsxn=jdpn+nsoil_)  ! maximum sm at soil node (m3/m3)
      parameter(jepn=jsxn+nsoil_)  ! Para the vari of Ksat at soil node (N/A)
      parameter(jbbn=jepn+nsoil_)  ! bubbling pressure at soil node (cm)
      parameter(japn=jbbn+nsoil_)  ! para alpha at soil node
      parameter(jbtn=japn+nsoil_)  ! para beta at soil node
      parameter(jgmn=jbtn+nsoil_)  ! para gamma at soil node
      parameter(jlai=jgmn+nsoil_)  ! leaf area index
      parameter(jslz=jlai+1)       ! surface roughness of bare soil (m)
      parameter(jsnz=jslz+1)       ! surface roughness of snow pack (m)
      parameter(jsic=jsnz+1)       ! soil ice content (mm)
      parameter(jcsn=jsic+lsoil_)  ! canopy snow  (mm h2o)
      parameter(jrsn=jcsn+1)       ! snow density (kg/m^3)
      parameter(jtsn=jrsn+1)       ! snow surface temperature (K)
      parameter(jtpk=jtsn+1)       ! snow pack temperature (K)
      parameter(jsfw=jtpk+1)       ! surface snow water equivalent (mm h2o)
      parameter(jpkw=jsfw+1)       ! snow pack snow water equivalent (mm h2o)
      parameter(jlst=jpkw+1)       ! time step since last snow fall
#ifdef ADJSOILW
      parameter(jobp=jlst+1)       ! observed precipitation
#endif
#endif
c
#ifdef VICLSM2
c  viclsm2
c
      integer kslmb, nslmb
      parameter (kslmb=lsoil_*msub_, nslmb=nsoil_*msub_)
#ifdef ADJSOILW
      parameter(numsfcs=25+lsoil_*12+nsoil_*8+msub_*11+
     &          kslmb*4+nslmb+lalbd_+2+1)
      parameter(numsfcv=64)
#else
      parameter(numsfcs=25+lsoil_*12+nsoil_*8+msub_*11+
     &          kslmb*4+nslmb+lalbd_+2)
      parameter(numsfcv=63)
#endif
      integer jtsf,jsmc,jsno,jstc,jtg3,jzor,jcv ,jcvb,jcvt,jalb,
     2        jsli,jveg,jcpy,jf10,jvet,jrot,jalf,just,jffm,jffh,
     3        jprc,jsrf,jbif,jds ,jdsm,jws ,jcef,jexp,jkst,jdph,
     4        jbub,jqrt,jbkd,jsld,jwcr,jwpw,jsmr,jsmx,jdpn,jsxn,
     5        jepn,jbbn,japn,jbtn,jgmn,jslz,jsnz,jnve,jmfr,jmcp,
     6        jmvt,jlai,jmsn,jmsm,jmsi,jmst,jcsn,jrsn,jtsn,jtpk,
     7        jsfw,jpkw,jlst
#ifdef ADJSOILW
      integer jobp
#endif

      parameter(jtsf=     1     )  ! sst and land surface temperature
      parameter(jsmc=jtsf+1     )  ! soil moisture
      parameter(jsno=jsmc+lsoil_)  ! snow depth
      parameter(jstc=jsno+1     )  ! soil temperature
      parameter(jtg3=jstc+nsoil_)  ! deep soil temperature
      parameter(jzor=jtg3+1     )  ! surface roughness
      parameter(jcv =jzor+1     )  ! convective cloud cover
      parameter(jcvb=jcv +1     )  ! convective cloud base
      parameter(jcvt=jcvb+1     )  ! convective cloud top
      parameter(jalb=jcvt+1     )  ! albedo
      parameter(jsli=jalb+lalbd_)  ! land/sea ice/snow mask
      parameter(jveg=jsli+1     )  ! fraction of vegetation cover
      parameter(jcpy=jveg+1     )  ! canopy water content
      parameter(jf10=jcpy+1     )  ! 10m and sig1 conv factor
      parameter(jvet=jf10+1     )  ! vegetation type
      parameter(jrot=jvet+1     )  ! root content
      parameter(jalf=jrot+kslmb )  ! albedo fraction
      parameter(just=jalf+2     )  ! frictional velocity
      parameter(jffm=just+1     )  ! coefficient for momentun exchange
      parameter(jffh=jffm+1     )  ! coefficient for heat exchange
      parameter(jprc=jffh+1     )  ! precip for vic
      parameter(jsrf=jprc+1     )  ! sr flag for vic
      parameter(jbif=jsrf+1     )  ! Variable infil curve parameter (N/A)
      parameter(jds =jbif+1     )  ! Fract of Dsm nonlinear baseflow begins
      parameter(jdsm=jds +1     )  ! Maximum velocity of baseflow (mm/day)
      parameter(jws =jdsm+1     )  ! Fract maxi sm nonlinear baseflow occurs
      parameter(jcef=jws +1     )  ! c
      parameter(jexp=jcef+1     )  ! Para the vari of Ksat with sm
      parameter(jkst=jexp+lsoil_)  ! Saturated hydrologic conductivity (mm/day)
      parameter(jdph=jkst+lsoil_)  ! thickness of soil layer (m)
      parameter(jbub=jdph+lsoil_)  ! Bubbling pressure of soil layer (cm)
      parameter(jqrt=jbub+lsoil_)  ! Quartz content of soil layer (fraction)
      parameter(jbkd=jqrt+lsoil_)  ! Bulk density of soil layer (kg/m3)
      parameter(jsld=jbkd+lsoil_)  ! Soil density of soil layer (kg/m3)
      parameter(jwcr=jsld+lsoil_)  ! sm content at the critical point (mm)
      parameter(jwpw=jwcr+lsoil_)  ! sm content wilting point (mm)
      parameter(jsmr=jwpw+lsoil_)  ! Soil moisture residual moisture (mm)
      parameter(jsmx=jsmr+lsoil_)  ! maximum soil moisture (mm)
      parameter(jdpn=jsmx+lsoil_)  ! thickness of soil node (m)
      parameter(jsxn=jdpn+nsoil_)  ! maximum sm at soil node (m3/m3)
      parameter(jepn=jsxn+nsoil_)  ! Para the vari of Ksat at soil node (N/A)
      parameter(jbbn=jepn+nsoil_)  ! bubbling pressure at soil node (cm)
      parameter(japn=jbbn+nsoil_)  ! para alpha at soil node
      parameter(jbtn=japn+nsoil_)  ! para beta at soil node
      parameter(jgmn=jbtn+nsoil_)  ! para gamma at soil node
      parameter(jslz=jgmn+nsoil_)  ! surface roughness of bare soil (m)
      parameter(jsnz=jslz+1     )  ! surface roughness of snow pack (m)
      parameter(jnve=jsnz+1     )  ! number of tiles in a grid (N/A)
      parameter(jmfr=jnve+1     )  ! vegetation fraction for tiles
      parameter(jmcp=jmfr+msub_ )  ! canopy water for tiles
      parameter(jmvt=jmcp+msub_ )  ! vegetation type for tiles
      parameter(jlai=jmvt+msub_ )  ! leaf area index
      parameter(jmsn=jlai+msub_ )  ! snow depth for tiles
      parameter(jmsm=jmsn+msub_ )  ! snow moisture for tiles (mm)
      parameter(jmsi=jmsm+kslmb )  ! soil ice content for tiles (mm)
      parameter(jmst=jmsi+kslmb )  ! soil temperature for tiles (K)
      parameter(jcsn=jmst+nslmb )  ! canopy snow  (mm h2o)
      parameter(jrsn=jcsn+msub_ )  ! snow density (kg/m^3)
      parameter(jtsn=jrsn+msub_ )  ! snow surface temperature (K)
      parameter(jtpk=jtsn+msub_ )  ! snow pack temperature (K)
      parameter(jsfw=jtpk+msub_ )  ! surface snow water equivalent (mm h2o)
      parameter(jpkw=jsfw+msub_ )  ! snow pack snow water equivalent (mm h2o)
      parameter(jlst=jpkw+msub_ )  ! time step since last snow fall
#ifdef ADJSOILW
      parameter(jobp=jlst+msub_ )  ! observed precipitation
#endif
#endif
c
c  2. number of layers
c
c  define sfc file variables and their properties.
c
#ifdef OSULSM1
c osulsm1
c
      integer ksfc(numsfcs)
      data ksfc(jtsf)/1/
      data ksfc(jsmc)/lsoil_/
      data ksfc(jsno)/1/
      data ksfc(jstc)/lsoil_/
      data ksfc(jtg3)/1/
      data ksfc(jzor)/1/
      data ksfc(jcv )/1/
      data ksfc(jcvb)/1/
      data ksfc(jcvt)/1/
      data ksfc(jalb)/lalbd_/
      data ksfc(jsli)/1/
      data ksfc(jplr)/1/
      data ksfc(jcpy)/1/
      data ksfc(jf10)/1/
#ifdef ADJSOILW
      data ksfc(jobp)/1/
#endif
#endif
c
#ifdef OSULSM2
c osulsm2
c
      integer ksfc(numsfcs)
      data ksfc(jtsf)/1/
      data ksfc(jsmc)/lsoil_/
      data ksfc(jsno)/1/
      data ksfc(jstc)/lsoil_/
      data ksfc(jtg3)/1/
      data ksfc(jzor)/1/
      data ksfc(jcv )/1/
      data ksfc(jcvb)/1/
      data ksfc(jcvt)/1/
      data ksfc(jalb)/lalbd_/
      data ksfc(jsli)/1/
      data ksfc(jveg)/1/
      data ksfc(jcpy)/1/
      data ksfc(jf10)/1/
      data ksfc(jvet)/1/
      data ksfc(jsot)/1/
      data ksfc(jalf)/2/
      data ksfc(just)/1/
      data ksfc(jffm)/1/
      data ksfc(jffh)/1/
#ifdef ADJSOILW
      data ksfc(jobp)/1/
#endif
#endif
c
#ifdef NOALSM1
c
      integer ksfc(numsfcs)
      data ksfc(jtsf)/1/
      data ksfc(jsmc)/lsoil_/
      data ksfc(jsno)/1/
      data ksfc(jstc)/lsoil_/
      data ksfc(jtg3)/1/
      data ksfc(jzor)/1/
      data ksfc(jcv )/1/
      data ksfc(jcvb)/1/
      data ksfc(jcvt)/1/
      data ksfc(jalb)/lalbd_/
      data ksfc(jsli)/1/
      data ksfc(jveg)/1/
      data ksfc(jcpy)/1/
      data ksfc(jf10)/1/
      data ksfc(jvet)/1/
      data ksfc(jsot)/1/
      data ksfc(jalf)/2/
      data ksfc(just)/1/
      data ksfc(jffm)/1/
      data ksfc(jffh)/1/
      data ksfc(jprc)/1/
      data ksfc(jsrf)/1/
      data ksfc(jsnd)/1/
      data ksfc(jslc)/lsoil_/
      data ksfc(jsmn)/1/
      data ksfc(jsmx)/1/
      data ksfc(jslo)/1/
      data ksfc(jsna)/1/
#ifdef ADJSOILW
      data ksfc(jobp)/1/
#endif
#endif
c
#ifdef VICLSM1
c viclsm1
c
      integer ksfc(numsfcs)
      data ksfc(jtsf)/1/
      data ksfc(jsmc)/lsoil_/
      data ksfc(jsno)/1/
      data ksfc(jstc)/nsoil_/
      data ksfc(jtg3)/1/
      data ksfc(jzor)/1/
      data ksfc(jcv )/1/
      data ksfc(jcvb)/1/
      data ksfc(jcvt)/1/
      data ksfc(jalb)/lalbd_/
      data ksfc(jsli)/1/
      data ksfc(jveg)/1/
      data ksfc(jcpy)/1/
      data ksfc(jf10)/1/
      data ksfc(jvet)/1/
      data ksfc(jrot)/lsoil_/
      data ksfc(jalf)/2/
      data ksfc(just)/1/
      data ksfc(jffm)/1/
      data ksfc(jffh)/1/
      data ksfc(jprc)/1/
      data ksfc(jsrf)/1/
      data ksfc(jbif)/1/
      data ksfc(jds )/1/
      data ksfc(jdsm)/1/
      data ksfc(jws )/1/
      data ksfc(jcef)/1/
      data ksfc(jexp)/lsoil_/
      data ksfc(jkst)/lsoil_/
      data ksfc(jdph)/lsoil_/
      data ksfc(jbub)/lsoil_/
      data ksfc(jqrt)/lsoil_/
      data ksfc(jbkd)/lsoil_/
      data ksfc(jsld)/lsoil_/
      data ksfc(jwcr)/lsoil_/
      data ksfc(jwpw)/lsoil_/
      data ksfc(jsmr)/lsoil_/
      data ksfc(jsmx)/lsoil_/
      data ksfc(jdpn)/nsoil_/
      data ksfc(jsxn)/nsoil_/
      data ksfc(jepn)/nsoil_/
      data ksfc(jbbn)/nsoil_/
      data ksfc(japn)/nsoil_/
      data ksfc(jbtn)/nsoil_/
      data ksfc(jgmn)/nsoil_/
      data ksfc(jlai)/1/
      data ksfc(jslz)/1/
      data ksfc(jsnz)/1/
      data ksfc(jsic)/lsoil_/
      data ksfc(jcsn)/1/
      data ksfc(jrsn)/1/
      data ksfc(jtsn)/1/
      data ksfc(jtpk)/1/
      data ksfc(jsfw)/1/
      data ksfc(jpkw)/1/
      data ksfc(jlst)/1/
#ifdef ADJSOILW
      data ksfc(jobp)/1/
#endif
#endif
c
#ifdef VICLSM2
c viclsm2
c
      integer ksfc(numsfcs)
      data ksfc(jtsf)/1/
      data ksfc(jsmc)/lsoil_/
      data ksfc(jsno)/1/
      data ksfc(jstc)/nsoil_/
      data ksfc(jtg3)/1/
      data ksfc(jzor)/1/
      data ksfc(jcv )/1/
      data ksfc(jcvb)/1/
      data ksfc(jcvt)/1/
      data ksfc(jalb)/lalbd_/
      data ksfc(jsli)/1/
      data ksfc(jveg)/1/
      data ksfc(jcpy)/1/
      data ksfc(jf10)/1/
      data ksfc(jvet)/1/
      data ksfc(jrot)/kslmb /
      data ksfc(jalf)/2/
      data ksfc(just)/1/
      data ksfc(jffm)/1/
      data ksfc(jffh)/1/
      data ksfc(jprc)/1/
      data ksfc(jsrf)/1/
      data ksfc(jbif)/1/
      data ksfc(jds )/1/
      data ksfc(jdsm)/1/
      data ksfc(jws )/1/
      data ksfc(jcef)/1/
      data ksfc(jexp)/lsoil_/
      data ksfc(jkst)/lsoil_/
      data ksfc(jdph)/lsoil_/
      data ksfc(jbub)/lsoil_/
      data ksfc(jqrt)/lsoil_/
      data ksfc(jbkd)/lsoil_/
      data ksfc(jsld)/lsoil_/
      data ksfc(jwcr)/lsoil_/
      data ksfc(jwpw)/lsoil_/
      data ksfc(jsmr)/lsoil_/
      data ksfc(jsmx)/lsoil_/
      data ksfc(jdpn)/nsoil_/
      data ksfc(jsxn)/nsoil_/
      data ksfc(jepn)/nsoil_/
      data ksfc(jbbn)/nsoil_/
      data ksfc(japn)/nsoil_/
      data ksfc(jbtn)/nsoil_/
      data ksfc(jgmn)/nsoil_/
      data ksfc(jslz)/1/
      data ksfc(jsnz)/1/
      data ksfc(jnve)/1/
      data ksfc(jmfr)/msub_/
      data ksfc(jmcp)/msub_/
      data ksfc(jmvt)/msub_/
      data ksfc(jlai)/msub_/
      data ksfc(jmsn)/msub_/
      data ksfc(jmsm)/kslmb/
      data ksfc(jmsi)/kslmb/
      data ksfc(jmst)/nslmb/
      data ksfc(jcsn)/msub_/
      data ksfc(jrsn)/msub_/
      data ksfc(jtsn)/msub_/
      data ksfc(jtpk)/msub_/
      data ksfc(jsfw)/msub_/
      data ksfc(jpkw)/msub_/
      data ksfc(jlst)/msub_/
#ifdef ADJSOILW
      data ksfc(jobp)/1/
#endif
#endif
c
c  define sfc file variables and their properties.
c
c  3.  qc parameters (vmaxmin)
c    max over ocean   without sea ice, min over ocean   without sea ice
c    max over  land   without snow   , min over land    without snow
c    max over sea ice without snow   , min over sea ice without snow
c    max over  land   with    snow   , min over land    with    snow
c    max over sea ice with    snow   , min over sea ice with    snow
c     does not check if vmaxmin(1,1,:).lt.vmaxmin(2,1,:)
c     these checks are applied to the value after scaling
c
#ifdef OSULSM1
c osulsm1
c
      real vmaxmin(2,5,numsfcs)
      integer isf,jsf
      data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)   ! sst and land temp
     &       / 313.00, 271.21,
     &         400.00, 100.00,
     &         271.21, 100.00,
     &         273.16, 100.00,
     &         273.16, 100.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /  -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0./
#else
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /   0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55/
#endif
      data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)   ! snow depth
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)   ! soilt 
     &       / 313.00, 200.00,
     &         353.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)   ! deep soil temp
     &       / 310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)   ! surface roughness
     &       /   1.00, 1.e-05,
     &         300.00,   0.05,
     &           1.00,   0.05,
     &         300.00,   0.05,
     &           1.00,   0.05/
      data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)   ! convective cloud cover
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)   ! 1-type albedo
     &       /   0.06,   0.06,
     &           0.80,   0.06,
     &           0.80,   0.80,
     &           0.80,   0.06,
     &           0.80,   0.80/
      data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)   ! land-sea-seaice mask
     &       /     0.,     0.,
     &             1.,     1.,
     &             2.,     2.,
     &             1.,     1.,
     &             2.,     2./
      data ((vmaxmin(isf,jsf,jplr),isf=1,2),jsf=1,5)   ! plant resis
     &       / 1000.0,    0.0,
     &         1000.0,    0.0,
     &         1000.0,    0.0,
     &         1000.0,    0.0,
     &         1000.0,    0.0/
      data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)   ! canopy water
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)   ! 10m and sig1 conv factor
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jobp),isf=1,2),jsf=1,5)   ! observed precipitation
     &       /    -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1./
#endif
#endif
c
#ifdef OSULSM2
c osulsm2
c
      real vmaxmin(2,5,numsfcs)
      integer isf,jsf
      data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)   ! sst and land temp
     &       / 313.00, 271.21,
     &         400.00, 100.00,
     &         271.21, 100.00,
     &         273.16, 100.00,
     &         273.16, 100.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /  -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0./
#else
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /   0.55,   0.55,
     &           0.55,   0.00,
     &           0.55,   0.55,
     &           0.55,   0.00,
     &           0.55,   0.55/
#endif
      data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)   ! snow depth
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)   ! soilt 
     &       / 313.00, 200.00,
     &         353.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)   ! deep soil temp
     &       / 310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)   ! surface roughness
     &       /   1.00, 1.e-05,
     &         300.00,   0.05,
     &           1.00,   0.05,
     &         300.00,   0.05,
     &           1.00,   0.05/
      data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)   ! convective cloud cover
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)   ! 4-type albedo
     &       /   0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01/
      data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)   ! land-sea-seaice mask
     &       /     0.,     0.,
     &             1.,     1.,
     &             2.,     2.,
     &             1.,     1.,
     &             2.,     2./
      data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)   ! vegetation cover
     &       /     0.,     0.,
     &             1.,     0.,
     &             1.,     0.,
     &             1.,     0.,
     &             1.,     0./
      data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)   ! canopy water
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)   ! 10m and sig1 conv factor
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
#ifdef USGS
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          12.00,   1.00,
     &          12.00,  12.00,
     &          12.00,   1.00,
     &          12.00,  12.00/
#else
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          13.00,   1.00,
     &          13.00,  13.00,
     &          13.00,   1.00,
     &          13.00,  13.00/
#endif
#ifdef USGS
      data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)   ! soil type
     &       /  14.00,  14.00,
     &          15.00,   0.00,
     &           0.00,   0.00,
     &          15.00,   0.00,
     &           0.00,   0.00/
#else
      data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)   ! soil type
     &       /   0.00,   0.00,
     &           9.00,   1.00,
     &           9.00,   9.00,
     &           9.00,   1.00,
     &           9.00,   9.00/
#endif
      data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)   ! albedo frac
     &       /   0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00/
      data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)   ! ustar
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)   ! drag coeff for m
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)   ! drag coeff for h
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jobp),isf=1,2),jsf=1,5)   ! observed precipitation
     &       /    -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1./
#endif
#endif
c
#ifdef NOALSM1
c noalsm1
c
      real vmaxmin(2,5,numsfcs)
      integer isf,jsf
      data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)   ! sst and land temp
     &       / 313.00, 271.21,
     &         400.00, 100.00,
     &         271.21, 100.00,
     &         273.16, 100.00,
     &         273.16, 100.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /  -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0./
#else
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /   0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55/
#endif
      data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)   ! snow depth
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)   ! soilt 
     &       / 313.00, 200.00,
     &         353.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)   ! deep soil temp
     &       / 310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)   ! surface roughness
     &       /   1.00, 1.e-05,
     &         300.00,   0.05,
     &           1.00,   0.05,
     &         300.00,   0.05,
     &           1.00,   0.05/
      data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)   ! convective cloud cover
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)   ! 4-type albedo
     &       /   0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01/
      data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)   ! land-sea-seaice mask
     &       /     0.,     0.,
     &             1.,     1.,
     &             2.,     2.,
     &             1.,     1.,
     &             2.,     2./
      data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)   ! vegetation cover
     &       /     0.,     0.,
     &             1.,     0.,
     &             0.,     0.,
     &             1.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)   ! canopy water
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)   ! 10m and sig1 conv factor
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
#ifdef USGS
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          12.00,   1.00,
     &          12.00,  12.00,
     &          12.00,   1.00,
     &          12.00,  12.00/
#else
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          13.00,   1.00,
     &          13.00,  13.00,
     &          13.00,   1.00,
     &          13.00,  13.00/
#endif
#ifdef USGS
      data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)   ! soil type
     &       /  14.00,  14.00,
     &          15.00,   1.00,
     &           0.00,   0.00,
     &          15.00,   1.00,
     &           0.00,   0.00/
#else
      data ((vmaxmin(isf,jsf,jsot),isf=1,2),jsf=1,5)   ! soil type
     &       /   0.00,   0.00,
     &           9.00,   1.00,
     &           9.00,   9.00,
     &           9.00,   1.00,
     &           9.00,   9.00/
#endif
      data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)   ! albedo frac
     &       /   0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00/
      data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)   ! ustar
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)   ! drag coeff for m
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)   ! drag coeff for h
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jprc),isf=1,2),jsf=1,5)   ! precip
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jsrf),isf=1,2),jsf=1,5)   ! sr flag
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jsnd),isf=1,2),jsf=1,5)   ! snow depth
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jslc),isf=1,2),jsf=1,5)   ! slc
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jsmn),isf=1,2),jsf=1,5)   ! min veg cov
     &       /   1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00/
      data ((vmaxmin(isf,jsf,jsmx),isf=1,2),jsf=1,5)   ! max veg cov
     &       /   1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00/
      data ((vmaxmin(isf,jsf,jslo),isf=1,2),jsf=1,5)   ! slope type
     &       /   0.00,   0.00,
     &          12.00,   1.00,
     &          13.00,  13.00,
     &          12.00,   1.00,
     &          13.00,  13.00/
      data ((vmaxmin(isf,jsf,jsna),isf=1,2),jsf=1,5)   ! snow alb
     &       /   1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jobp),isf=1,2),jsf=1,5)   ! observed precipitation
     &       /    -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1./
#endif
#endif
c
#ifdef VICLSM1
c viclsm1
c
      real vmaxmin(2,5,numsfcs)
      integer isf,jsf
      data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)   ! sst and land temp
     &       / 313.00, 271.21,
     &         400.00, 100.00,
     &         271.21, 100.00,
     &         273.16, 100.00,
     &         273.16, 100.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /  -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0./
#else
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /   0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55/
#endif
      data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)   ! snow depth
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)   ! soilt 
     &       / 313.00, 200.00,
     &         353.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)   ! deep soil temp
     &       / 310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)   ! surface roughness
     &       /   1.00, 1.e-05,
     &         300.00,   0.05,
     &           1.00,   0.05,
     &         300.00,   0.05,
     &           1.00,   0.05/
      data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)   ! convective cloud cover
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)   ! 4-type albedo
     &       /   0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01/
      data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)   ! land-sea-seaice mask
     &       /     0.,     0.,
     &             1.,     1.,
     &             2.,     2.,
     &             1.,     1.,
     &             2.,     2./
      data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)   ! vegetation cover
     &       /     0.,     0.,
     &             1.,     0.,
     &             0.,     0.,
     &             1.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)   ! canopy water
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)   ! 10m and sig1 conv factor
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
#ifdef USGS
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          12.00,   1.00,
     &          12.00,  12.00,
     &          12.00,   1.00,
     &          12.00,  12.00/
#else
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          13.00,   1.00,
     &          13.00,  13.00,
     &          13.00,   1.00,
     &          13.00,  13.00/
#endif
      data ((vmaxmin(isf,jsf,jrot),isf=1,2),jsf=1,5)   ! root content
     &       /   1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00/
      data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)   ! albedo frac
     &       /   0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00/
      data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)   ! ustar
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)   ! drag coeff for m
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)   ! drag coeff for h
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jprc),isf=1,2),jsf=1,5)   ! precip
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jsrf),isf=1,2),jsf=1,5)   ! sr flag
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jbif),isf=1,2),jsf=1,5)   ! bif
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jds ),isf=1,2),jsf=1,5)   ! ds
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jdsm),isf=1,2),jsf=1,5)   ! dsm
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jws ),isf=1,2),jsf=1,5)   ! ws
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jcef),isf=1,2),jsf=1,5)   ! cef
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jexp),isf=1,2),jsf=1,5)   ! expt
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jkst),isf=1,2),jsf=1,5)   ! kst
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jdph),isf=1,2),jsf=1,5)   ! dph
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbub),isf=1,2),jsf=1,5)   ! bub
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jqrt),isf=1,2),jsf=1,5)   ! qrt
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbkd),isf=1,2),jsf=1,5)   ! bkd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsld),isf=1,2),jsf=1,5)   ! sld
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jwcr),isf=1,2),jsf=1,5)   ! wcr
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jwpw),isf=1,2),jsf=1,5)   ! wpw
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsmr),isf=1,2),jsf=1,5)   ! smr
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsmx),isf=1,2),jsf=1,5)   ! smx
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jdpn),isf=1,2),jsf=1,5)   ! dphnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsxn),isf=1,2),jsf=1,5)   ! smxnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jepn),isf=1,2),jsf=1,5)   ! exptnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbbn),isf=1,2),jsf=1,5)   ! bubnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,japn),isf=1,2),jsf=1,5)   ! alpnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbtn),isf=1,2),jsf=1,5)   ! betnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jgmn),isf=1,2),jsf=1,5)   ! gamnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jlai),isf=1,2),jsf=1,5)   ! flai
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jslz),isf=1,2),jsf=1,5)   ! silz
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsnz),isf=1,2),jsf=1,5)   ! snwz
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsic),isf=1,2),jsf=1,5)   ! sic
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jcsn),isf=1,2),jsf=1,5)   ! csn
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jrsn),isf=1,2),jsf=1,5)   ! rsn
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jtsn),isf=1,2),jsf=1,5)   ! tsf
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jtpk),isf=1,2),jsf=1,5)   ! tpk
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsfw),isf=1,2),jsf=1,5)   ! sfw
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jpkw),isf=1,2),jsf=1,5)   ! pkw
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jlst),isf=1,2),jsf=1,5)   ! last snowfall
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jobp),isf=1,2),jsf=1,5)   ! observed precipitation
     &       /    -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1./
#endif
#endif
c
#ifdef VICLSM2
c viclsm2
c
      real vmaxmin(2,5,numsfcs)
      integer isf,jsf
      data ((vmaxmin(isf,jsf,jtsf),isf=1,2),jsf=1,5)   ! sst and land temp
     &       / 313.00, 271.21,
     &         400.00, 100.00,
     &         271.21, 100.00,
     &         273.16, 100.00,
     &         273.16, 100.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /  -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0.,
     &          -1.,   0./
#else
      data ((vmaxmin(isf,jsf,jsmc),isf=1,2),jsf=1,5)   ! soilm 
     &       /   0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55,
     &           0.55,   0.05,
     &           0.55,   0.55/
#endif
      data ((vmaxmin(isf,jsf,jsno),isf=1,2),jsf=1,5)   ! snow depth
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jstc),isf=1,2),jsf=1,5)   ! soilt 
     &       / 313.00, 200.00,
     &         353.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jtg3),isf=1,2),jsf=1,5)   ! deep soil temp
     &       / 310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00,
     &         310.00, 200.00/
      data ((vmaxmin(isf,jsf,jzor),isf=1,2),jsf=1,5)   ! surface roughness
     &       /   1.00, 1.e-05,
     &         300.00,   0.05,
     &           1.00,   0.05,
     &         300.00,   0.05,
     &           1.00,   0.05/
      data ((vmaxmin(isf,jsf,jcv ),isf=1,2),jsf=1,5)   ! convective cloud cover
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvb),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcvt),isf=1,2),jsf=1,5)   ! convective cloud base
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jalb),isf=1,2),jsf=1,5)   ! 4-type albedo
     &       /   0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01,
     &           0.80,   0.01,
     &           0.01,   0.01/
      data ((vmaxmin(isf,jsf,jsli),isf=1,2),jsf=1,5)   ! land-sea-seaice mask
     &       /     0.,     0.,
     &             1.,     1.,
     &             2.,     2.,
     &             1.,     1.,
     &             2.,     2./
      data ((vmaxmin(isf,jsf,jveg),isf=1,2),jsf=1,5)   ! vegetation cover
     &       /     0.,     0.,
     &             1.,     0.,
     &             0.,     0.,
     &             1.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jcpy),isf=1,2),jsf=1,5)   ! canopy water
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jf10),isf=1,2),jsf=1,5)   ! 10m and sig1 conv factor
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
#ifdef USGS
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          12.00,   1.00,
     &          12.00,  12.00,
     &          12.00,   1.00,
     &          12.00,  12.00/
#else
      data ((vmaxmin(isf,jsf,jvet),isf=1,2),jsf=1,5)   ! veg type
     &       /   0.00,   0.00,
     &          13.00,   1.00,
     &          13.00,  13.00,
     &          13.00,   1.00,
     &          13.00,  13.00/
#endif
      data ((vmaxmin(isf,jsf,jrot),isf=1,2),jsf=1,5)   ! root content
     &       /   1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00,
     &           1.00,   0.00/
      data ((vmaxmin(isf,jsf,jalf),isf=1,2),jsf=1,5)   ! albedo frac
     &       /   0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00,
     &           1.00,   0.00,
     &           0.00,   0.00/
      data ((vmaxmin(isf,jsf,just),isf=1,2),jsf=1,5)   ! ustar
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffm),isf=1,2),jsf=1,5)   ! drag coeff for m
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jffh),isf=1,2),jsf=1,5)   ! drag coeff for h
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jprc),isf=1,2),jsf=1,5)   ! precip
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jsrf),isf=1,2),jsf=1,5)   ! sr flag
     &       /    -1.,     1.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0.,
     &             0.,     0./
      data ((vmaxmin(isf,jsf,jbif),isf=1,2),jsf=1,5)   ! bif
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jds ),isf=1,2),jsf=1,5)   ! ds
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jdsm),isf=1,2),jsf=1,5)   ! dsm
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jws ),isf=1,2),jsf=1,5)   ! ws
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jcef),isf=1,2),jsf=1,5)   ! cef
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jexp),isf=1,2),jsf=1,5)   ! expt
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jkst),isf=1,2),jsf=1,5)   ! kst
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jdph),isf=1,2),jsf=1,5)   ! dph
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbub),isf=1,2),jsf=1,5)   ! bub
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jqrt),isf=1,2),jsf=1,5)   ! qrt
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbkd),isf=1,2),jsf=1,5)   ! bkd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsld),isf=1,2),jsf=1,5)   ! sld
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jwcr),isf=1,2),jsf=1,5)   ! wcr
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jwpw),isf=1,2),jsf=1,5)   ! wpw
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsmr),isf=1,2),jsf=1,5)   ! smr
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsmx),isf=1,2),jsf=1,5)   ! smx
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jdpn),isf=1,2),jsf=1,5)   ! dphnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsxn),isf=1,2),jsf=1,5)   ! smxnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jepn),isf=1,2),jsf=1,5)   ! exptnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbbn),isf=1,2),jsf=1,5)   ! bubnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,japn),isf=1,2),jsf=1,5)   ! alpnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jbtn),isf=1,2),jsf=1,5)   ! betnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jgmn),isf=1,2),jsf=1,5)   ! gamnd
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jslz),isf=1,2),jsf=1,5)   ! silz
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsnz),isf=1,2),jsf=1,5)   ! snwz
     &       /  0.0000,  1.00,
     &          0.0005,  1.00,
     &          0.0005,  1.00,
     &          0.0005,  1.00,
     &          0.0005,  1.00/
      data ((vmaxmin(isf,jsf,jnve),isf=1,2),jsf=1,5)   ! nve
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmfr),isf=1,2),jsf=1,5)   ! mfr
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmcp),isf=1,2),jsf=1,5)   ! mcp
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmvt),isf=1,2),jsf=1,5)   ! mvt
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jlai),isf=1,2),jsf=1,5)   ! flai
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmsn),isf=1,2),jsf=1,5)   ! msn
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmsm),isf=1,2),jsf=1,5)   ! msm
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmsi),isf=1,2),jsf=1,5)   ! msic
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jmst),isf=1,2),jsf=1,5)   ! mst
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jcsn),isf=1,2),jsf=1,5)   ! csn
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jrsn),isf=1,2),jsf=1,5)   ! rsn
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jtsn),isf=1,2),jsf=1,5)   ! tsf
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jtpk),isf=1,2),jsf=1,5)   ! tpk
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jsfw),isf=1,2),jsf=1,5)   ! sfw
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jpkw),isf=1,2),jsf=1,5)   ! pkw
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
      data ((vmaxmin(isf,jsf,jlst),isf=1,2),jsf=1,5)   ! last snowfall
     &       /  -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00,
     &          -1.00,  1.00/
#ifdef ADJSOILW
      data ((vmaxmin(isf,jsf,jobp),isf=1,2),jsf=1,5)   ! observed precipitation
     &       /    -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1.,
     &            -1.,     1./
#endif
#endif
c
c
c end of varsfc.h
c
