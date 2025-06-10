c
c  define grib surface variables and their properties.
c
c  if you are introducing new grib fields, you need to increase 
c  numgrbs and add the subsequent 9 parameters for each new field.
c
      integer numgrbs
      parameter(numgrbs=54)
      integer igb
c
      integer igla,imxi,imsk,iais,iscv,isno,itsf,iab4,iab1,ialf,
     1        iso2,isn1,izor,iplr,itg3,ito2,itn1,iveg,ivet,isot,islo,
     2        isna,ismx,ismn
      integer ivgv,ivtv,ivrt,ivaf,ibif,ids ,idsm,iws ,icef,iexp,ikst,
     1        idph,ibub,iqrt,ibkd,isld,iwcr,iwpw,ismr,islz,isnz,ilai,
     2        ivg2,irt1,irt2,irt3,isv1,itv1,isv2,itv2,iobp
c
      parameter(igla= 1)  ! glacier points
      parameter(imxi= 2)  ! max ice extent
      parameter(imsk= 3)  ! land sea mask
      parameter(iais= 4)  ! sea ice mask/concentration
      parameter(iscv= 5)  ! snow cover
      parameter(isno= 6)  ! snow depth
      parameter(itsf= 7)  ! sst and land skin temp
      parameter(iab4= 8)  ! 4-type albedo
      parameter(iab1= 9)  ! 1-type albedo
      parameter(ialf=10)  ! albedo fraction
      parameter(iso2=11)  ! soil moisture osu2
      parameter(isn1=12)  ! soil moisture noa1
      parameter(izor=13)  ! surface roughness
      parameter(iplr=14)  ! stomato resistance
      parameter(itg3=15)  ! deep layer soil temp
      parameter(ito2=16)  ! soil temp osu2
      parameter(itn1=17)  ! soil temp noa1
      parameter(iveg=18)  ! vegetation cover
      parameter(ivet=19)  ! vegetation type
      parameter(isot=20)  ! soil type
      parameter(islo=21)  ! slope type
      parameter(isna=22)  ! snow albedo
      parameter(ismx=23)  ! max vegetation cover
      parameter(ismn=24)  ! min vegetation cover
c-vic addition
      parameter(ivgv=25)  ! vic vegetation cover
      parameter(ivtv=26)  ! vic vegetation type
      parameter(ivrt=27)  ! root content for one vegetation type
      parameter(ibif=28)  ! vic curve parameter (N/A)
      parameter(ids =29)  ! Fract of dsm nonlinear baseflow begins
      parameter(idsm=30)  ! Maximum velocity of baseflow (mm/day)
      parameter(iws =31)  ! Fract maxi sm nonlinear baseflow occurs
      parameter(icef=32)  ! exponent used in infiltration curve (N/A)
      parameter(iexp=33)  ! Para the vari of Ksat with sm (N/A)
      parameter(ikst=34)  ! Saturated hydrologic conductivity (mm/day)
      parameter(idph=35)  ! Soil layer thickness (m)
      parameter(ibub=36)  ! Bubbling pressure of soil layer (cm)
      parameter(iqrt=37)  ! Quartz content of soil layer (fraction)
      parameter(ibkd=38)  ! Bulk density of soil layer (kg/m3)
      parameter(isld=39)  ! Soil density of soil layer (kg/m3)
      parameter(iwcr=40)  ! sm content at the critical point (mm)
      parameter(iwpw=41)  ! sm content wilting point (mm)
      parameter(ismr=42)  ! sm residual moisture (mm)
      parameter(islz=43)  ! surface roughness of bare soil (m)
      parameter(isnz=44)  ! surface roughness of snow pack (m)
      parameter(ilai=45)  ! leaf area index
      parameter(ivg2=46)  ! vic2 vegetation cover
      parameter(irt1=47)  ! root content for layer 1
      parameter(irt2=48)  ! root content for layer 2
      parameter(irt3=49)  ! root content for layer 3
      parameter(isv1=50)  ! soil moisture vic1
      parameter(itv1=51)  ! soil temperature vic1
      parameter(isv2=52)  ! soil moisture vic2
      parameter(itv2=53)  ! soil temperature vic2
      parameter(iobp=54)  ! observed precipitation
c
c  define various parameters for each field
c
c  1. number of layers (kgrb)
c
      integer kgrb(numgrbs)
      data kgrb(igla)/1/  ! glacier points
      data kgrb(imxi)/1/  ! max ice extent
      data kgrb(imsk)/1/  ! land sea mask
      data kgrb(iais)/1/  ! sea ice mask/concentration
      data kgrb(iscv)/1/  ! snow cover
      data kgrb(isno)/1/  ! snow depth
      data kgrb(itsf)/1/  ! sst and land skin temp
      data kgrb(iab4)/4/  ! 4-type albedo
      data kgrb(iab1)/1/  ! 1-type albedo
      data kgrb(ialf)/2/  ! albedo fraction
      data kgrb(iso2)/2/  ! soil moisture osu2
      data kgrb(isn1)/4/  ! soil moisture noa1
      data kgrb(izor)/1/  ! surface roughness
      data kgrb(iplr)/1/  ! stomata resistance
      data kgrb(itg3)/1/  ! deep layer soil temp
      data kgrb(ito2)/2/  ! soil temp osu2
      data kgrb(itn1)/4/  ! soil temp noa1
      data kgrb(iveg)/1/  ! vegetation cover
      data kgrb(ivet)/1/  ! vegetation type
      data kgrb(isot)/1/  ! soil type
      data kgrb(islo)/1/  ! slope type
      data kgrb(isna)/1/  ! snow albedo
      data kgrb(ismx)/1/  ! max vegetation cover
      data kgrb(ismn)/1/  ! min vegetation cover
c-vic addition
      data kgrb(ivgv)/1/  ! vic1 vegetation cover
      data kgrb(ivtv)/1/  ! vic1 vegetation type
      data kgrb(ivrt)/3/  ! root content for one vegetation type
      data kgrb(ibif)/1/  ! vic curve data kgrb (N/A)
      data kgrb(ids )/1/  ! Fract of dsm nonlinear baseflow begins
      data kgrb(idsm)/1/  ! Maximum velocity of baseflow (mm/day)
      data kgrb(iws )/1/  ! Fract maxi sm nonlinear baseflow occurs
      data kgrb(icef)/1/  ! exponent used in infiltration curve (N/A)
      data kgrb(iexp)/3/  ! Para the vari of Ksat with sm (N/A)
      data kgrb(ikst)/3/  ! Saturated hydrologic conductivity (mm/day)
      data kgrb(idph)/3/  ! Soil layer thickness (m)
      data kgrb(ibub)/3/  ! Bubbling pressure of soil layer (cm)
      data kgrb(iqrt)/3/  ! Quartz content of soil layer (fraction)
      data kgrb(ibkd)/3/  ! Bulk density of soil layer (kg/m3)
      data kgrb(isld)/3/  ! Soil density of soil layer (kg/m3)
      data kgrb(iwcr)/3/  ! sm content at the critical point (mm)
      data kgrb(iwpw)/3/  ! sm content wilting point (mm)
      data kgrb(ismr)/3/  ! sm residual moisture (mm)
      data kgrb(islz)/1/  ! surface roughness of bare soil (m)
      data kgrb(isnz)/1/  ! surface roughness of snow pack (m)
      data kgrb(ilai)/11/ ! leaf area index
      data kgrb(ivg2)/12/ ! vic2 vegetation cover
      data kgrb(irt1)/11/ ! root content for layer 1
      data kgrb(irt2)/11/ ! root content for layer 2
      data kgrb(irt3)/11/ ! root content for layer 3
      data kgrb(isv1)/3/  ! soil moisture vic1
      data kgrb(itv1)/5/  ! soil temperature vic1
      data kgrb(isv2)/3/  ! soil moisture vic2
      data kgrb(itv2)/5/  ! soil temperature vic2
      data kgrb(iobp)/1/  ! observed precipitation
c
c  2.  grib header record that identifies the field (kpd)
c       number of soil-layers/field-types should not exceed kgrbmax
c
      integer kgrbmax
      parameter(kgrbmax=100)
      integer kpd(25,kgrbmax,numgrbs)
      data (kpd(igb,1,igla),igb=1,25)/4*-1,238, -1,   -1,18*-1/ ! glcier points
      data (kpd(igb,1,imxi),igb=1,25)/4*-1, 91, -1,   -1,18*-1/ ! max ice extnt
      data (kpd(igb,1,imsk),igb=1,25)/4*-1, 91, -1,   -1,18*-1/ ! land sea mask
      data (kpd(igb,1,iais),igb=1,25)/4*-1, 91, -1,   -1,18*-1/ ! slimsk/icecon
      data (kpd(igb,1,iscv),igb=1,25)/4*-1,238, -1,   -1,18*-1/ ! snow cover
      data (kpd(igb,1,isno),igb=1,25)/4*-1, 65, -1,   -1,18*-1/ ! snow depth
      data (kpd(igb,1,itsf),igb=1,25)/4*-1, 11, -1,    0,18*-1/ ! sst and sknt
      data (kpd(igb,1,iab4),igb=1,25)/4*-1,212, -1,   -1,18*-1/ ! 4-type albd 1
      data (kpd(igb,2,iab4),igb=1,25)/4*-1,215, -1,   -1,18*-1/ ! 4-type albd 2
      data (kpd(igb,3,iab4),igb=1,25)/4*-1,213, -1,   -1,18*-1/ ! 4-type albe 3
      data (kpd(igb,4,iab4),igb=1,25)/4*-1,216, -1,   -1,18*-1/ ! 4-type albe 4
      data (kpd(igb,1,iab1),igb=1,25)/4*-1, 84, -1,   -1,18*-1/ ! 1-type albedo
      data (kpd(igb,1,ialf),igb=1,25)/4*-1,214, -1,   -1,18*-1/ ! albd frac 1
      data (kpd(igb,2,ialf),igb=1,25)/4*-1,217, -1,   -1,18*-1/ ! albd frac 2
      data (kpd(igb,1,iso2),igb=1,25)/4*-1,144, -1,   10,18*-1/ ! soilm osu2 1
      data (kpd(igb,2,iso2),igb=1,25)/4*-1,144, -1, 2760,18*-1/ ! soilm osu2 2
      data (kpd(igb,1,isn1),igb=1,25)/4*-1,144, -1,   10,18*-1/ ! soilm noa1 1
      data (kpd(igb,2,isn1),igb=1,25)/4*-1,144, -1, 2600,18*-1/ ! soilm noa1 2
      data (kpd(igb,3,isn1),igb=1,25)/4*-1,144, -1,10340,18*-1/ ! soilm noa1 3
      data (kpd(igb,4,isn1),igb=1,25)/4*-1,144, -1,25800,18*-1/ ! soilm noa1 4
      data (kpd(igb,1,izor),igb=1,25)/4*-1, 83, -1,   -1,18*-1/ ! sfc rough
      data (kpd(igb,1,iplr),igb=1,25)/4*-1,224, -1,   -1,18*-1/ ! plant resis
      data (kpd(igb,1,itg3),igb=1,25)/4*-1, 11,111,  500,18*-1/ ! deep soil t
      data (kpd(igb,1,ito2),igb=1,25)/4*-1, 11, -1,   10,18*-1/ ! soilt osu2 1
      data (kpd(igb,2,ito2),igb=1,25)/4*-1, 11, -1, 2760,18*-1/ ! soilt osu2 2
      data (kpd(igb,1,itn1),igb=1,25)/4*-1, 11, -1,   10,18*-1/ ! soilt noa1 1
      data (kpd(igb,2,itn1),igb=1,25)/4*-1, 11, -1, 2600,18*-1/ ! soilt noa1 2
      data (kpd(igb,3,itn1),igb=1,25)/4*-1, 11, -1,10340,18*-1/ ! soilt noa1 3
      data (kpd(igb,4,itn1),igb=1,25)/4*-1, 11, -1,25800,18*-1/ ! soilt noa1 4
      data (kpd(igb,1,iveg),igb=1,25)/4*-1, 87, -1,   -1,18*-1/ ! veg cover
      data (kpd(igb,1,ivet),igb=1,25)/4*-1,225, -1,   -1,18*-1/ ! veg type
      data (kpd(igb,1,isot),igb=1,25)/4*-1,230, -1,   -1,18*-1/ ! soil type
      data (kpd(igb,1,islo),igb=1,25)/4*-1,236, -1,   -1,18*-1/ ! slope type
      data (kpd(igb,1,isna),igb=1,25)/4*-1,255, -1,   -1,18*-1/ ! snow albedo
      data (kpd(igb,1,ismx),igb=1,25)/4*-1,255, -1,   -1,18*-1/ ! max veg cov
      data (kpd(igb,1,ismn),igb=1,25)/4*-1,255, -1,   -1,18*-1/ ! min veg cov
      data (kpd(igb,1,ivgv),igb=1,25)/4*-1, 87, -1,   -1,18*-1/ ! vic veg cover
      data (kpd(igb,1,ivtv),igb=1,25)/4*-1,225, -1,   -1,18*-1/ ! veg type
      data (kpd(igb,1,ivrt),igb=1,25)/4*-1,118,112,   10,18*-1/ ! root viclsm1 1
      data (kpd(igb,2,ivrt),igb=1,25)/4*-1,118,112, 2600,18*-1/ ! root viclsm1 2
      data (kpd(igb,3,ivrt),igb=1,25)/4*-1,118,112,10380,18*-1/ ! root viclsm1 3
      data (kpd(igb,1,ibif),igb=1,25)/4*-1,119, -1,   -1,18*-1/ ! binf
      data (kpd(igb,1,ids ),igb=1,25)/4*-1,140, -1,   -1,18*-1/ ! Fract of Dsm
      data (kpd(igb,1,idsm),igb=1,25)/4*-1,141, -1,   -1,18*-1/ ! Max baseflow
      data (kpd(igb,1,iws ),igb=1,25)/4*-1,142, -1,   -1,18*-1/ ! Fract maxi sm
      data (kpd(igb,1,icef),igb=1,25)/4*-1,151, -1,   -1,18*-1/ ! c
      data (kpd(igb,1,iexp),igb=1,25)/4*-1,143,112,   10,18*-1/ ! expt viclsm1 1
      data (kpd(igb,2,iexp),igb=1,25)/4*-1,143,112, 2600,18*-1/ ! expt viclsm1 2
      data (kpd(igb,3,iexp),igb=1,25)/4*-1,143,112,10380,18*-1/ ! expt viclsm1 3
      data (kpd(igb,1,ikst),igb=1,25)/4*-1,150,112,   10,18*-1/ ! ksat viclsm1 1
      data (kpd(igb,2,ikst),igb=1,25)/4*-1,150,112, 2600,18*-1/ ! ksat viclsm1 2
      data (kpd(igb,3,ikst),igb=1,25)/4*-1,150,112,10380,18*-1/ ! ksat viclsm1 3
      data (kpd(igb,1,idph),igb=1,25)/4*-1, 24,112,   10,18*-1/ ! laythick 1
      data (kpd(igb,2,idph),igb=1,25)/4*-1, 24,112, 2600,18*-1/ ! laythick 2
      data (kpd(igb,3,idph),igb=1,25)/4*-1, 24,112,10380,18*-1/ ! laythick 3
      data (kpd(igb,1,ibub),igb=1,25)/4*-1,152,112,   10,18*-1/ ! bub viclsm1 1
      data (kpd(igb,2,ibub),igb=1,25)/4*-1,152,112, 2600,18*-1/ ! bub viclsm1 2
      data (kpd(igb,3,ibub),igb=1,25)/4*-1,152,112,10380,18*-1/ ! bub viclsm1 3
      data (kpd(igb,1,iqrt),igb=1,25)/4*-1,154,112,   10,18*-1/ ! qrtz viclsm1 1
      data (kpd(igb,2,iqrt),igb=1,25)/4*-1,154,112, 2600,18*-1/ ! qrtz viclsm1 2
      data (kpd(igb,3,iqrt),igb=1,25)/4*-1,154,112,10380,18*-1/ ! qrtz viclsm1 3
      data (kpd(igb,1,ibkd),igb=1,25)/4*-1,173,112,   10,18*-1/ ! bulk viclsm1 1
      data (kpd(igb,2,ibkd),igb=1,25)/4*-1,173,112, 2600,18*-1/ ! bulk viclsm1 2
      data (kpd(igb,3,ibkd),igb=1,25)/4*-1,173,112,10380,18*-1/ ! bulk viclsm1 3
      data (kpd(igb,1,isld),igb=1,25)/4*-1,174,112,   10,18*-1/ ! sden  viclsm1 1
      data (kpd(igb,2,isld),igb=1,25)/4*-1,174,112, 2600,18*-1/ ! sden  viclsm1 2
      data (kpd(igb,3,isld),igb=1,25)/4*-1,174,112,10380,18*-1/ ! sden  viclsm1 3
      data (kpd(igb,1,iwcr),igb=1,25)/4*-1,180,112,   10,18*-1/ ! wcr viclsm1 1
      data (kpd(igb,2,iwcr),igb=1,25)/4*-1,180,112, 2600,18*-1/ ! wcr viclsm1 2
      data (kpd(igb,3,iwcr),igb=1,25)/4*-1,180,112,10380,18*-1/ ! wcr viclsm1 3
      data (kpd(igb,1,iwpw),igb=1,25)/4*-1,201,112,   10,18*-1/ ! wpw viclsm1 1
      data (kpd(igb,2,iwpw),igb=1,25)/4*-1,201,112, 2600,18*-1/ ! wpw viclsm1 2
      data (kpd(igb,3,iwpw),igb=1,25)/4*-1,201,112,10380,18*-1/ ! wpw viclsm1 3
      data (kpd(igb,1,ismr),igb=1,25)/4*-1,120,112,   10,18*-1/ ! smr 1
      data (kpd(igb,2,ismr),igb=1,25)/4*-1,120,112, 2600,18*-1/ ! smr 2
      data (kpd(igb,3,ismr),igb=1,25)/4*-1,120,112,10380,18*-1/ ! smr 3
      data (kpd(igb,1,islz),igb=1,25)/4*-1, 83, -1,   -1,18*-1/ ! silz
      data (kpd(igb,1,isnz),igb=1,25)/4*-1, 77, -1,   -1,18*-1/ ! snwz
      data (kpd(igb, 1,ilai),igb=1,25)/4*-1, 24,112,  256,18*-1/ ! lai
      data (kpd(igb, 2,ilai),igb=1,25)/4*-1, 24,112,  512,18*-1/ ! lai
      data (kpd(igb, 3,ilai),igb=1,25)/4*-1, 24,112,  768,18*-1/ ! lai
      data (kpd(igb, 4,ilai),igb=1,25)/4*-1, 24,112, 1024,18*-1/ ! lai
      data (kpd(igb, 5,ilai),igb=1,25)/4*-1, 24,112, 1280,18*-1/ ! lai
      data (kpd(igb, 6,ilai),igb=1,25)/4*-1, 24,112, 1536,18*-1/ ! lai
      data (kpd(igb, 7,ilai),igb=1,25)/4*-1, 24,112, 1792,18*-1/ ! lai
      data (kpd(igb, 8,ilai),igb=1,25)/4*-1, 24,112, 2048,18*-1/ ! lai
      data (kpd(igb, 9,ilai),igb=1,25)/4*-1, 24,112, 2304,18*-1/ ! lai
      data (kpd(igb,10,ilai),igb=1,25)/4*-1, 24,112, 2560,18*-1/ ! lai
      data (kpd(igb,11,ilai),igb=1,25)/4*-1, 24,112, 2816,18*-1/ ! lai
      data (kpd(igb, 1,ivg2),igb=1,25)/4*-1,  4,112,  256,18*-1/ ! veg cover
      data (kpd(igb, 2,ivg2),igb=1,25)/4*-1,  4,112,  512,18*-1/ ! veg cover
      data (kpd(igb, 3,ivg2),igb=1,25)/4*-1,  4,112,  768,18*-1/ ! veg cover
      data (kpd(igb, 4,ivg2),igb=1,25)/4*-1,  4,112, 1024,18*-1/ ! veg cover
      data (kpd(igb, 5,ivg2),igb=1,25)/4*-1,  4,112, 1280,18*-1/ ! veg cover
      data (kpd(igb, 6,ivg2),igb=1,25)/4*-1,  4,112, 1536,18*-1/ ! veg cover
      data (kpd(igb, 7,ivg2),igb=1,25)/4*-1,  4,112, 1792,18*-1/ ! veg cover
      data (kpd(igb, 8,ivg2),igb=1,25)/4*-1,  4,112, 2048,18*-1/ ! veg cover
      data (kpd(igb, 9,ivg2),igb=1,25)/4*-1,  4,112, 2304,18*-1/ ! veg cover
      data (kpd(igb,10,ivg2),igb=1,25)/4*-1,  4,112, 2560,18*-1/ ! veg cover
      data (kpd(igb,11,ivg2),igb=1,25)/4*-1,  4,112, 2816,18*-1/ ! veg cover
      data (kpd(igb,12,ivg2),igb=1,25)/4*-1,  4,112, 3072,18*-1/ ! veg cover
      data (kpd(igb, 1,irt1),igb=1,25)/4*-1,118,112,  256,18*-1/ ! rt1
      data (kpd(igb, 2,irt1),igb=1,25)/4*-1,118,112,  512,18*-1/ ! rt1
      data (kpd(igb, 3,irt1),igb=1,25)/4*-1,118,112,  768,18*-1/ ! rt1
      data (kpd(igb, 4,irt1),igb=1,25)/4*-1,118,112, 1024,18*-1/ ! rt1
      data (kpd(igb, 5,irt1),igb=1,25)/4*-1,118,112, 1280,18*-1/ ! rt1
      data (kpd(igb, 6,irt1),igb=1,25)/4*-1,118,112, 1536,18*-1/ ! rt1
      data (kpd(igb, 7,irt1),igb=1,25)/4*-1,118,112, 1792,18*-1/ ! rt1
      data (kpd(igb, 8,irt1),igb=1,25)/4*-1,118,112, 2048,18*-1/ ! rt1
      data (kpd(igb, 9,irt1),igb=1,25)/4*-1,118,112, 2304,18*-1/ ! rt1
      data (kpd(igb,10,irt1),igb=1,25)/4*-1,118,112, 2560,18*-1/ ! rt1
      data (kpd(igb,11,irt1),igb=1,25)/4*-1,118,112, 2816,18*-1/ ! rt1
      data (kpd(igb, 1,irt2),igb=1,25)/4*-1,119,112,  256,18*-1/ ! rt2
      data (kpd(igb, 2,irt2),igb=1,25)/4*-1,119,112,  512,18*-1/ ! rt2
      data (kpd(igb, 3,irt2),igb=1,25)/4*-1,119,112,  768,18*-1/ ! rt2
      data (kpd(igb, 4,irt2),igb=1,25)/4*-1,119,112, 1024,18*-1/ ! rt2
      data (kpd(igb, 5,irt2),igb=1,25)/4*-1,119,112, 1280,18*-1/ ! rt2
      data (kpd(igb, 6,irt2),igb=1,25)/4*-1,119,112, 1536,18*-1/ ! rt2
      data (kpd(igb, 7,irt2),igb=1,25)/4*-1,119,112, 1792,18*-1/ ! rt2
      data (kpd(igb, 8,irt2),igb=1,25)/4*-1,119,112, 2048,18*-1/ ! rt2
      data (kpd(igb, 9,irt2),igb=1,25)/4*-1,119,112, 2304,18*-1/ ! rt2
      data (kpd(igb,10,irt2),igb=1,25)/4*-1,119,112, 2560,18*-1/ ! rt2
      data (kpd(igb,11,irt2),igb=1,25)/4*-1,119,112, 2816,18*-1/ ! rt2
      data (kpd(igb, 1,irt3),igb=1,25)/4*-1,120,112,  256,18*-1/ ! rt3
      data (kpd(igb, 2,irt3),igb=1,25)/4*-1,120,112,  512,18*-1/ ! rt3
      data (kpd(igb, 3,irt3),igb=1,25)/4*-1,120,112,  768,18*-1/ ! rt3
      data (kpd(igb, 4,irt3),igb=1,25)/4*-1,120,112, 1024,18*-1/ ! rt3
      data (kpd(igb, 5,irt3),igb=1,25)/4*-1,120,112, 1280,18*-1/ ! rt3
      data (kpd(igb, 6,irt3),igb=1,25)/4*-1,120,112, 1536,18*-1/ ! rt3
      data (kpd(igb, 7,irt3),igb=1,25)/4*-1,120,112, 1792,18*-1/ ! rt3
      data (kpd(igb, 8,irt3),igb=1,25)/4*-1,120,112, 2048,18*-1/ ! rt3
      data (kpd(igb, 9,irt3),igb=1,25)/4*-1,120,112, 2304,18*-1/ ! rt3
      data (kpd(igb,10,irt3),igb=1,25)/4*-1,120,112, 2560,18*-1/ ! rt3
      data (kpd(igb,11,irt3),igb=1,25)/4*-1,120,112, 2816,18*-1/ ! rt3
      data (kpd(igb,1,isv1),igb=1,25)/4*-1,144, -1,   10,18*-1/ ! soilm vic1 1
      data (kpd(igb,2,isv1),igb=1,25)/4*-1,144, -1, 2600,18*-1/ ! soilm vic1 2
      data (kpd(igb,3,isv1),igb=1,25)/4*-1,144, -1,10380,18*-1/ ! soilm vic1 3
      data (kpd(igb,1,itv1),igb=1,25)/4*-1, 11, -1,   10,18*-1/ ! soilt vic1 1
      data (kpd(igb,2,itv1),igb=1,25)/4*-1, 11, -1, 2600,18*-1/ ! soilt vic1 2
      data (kpd(igb,3,itv1),igb=1,25)/4*-1, 11, -1,10380,18*-1/ ! soilt vic1 3
      data (kpd(igb,4,itv1),igb=1,25)/4*-1, 11, -1, 2601,18*-1/ ! soilt vic1 2.1
      data (kpd(igb,5,itv1),igb=1,25)/4*-1, 11, -1,10381,18*-1/ ! soilt vic1 3.1
      data (kpd(igb,1,isv2),igb=1,25)/4*-1,144, -1,   10,18*-1/ ! soilm vic2 1
      data (kpd(igb,2,isv2),igb=1,25)/4*-1,144, -1, 2600,18*-1/ ! soilm vic2 2
      data (kpd(igb,3,isv2),igb=1,25)/4*-1,144, -1,10380,18*-1/ ! soilm vic2 3
      data (kpd(igb,1,itv2),igb=1,25)/4*-1, 11, -1,   10,18*-1/ ! soilt vic2 1
      data (kpd(igb,2,itv2),igb=1,25)/4*-1, 11, -1, 2600,18*-1/ ! soilt vic2 2
      data (kpd(igb,3,itv2),igb=1,25)/4*-1, 11, -1,10380,18*-1/ ! soilt vic2 3
      data (kpd(igb,4,itv2),igb=1,25)/4*-1, 11, -1, 2601,18*-1/ ! soilt vic2 2.1
      data (kpd(igb,5,itv2),igb=1,25)/4*-1, 11, -1,10381,18*-1/ ! soilt vic2 3.1
      data (kpd(igb,1,iobp),igb=1,25)/4*-1, 61, -1,   -1,18*-1/ ! observed prcp
c
c  3. spacial interpolation type
c     inttyp=0  bilinear
c           =1  nearest neighbor
c           =2  dominant type
c
      integer inttyp(numgrbs)
      data inttyp(igla)/1/  ! glacier points
      data inttyp(imxi)/1/  ! max ice extent
      data inttyp(imsk)/0/  ! land sea mask (non grib on model grid)
      data inttyp(iais)/1/  ! sea ice mask/concentration
      data inttyp(iscv)/1/  ! snow cover
      data inttyp(isno)/0/  ! snow depth
      data inttyp(itsf)/0/  ! sst and land skin temp
      data inttyp(iab4)/0/  ! 4-type albedo
      data inttyp(iab1)/0/  ! 1-type albedo
      data inttyp(ialf)/0/  ! albedo fraction
      data inttyp(iso2)/0/  ! soil moisture osu2
      data inttyp(isn1)/0/  ! soil moisture noa1
      data inttyp(izor)/0/  ! surface roughness
      data inttyp(iplr)/0/  ! stomata resistance
      data inttyp(itg3)/0/  ! deep layer soil temp
      data inttyp(ito2)/0/  ! soil temp osu2
      data inttyp(itn1)/0/  ! soil temp noa1
      data inttyp(iveg)/0/  ! vegetation cover
      data inttyp(ivet)/2/  ! vegetation type
      data inttyp(isot)/2/  ! soil type
      data inttyp(islo)/2/  ! slope type
      data inttyp(isna)/0/  ! snow albedo
      data inttyp(ismx)/0/  ! max vegetation cover
      data inttyp(ismn)/0/  ! min vegetation cover
c-vic addition
      data inttyp(ivgv)/0/  ! vic vegetation cover
      data inttyp(ivtv)/2/  ! vic vegetation type
      data inttyp(ivrt)/1/  ! root content for one vegetation type
      data inttyp(ibif)/1/  ! vic curve parameter (N/A)
      data inttyp(ids )/1/  ! Fract of dsm nonlinear baseflow begins
      data inttyp(idsm)/1/  ! Maximum velocity of baseflow (mm/day)
      data inttyp(iws )/1/  ! Fract maxi sm nonlinear baseflow occurs
      data inttyp(icef)/1/  ! exponent used in infiltration curve (N/A)
      data inttyp(iexp)/1/  ! Para the vari of Ksat with sm (N/A)
      data inttyp(ikst)/1/  ! Saturated hydrologic conductivity (mm/day)
      data inttyp(idph)/1/  ! Soil layer thickness (m)
      data inttyp(ibub)/1/  ! Bubbling pressure of soil layer (cm)
      data inttyp(iqrt)/1/  ! Quartz content of soil layer (fraction)
      data inttyp(ibkd)/1/  ! Bulk density of soil layer (kg/m3)
      data inttyp(isld)/1/  ! Soil density of soil layer (kg/m3)
      data inttyp(iwcr)/1/  ! sm content at the critical point (mm)
      data inttyp(iwpw)/1/  ! sm content wilting point (mm)
      data inttyp(ismr)/1/  ! sm residual moisture (mm)
      data inttyp(islz)/1/  ! surface roughness of bare soil (m)
      data inttyp(isnz)/1/  ! surface roughness of snow pack (m)
      data inttyp(ilai)/1/  ! leaf area index
      data inttyp(ivg2)/1/  ! vic2 vegetation cover
      data inttyp(irt1)/1/  ! root content for layer 1
      data inttyp(irt2)/1/  ! root content for layer 2
      data inttyp(irt3)/1/  ! root content for layer 3
      data inttyp(isv1)/0/  ! soil moisture vic1
      data inttyp(itv1)/0/  ! soil temperature vic1
      data inttyp(isv2)/0/  ! soil moisture vic2
      data inttyp(itv2)/0/  ! soil temperature vic2
      data inttyp(iobp)/1/  ! observed precipitation
c
c  4. maxmum possible dimension size of the grib input file (including
c    climatology and analysis files (imdata,jmdata)
c
      integer imdata(numgrbs),jmdata(numgrbs)
      data imdata(igla),jmdata(igla)/ 180,  91/ ! glacier points
      data imdata(imxi),jmdata(imxi)/ 180,  91/ ! max ice extent
      data imdata(imsk),jmdata(imsk)/ 180,  91/ ! land sea mask (non grib)
      data imdata(iais),jmdata(iais)/ 721, 361/ ! sea ice mask/concentration
      data imdata(iscv),jmdata(iscv)/ 721, 361/ ! snow cover
      data imdata(isno),jmdata(isno)/ 721, 361/ ! snow depth
      data imdata(itsf),jmdata(itsf)/ 721, 361/ ! sst and land skin temp
      data imdata(iab4),jmdata(iab4)/ 361, 181/ ! 4-type albedo
      data imdata(iab1),jmdata(iab1)/ 361, 181/ ! 1-type albedo
      data imdata(ialf),jmdata(ialf)/ 361, 181/ ! albedo fraction
      data imdata(iso2),jmdata(iso2)/ 192,  94/ ! soil moisture osu2
      data imdata(isn1),jmdata(isn1)/ 361, 181/ ! soil moisture noa1
      data imdata(izor),jmdata(izor)/ 361, 181/ ! surface roughness
      data imdata(iplr),jmdata(iplr)/ 361, 181/ ! stomata resistance
      data imdata(itg3),jmdata(itg3)/ 138, 116/ ! deep layer soil temp
      data imdata(ito2),jmdata(ito2)/ 192,  94/ ! soil temp osu2
      data imdata(itn1),jmdata(itn1)/ 361, 181/ ! soil temp noa1
      data imdata(iveg),jmdata(iveg)/2500,1250/ ! vegetation cover
      data imdata(ivet),jmdata(ivet)/2500,1250/ ! vegetation type
      data imdata(isot),jmdata(isot)/2500,1250/ ! soil type
      data imdata(islo),jmdata(islo)/ 361, 181/ ! slope type
      data imdata(isna),jmdata(isna)/ 361, 181/ ! snow albedo
      data imdata(ismx),jmdata(ismx)/2500,1250/ ! max vegetation cover
      data imdata(ismn),jmdata(ismn)/2500,1250/ ! min vegetation cover
c-vic addition
      data imdata(ivgv),jmdata(ivgv)/ 720, 280/ ! vegetation cover
      data imdata(ivtv),jmdata(ivtv)/ 720, 280/ ! vegetation type
      data imdata(ivrt),jmdata(ivrt)/ 720, 280/ ! soil type
      data imdata(ibif),jmdata(ibif)/ 720, 280/ ! binf
      data imdata(ids ),jmdata(ids )/ 720, 280/ ! Fract of Dsm VICLSM1
      data imdata(idsm),jmdata(idsm)/ 720, 280/ ! Dsm VICLSM1
      data imdata(iws ),jmdata(iws )/ 720, 280/ ! Ws VICLSM1
      data imdata(icef),jmdata(icef)/ 720, 280/ ! c
      data imdata(iexp),jmdata(iexp)/ 720, 280/ ! expt VICLSM1
      data imdata(ikst),jmdata(ikst)/ 720, 280/ ! ksat VICLSM1
      data imdata(idph),jmdata(idph)/ 720, 280/ ! layer thickness
      data imdata(ibub),jmdata(ibub)/ 720, 280/ ! bub VICLSM1
      data imdata(iqrt),jmdata(iqrt)/ 720, 280/ ! quartz VICLSM1
      data imdata(ibkd),jmdata(ibkd)/ 720, 280/ ! bkd VICLSM1
      data imdata(isld),jmdata(isld)/ 720, 280/ ! sld VICLSM1
      data imdata(iwcr),jmdata(iwcr)/ 720, 280/ ! wcr VICLSM1
      data imdata(iwpw),jmdata(iwpw)/ 720, 280/ ! wpw VICLSM1
      data imdata(ismr),jmdata(ismr)/ 720, 280/ ! smr
      data imdata(islz),jmdata(islz)/ 720, 280/ ! silz
      data imdata(isnz),jmdata(isnz)/ 720, 280/ ! snwz
      data imdata(ilai),jmdata(ilai)/ 720, 280/ ! lai
      data imdata(ivg2),jmdata(ivg2)/ 720, 280/ ! vic2 vegetation cover
      data imdata(irt1),jmdata(irt1)/ 720, 280/ ! rt1
      data imdata(irt2),jmdata(irt2)/ 720, 280/ ! rt2
      data imdata(irt3),jmdata(irt3)/ 720, 280/ ! rt3
      data imdata(isv1),jmdata(isv1)/ 720, 280/ ! soil moisture vic1
      data imdata(itv1),jmdata(itv1)/ 720, 280/ ! soil temperature vic1
      data imdata(isv2),jmdata(isv2)/ 720, 280/ ! soil moisture vic2
      data imdata(itv2),jmdata(itv2)/ 720, 280/ ! soil temperature vic2
      data imdata(iobp),jmdata(iobp)/2500,1250/ ! observed precipitation
c
c  5. flag whether to use land sea mask for interpolation
c      mask=0  .. do not use mask
c          =1  .. use land mask
c          =2  .. use sea mask
c     if mask >0;
c       land sea mask definition of the input grib file grid (limask)
c       limask=.true.  .. use grib file field itself to define the mask
c              .false. .. use ready made 1x1 resolution land sea mask
c       if limask is .true.; 
c         additional variables, 'cvalin' and 'condin' 
c         need to be defined to obtain the land ocean mask.  cvalin is the 
c         criteria value and condin is the condition (lt,le,gt,ge,eq,ne)
c
      integer mask(numgrbs)
      logical limask(numgrbs)
      real cvalin(numgrbs)
      character*2 condin(numgrbs)
      data mask(igla)/0/ ! glacier points
      data mask(imxi)/0/ ! max ice extent
      data mask(imsk)/0/ ! land sea mask
      data mask(iais)/0/ ! sea ice mask/concentration
      data mask(iscv)/0/ ! snow cover
      data mask(isno)/0/ ! snow depth
      data mask(itsf)/1/ ! sst and land skin temp
        data limask(itsf)/.true./
          data cvalin(itsf),condin(itsf)/273.15,'lt'/
      data mask(iab4)/0/
      data mask(iab1)/0/ ! 1-type albedo
#ifndef VICLSM1
      data mask(ialf)/0/ ! albedo fraction
#endif
      data mask(iso2)/1/ ! soil moisture osu2
        data limask(iso2)/.true./
          data cvalin(iso2),condin(iso2)/0.469,'lt'/
      data mask(isn1)/0/ ! soil moisture noa1
      data mask(izor)/1/ ! surface roughness
        data limask(izor)/.true./
          data cvalin(izor),condin(izor)/9.9,'lt'/
      data mask(iplr)/0/ ! stomata resistance
      data mask(itg3)/0/ ! deep layer soil temp
      data mask(ito2)/0/ ! soil temp osu2
      data mask(itn1)/0/ ! soil temp noa1
      data mask(iveg)/0/ ! vegetation cover
      data mask(ivet)/0/ ! vegetation type
      data mask(isot)/0/ ! soil type
      data mask(islo)/0/ ! slope type
      data mask(isna)/0/ ! snow albedo
      data mask(ismx)/0/ ! max vegetation cover
      data mask(ismn)/0/ ! min vegetation cover
c-vic addition
      data mask(ivgv)/0/ ! vegetation cover vic1
      data mask(ivtv)/0/ ! vegetation type vic1
      data mask(ivrt)/1/ ! vegetation root
        data limask(ivrt)/.false./
#ifdef VICLSM1
      data mask(ialf)/1/ ! albedo fraction
        data limask(ialf)/.false./
#endif
      data mask(ibif)/1/ ! fract of Dsm
        data limask(ibif)/.false./
      data mask(ids )/1/ ! fract of Dsm
        data limask(ids )/.false./
      data mask(idsm)/1/ ! Dsm
        data limask(idsm)/.false./
      data mask(iws )/1/ ! Ws
        data limask(iws )/.false./
      data mask(icef)/1/ ! Ws
        data limask(icef)/.false./
      data mask(iexp)/1/ ! expt
        data limask(iexp)/.false./
      data mask(ikst)/1/ ! ksat
        data limask(ikst)/.false./
      data mask(idph)/1/ ! dp
        data limask(idph)/.false./
      data mask(ibub)/1/ ! bubble
        data limask(ibub)/.false./
      data mask(iqrt)/1/ ! quartz
        data limask(iqrt)/.false./
      data mask(ibkd)/1/ ! bulk density
        data limask(ibkd)/.false./
      data mask(isld)/1/ ! soil density
        data limask(isld)/.false./
      data mask(iwcr)/1/ ! wcr
        data limask(iwcr)/.false./
      data mask(iwpw)/1/ ! wpw
        data limask(iwpw)/.false./
      data mask(ismr)/1/ ! wpw
        data limask(ismr)/.false./
      data mask(islz)/1/ ! wpw
        data limask(islz)/.false./
      data mask(isnz)/1/ ! wpw
        data limask(isnz)/.false./
      data mask(ilai)/0/ ! lai
      data mask(ivg2)/0/ ! vic2 vegetation cover
      data mask(irt1)/1/ ! rt1
        data limask(irt1)/.false./
      data mask(irt2)/1/ ! rt2
        data limask(irt2)/.false./
      data mask(irt3)/1/ ! rt3
        data limask(irt3)/.false./
      data mask(isv1)/0/ ! soil moisture vic1
      data mask(itv1)/0/ ! soil tempeature vic1
      data mask(isv2)/0/ ! soil moisture vic2
      data mask(itv2)/0/ ! soil tempeature vic2
      data mask(iobp)/0/ ! observed precipitation
c
c  6. scaling of the result by multiplying scale
c
      real scale(numgrbs)
      data scale(igla)/  1./ ! glacier points
      data scale(imxi)/  1./ ! max ice extent
      data scale(imsk)/  1./ ! lad sea mask
      data scale(iais)/  1./ ! sea ice mask/concentration
      data scale(iscv)/  1./ ! snow cover
      data scale(isno)/  1./ ! snow depth
      data scale(itsf)/  1./ ! sst and land skin temp
      data scale(iab4)/0.01/ ! 4-type albedo
      data scale(iab1)/0.01/ ! 1-type albedo
      data scale(ialf)/0.01/ ! albedo fraction
      data scale(iso2)/  1./ ! soil moisture osu2
      data scale(isn1)/  1./ ! soil moisture noa1
      data scale(izor)/100./ ! surface roughness
      data scale(iplr)/  1./ ! stomata resistance
      data scale(itg3)/  1./ ! deep layer soil temp
      data scale(ito2)/  1./ ! soil temp osu2
      data scale(itn1)/  1./ ! soil temp noa1
#ifdef USGS
      data scale(iveg)/0.01/ ! vegetation cover
#else
      data scale(iveg)/  1./ ! vegetation cover
#endif
      data scale(ivet)/  1./ ! vegetation type
      data scale(isot)/  1./ ! soil type
      data scale(islo)/  1./ ! slope type
      data scale(isna)/0.01/ ! snow albedo
      data scale(ismx)/0.01/ ! max vegetation cover
      data scale(ismn)/0.01/ ! min vegetation cover
c- vic addition
      data scale(ivgv)/0.01/ ! vic vegetation cover
      data scale(ivtv)/  1./ ! vic vegetation type
      data scale(ivrt)/  1./ ! root content for one vegetation type
c     data scale(ialf)/0.01/ ! albedo fraction
      data scale(ibif)/  1./ ! binf
      data scale(ids )/  1./ ! Fract of Dsm
      data scale(idsm)/  1./ ! Maximum velocity of baseflow (mm/day)
      data scale(iws )/  1./ ! Fract maxi sm nonlinear baseflow occurs
      data scale(icef)/  1./ ! c
      data scale(iexp)/  1./ ! Para the vari of Ksat with sm (N/A)
      data scale(ikst)/  1./ ! Saturated hydrologic conductivity (mm/day)
      data scale(idph)/  1./ ! layer thickness (m)
      data scale(ibub)/  1./ ! Bubbling pressure of soil layer (cm)
      data scale(iqrt)/  1./ ! Quartz content of soil layer (fraction)
      data scale(ibkd)/  1./ ! Bulk density of soil layer (kg/m3)
      data scale(isld)/  1./ ! Soil density of soil layer (kg/m3)
      data scale(iwcr)/  1./ ! sm content at the critical point (mm)
      data scale(iwpw)/  1./ ! sm content wilting point (mm)
      data scale(ismr)/  1./ ! smr (mm)
      data scale(islz)/  1./ ! silz (m)
      data scale(isnz)/  1./ ! snwz (m)
      data scale(ilai)/  1./ ! lai
      data scale(ivg2)/  1./ ! vegetation cover
      data scale(irt1)/  1./ ! root content in layer 1
      data scale(irt2)/  1./ ! root content in layer 2
      data scale(irt3)/  1./ ! root content in layer 3
      data scale(isv1)/  1./ ! soil moisture vic1
      data scale(itv1)/  1./ ! soil temperature vic1
      data scale(isv2)/  1./ ! soil moisture vic2
      data scale(itv2)/  1./ ! soil temperature vic2
      data scale(iobp)/  1.0e-3/ ! observed precipitation (from mm to m)
c
c  7. flag for ouputting mask
c       lomask=.true.  .. output 1 or 0 mask
c              .false. .. regular output 
c       if lomask is .true.; 
c         additional variables, 'cvalout' and 'condout' 
c         need to be defined to obtain the output mask.  cvalout is the 
c         criteria value and condout is the condition (lt,le,gt,ge,eq,ne)
c
      logical lomask(numgrbs)
      real cvalout(numgrbs)
      character*2 condout(numgrbs)
      data lomask(igla)/.true. / ! glacier points
          data cvalout(igla),condout(igla)/0.5,'ge'/
      data lomask(imxi)/.true. / ! max ice extent
          data cvalout(imxi),condout(imxi)/0.5,'ge'/
      data lomask(imsk)/.true. / ! land sea mask (non grib)
          data cvalout(imsk),condout(imsk)/0.5,'ge'/
      data lomask(iais)/.true. / ! sea ice mask/concentration
          data cvalout(iais),condout(iais)/0.5,'ge'/
      data lomask(iscv)/.true. / ! snow cover
          data cvalout(iscv),condout(iscv)/0.5,'ge'/
      data lomask(isno)/.false./ ! snow depth
      data lomask(itsf)/.false./ ! sst and land skin temp
      data lomask(iab4)/.false./ ! 4-type albedo
      data lomask(iab1)/.false./ ! 1-type albedo
      data lomask(ialf)/.false./ ! albedo fraction
      data lomask(iso2)/.false./ ! soil moisture osu2
      data lomask(isn1)/.false./ ! soil moisture noa1
      data lomask(izor)/.false./ ! surface roughness
      data lomask(iplr)/.false./ ! stomata resistance
      data lomask(itg3)/.false./ ! deep layer soil temp
      data lomask(ito2)/.false./ ! soil temp osu2
      data lomask(itn1)/.false./ ! soil temp noa1
      data lomask(iveg)/.false./ ! vegetation cover
      data lomask(ivet)/.false./ ! vegetation type
      data lomask(isot)/.false./ ! soil type
      data lomask(islo)/.false./ ! slope type
      data lomask(isna)/.false./ ! snow albedo
      data lomask(ismx)/.false./ ! max vegetation cover
      data lomask(ismn)/.false./ ! min vegetation cover
c-vic addition
      data lomask(ivgv)/.false./ ! vegetation cover vic
      data lomask(ivtv)/.false./ ! vegetation type vic
      data lomask(ivrt)/.false./ ! root content for one vegetation type
c     data lomask(ialf)/.false./ ! albedo fraction
      data lomask(ibif)/.false./ ! binf
      data lomask(ids )/.false./ ! Fract of Dsm nonlinear baseflow begins
      data lomask(idsm)/.false./ ! Maximum velocity of baseflow (mm/day)
      data lomask(iws )/.false./ ! Fract maxi sm nonlinear baseflow occurs
      data lomask(icef)/.false./ ! c
      data lomask(iexp)/.false./ ! Para the vari of Ksat with sm (N/A)
      data lomask(ikst)/.false./ ! Saturated hydrologic conductivity (mm/day)
      data lomask(idph)/.false./ ! soil layer thickness (m)
      data lomask(ibub)/.false./ ! Bubbling pressure of soil layer (cm)
      data lomask(iqrt)/.false./ ! Quartz content of soil layer (fraction)
      data lomask(ibkd)/.false./ ! Bulk density of soil layer (kg/m3)
      data lomask(isld)/.false./ ! Soil density of soil layer (kg/m3)
      data lomask(iwcr)/.false./ ! sm content at the critical point (mm)
      data lomask(iwpw)/.false./ ! sm content wilting point (mm)
      data lomask(ismr)/.false./ ! smr (mm)
      data lomask(islz)/.false./ ! silz (m)
      data lomask(isnz)/.false./ ! snwz (m)
      data lomask(ilai)/.false./ ! lai
      data lomask(ivg2)/.false./ ! vegetation cover
      data lomask(irt1)/.false./ ! root content for layer 1
      data lomask(irt2)/.false./ ! root content for layer 2
      data lomask(irt3)/.false./ ! root content for layer 3
      data lomask(isv1)/.false./ ! soil moisture vic1
      data lomask(itv1)/.false./ ! soil temperature vic1
      data lomask(isv2)/.false./ ! soil moisture vic2
      data lomask(itv2)/.false./ ! soil temperature vic2
      data lomask(iobp)/.false./ ! observed precipitation
c
c  8.  flag for using the grib field over different land surface
c      characteristics.
c
c     lflgrb(1,:) .. open ocean 
c     lflgrb(2,:) .. land without snow
c     lflgrb(3,:) .. sea ice without snow
c     lflgrb(4,:) .. land with snow cover
c     lflgrb(5,:) .. sea ice with snow cover
c  
c   =.t. overwrites forecast sfc rec
c   =.f. keep forecast sfc rec
c     
      logical lflgrb(5,numgrbs)
      data (lflgrb(igb,igla),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! glacier points (nver mrgd)
      data (lflgrb(igb,imxi),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! maxice (nver mrgd)
      data (lflgrb(igb,imsk),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! land sea mask (nver mrgd)
      data (lflgrb(igb,iais),igb=1,5)/
     1    .true.,.false.,.true.,.false.,.true./ ! sea ice mask/conc.
      data (lflgrb(igb,iscv),igb=1,5)/
     1    .false.,.true.,.true.,.true.,.true./ ! snow cover
      data (lflgrb(igb,isno),igb=1,5)/
     1    .false.,.true.,.true.,.true.,.true./ ! snow depth
      data (lflgrb(igb,itsf),igb=1,5)/
     1    .true.,.false.,.false.,.false.,.false./ ! sst and skin temp
      data (lflgrb(igb,iab4),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! 4-type albedo
      data (lflgrb(igb,iab1),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! 1-type albedo
      data (lflgrb(igb,ialf),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! albedo fraction
      data (lflgrb(igb,iso2),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soilm osu2
      data (lflgrb(igb,isn1),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soilm noa1
      data (lflgrb(igb,izor),igb=1,5)/
     1    .false.,.true.,.true.,.true.,.true./ ! sfc roughness
      data (lflgrb(igb,iplr),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! stomato resistance
      data (lflgrb(igb,itg3),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! tg3
      data (lflgrb(igb,ito2),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soil temp osu2
      data (lflgrb(igb,itn1),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soil temp noa1
      data (lflgrb(igb,iveg),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! veg cover
      data (lflgrb(igb,ivet),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! veg type
      data (lflgrb(igb,isot),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! soil type
      data (lflgrb(igb,islo),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! slope type
      data (lflgrb(igb,isna),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! snow albedo
      data (lflgrb(igb,ismx),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! max veg cov
      data (lflgrb(igb,ismn),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! min veg cov
c-vic addition
      data (lflgrb(igb,ivtv),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! vic veg cover
      data (lflgrb(igb,ivgv),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! vic veg type
      data (lflgrb(igb,ivrt),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! vrt
c     data (lflgrb(igb,ialf),igb=1,5)/
c    1    .true.,.true.,.true.,.true.,.true./ ! albedo fraction
      data (lflgrb(igb,ibif),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! binf
      data (lflgrb(igb,ids ),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! Fract of Dsm
      data (lflgrb(igb,idsm),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! Dsm
      data (lflgrb(igb,iws ),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! Ws
      data (lflgrb(igb,icef),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! c
      data (lflgrb(igb,iexp),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! expt
      data (lflgrb(igb,ikst),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! ksat
      data (lflgrb(igb,idph),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! dph
      data (lflgrb(igb,ibub),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! bubble
      data (lflgrb(igb,iqrt),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! quartz
      data (lflgrb(igb,ibkd),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! bkd
      data (lflgrb(igb,isld),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! sld
      data (lflgrb(igb,iwcr),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! wcr
      data (lflgrb(igb,iwpw),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! wpw
      data (lflgrb(igb,ismr),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! smr
      data (lflgrb(igb,islz),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! silz
      data (lflgrb(igb,isnz),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! snwz
      data (lflgrb(igb,ilai),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! lai
      data (lflgrb(igb,ivg2),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! veg cover
      data (lflgrb(igb,irt1),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! rt1
      data (lflgrb(igb,irt2),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! rt2
      data (lflgrb(igb,irt3),igb=1,5)/
     1    .true.,.true.,.true.,.true.,.true./ ! rt3
      data (lflgrb(igb,isv1),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soilt vic1
      data (lflgrb(igb,itv1),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soilt vic1
      data (lflgrb(igb,isv2),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soilt vic2
      data (lflgrb(igb,itv2),igb=1,5)/
     1    .false.,.true.,.false.,.true.,.false./ ! soilt vic2
      data (lflgrb(igb,iobp),igb=1,5)/
     1    .false.,.false.,.false.,.false.,.false./ ! observed precipitation
c
c use observed for soil moisture adjustment only,  not to supplant model
c
c
c  9.  full variable name (for print out purpose)
c
      character*8 gvar(numgrbs)
      data gvar(igla)/'glacier '/
      data gvar(imxi)/'maxisice'/
      data gvar(imsk)/'lndseamk'/
      data gvar(iais)/'seaice  '/
      data gvar(iscv)/'snowcovr'/
      data gvar(isno)/'snowdept'/
      data gvar(itsf)/'sstskint'/
      data gvar(iab4)/'4-typalb'/
      data gvar(iab1)/'1-typalb'/
      data gvar(ialf)/'albdfrac'/
      data gvar(iso2)/'soilosu2'/
      data gvar(isn1)/'soilnoa1'/
      data gvar(izor)/'srfrough'/
      data gvar(iplr)/'stmaeres'/
      data gvar(itg3)/'depsoilt'/
      data gvar(ito2)/'soltosu2'/
      data gvar(itn1)/'soltnoa1'/
      data gvar(iveg)/'vegcover'/
      data gvar(ivet)/'veg type'/
      data gvar(isot)/'soiltype'/
      data gvar(islo)/'sloptype'/
      data gvar(isna)/'snwalbdo'/
      data gvar(ismx)/'mxvegcov'/
      data gvar(ismn)/'mnvegcov'/
c-vic addtion
      data gvar(ivgv)/'vegcvic1'/
      data gvar(ivtv)/'vegtvic1'/
      data gvar(ivrt)/'vegroot '/
c     data gvar(ialf)/'albdfrac'/
			data gvar(ibif)/'binfvic '/
      data gvar(ids )/'dsviclsm'/
      data gvar(idsm)/'dsmvic  '/
      data gvar(iws )/'wsviclsm'/
      data gvar(icef)/'cefvic  '/
      data gvar(iexp)/'exptvic '/
      data gvar(ikst)/'ksatvic '/
      data gvar(idph)/'dphvic  '/
      data gvar(ibub)/'bubblepr'/
      data gvar(iqrt)/'quartz  '/
      data gvar(ibkd)/'bulkdens'/
      data gvar(isld)/'soildens'/
      data gvar(iwcr)/'wcrvic  '/
      data gvar(iwpw)/'wpwvic  '/
      data gvar(ismr)/'smrvic  '/
      data gvar(islz)/'silzvic '/
      data gvar(isnz)/'snwzvic '/
      data gvar(ilai)/'laivic  '/
      data gvar(ivg2)/'vgcvvic2'/
      data gvar(irt1)/'vegroot1'/
      data gvar(irt2)/'vegroot2'/
      data gvar(irt3)/'vegroot3'/
      data gvar(isv1)/'solmvic1'/
      data gvar(itv1)/'soltvic1'/
      data gvar(isv2)/'solmvic2'/
      data gvar(itv2)/'soltvic2'/
      data gvar(iobp)/'obsprecp'/
c
c  end definition of surface variable properties

