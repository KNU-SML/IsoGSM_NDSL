      subroutine idsdef(iptv,ids)                                               
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram: idsdef         sets default decimal scalings                      
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: sets decimal scalings defaults for various parameters.              
c   a decimal scaling of -3 means data is packed in kilo-si units.              
c                                                                               
c program history log:                                                          
c   92-10-31  iredell                                                           
c                                                                               
c usage:    call idsdef(iptv,ids)                                               
c   input arguments:                                                            
c     iptv         paramter table version (only 1 or 2 is recognized)           
c   output arguments:                                                           
c     ids          integer (255) decimal scalings                               
c                  (unknown decimal scalings will not be set)                   
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
#include <define.h>
      dimension ids(255)                                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do k=1,255
        ids(k)=0
      enddo

      if(iptv.eq.1.or.iptv.eq.2) then                                           
        ids(001)=0      ! pressure (pa)                                         
        ids(002)=0      ! sea-level pressure (pa)                               
        ids(003)=4      ! pressure tendency (pa/s)                              
                        !                                                       
                        !                                                       
        ids(006)=0      ! geopotential (m2/s2)  
        ids(007)=1      ! geopotential height (m)                               
        ids(008)=1      ! geometric height (m)                                  
        ids(009)=1      ! standard deviation of height (m)                      
                        !                                                       
        ids(011)=2      ! temperature (k)                                       
        ids(012)=2      ! virtual temperature (k)                               
        ids(013)=2      ! potential temperature (k)                             
        ids(014)=2      ! pseudo-adiabatic potential temperature (k)            
        ids(015)=2      ! maximum temperature (k)                               
        ids(016)=2      ! minimum temperature (k)                               
        ids(017)=2      ! dewpoint temperature (k)                              
        ids(018)=2      ! dewpoint depression (k)                               
        ids(019)=5      ! temperature lapse rate (k/m)                          
        ids(020)=0      ! visibility (m)                                        
        ids(021)=5      ! river storage (log10,mm)
        ids(022)=5      ! river flow (log10,kg/s)
        ids(023)=0      ! river direction (0-9)
        ids(024)=9      ! total runoff to river channel (kg/m2/s)
        ids(025)=6      ! river storage isotope ratio 1 (ratio)
        ids(026)=6      ! river flow isotope ratio 1 (ratio)
        ids(027)=9      ! total runoff to river channel (kg/m2/s) 1
        ids(028)=6      ! river storage isotope ratio 2 (ratio)
        ids(029)=6      ! river flow isotope ratio 2 (ratio)
        ids(030)=9      ! total runoff to river channel (kg/m2/s) 2
        ids(031)=1      ! wind direction (degrees)                              
        ids(032)=2      ! wind speed (m/s)                                      
        ids(033)=2      ! zonal wind (m/s)                                      
        ids(034)=2      ! meridional wind (m/s)                                 
        ids(035)=-4     ! streamfunction (m2/s)                                 
        ids(036)=-4     ! velocity potential (m2/s)                             
        ids(037)=0      ! montgomery stream function (m2/s2)                    
        ids(038)=9      ! sigma vertical velocity (1/s)                         
        ids(039)=2      ! pressure vertical velocity (pa/s)                     
        ids(040)=5      ! geometric vertical velocity (m/s)                     
        ids(041)=7      ! absolute vorticity (1/s)                              
        ids(042)=7      ! absolute divergence (1/s)                             
        ids(043)=7      ! relative vorticity (1/s)                              
        ids(044)=7      ! relative divergence (1/s)                             
        ids(045)=5      ! vertical u shear (1/s)                                
        ids(046)=5      ! vertical v shear (1/s)                                
        ids(047)=0      ! direction of current (degrees)                        
                        ! speed of current (m/s)                                
                        ! u of current (m/s)                                    
                        ! v of current (m/s)                                    
        ids(051)=9      ! specific humidity (kg/kg)                             
        ids(052)=3      ! relative humidity (percent)                           
        ids(053)=9      ! humidity mixing ratio (kg/kg)                         
        ids(054)=5      ! precipitable water (kg/m2)                            
        ids(055)=0      ! vapor pressure (pa)                                   
        ids(056)=0      ! saturation deficit (pa)                               
        ids(057)=9      ! evaporation (kg/m2/s)                                   
        ids(058)=2      ! cloud ice (kg/m2)                                     
        ids(059)=9      ! precipitation rate (kg/m2/s)                          
        ids(060)=0      ! thunderstorm probability (percent)                    
        ids(061)=2      ! total precipitation (kg/m2)                           
        ids(062)=2      ! large-scale precipitation (kg/m2)                     
        ids(063)=2      ! convective precipitation (kg/m2)                      
#ifdef NOAHYDRO
        ids(064)=8      ! water equivalent snowfall rate (kg/m2/s)              
#ifdef LINUX
        ids(065)=2      ! water equivalent of snow depth (kg/m2)                
#else
        ids(065)=4      ! water equivalent of snow depth (kg/m2)                
#endif
#else
        ids(064)=7      ! water equivalent snowfall rate (kg/m2/s)              
        ids(065)=1      ! water equivalent of snow depth (kg/m2)                
#endif
        ids(066)=3      ! snow depth (m)                                        
#ifdef ISOTOPE
        ids(067)=9      ! precipitation rate for isotope1 (kg/m2/s*SMOW)        
        ids(068)=9      ! precipitation rate for isotope2 (kg/m2/s*SMOW)        
        ids(069)=9      ! conv-precip rate for isotope1 (kg/m2/s*SMOW)          
        ids(070)=9      ! conv-precip rate for isotope2 (kg/m2/s*SMOW)     
#else
                        ! mixed-layer depth (m)                                 
                        ! transient thermocline depth (m)                       
                        ! main thermocline depth (m)                            
                        ! main thermocline anomaly (m)                          
#endif
        ids(071)=1      ! total cloud cover (percent)                           
        ids(072)=1      ! convective cloud cover (percent)                      
        ids(073)=1      ! low cloud cover (percent)                             
        ids(074)=1      ! middle cloud cover (percent)                          
        ids(075)=1      ! high cloud cover (percent)                            
        ids(076)=5      ! cloud water (kg/m2)                                   
        ids(077)=5      ! cloud rain  (kg/m2)  
        ids(078)=2      ! convective snow (kg/m2)                               
        ids(079)=2      ! large scale snow (kg/m2)                              
        ids(080)=2      ! water temperature (k)                                 
        ids(081)=0      ! sea-land mask ()                                      
                        ! deviation of sea level from mean (m)                  
        ids(083)=6      ! roughness (m)                                         
        ids(084)=2      ! albedo (percent)                                      
        ids(085)=2      ! soil temperature (k)                                  
        ids(086)=1      ! soil wetness (kg/m2)                                  
        ids(087)=1      ! vegetation (percent)                                  
                        ! salinity (kg/kg)                                      
        ids(089)=5      ! density (kg/m3)                                       
        ids(090)=9      ! runoff rate (kg/m2/s)
        ids(091)=1      ! ice concentration ()                                  
        ids(092)=9      ! runoff rate isotope 1
        ids(093)=9      ! runoff rate isotope 2
#ifdef ISOTOPE
        ids(094)=5      ! surface reservoir isotope ratio 1 (ratio)
        ids(095)=5      ! surface reservoir isotope ratio 2 (ratio)
#else
                        ! speed of ice drift (m/s)                              
                        ! u of ice drift (m/s)                                  
#endif
        ids(096)=6      ! 24 hour accumulated precip
        ids(097)=2      ! 24 hour max temp.
        ids(098)=2      ! 24 hour min temp.
#ifndef NOAHYDRO
        ids(099)=7      ! snow melt rate (kg/m2/s)
#else
        ids(099)=8      ! snow melt rate (kg/m2/s)
#endif
#ifndef NFDR
                        ! sig height of waves and swell (m)                     
        ids(101)=0      ! direction of wind waves (degrees)                     
                        ! sig height of wind waves (m)                          
                        ! mean period of wind waves (s)                         
        ids(104)=0      ! direction of swell waves (degrees)                    
                        ! sig height of swell waves (m)                         
                        ! mean period of swell waves (s)                        
        ids(107)=0      ! primary wave direction (degrees)                      
                        ! primary wave mean period (s)                          
        ids(109)=0      ! secondary wave direction (degrees)                    
                        ! secondary wave mean period (s)                        
#else
        ids(100)=2      ! 2-meter temperature at 2pm local (degrees)
        ids(101)=1      ! 2-meter relative humidity at 2pm
        ids(102)=1      ! total cloudiness at 2pm
        ids(103)=2      ! 10-meter wind speed at 2pm
        ids(104)=0      ! Precipitation duration (sec)
        ids(105)=1      ! Maximum relative humidity
        ids(106)=1      ! Minimum relative humidity
        ids(107)=0      ! Fire Danger Spread Component
        ids(108)=0      ! Fire Danger Burning Index
        ids(109)=0      ! Fire Danger Energy Release
        ids(110)=0      ! Fire Danger Ignition Component
#endif
#ifdef ISOTOPE
        ids(111)=9      ! total condensation (kg/m2/s)
        ids(112)=9      ! total condensation isotope 1 (kg/m2/s*SMOW)
        ids(113)=9      ! total condensation isotope 2 (kg/m2/s*SMOW)
        ids(114)=9      ! convective condensation (kg/m2/s)
        ids(115)=9      ! convective condensation isotope 1 (kg/m2/s*SMOW)
        ids(116)=9      ! convective condensation isotope 2 (kg/m2/s*SMOW)
#else
        ids(111)=1      ! net solar radiative flux at surface (w/m2)            
        ids(112)=1      ! net longwave radiative flux at surface (w/m2)         
        ids(113)=1      ! net solar radiative flux at top (w/m2)                
        ids(114)=1      ! net longwave radiative flux at top (w/m2)             
        ids(115)=1      ! net longwave radiative flux (w/m2)                    
        ids(116)=1      ! net solar radiative flux (w/m2)                       
#endif
        ids(117)=1      ! total radiative flux (w/m2)                           
                        !                                                       
#ifdef CO2PROG
        ids(119)=9      ! co2 flux (kgCO2/m2/s)
#endif
#ifdef ISOTOPE
        ids(119)=5      ! latent heat flux isotope 1 (w/m2)
        ids(120)=5      ! latent heat flux isotope 2 (w/m2)
#endif          
        ids(121)=5      ! latent heat flux (w/m2)                               
        ids(122)=2      ! sensible heat flux (w/m2)                             
        ids(123)=1      ! boundary layer dissipation (w/m2)                     
        ids(124)=4      ! u wind stress (n/m2)                                  
        ids(125)=4      ! v wind stress (n/m2)                                  
                        ! wind mixing energy (j)                                
#ifndef NFDR
                        ! image data ()                                         
#else
        ids(127)=0      ! Keeth-Bryan Drought index
#endif
        ids(128)=0      ! mean sea-level pressure (stdatm) (pa)                 
        ids(129)=0      ! mean sea-level pressure (maps) (pa)                   
        ids(130)=0      ! mean sea-level pressure (eta) (pa)                    
        ids(131)=2      ! surface lifted index (k)                              
        ids(132)=2      ! best lifted index (k)                                 
        ids(133)=2      ! k index (k)                                           
        ids(134)=2      ! sweat index (k)                                       
        ids(135)=11     ! horizontal moisture divergence (kg/kg/s)              
        ids(136)=5      ! speed shear (1/s)                                     
        ids(137)=6      ! 3-hr pressure tendency (pa/s)                         
        ids(138)=7      ! brunt-vaisala frequency squared (1/s2)                
        ids(139)=12     ! potential vorticity (mass-weighted) (1/s/m)           
        ids(140)=0      ! rain mask ()                                          
        ids(141)=0      ! freezing rain mask ()                                 
        ids(142)=0      ! ice pellets mask ()                                   
        ids(143)=0      ! snow mask ()                                          
        ids(144)=8      ! volumetric soil moisture content (fraction)           
        ids(145)=1      ! potential evaporation rate (w/m2)                     
        ids(146)=1      ! cloud workfunction (j/kg)                             
        ids(147)=4      ! u gravity wave stress (n/m2)                          
        ids(148)=4      ! v gravity wave stress (n/m2)                          
        ids(149)=11     ! potential vorticity (m2/s/kg)                         
                        ! covariance between v and u (m2/s2)                    
                        ! covariance between u and t (k*m/s)                    
                        ! covariance between v and t (k*m/s)                    
#ifdef ISOTOPE
        ids(150)=5      ! precipitable water isotope 1 (kg/m/m)                 
        ids(151)=5      ! precipitable water isotope 2 (kg/m/m)                 
        ids(152)=9      ! precipitatation water (kg/kg)
        ids(153)=9      ! cloud water water (kg/kg)
#else
                        !                                                       
                        !                                                       
        ids(152)=6      ! precipitatation water (kg/kg)
        ids(153)=6      ! cloud water water (kg/kg)
#endif
#ifdef NOAHYDRO
        ids(154)=4      ! soil type (integer)
#else
        ids(154)=2      ! ozone mixing ratio (kg/kg)
#endif
        ids(155)=2      ! ground heat flux (w/m2)                               
        ids(156)=1      ! convective inhibition (w/m2)                          
        ids(157)=1      ! convective ape (j/kg)                                 
        ids(158)=1      ! turbulent ke (j/kg)                                   
        ids(159)=0      ! condensation pressure of lifted parcel (pa)           
        ids(160)=1      ! clear sky upward solar flux (w/m2)                    
        ids(161)=1      ! clear sky downward solar flux (w/m2)                  
        ids(162)=1      ! clear sky upward longwave flux (w/m2)                 
        ids(163)=1      ! clear sky downward longwave flux (w/m2)               
        ids(164)=1      ! cloud forcing net solar flux (w/m2)                   
        ids(165)=1      ! cloud forcing net longwave flux (w/m2)                
        ids(166)=1      ! visible beam downward solar flux (w/m2)               
        ids(167)=1      ! visible diffuse downward solar flux (w/m2)            
        ids(168)=1      ! near ir beam downward solar flux (w/m2)               
        ids(169)=1      ! near ir diffuse downward solar flux (w/m2)            
                        !                                                       
                        !                                                       
        ids(172)=4      ! momentum flux (n/m2)                                  
        ids(173)=0      ! mass point model surface ()                           
        ids(174)=0      ! velocity point model surface ()                       
        ids(175)=0      ! sigma layer number ()                                 
        ids(176)=2      ! latitude (degrees)                                    
        ids(177)=2      ! east longitude (degrees)                              
                        !                                                       
                        !                                                       
#ifndef NOAHYDRO
        ids(180)=4      ! canopy re-evaporation (w/m2)
#else
        ids(180)=3      ! canopy re-evaporation (w/m2)
#endif
        ids(181)=10     ! x-gradient log pressure (1/m)                         
        ids(182)=10     ! y-gradient log pressure (1/m)                         
        ids(183)=6      ! x-gradient height (m/m)                               
        ids(184)=6      ! y-gradient height (m/m)                               
        ids(185)=10     ! flux divergence of temperature (K/s)
        ids(186)=11     ! flux divergence of specific humidity (/s)
        ids(187)=0      ! vertically integrated enthalpy (g/s2)
        ids(188)=0      ! vertically integrated latent energy (g/s2)
        ids(189)=8      ! vertically integrated enthalpy flux div (g/s3)
        ids(190)=8      ! vertically integrated latent energy flux div (g/s3)
        ids(191)=9      ! x-derivative of temperature (K/m)
        ids(192)=9      ! y-derivative of temperature (K/m)
        ids(193)=12     ! x-derivative of specific humidity (1/m)
        ids(194)=12     ! y-derivative of specific humidity (1/m)
        ids(195)=6      ! x-derivative of geopotential (non-dim)
        ids(196)=6      ! y-derivative of geopotential (non-dim)
        ids(197)=3      ! flux divergence of height (1/m)
        ids(198)=0      ! vertically ingegrated gz flux div
                        !
                        !
        ids(201)=0      ! ice-free water surcace (percent)                      
        ids(202)=2      ! Integrated moisture u-flux (m/s)
        ids(203)=2      ! Integrated moisture v-flux (m/s)
        ids(204)=1      ! downward solar radiative flux (w/m2)                  
        ids(205)=1      ! downward longwave radiative flux (w/m2)               
                        !                                                       
        ids(207)=0      ! moisture availability (percent)                       
                        ! exchange coefficient (kg/m2/s)                        
        ids(209)=0      ! fire weather index(non-dim)
                        !                                                       
        ids(211)=1      ! upward solar radiative flux (w/m2)                    
        ids(212)=1      ! upward longwave radiative flux (w/m2)                 
        ids(213)=0      ! non-convective cloud cover (percent)                  
        ids(214)=9      ! convective precipitation rate (kg/m2/s)               
        ids(215)=8      ! total diabatic heating rate (k/s)                     
        ids(216)=8      ! total radiative heating rate (k/s)                    
        ids(217)=8      ! total diabatic nonradiative heating rate (k/s)        
        ids(218)=2      ! precipitation index (fraction)                        
        ids(219)=2      ! std dev of ir t over 1x1 deg area (k)                 
        ids(220)=5      ! natural log of surface pressure over 1 kpa ()         
        ids(221)=1      ! planetary boundary layer height (m)                   
        ids(222)=1      ! 5-wave geopotential height (m)                        
#ifndef NOAHYDRO
        ids(223)=6      ! plant canopy surface water (kg/m2)                    
#else
        ids(223)=4      ! plant canopy surface water (kg/m2)                    
#endif
        ids(224)=8      ! baseflow-groundwater runoff rate (kg/m2/s)
        ids(225)=0      ! vegetation type (no dimension) 
                        ! storm surface runoff rate (kg/m2/s)   
        ids(226)=9      ! base flow isotope 1
        ids(227)=9      ! base flow isotope 2
                        ! asymptotic mixing length (m)                          
        ids(228)=2      ! potential evaporation (kg/m2)                         
        ids(229)=1      ! snow phase-change heat flux (w/m2)                    
#ifndef NOAHYDRO
        ids(230)=1      ! snow sublimation heat flux (w/m2)
#else
        ids(230)=4      ! snow sublimation heat flux (w/m2)
#endif
        ids(231)=4      ! convective cloud mass flux (pa/s)                     
        ids(232)=1      ! downward total radiation flux (w/m2)                  
        ids(233)=1      ! upward total radiation flux (w/m2)                    
        ids(234)=9      ! base flow run off rate (kg/m2/s)
                        !
#ifndef NOAHYDRO
        ids(236)=1      ! another latent heat flux (w/m2)                          
#else
        ids(236)=4      ! another latent heat flux (w/m2)                          
#endif
                        !
        ids(238)=1      ! snow cover (percent)                                  
        ids(239)=2      ! snow temperature (k)                                  
                        !                                                       
        ids(241)=8      ! large scale condensation heating rate (k/s)           
        ids(242)=8      ! deep convective heating rate (k/s)                    
        ids(243)=11     ! deep convective moistening rate (kg/kg/s)             
        ids(244)=8      ! shallow convective heating rate (k/s)                 
        ids(245)=11     ! shallow convective moistening rate (kg/kg/s)          
        ids(246)=8      ! vertical diffusion heating rate (kg/kg/s)             
        ids(247)=8      ! vertical diffusion zonal acceleration (m/s/s)         
        ids(248)=8      ! vertical diffusion merid acceleration (m/s/s)         
        ids(249)=11     ! vertical diffusion moistening rate (kg/kg/s)          
        ids(250)=8      ! solar radiative heating rate (k/s)                    
        ids(251)=8      ! longwave radiative heating rate (k/s)                 
                        ! drag coefficient () 
                        ! friction velocity (m/s)                               
                        ! richardson number ()                                  
#ifndef NOAHYDRO
                        !                           
#else
        ids(255)=8      ! additional precipitation rate (kg/m2/s)   
#endif
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
