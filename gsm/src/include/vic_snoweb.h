* ------------------------ code history --------------------------------
* source file:       vic_snoweb.h
* purpose:           variable common block for snow energy balance
* PRGMMR:            Ji Chen
* ORG:               ECPC/CRD/SIO/UCSD 
* DATE:              November 2003
* 
* Notes:
*    Upcase "E" before variables means computing Energy balance
* ----------------------------------------------------------------------

* ------------------- input variables ----------------------------------
 
* -- model basic parameter
      real Edtime        ! time step (second)
      integer Emonth     ! month of current time step

* -- atmosphere data

      real EfrTF         ! Rain through fall (m/time step)
      real ETair         ! Air temperature (c)
      real Erhoair       ! air density  (kg/m3)
      real Evpair        ! Actual vapor pressure of air (Pa)
      real Evpd          ! Vapor pressure deficit (Pa)
      real ERa           ! aerodynamic resistance (s/m)
      real Ewind         ! Wind speed (m/s)
      real Epgcm         ! Air pressure (Pa)

      real Enetsht       ! net shortwave radiation (w/m2)
      real Elongw        ! Longwave radiation (w/m2)

* -- vegetation parameters

      integer Evtype     ! vegetation type

* -- land surface parameters

      real EoldTsf       ! Snow Surface temp for previous time (c)
      real ETgrnd        ! ground surface temperature (c)
      real Efrsn         ! snow density (kg/m3)
      real Esnowd        ! snow depth (m)
      real Esfswq        ! Snow water equivalent on snow surface (m)
      real Esfwater      ! Liquid water in snow surface (m)
* ----------------------------------------------------------------------

      common /vicsnoweb_inp/Edtime,  Emonth,    EfrTF,  ETair,
     2                     Erhoair,  Evpair,     Evpd,    ERa,  
     3                       Ewind,   Epgcm,  Enetsht, Elongw,
     4                      Evtype, EoldTsf,   ETgrnd,  Efrsn,
     5                      Esnowd,  Esfswq, Esfwater

* ------------------- output variables ---------------------------------
      real Orefreeze     ! Refreeze energy (W/m2)
      real OlatHs        ! Latent heat exchange at surface (W/m2)
      real OsenHs        ! Sensible heat exchange at surface (W/m2)
* ----------------------------------------------------------------------

      common /vicsnoweb_out/Orefreeze, OlatHs, OsenHs
