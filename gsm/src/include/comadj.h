c
c $Id: comadj.h,v 1.3 2004/05/21 23:38:52 kana Exp $
c $Author: kana $
c
C
C Convective adjustment
C
      common/comadj/cappa   ,cpair   ,epsilo  ,gravit  ,latvap  ,
     $              rhoh2o  ,cldcp   ,clrh2o  ,latice  ,rair    ,
     $              rh2o    ,zvir    ,t0      ,nlvdry
C
      real cappa            ! R/cp
      real cpair            ! Specific heat of dry air
      real epsilo           ! Ratio of h2o to dry air molecular weights 
      real gravit           ! Gravitational acceleration 
      real latvap           ! Latent heat of vaporization
      real rhoh2o           ! Density of liquid water (STP)
      real cldcp            ! Latvap/cpair (L/cp)
      real clrh2o           ! Ratio of latvap to water vapor gas const
      real latice      ! Latent heat of fusion
      real rair        ! Gas constant for dry air
      real rh2o        ! Gas constant for water vapor
      real zvir        ! rh2o/rair - 1
      real t0          ! Reference temperature for t-prime computations
C
      integer nlvdry        ! Number of levels to apply dry adjustment

 
