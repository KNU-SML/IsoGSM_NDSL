#define NSTYPE nstype
c
c  grddf sets up moisture diffusivity and hydrolic conductivity
c  for all soil types
c  grddfs sets up thermal diffusivity for all soil types
c
      block data dfkt
#include <paramodel.h>
      parameter(ntype=nstype_)
      common /comgdfkt/ b(ntype),satpsi(ntype),satkt(ntype),tsat(ntype)
c
c  the nine soil types are:
c    1  ... loamy sand (coarse)
c    2  ... silty clay loam (medium)
c    3  ... light clay (fine)
c    4  ... sandy loam (coarse-medium)
c    5  ... sandy clay (coarse-fine)
c    6  ... clay loam  (medium-fine)
c    7  ... sandy clay loam (coarse-med-fine)
c    8  ... loam  (organic)
c    9  ... ice (use loamy sand property)
c
c    nstype = 16 from stagsgo soil data
c *************************************************************************
c
c     soil types    statsgo (miller ??, 199?)  cosby et al (1984)
c             1          sand                  sand
c             2          loamy sand            loamy sand
c             3          sandy loam            sandy loam
c             4          silt loam             silty loam
c             5          silt                  silty loam
c             6          loam                  loam
c             7          sandy clay loam       sandy clay loam
c             8          silty clay loam       silty clay loam
c             9          clay loam             clay loam
c            10          sandy clay            sandy clay
c            11          silty clay            silty clay
c            12          clay                  light clay
c            13          organic materials     loam
c            14          water
c            15          bedrock
c            16          other (land-ice)
c                        the value of this class is the same as in classe-2
c
#ifndef OSULSM1
#ifdef STATSGO_SOIL
      data b    /2.79,  4.26,  4.74,  5.33,  5.33,  5.25,  6.66,
     &  8.72,  8.17,  10.73, 10.39, 11.55, 5.25,  0.0, 2.79,  4.26/
      data satpsi/0.069, 0.036, 0.141, 0.759, 0.759, 0.355, 0.135,
     &  0.617, 0.263, 0.098, 0.324, 0.468, 0.355, 0.0, 0.069, 0.036/
      data satkt /1.07e-6, 1.41e-5, 5.23e-6, 2.81e-6, 2.81e-6,
     &  3.38e-6, 4.45e-6, 2.04e-6, 2.45e-6, 7.22e-6, 1.34e-6,
     &  9.74e-7, 3.38e-6, 0.0, 1.41e-4, 1.41e-5/
      data tsat/0.339, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404,
     &  0.464, 0.465, 0.406, 0.468, 0.468, 0.439, 1.0, 0.20, 0.421/
#else
      data b/4.26,8.72,11.55,4.74,10.73,8.17,6.77,5.25,4.26/
      data satpsi/.04,.62,.47,.14,.10,.26,.14,.36,.04/
      data satkt/1.41e-5,.20e-5,.10e-5,.52e-5,.72e-5,
     &           .25e-5,.45e-5,.34e-5,1.41e-5/
      data tsat/.421,.464,.468,.434,.406,.465,.404,.439,.421/
#endif
#else
      data b/4.05,4.38,4.9,5.3,5.39,7.12,7.75,8.52,10.4/
c    &       ,10.4,11.4/
      data satpsi/.121,.09,.218,.786,.478,.299,.356,.63,.153/
c    &            ,.49,.405/
      data satkt/1.76e-4,1.5633e-4,3.467e-5,7.2e-6,6.95e-6
     &           ,6.3e-6,1.7e-6,2.45e-6,2.167e-6/
c    &           ,1.033e-6,1.283e-6/
      data tsat/.395,.41,.435,.485,.451,.42,.477,.476,.426/
c    &          ,.492,.482/
#endif
      end
