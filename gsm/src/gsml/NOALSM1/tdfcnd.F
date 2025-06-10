      subroutine tdfcnd ( df, smc, qz,  smcmax, sh2o)
c
      implicit none
c
c ----------------------------------------------------------------------
c subroutine tdfcnd
c ----------------------------------------------------------------------
c calculate thermal diffusivity and conductivity of the soil for a given
c point and time.
c ----------------------------------------------------------------------
c peters-lidard approach (peters-lidard et al., 1998)
c june 2001 changes: frozen soil condition.
c ----------------------------------------------------------------------
       real df
       real gammd
       real thkdry
       real ake
       real thkice
       real thko
       real thkqtz
       real thksat
       real thks
       real thkw
       real qz
       real satratio
       real sh2o
       real smc
       real smcmax
       real xu
       real xunfroz
c
c
c ----------------------------------------------------------------------
c we now get quartz as an input argument (set in routine redprm):
c      data quartz /0.82, 0.10, 0.25, 0.60, 0.52, 
c     &             0.35, 0.60, 0.40, 0.82/
c ----------------------------------------------------------------------
c if the soil has any moisture content compute a partial sum/product
c otherwise use a constant value which works well with most soils
c ----------------------------------------------------------------------
c  thkw ......water thermal conductivity
c  thkqtz ....thermal conductivity for quartz
c  thko ......thermal conductivity for other soil components
c  thks ......thermal conductivity for the solids combined(quartz+other)
c  thkice ....ice thermal conductivity
c  smcmax ....porosity (= smcmax)
c  qz .........quartz content (soil type dependent)
c ----------------------------------------------------------------------
c use as in peters-lidard, 1998 (modif. from johansen, 1975).
c
c                                  pablo grunmann, 08/17/98
c refs.:
c      farouki, o.t.,1986: thermal properties of soils. series on rock 
c              and soil mechanics, vol. 11, trans tech, 136 pp.
c      johansen, o., 1975: thermal conductivity of soils. ph.d. thesis,
c              university of trondheim,
c      peters-lidard, c. d., et al., 1998: the effect of soil thermal 
c              conductivity parameterization on surface energy fluxes
c              and temperatures. journal of the atmospheric sciences,
c              vol. 55, pp. 1209-1224.
c ----------------------------------------------------------------------
c needs parameters
c porosity(soil type):
c      poros = smcmax
c saturation ratio:
      satratio = smc/smcmax
c parameters  w/(m.k)
      thkice = 2.2
      thkw = 0.57
      thko = 2.0
c      if (qz .le. 0.2) thko = 3.0
      thkqtz = 7.7
c solids' conductivity      
      thks = (thkqtz**qz)*(thko**(1.- qz))
c unfrozen fraction (from 1., i.e., 100%liquid, to 0. (100% frozen))

      if(smc.gt.0) then
      xunfroz = sh2o /smc   !  ? both tend to zero simultaneously
      else
      xunfroz=0.
      endif

c unfrozen volume for saturation (porosity*xunfroz)
      xu=xunfroz*smcmax 
c saturated thermal conductivity
      thksat = thks**(1.-smcmax)*thkice**(smcmax-xu)*thkw**(xu)
c dry density in kg/m3
      gammd = (1. - smcmax)*2700.
c dry thermal conductivity in w.m-1.k-1
      thkdry = (0.135*gammd + 64.7)/(2700. - 0.947*gammd)
      if ( (sh2o + 0.0005) .lt. smc ) then
c frozen
              ake = satratio
      else
c unfrozen
c range of validity for the kersten number (ake)
          if ( satratio .gt. 0.1 ) then
c kersten number (using "fine" formula, valid for soils containing at 
c least 5% of particles with diameter less than 2.e-6 meters.)
c (for "coarse" formula, see peters-lidard et al., 1998).
              ake = log10(satratio) + 1.0
          else
c use k = kdry
              ake = 0.0
          endif
      endif
c  thermal conductivity
       df = ake*(thksat - thkdry) + thkdry
c ----------------------------------------------------------------------
c end subroutine tdfcnd
c ----------------------------------------------------------------------
      return
      end
