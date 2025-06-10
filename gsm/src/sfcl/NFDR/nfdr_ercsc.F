      subroutine nfdr_ercsc(w1dp,w10d,w100d,w1000d,wwood,wherbc,
     *                      sig1d,sig10d,sg100d,s1000d,sgherb,
     *                      sgwood,depth,fm1,fm10,fm100,
     *                      fm1000,fmwood,fherbc,extmoi,hd,wspix,
     *                      wndftr,slppix,isc,ierc,rire,
     *                      sgbrt,ros)
c
      implicit none
c 
c  This subroutine calculates the weighting factors for both the isc
c  and ierc and solves the various equations of the Rothermel spread
c  model. Fuels are weighted by surface area for the isc and by load-
c  ing for the ierc. This set of equations is current as of Jan 76.
c
      real w1dp,w10d,w100d,w1000d,wwood,wherbc,
     *     sig1d,sig10d,sg100d,s1000d,sgherb,
     *     sgwood,depth,fm1,fm10,fm100,
     *     fm1000,fmwood,fherbc,extmoi,hd,wspix,
     *     wndftr,slppix,rire,ros
      real fixrho,fct1e,fct10e,f1000e,sgbrl,wtotdn,wtfmd,phislp
      integer ierc,isc
c
c  integer islop
c
      real sd,sl,rhod,rhol,std,stl
      data sd,sl,rhod,rhol,std,stl /0.01,0.01,32.,32.,.0555,.0555/
      real slpfct(5)
      data slpfct /.267,.533,1.068,2.134,4.273/
c
      real h1,fixho,wtotd,wtotl,sa1,sa10,sa100,sherbc,sawood,
     *     sadead,salive,fct1,fct10,fct100,fcherb,fcwood,fwoode,
     *     fhrbce,extliv,fcded,fcliv,wtot,
     *     fctle,fc10e,fc100e,fc1000e,fcdede,fclive,
     *     sgbrd,sgbr1,sgbrde,sgbrle,sgbrt,sgbrte,
     *     w1n,w10n,w100n,whernc,wwoodn,wtotln,wbar,wliv,
     *     rhobar,betbar,zeta,betop,betope,ad,ade,as,ase,
     *     gmapm,gmapme,etamd,etaml,c,b,e,ufact,
     *     wtrmd,wtfml,wtfmde,wtfmle,hnu1,hnu10,hnu100,wrat,
     *     fmff,etatd,etatl,etatde,etatle,
     *     rir,xxws,phiwnd,htsink,tau,xint,hl
c
c  islop = slppix
c
      hl=hd
      fixrho = (w1dp+w10d+w100d+wwood+wherbc)/depth
c 
c  Total loadings of live and dead fuel categories.
c
      wtotd=w1dp+w10d+w100d+w1000d
      wtotl=wwood+wherbc
c 
c  Weighting factors for rate-of-spread. (by sfc area)
c  SFC area of dead and live fuel classes.
c
      sa1=w1dp*sig1d/rhod
      sa10=w10d*sig10d/rhod
      sa100=w100d*sg100d/rhod
      sherbc=wherbc*sgherb/rhol
      sawood=wwood*sgwood/rhol
c 
c  Total sfc areas of dead and live fuel categories (by sfc area).
c
      sadead=sa1+sa10+sa100
      salive=sawood+sherbc
c 
c  Weighting factors for dead and live fuel classes (by sfc area).
c
      fct1=sa1/sadead
      fct10=sa10/sadead
      fct100=sa100/sadead
      fcherb=0.
      fcwood = 0
      fwoode = 0
      fhrbce = 0
      extliv = 0
c 
c  If there are no live fuels present
c
      if(wtotl.le.0.) go to 2
      fcherb=sherbc/salive
      fcwood=sawood/salive
c 
c  Weighting factors for dead and live fuel categories.
c
    2 continue
      fcded=sadead/(sadead+salive)
      fcliv=salive/(sadead+salive)
c 
c  Weighting factors for energy release. (by loading)
c
c  Total fuel loading
c
      wtot=wtotd+wtotl
c 
c  Weighting factors for dead and live fuel classes (by load).
c
      fct1e =w1dp/wtotd
      fct10e=w10d/wtotd
      fc100e=w100d/wtotd
      f1000e=w1000d/wtotd
      if(wtotl.le.0.) go to 3
      fwoode=wwood/wtotl
      fhrbce=wherbc/wtotl
c 
c  Weighting factors for dead and live fuel categories (by load).
c
    3 continue
      fcdede=wtotd/wtot
      fclive=wtotl/wtot
c 
c  Following are basic equations of the rate of spread and energy
c  release models. They are taken from Rothermels work 1970-1976.
c  variables with names ending in /e/ are for use in the ierc calc.
c
c  Weighted sfc-vol ratios of dead and live fuel categories
c
      sgbrd=fct1*sig1d+fct10*sig10d+fct100*sg100d
      sgbrl=fcwood*sgwood+fcherb*sgherb
      sgbrde=fct1e*sig1d+fct10e*sig10d+fc100e*sg100d+f1000e*s1000d
      sgbrle=fwoode*sgwood+fhrbce*sgherb
c 
c  Characteristic sfc-vol ratios
c
      sgbrt=sgbrd*fcded+sgbrl*fcliv
      sgbrte=sgbrde*fcdede+sgbrle*fclive
c 
c  Compute net fuel loading
c
      w1n=w1dp*(1.-std)
      w10n=w10d*(1.-std)
      w100n=w100d*(1.-std)
      whernc=wherbc*(1.-stl)
      wwoodn=wwood*(1.-stl)
      wtotln=wtotl*(1.-stl)
      wtotdn=wtotd*(1.-std)
c 
c  Characteristic loading of dead and live fuel categories
c
      wbar=fct1*w1n+fct10*w10n+fct100*w100n
      wliv=fcwood*wwoodn+fcherb*whernc
      if(sgwood.ge.1200..and.sgherb.ge.1200.)wliv=wtotln
c 
c  ******* new code for increasing load as a function of kdi *******
c  Calculate a new depth to preserve packing ratio of original model
c
      depth=(wtot-w1000d)/fixrho
c
c  ********************** end new code ****************************
c 
c  Bulk density
c
      rhobar=(wtot-w1000d)/depth
c 
c  Packing ratio
c
      betbar=rhobar/rhod
c 
c  Propogating flux ratio
c
      zeta=(192.0+0.2595*sgbrt)**(-1)*exp((0.792+0.681*sgbrt**0.5)*
     &        (betbar+0.1))
c 
c  Optimum packing ratios
c
      betop=3.348*(sgbrt**(-0.8189))
      betope=3.348*(sgbrte**(-0.8189))
c 
c  Optimum reaction velocities
c
      ad=133.*sgbrt**(-0.7913)
      ade=133.*sgbrte**(-0.7913)
      as=sgbrt**1.5
      ase=sgbrte**1.5
c 
      gmapm=as/(495.0+0.0594*as)*(betbar/betop)**ad*exp(ad*(1.0-
     &              betbar/betop))
      gmapme=ase/(495.0+0.0594*ase)*(betbar/betope)**ade*exp(ade*(
     &         1.0-betbar/betope))
c 
c  Mineral damping coefficients.
c
      etamd=0.174*sd**(-0.19)
      if(etamd.gt.1.)etamd=1.
      etaml=0.174*sl**(-0.19)
      if(etaml.gt.1.)etaml=1.
 
c     Fuel bed wind controlling factor.
      c=7.47*exp(-0.133*(sgbrt**0.55))
      b=0.02526*(sgbrt**0.54)
      e=0.715*exp(-3.59e-4*sgbrt)
      ufact=c*((betbar/betop)**(-e))
c 
c  Weighted moisture content of dead and live fuel categories
c
      wtfmd=fct1*fm1+fct10*fm10+fct100*fm100
      wtfml=fcwood*fmwood+fcherb*fherbc
      wtfmde=fct1e*fm1+fct10e*fm10+fc100e*fm100+f1000e*fm1000
      wtfmle=fwoode*fmwood+fhrbce*fherbc
c 
c  If there are no live fuels present
c
      if(wtotln.le.0.) go to 11
c 
c  Heating number for extliv.
c
      hnu1=w1n*exp(-138./sig1d)
      hnu10=w10n*exp(-138./sig10d)
      hnu100=w100n*exp(-138./sg100d)
c 
c  Dead-live loading ratio for extliv.
c
      wrat=(hnu1+hnu10+hnu100)/(wwoodn*exp(-500./sgwood)
     &     +whernc*exp(-500./sgherb))
c 
c  Fine dead fuel moisture content for extliv.
c
      fmff=(fm1*hnu1+fm10*hnu10+fm100*hnu100)/(hnu1+hnu10+hnu100)
c 
c  Moisture of extinction for live fuels (extliv).
c
      extliv=(2.9*wrat*(1.-fmff/extmoi)-0.226)*100.
   11 continue
      if(extliv.lt.extmoi)extliv=extmoi
c 
c  Moisture damping coefficients.
c
      etatd=1.-2.59*(wtfmd/extmoi)+5.11*(wtfmd/extmoi)
     &      **2-3.52*(wtfmd/extmoi)**3
      etatl=1.-2.59*(wtfml/extliv)+5.11*(wtfml/extliv)
     &      **2-3.52*(wtfml/extliv)**3
      etatde=1.-2.*(wtfmde/extmoi)+1.5*(wtfmde/extmoi)
     &      **2-0.5*(wtfmde/extmoi)**3
      etatle=1.-2.*(wtfmle/extliv)+1.5*(wtfmle/extliv)
     &      **2-0.5*(wtfmle/extliv)**3
c 
      etatd=amin1(1.,amax1(0.,etatd))
      etatl=amin1(1.,amax1(0.,etatl))
      etatde=amin1(1.,amax1(0.,etatde))
      etatle=amin1(1.,amax1(0.,etatle))
c 
c  Reaction intensities.
c
      rir=gmapm*((wbar*hd*etamd*etatd)+(wliv*hl*etaml*etatl))
      rire=gmapme*((fcdede*wtotdn*hd*etamd*etatde)+(fclive*wtotln*
     &      hl*etaml*etatle))
c 
c  Wind coefficient
c
      xxws = wspix*wndftr*88.
c
c  Limit the effective wind speed to 9tenths the reaction
c  intensity to approximate the dropping off of spread when
c  the wind speed acts as a heat sink.
c
      if(xxws.gt.(.9*rir)) xxws = .9*rir
      phiwnd = ufact*(xxws**b)
c 
c  Slope coefficient
c
      phislp=5.275*betbar**(-0.3)*(slppix*0.01)**2
c 
c  Heat sink
c
      htsink=rhobar*(fcded*(fct1*exp(-138./sig1d)*(250.+11.16*fm1)+fct10
     &*exp(-138./sig10d)*(250.+11.16*fm10)+fct100*exp(-138./sg100d)*(250
     &.+11.16*fm100))+fcliv*(fcherb*exp(-138./sgherb)*(250.+11.16*fherbc
     &)+fcwood*exp(-138./sgwood)*(250.+11.16*fmwood)))
 
c 
c  Spread component
c
      ros=rir*zeta*(1.+phiwnd+phislp)/htsink
      isc=nint(ros)
c 
c  Calculations for the energy release component as modified by
c  Deeming
c 
c  Residence time (tau) as modified  in 1976 by rothermel.
c
      tau=384./sgbrt
      xint=tau*rire
c 
c  Energy release component
c
      ierc=nint(xint/25.)
 
      return
      end
