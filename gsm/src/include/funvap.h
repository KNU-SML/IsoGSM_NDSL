c
c  water vapor functions
c
c  virtual temperature from temperature and specific humidity
c
      ftv(xt,xq)=xt*
     1   (1.e0+(rv_/rd_-1.e0)*xq)
c
c  temperature from virtual temperature and specific humidity
c
      ft(xtv,xq)=xtv/
     1   (1.e0+(rv_/rd_-1.e0)*xq)
c
c  specific humidity from vapor pressure and air pressure
c
      fq(xpv,xp)=(rd_/rv_)*xpv/
     1    (xp+(rd_/rv_-1.e0)*xpv)
c
c  vapor pressure from specific humidity and air pressure
c
      fpv(xq,xp)=xp*xq/
     1    ((rd_/rv_)-
     2     (rd_/rv_-1.e0)*xq)
c
c  saturation vapor pressure (kpa) from temperature
c
      fpvs(xt)=(_psat_*1.e-3)*
     1     exp((hvap_/rv_)/ttp_
     2         -(hvap_/rv_)/xt)
c
c  saturation specific humidity from temperature and pressure (kpa)
c
      fqs(xt,xp)=fq(fpvs(xt),xp)
