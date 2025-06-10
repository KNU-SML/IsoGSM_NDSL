      subroutine pder(m,enn1,elonn1,eon,eontop,pln,plntop,
     &                plndx,plndy)
c$$$  subprogram documentation block
c
c subprogram:    pder        compute legendre polynomials derivatives
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: computes the zonal and meridional derivatives
c           of the normalized associated legendre polynomials
c           in the spectral triangle at a given latitude.
c           subprogram gspc and pleg should be called already.
c
c program history log:
c   91-10-31  mark iredell
c
c usage:    call pder(m,enn1,elonn1,eon,eontop,pln,plntop,
c    &                plndx,plndy)
c
c   input argument list:
c     m        - integer spectral truncation
c     enn1     - real ((m+1)*(m+2)/2) n*(n+1)/a**2
c     elonn1   - real ((m+1)*(m+2)/2) l/(n*(n+1))*a
c     eon      - real ((m+1)*(m+2)/2) epsilon/n*a
c     eontop   - real (m+1) epsilon/n*a over top
c     pln      - real ((m+1)*(m+2)/2) legendre polynomial
c     plntop   - real (m+1) legendre polynomial over top
c
c   output argument list:
c     plndx    - real ((m+1)*(m+2)/2) zonal derivatives (no i)
c     plndy    - real ((m+1)*(m+2)/2) meridional derivatives
c
c attributes:
c   language: cray fortran
c
c$$$
cfpp$ noconcur r
      real enn1((m+1)*(m+2)/2),elonn1((m+1)*(m+2)/2)
      real eon((m+1)*(m+2)/2),eontop(m+1)
      real pln((m+1)*(m+2)/2),plntop(m+1)
      real plndx((m+1)*(m+2)/2),plndy((m+1)*(m+2)/2)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute polynomial derivatives for spectral analysis
      i=1
      plndy(i)=eon(i+1)*pln(i+1)
      do i=2,(m+1)*(m+2)/2-1
        plndy(i)=eon(i+1)*pln(i+1)-eon(i)*pln(i-1)
      enddo
      i=(m+1)*(m+2)/2
      plndy(i)=-eon(i)*pln(i-1)
      do l=0,m
        nml=m-l
        i=l*(2*m+3-l)/2+(nml+1)
        plndy(i)=plndy(i)+eontop(l+1)*plntop(l+1)
      enddo
      do i=1,(m+1)*(m+2)/2
        plndx(i)=enn1(i)*elonn1(i)*pln(i)
        plndy(i)=enn1(i)*plndy(i)
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
