#define IGRD1 igrd1
#define IGRD igrd
      subroutine liftix(im,ix,kt,pt,ps,t,q,tm,sli,bli)
c$$$  subprogram documentation block
c
c  subprogram:    liftix      compute lifted indices from sigma
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: computes both the surface lifted index and best lifted index
c   from profiles in constant pressure thickness layers above ground.
c   the surface lifted index is computed by raising the lowest layer
c   to 500 mb and subtracting its parcel temperature
c   from the environment temperature.
c   the best lifted index is computed by finding the parcel
c   with the warmest equivalent potential temperature,
c   then raising it to 500 mb and subtracting its parcel temperature
c   from the environment temperature.
c
c program history log:
c   92-10-31  iredell
c   94-04-28  iredell   fixed parameters
c   94-06-03  iredell   returned two indices
c   94-07-29  iredell   used constant pressure thickness profiles
c
c usage:    call liftix(im,ix,kt,pt,ps,t,q,tm,sli,bli)
c
c   input argument list:
c     im       - integer number of points
c     ix       - integer first dimension of upper air data
c     kt       - integer number of layers in profile
c     pt       - real pressure thickness in kpa
c     ps       - real (im) surface pressure in kpa
c     t        - real (ix,kt) temperature in k
c     q        - real (ix,kt) specific humidity in kg/kg
c     tm       - real (im) 500 mb temperature in k
c
c   output argument list:
c     sli      - real (ix) surface lifted index in k
c     bli      - real (ix) best lifted index in k
c
c subprograms called:
c   (fpkap)   - function to compute pressure to the kappa
c   (ftdp)    - function to compute dewpoint temperature
c   (ftlcl)   - function to compute lifting condensation level
c   (fthe)    - function to compute equivalent potential temperature
c   (ftma)    - function to compute moist adiabat temperature
c
c attributes:
c   language: cray fortran
c
c$$$
cfppc expand(fpkap,ftdp,ftlcl,fthe,ftma)
#include <paramodel.h>
      dimension ps(im),t(im,kt),q(im,kt),tm(im),sli(im),bli(im)
      parameter(cp= 1.0046e+3 ,rd= 2.8705e+2 ,rv= 4.6150e+2 )
      parameter(rk=rd/cp,eps=rd/rv,epsm1=rd/rv-1.)
      parameter(plift=50.)
      dimension p2kmas(im),themas(im),p2kmab(im),themab(im)
c
c  select the warmest equivalent potential temperature
c
      do k=1,kt
        do i=1,im
          p=ps(i)-(k-0.5)*pt
          pv=p*q(i,k)/(eps-epsm1*q(i,k))
          tdpd=max(t(i,k)-ftdp(pv),0.)
          tlcl=ftlcl(t(i,k),tdpd)
          p2klcl=fpkap(p)*tlcl/t(i,k)
          thelcl=fthe(tlcl,p2klcl)
          if(k.eq.1) then
            p2kmas(i)=p2klcl
            themas(i)=thelcl
            p2kmab(i)=p2klcl
            themab(i)=thelcl
          elseif(thelcl.gt.themab(i)) then
            p2kmab(i)=p2klcl
            themab(i)=thelcl
          endif
        enddo
      enddo
c 
c  lift the parcel to 500 mb along a dry adiabat below the lcl
c  or along a moist adiabat above the lcl.
c  the lifted index is the environment minus parcel temperature.
c
      pliftk=(plift/100.)**rk
      do i=1,im
        if(ps(i).gt.plift) then
          p2ks=min(pliftk,p2kmas(i))
          sli(i)=tm(i)-pliftk/p2ks*ftma(themas(i),p2ks,qma)
          p2kb=min(pliftk,p2kmab(i))
          bli(i)=tm(i)-pliftk/p2kb*ftma(themab(i),p2kb,qma)
        else
          sli(i)=0.
          bli(i)=0.
        endif
      enddo
c
      return
      end
