      subroutine dtfrznew(tu,p,thteu,qu,qfrz,qice,aliq,bliq,cliq,dliq)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      real,         intent(in   )   :: p,qfrz,aliq,bliq,cliq,dliq
      real,         intent(inout)   :: tu,thteu,qu,qice
      real    ::    rlc,rls,rlf,cpp,a,dtfrz,es,qs,dqevap,pii
!-----------------------------------------------------------------------
!
!..allow the freezing of liquid water in the updraft to proceed as an 
!..approximately linear function of temperature in the temperature range
!..ttfrz to tbfrz...
!..for colder termperatures, freeze all liquid water...
!..thermodynamic properties are still calculated with respect to liquid water
!..to allow the use of lookup table to extract tmp from thetae...
!
      rlc=2.5e6-2369.276*(tu-273.16)
      rls=2833922.-259.532*(tu-273.16)
      rlf=rls-rlc
      cpp=1005.7*(1.+0.89*qu)
!
!  a = d(es)/dt is that calculated from buck (1981) emperical formulas
!  for saturation vapor pressure...
!
      a=(cliq-bliq*dliq)/((tu-dliq)*(tu-dliq))
      dtfrz = rlf*qfrz/(cpp+rls*qu*a)
      tu = tu+dtfrz
      
      es = aliq*exp((bliq*tu-cliq)/(tu-dliq))
      qs = es*0.622/(p-es)
!
!...freezing warms the air and it becomes unsaturated...assume that some of the 
!...liquid water that is available for freezing evaporates to maintain satura-
!...tion...since this water has already been transferred to the ice category,
!...subtract it from ice concentration, then set updraft mixing ratio at the new
!...temperature to the saturation value...
!
      dqevap = qs-qu
      qice = qice-dqevap
      qu = qu+dqevap
      pii=(1.e5/p)**(0.2854*(1.-0.28*qu))
      thteu=tu*pii*exp((3374.6525/tu-2.5403)*qu*(1.+0.81*qu))
!
      end subroutine dtfrznew
