      subroutine nfdr_oneten(dbpix,rhpix,cldpix,fm100,tf,fm1,fm10)
c
      implicit none
c
c     This subroutine calculates the 1-hr-tl fuel moisture (fm1) and the
c     10-hr-tl fuel moisture (fm10).  FM10 is not dependent on fm100 value
c 
c  input 
c 
      real dbpix,rhpix,cldpix,fm100
c
c  output
c
      real fm1,fm10
c
      integer isw
      real tf,hf,emc
c
      real tfct(4),hfct(4)
      data tfct/25.,19.,12.,5./
      data hfct/.75,.83,.92,1.00/
c 
c  calculate temperature and relative humidity at the fuel level
c
      isw=0
      if(cldpix .gt. 0.10) isw=1
      if(cldpix .gt. 0.50) isw=2
      if(cldpix .gt. 0.90) isw=3
      if(isw .gt. 3) isw = 3
      tf=tfct(isw+1)+dbpix
      hf=hfct(isw+1)*rhpix
      call nfdr_eqmc(tf,hf,emc)
c
c  Calculate the 10 hour and 1 hour fuel moistures 
c
      fm10=1.28*emc*0.8+0.2*fm100
      if (fm10 .lt. 2.0) fm10=2.0
      fm1=(4.*1.03*emc+fm10)/5.
      if(fm1 .lt. 2.0)fm1=2.0
c
      return
      end
