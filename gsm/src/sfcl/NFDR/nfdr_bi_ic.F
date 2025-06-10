      subroutine nfdr_bi_ic(sgbrt,ire,extmoi,fm1,tf,ros,scm,bi,ic)
c
      implicit none
c
      real sgbrt,ire,extmoi,fm1,tf,ros,scm
      integer bi,ic
c
      real qigmx,pwr,const,pnorm1,pnorm2
      data qigmx/344./,pwr/3.66/,const/.000923/,pnorm1/.00232/
      data pnorm2/.99767/
c
      real tau,fl,qign,chi,tmpcls,test,pi,scn
c 
c  Residence time (tau). Modified 1976 by Rothermel
c
      tau=384./sgbrt
c 
c  Byrams flame length.
c
      fl=.45*(ros/60.*ire*tau)**.46
c 
c  Burning index
c
      bi = nint(fl*10.0)
c
c  Ignition component calculation.
c
      qign=0.
      chi=0.
c
      if(fm1.ge.extmoi)go to 10
      tmpcls=(tf-32.)/1.8
      qign=144.5-0.266*tmpcls-5.8e-4*tmpcls*tmpcls-tmpcls*fm1/100.
     &        +18.54*(1.-exp(-.151*fm1))+6.4*fm1
c
c  If the ignition energy required is too high leave.
c
      if(qign.gt.qigmx) goto 10
      chi=(qigmx-qign)/10.
      test=chi**pwr*const/50.
      if(test.le.pnorm1) go to 10
      pi=min(100.,max(0.,(test-pnorm1)*100./pnorm2))
      scn = min(100.0,(100.0*ros/scm))
      ic=int(min(100.0,max(0.0,(0.1*pi*scn**0.5))))
      return
   10 continue
      ic=0
      return
      end
