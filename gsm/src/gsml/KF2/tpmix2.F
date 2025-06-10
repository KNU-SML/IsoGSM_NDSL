      subroutine tpmix2(p,thes,tu,qu,qliq,qice,qnewlq,qnewic,xlv1,xlv0)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!
! lookup table variables:
      integer, parameter :: kfnt=250,kfnp=220
      real, dimension(1:kfnt,1:kfnp) :: ttab,qstab
      real, dimension(1:kfnp) :: the0k
      real, dimension(1:200) :: alu
      real :: rdpr,rdthk,plutop
! end of lookup table variables:
      real,         intent(in   )   :: p,thes,xlv1,xlv0
      real,         intent(out  )   :: qnewlq,qnewic
      real,         intent(inout)   :: tu,qu,qliq,qice
      real    ::    tp,qq,bth,tth,pp,t00,t10,t01,t11,q00,q10,q01,q11,  
     1           temp,qs,qnew,dq,qtot,rll,cpp
      integer ::    iptb,ithtb
!-----------------------------------------------------------------------

!c******** lookup table variables... ****************************
      common/kflut/ ttab,qstab,the0k,
     *              alu,rdpr,rdthk,plutop 
!c*************************************************************** 
!c
!c**********************************************************************
!c     scaling pressure and tt table index                         
!c**********************************************************************
!c
      tp=(p-plutop)*rdpr
      qq=tp-aint(tp)
      iptb=int(tp)+1

!
!***********************************************************************
!              base and scaling factor for the                          
!***********************************************************************
!
!  scaling the and tt table index                                      
      bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
      tth=(thes-bth)*rdthk
      pp   =tth-aint(tth)
      ithtb=int(tth)+1
       if(iptb.ge.220 .or. iptb.le.1 .or. ithtb.ge.250
     1      .or. ithtb.le.1)then
         write(98,*)'**** out of bounds *********'
!        call flush(98)
       endif
!
      t00=ttab(ithtb  ,iptb  )
      t10=ttab(ithtb+1,iptb  )
      t01=ttab(ithtb  ,iptb+1)
      t11=ttab(ithtb+1,iptb+1)
!
      q00=qstab(ithtb  ,iptb  )
      q10=qstab(ithtb+1,iptb  )
      q01=qstab(ithtb  ,iptb+1)
      q11=qstab(ithtb+1,iptb+1)
!
!***********************************************************************
!            parcel temperature                                        
!***********************************************************************
!
      temp=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
!
      qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)
!
      dq=qs-qu
      if(dq.le.0.)then
        qnew=qu-qs
        qu=qs
      else 
!
!   if the parcel is subsaturated, temperature and mixing ratio must be
!   adjusted...if liquid water is present, it is allowed to evaporate
! 
        qnew=0.
        qtot=qliq+qice
!
!   if there is enough liquid or ice to saturate the parcel, temp stays at its
!   wet bulb value, vapor mixing ratio is at saturated level, and the mixing
!   ratios of liquid and ice are adjusted to make up the original saturation
!   deficit... otherwise, any available liq or ice vaporizes and appropriate
!   adjustments to parcel temp; vapor, liquid, and ice mixing ratios are made.
!
!...subsaturated values only occur in calculations involving various mixtures of
!...updraft and environmental air for estimation of entrainment and detrainment.
!...for these purposes, assume that reasonable estimates can be given using 
!...liquid water saturation calculations only - i.e., ignore the effect of the
!...ice phase in this process only...will not affect conservative properties...
!
        if(qtot.ge.dq)then
          qliq=qliq-dq*qliq/(qtot+1.e-10)
          qice=qice-dq*qice/(qtot+1.e-10)
          qu=qs
        else
          rll=xlv0-xlv1*temp
          cpp=1005.7*(1.+0.89*qu)
          if(qtot.lt.1.e-10)then
!
!...if no liquid water or ice is available, temperature is given by:
            temp=temp+rll*(dq/(1.+dq))/cpp
          else
!
!...if some liq water/ice is available, but not enough to achieve saturation,
!   the temperature is given by:
!
            temp=temp+rll*((dq-qtot)/(1+dq-qtot))/cpp
            qu=qu+qtot
            qtot=0.
            qliq=0.
            qice=0.
          endif
        endif
      endif
      tu=temp
      qnewlq=qnew
      qnewic=0.
!
      end subroutine tpmix2
