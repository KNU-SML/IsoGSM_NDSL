      subroutine envirtht(p1,t1,q1,tht1,aliq,bliq,cliq,dliq) 
!
!--------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      real,    intent(in   )   :: p1,t1,q1,aliq,bliq,cliq,dliq
      real,    intent(inout)   :: tht1
      real  :: ee,tlog,astrt,ainc,a1,tp,value,aintrp,tdpt,tsat,tht, 
     1           t00,p00,c1,c2,c3,c4,c5
      integer ::    indlu
! lookup table variables:
      integer, parameter :: kfnt=250,kfnp=220
      real, dimension(1:kfnt,1:kfnp) :: ttab,qstab
      real, dimension(1:kfnp) :: the0k
      real, dimension(1:200) :: alu
      real :: rdpr,rdthk,plutop
! end of lookup table variables:
!-----------------------------------------------------------------------

!
!*** lookup table variables (f77 format)... ****************************
      common/kflut/ ttab,qstab,the0k,alu,rdpr,rdthk,plutop
!***************************************************************
!-----------------------------------------------------------------------
      data t00,p00,c1,c2,c3,c4,c5/273.16,1.e5,3374.6525,2.5403,3114.834,
     1     0.278296,1.0723e-3/  
!                                                                       
!  calculate environmental equivalent potential temperature...          
!                                                                       
! note: calculations for mixed/ice phase no longer used...jsk 8/00
!
      ee=q1*p1/(0.622+q1)                                             
!     tlog=alog(ee/aliq)                                              
! ...calculate log term using lookup table...
!
      astrt=1.e-3
      ainc=0.075
      a1=ee/aliq
      tp=(a1-astrt)/ainc
      indlu=int(tp)+1
      value=(indlu-1)*ainc+astrt
      aintrp=(a1-value)/ainc
      tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
!
      tdpt=(cliq-dliq*tlog)/(bliq-tlog)                               
      tsat=tdpt-(.212+1.571e-3*(tdpt-t00)-4.36e-4*(t1-t00))*(t1-tdpt) 
      tht=t1*(p00/p1)**(0.2854*(1.-0.28*q1))                          
      tht1=tht*exp((c1/tsat-c2)*q1*(1.+0.81*q1))                      
!
      end subroutine envirtht                                     
