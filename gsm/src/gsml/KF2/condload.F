      subroutine condload(qliq,qice,wtw,dz,boterm,enterm,rate,qnewlq,
     1                    qnewic,qlqout,qicout,g)

!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!  9/18/88...this precipitation fallout scheme is based on the scheme us
!  by ogura and cho (1973).  liquid water fallout from a parcel is cal-
!  culated using the equation dq=-rate*q*dt, but to simulate a quasi-
!  continuous process, and to eliminate a dependency on vertical
!  resolution this is expressed as q=q*exp(-rate*dz).

      real, intent(in   )   :: g
      real, intent(in   )   :: dz,boterm,enterm,rate
      real, intent(inout)   :: qlqout,qicout,wtw,qliq,qice,qnewlq,qnewic
      real :: qtot,qnew,qest,g1,wavg,conv,ratio3,oldq,ratio4,dq,pptdrg

!
!  9/18/88...this precipitation fallout scheme is based on the scheme us
!  by ogura and cho (1973).  liquid water fallout from a parcel is cal- 
!  culated using the equation dq=-rate*q*dt, but to simulate a quasi-   
!  continuous process, and to eliminate a dependency on vertical        
!  resolution this is expressed as q=q*exp(-rate*dz).                   
      qtot=qliq+qice                                                    
      qnew=qnewlq+qnewic                                                
!                                                                       
!  estimate the vertical velocity so that an average vertical velocity 
!  be calculated to estimate the time required for ascent between model 
!  levels...                                                            
!                                                                       
      qest=0.5*(qtot+qnew)                                              
      g1=wtw+boterm-enterm-2.*g*dz*qest/1.5                             
      if(g1.lt.0.0)g1=0.                                                
      wavg=0.5*(sqrt(wtw)+sqrt(g1))                                    
      conv=rate*dz/wavg                                                 
!                                                                       
!  ratio3 is the fraction of liquid water in fresh condensate, ratio4 is
!  the fraction of liquid water in the total amount of condensate involv
!  in the precipitation process - note that only 60% of the fresh conden
!  sate is is allowed to participate in the conversion process...       
!                                                                       
      ratio3=qnewlq/(qnew+1.e-8)                                       
!     oldq=qtot                                                         
      qtot=qtot+0.6*qnew                                                
      oldq=qtot                                                         
      ratio4=(0.6*qnewlq+qliq)/(qtot+1.e-8)                            
      qtot=qtot*exp(-conv)                                              
!                                                                       
!  determine the amount of precipitation that falls out of the updraft  
!  parcel at this level...                                              
!                                                                       
      dq=oldq-qtot                                                      
      qlqout=ratio4*dq                                                  
      qicout=(1.-ratio4)*dq                                             
!                                                                       
!  estimate the mean load of condensate on the updraft in the layer, cal
!  late vertical velocity                                               
!                                                                       
      pptdrg=0.5*(oldq+qtot-0.2*qnew)                                   
      wtw=wtw+boterm-enterm-2.*g*dz*pptdrg/1.5                          
      if(abs(wtw).lt.1.e-4)wtw=1.e-4
!                                                                       
!  determine the new liquid water and ice concentrations including losse
!  due to precipitation and gains from condensation...                  
!                                                                       
      qliq=ratio4*qtot+ratio3*0.4*qnew                                  
      qice=(1.-ratio4)*qtot+(1.-ratio3)*0.4*qnew                        
      qnewlq=0.                                                         
      qnewic=0.                                                         

      end subroutine condload
