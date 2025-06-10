      subroutine tpmix2dd(p,thes,ts,qs)
!
      implicit none
!-----------------------------------------------------------------------
! lookup table variables:
      integer, parameter :: kfnt=250,kfnp=220
      real, dimension(1:kfnt,1:kfnp) :: ttab,qstab
      real, dimension(1:kfnp) :: the0k
      real, dimension(1:200) :: alu
      real :: rdpr,rdthk,plutop
! end of lookup table variables:
!-----------------------------------------------------------------------
      real,         intent(in   )   :: p,thes
      real,         intent(inout)   :: ts,qs
      real    ::    tp,qq,bth,tth,pp,t00,t10,t01,t11,q00,q10,q01,q11
      integer ::    iptb,ithtb
!-----------------------------------------------------------------------

!
!*** lookup table variables (f77 format)... ****************************
      common/kflut/ ttab,qstab,the0k,alu,rdpr,rdthk,plutop
!*************************************************************** 
!
!***********************************************************************
!     scaling pressure and tt table index                         
!***********************************************************************
!
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
!    parcel temperature and saturation mixing ratio      
!***********************************************************************
!
      ts=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
!
      qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)
!
      end subroutine tpmix2dd
