      subroutine kf_lutab(svp1,svp2,svp3,svpt0)
!
!  this subroutine is a lookup table.
!  given a series of series of saturation equivalent potential
!  temperatures, the temperature is calculated.
!
!--------------------------------------------------------------------
      implicit none
!--------------------------------------------------------------------
! lookup table variables
      integer, parameter :: kfnt=250,kfnp=220
      real, dimension(1:kfnt,1:kfnp) :: ttab,qstab
      real, dimension(1:kfnp) :: the0k
      real, dimension(1:200) :: alu
      real :: rdpr,rdthk,plutop
! end of lookup table variables

      integer :: kp,it,itcnt,i
      real :: dth,tmin,toler,pbot,dpr, 
     1      temp,p,es,qs,pi,thes,tgues,thgues,f0,t1,t0,thgs,f1,dt, 
     2      astrt,ainc,a1,thtgs
      real    :: aliq,bliq,cliq,dliq
      real, intent(in)    :: svp1,svp2,svp3,svpt0
!
      common/kflut/ ttab,qstab,the0k,
     *              alu,rdpr,rdthk,plutop
! equivalent potential temperature increment
      data dth/1./
! minimum starting temp
      data tmin/150./
! tolerance for accuracy of temperature
      data toler/0.001/
! top pressure (pascals)
      plutop=5000.0
! bottom pressure (pascals)
      pbot=110000.0

      aliq = svp1*1000.
      bliq = svp2
      cliq = svp2*svpt0
      dliq = svp3

!
! compute parameters
!
! 1._over_(sat. equiv. theta increment)
      rdthk=1./dth
! pressure increment
!
      dpr=(pbot-plutop)/real(kfnp-1)
!      dpr=(pbot-plutop)/real(kfnp-1)
! 1._over_(pressure increment)
      rdpr=1./dpr
! compute the spread of thes
!     thespd=dth*(kfnt-1)
!
! calculate the starting sat. equiv. theta
!
      temp=tmin
      p=plutop-dpr
      do kp=1,kfnp
        p=p+dpr
        es=aliq*exp((bliq*temp-cliq)/(temp-dliq))
        qs=0.622*es/(p-es)
        pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
        the0k(kp)=temp*pi*exp((3374.6525/temp-2.5403)*qs* 
     +         (1.+0.81*qs))
      enddo
!
! compute temperatures for each sat. equiv. potential temp.
!
      p=plutop-dpr
      do kp=1,kfnp
        thes=the0k(kp)-dth
        p=p+dpr
        do it=1,kfnt
! define sat. equiv. pot. temp.
          thes=thes+dth
! iterate to find temperature
! find initial guess
          if(it.eq.1) then
            tgues=tmin
          else
            tgues=ttab(it-1,kp)
          endif
          es=aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
          qs=0.622*es/(p-es)
          pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
          thgues=tgues*pi*exp((3374.6525/tgues-2.5403)*qs*  
     +         (1.+0.81*qs))
          f0=thgues-thes
          t1=tgues-0.5*f0
          t0=tgues
          itcnt=0
! iteration loop
          do itcnt=1,11
            es=aliq*exp((bliq*t1-cliq)/(t1-dliq))
            qs=0.622*es/(p-es)
            pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
            thtgs=t1*pi*exp((3374.6525/t1-2.5403)*qs*(1.+0.81*qs))
            f1=thtgs-thes
            if(abs(f1).lt.toler)then
              exit
            endif
!           itcnt=itcnt+1
            dt=f1*(t1-t0)/(f1-f0)
            t0=t1
            f0=f1
            t1=t1-dt
          enddo
          ttab(it,kp)=t1
          qstab(it,kp)=qs
        enddo
      enddo
!
! lookup table for tlog(emix/aliq)
!
! set up intial values for lookup tables
!
       astrt=1.e-3
       ainc=0.075
!
       a1=astrt-ainc
       do i=1,200
         a1=a1+ainc
         alu(i)=alog(a1)
       enddo
!
      end subroutine kf_lutab
