      module cldwat
!-------------------------------------------------------------
!-------------------------------------------------------------
        use ppgrid,        only: pcols, pver
        use physconst,     only: tmelt

        implicit none

        private
        save
        public cldwat_fice   ! Public interfaces
        real,public :: conke ! tunable constant for evaporation of precip

      contains

      subroutine cldwat_fice(ncol, t, fice, fsnow)
!-------------------------------------------------------------
!-------------------------------------------------------------
        implicit none

        integer,  intent(in)  :: ncol      ! number of active columns
        real, intent(in)      :: t(pcols,pver)     ! temperature
        real, intent(out)     :: fice(pcols,pver)  ! Fractional ice content within cloud
        real, intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

        integer :: i,k
        real    :: delt

        do k=1,pver
        do i=1,ncol
! copy from cldprpi.F          
          delt = 263.16 - t(i,k)
!--- fraction of ice
          fice(i,k) = amin1(1.0, amax1(0.0e0, delt/30.0e0 ))
          delt = 273.16 - t(i,k)
!--- fraction of snow in convection 
          fsnow(i,k) = amin1(1.0, amax1(0.0e0, delt/5.0e0 ))
        end do  
        end do  

      end subroutine cldwat_fice

      end module cldwat
