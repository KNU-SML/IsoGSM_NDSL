      module physics_types
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      use ppgrid,       only: pcols, pver

      implicit none
      private          ! Make default type private to the module

      public physics_state
      public physics_ptend

!-----------------------------------------------------------------
      type physics_state

        integer :: lchnk              ! chunk index
        integer :: ncol               ! number of active columns
        real    :: phis(pcols)        ! surface geopotential
        real    :: t(pcols,pver)      ! temperature (K)
        real    :: pmid(pcols,pver)   ! midpoint pressure (Pa)
        real    :: pdel(pcols,pver)   ! layer thickness (Pa)
        real    :: zm(pcols,pver)     ! geopotential height above surface at midpoints (m)
        real    :: q(pcols,pver,1)    ! constituent mixing ratio (kg/kg)
        real    :: pint(pcols,pver+1) ! interface pressure (Pa)
        real    :: zi(pcols,pver+1)   ! geopotential height above surface at interfaces (m)


      end type physics_state
!-----------------------------------------------------------------
      type physics_ptend

        character*24 :: name ! name of parameterization which produced tendencies.
        logical :: ls      ! true if dsdt is returned
!!        logical :: lu      ! true if dudt is returned
!!        logical :: lv      ! true if dvdt is returned
        logical :: lq(1)   ! true if dqdt() is returned

        real    :: s(pcols,pver)    ! heating rate (J/kg/s)
!!        real    :: u(pcols,pver)    ! u momentum tendency (m/s/s)
!!        real    :: v(pcols,pver)    ! v momentum tendency (m/s/s)
        real    :: q(pcols,pver,1)  ! consituent tendencies (kg/kg/s)

      end type physics_ptend

      end module physics_types

