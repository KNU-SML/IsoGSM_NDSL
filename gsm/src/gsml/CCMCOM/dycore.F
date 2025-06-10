      module dycore
!-------------------------------------------------------------

        implicit none

        private
        public :: dycore_is, get_resolution

      contains

!-------------------------------------------------------------
      logical function dycore_is (name)
 
        character(len=*), intent(in) :: name
!        dycore_is = .true.         
        dycore_is = .false.         
        return

      end function dycore_is

!-------------------------------------------------------------
      character(len=7) function get_resolution()

        get_resolution = 'T62'
        return

      end function get_resolution

      end module dycore
