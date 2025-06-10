      module shr_kind_mod
!---------------------------------------------------
!     reference
!       cam3.1/cam1/models/atm/cam/tools/cprnc/shr_kind_mod.F
!---------------------------------------------------
      public
#ifndef sun
      integer,parameter :: shr_kind_r8 = selected_real_kind(12) ! 8 byte real
#endif

      end module shr_kind_mod
