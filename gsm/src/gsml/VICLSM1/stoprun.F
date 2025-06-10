      subroutine stoprun
#include <paramodel.h>

* ------------------------ code history --------------------------------
* source file:       stoprun.F
* purpose:           terminate job
* date:              November 2003
* ----------------------------------------------------------------------

#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
      end
