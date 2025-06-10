#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define RD rd
#define RHOH2O rhoh2o
#define RV rv
      subroutine initcommons
#include <paramodel.h>
#include <comadj.h>
C
      rair   = rd_
      cpair  = cp_
      cappa  = rair/cpair
      epsilo = rd_/rv_
      gravit = g_
      latvap = hvap_
      latice = hfus_
      rhoh2o = rhoh2o_
      clrh2o = latvap/rhoh2o
      cldcp  = latvap/cpair
      rh2o   = rv_
      zvir   = rh2o/rair - 1.

      call mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhoh2o  )
      call esinti (epsilo  ,latvap  ,latice  ,rh2o  ,cpair   )
      return
      end
