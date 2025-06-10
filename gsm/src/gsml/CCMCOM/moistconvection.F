#define CP cp
#define G g
#define RD rd
#define LEVS levs

       module moistconvection
!--------------------------------------------------------------
!--------------------------------------------------------------
#include <paramodel.h>

        real,parameter :: cp     = cp_     ! = cpres = cpair
        real,parameter :: grav   = g_      ! = gravit
        real,parameter :: rgrav  = 1./grav
        real,parameter :: rgas   = rd_     ! gas constant for dry air
!!        limcnv : caluculate like physics/cam1/convect_deep.F9
!         Limit deep convection to regions below 40 mb
!!        integer limcnv     ! top interface level limit for convection
        integer,parameter :: limcnv = 2 ! top interface level limit for convection

       end module moistconvection
