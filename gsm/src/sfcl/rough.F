      subroutine rough(vegtyp,ijmax,z0)
#include <paramodel.h>
c
c  compute roughness from vegetation type (13 types)
c
c
c   Input:
c      vegtyp:  real array of ijmax.  vegetation type.
c      ijmax:  dimension of vegtype
c
c   Output:
c 
c      rough:  real array of ijmax.  surface roughness.
c     
      dimension vegtyp(ijmax),z0(ijmax)
      dimension ri(13)
c
      save
c
c  the clim values
c
c     data ri/0.05,0.5,0.2,0.5,0.5,0.1,0.5,2.7,1.1,0.8,
c    1        0.85,0.6,0.001/
c
c  mm5 values
c     data RI/10.,0.05,0.12,0.13,0.2,0.1,0.1,0.5,0.5,0.4,
c    1        0.5,0.5,0.005/
c
c   Fein Chen's values
      data RI/0.01,0.075,0.238,0.04,0.05,0.11,0.11,2.654,1.089
     1       ,0.81,0.832,0.826,0.11/
c
      do ij=1,ijmax
          iv=vegtyp(ij)+1
          if ( iv .gt. 13 ) then
            print *,'vegetation type greater than 13'
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
          z0(ij)=ri(iv)*100.
        enddo
c
      return
      end
             
