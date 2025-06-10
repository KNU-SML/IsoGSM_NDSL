#define CP cp
#define G g
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define RD rd
#define RVRDM1 rvrdm1
#define RV rv
#include <define.h>
      subroutine kf2cnv(imx,imx2,kmx,gu0,gv0,gt0,gq0,si,sl,del,psexp, 
     1              ncloud,clw,cli,
     1              vvel,deltim,dx,rain1,kbot,ktop,kuo,lat)
!
!------------------ --------------------------------------------------
!
!     the kain-fritsch convective scheme is adapted from wrfv2
!               implemeted by b.-j. jung, 2005 spring 
!                                                         
!---------------------------------------------------------------------
!
      implicit none
      integer, intent(in) :: lat, imx,imx2, kmx, ncloud
      integer, dimension(imx2), intent(out) :: kbot, ktop, kuo
      real, intent(in) :: deltim, dx
      real, dimension(imx2,kmx), intent(in) :: gu0, gv0, vvel
      real, dimension(imx2,kmx), intent(inout) :: gt0, gq0
      real, dimension(imx2,kmx), intent(out) :: clw,cli
      real, dimension(imx2,kmx), intent(in) :: sl, del
      real, dimension(imx2,kmx+1), intent(in) :: si
      real, dimension(imx2), intent(out) :: rain1
      real, dimension(imx2), intent(in) :: psexp
#include <paramodel.h>
c local array
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
!!!!!!!!!! local variables
      real, dimension(levs_) :: gu01d, gv01d, gt01d, gq01d, w01d, 
     1               p0, rho1d, tv, dzq
      real, dimension(levs_) :: dqdt, dqidt, dqcdt, dqrdt, dqsdt, dtdt
      real, dimension(levs_+1) :: zi 
      real :: dxsq, dlnsig, dt2
      real, dimension(ILOTS,1) :: raincv, nca
      integer :: i, k, p_qr, p_qi, p_qs, p_first_scalar, ntst
      integer :: cldtop, cldbot
      real :: xlv0,xlv1,xls0,xls1
      real :: cp,r,g,ep1,ep2
      real :: svp1,svp2,svp3,svpt0
      logical :: warm_rain
!!!!!!!!!! common variables
      integer, parameter :: kfnt=250,kfnp=220
      real, dimension(kfnt,kfnp) :: ttab,qstab
      real, dimension(kfnp) :: the0k
      real, dimension(200) :: alu
      real :: rdpr,rdthk,plutop
      common/kflut/ ttab,qstab,the0k,
     *              alu,rdpr,rdthk,plutop
!!!!!!!!!! end declare_variables
      parameter(cp=_cp_)
      parameter(r=_rd_)
      parameter(g=_g_)
      parameter(ep1=_rvrdm1_)
      parameter(ep2=_rd_/_rv_)
      svp1  =  0.6112
      svp2  =  17.67
      svp3  =  29.65
      svpt0 =  273.15
      xlv0  =  3.15e6
      xlv1  =  2370.
      xls0  =  2.905e6
      xls1  =  259.532
      p_qr  =  3
      p_qi  =  4
      p_qs  =  5
      if (ncloud.ge.4) then
        p_first_scalar  =  4
      else
        p_first_scalar  =  9
      endif
      ntst  =  1
#ifdef ICE
      warm_rain = .false. 
#else
      warm_rain = .true. 
#endif
      dxsq = dx * dx
      dt2 = deltim * 2.
!
      call kf_lutab(svp1,svp2,svp3,svpt0)
!
!!!!!!!!!! big-loop start
!
      do i = 1, imx
        do k = 1, kmx
          gu01d(k) = gu0(i,k)
          gv01d(k) = gv0(i,k)
          gt01d(k) = gt0(i,k)
          gq01d(k) = gq0(i,k)
          p0(k) = sl(i,k) * psexp(i) * 1000.
          tv(k) = gt0(i,k) * (1.+ep1*gq0(i,k))
          rho1d(k) = p0(k) / r / tv(k) 
          w01d(k) = - vvel(i,k) / g / rho1d(k) * 1000.
        enddo
!
        zi(1)=0.
        do k = 2, kmx
          dlnsig = log( si(i,k) / si(i,k-1) )
          zi(k) = zi(k-1) - dlnsig*r/g*tv(k-1)
        enddo
        zi(kmx+1) = zi(kmx) + tv(kmx)*r/g*del(i,kmx)/sl(i,kmx)
        do k = 1, kmx
          dzq(k) = zi(k+1) - zi(k)
        enddo
!
        do k = 1, kmx
          dqdt(k)  = 0.
          dqidt(k) = 0.
          dqcdt(k) = 0.
          dqrdt(k) = 0.
          dqsdt(k) = 0.
          dtdt(k)  = 0.
!     assumes initial value of cloud liquid/ice water to be zero
!     i.e. no environmental liquid/ice water used!!!
          clw(i,k)=0.0  
          cli(i,k)=0.0
        enddo
        rain1(i) = 0.
        raincv(i,1) = rain1(i)
!
!!!!!!!!!! call kf2-scheme with single column (~k) !!!!!!!!!!
!
        call kf_eta_para(i, lat,
     1           gu01d, gv01d, gt01d, gq01d, p0, dzq, w01d,
     2           dt2, dx, dxsq, rho1d,
     3           xlv0, xlv1, xls0, xls1, cp, r, g,
     4           ep2, svp1, svp2, svp3, svpt0,
     5           dqdt, dqidt, dqcdt, dqrdt, dqsdt, dtdt, raincv,
     6           nca, ntst,
     7           p_qi, p_qs, p_first_scalar, warm_rain,
     8           1, imx, lat, lat, 1, kmx,
     9           1, imx, lat, lat, 1, kmx,
     1           1, imx, lat, lat, 1, kmx,
     2           cldtop, cldbot)
!
!!!!!!!!!! post-process
!
        do k = 1, kmx
          gq0(i,k) = gq0(i,k) + dqdt(k) * dt2
          gt0(i,k) = gt0(i,k) + dtdt(k) * dt2
        enddo
        if (ncloud.ge.1) then  ! cloud water(/ice) 
          do k = 1, kmx
            clw(i,k) = clw(i,k) + dqcdt(k) * dt2
          enddo
        endif
        if (ncloud.ge.4) then  ! cloud ice
          do k = 1, kmx
            cli(i,k) = cli(i,k) + dqidt(k) * dt2
          enddo
        endif
        rain1(i) = raincv(i,1) / 1000.
        ktop(i) = cldtop
        kbot(i) = cldbot
        if (rain1(i) .gt. 0.) then
          kuo(i) = 1
        else
          kuo(i) = 0
        endif 
!
      enddo
!
      end subroutine kf2cnv
