	subroutine nodepara(msl,    nsl,    dph,   dphnd,
     &                      smx,    bub,   expt,   smxnd,
     &                    expnd,  bubnd,  alpha,    beta,
     &                    gamma)

	implicit none

************************************************************************
* source file:       nodepara.F
* PRGMMR:            Ji Chen   
* ORG:               ECPC/CRD/SIO/UCSD 
* DATE:              November 2003
* PRGM HISTORY:      modified from UW_VIC
*                    (set_node_parameters in soil_conduction.c)
*
* ABSTRACT: This subroutine sets the thermal node soil parameters to 
*           constant values based on those defined for the current grid 
*           cells soil type. Thermal node propertiers for the energy 
*           balance solution are also set (these constants are used to 
*           reduce the solution time required within each iteration).
************************************************************************

* ------------------------ input variables -----------------------------
	integer  msl	      ! number of soil moisture layers
	integer  nsl          ! number of soil thermal nodes
	real     dph(msl)     ! soil moisture layer thickness (m)
	real     dphnd(nsl)   ! thermal node thicknes (m)
	real     smx(msl)     ! maximum layer moisture content (mm)
	real     bub(msl)     ! soil layer bubbling pressure (cm)
	real     expt(msl)    ! soil moisture exponential (N/A)
* ----------------------------------------------------------------------

* ------------------------ return variables ----------------------------
	real     smxnd(nsl)   ! maximum moisture thermal node (m3/m3)
	real     expnd(nsl)   ! exponential at thermal node (N/A)
	real     bubnd(nsl)   ! bubbling pressure at thermal node (cm)
	real     alpha(nsl)   ! first thermal eqn term
	real     beta(nsl)    ! second thermal eqn term
	real     gamma(nsl)   ! third thermal eqn term
*-----------------------------------------------------------------------

* ------------------------ local variables ---------------------------
	integer nidx, lidx    ! loop indexes
	real Lsum, Zsum

	logical Past_Bottom   ! check whether pass soil column bottom
* --------------------------------------------------------------------

* --------------------------------------------------------------------
* initialize variables
* --------------------------------------------------------------------
	Past_Bottom = .FALSE.
	lidx = 1
	Lsum = 0.
	Zsum = 0.

	do nidx = 1, nsl	                ! set node parameters
	   Zsum = Zsum + dphnd(nidx)
	   if(Zsum.gt.Lsum.and..not.Past_Bottom) then
	      Lsum = Lsum + dph(lidx)
	      lidx = lidx + 1
	      if( lidx .ge. msl ) then
		 Past_Bottom = .TRUE.
		 lidx = msl
	      end if

	      do while(Zsum.gt.Lsum)
		 Lsum = Lsum + dph(lidx)
		 lidx = lidx + 1
		 if( lidx .ge. msl ) then
		    Past_Bottom = .TRUE.
		    lidx = msl
		 end if
	      end do
	   end if

	   if(Zsum.eq.Lsum.and.nidx.ne.1.and.      ! node on layer boundary
     &        lidx.ne.1) then
	      smxnd(nidx) = (smx(lidx-1)/dph(lidx-1) +
     &                       smx(lidx)/dph(lidx))/2.0/1000.0
	      expnd(nidx) = (expt(lidx-1)+expt(lidx))/2.0
	      bubnd(nidx) = (bub(lidx-1) + bub(lidx))/2.0
	   else                 !node completely in layer
	      smxnd(nidx) = smx(lidx)/dph(lidx)/1000.0
	      expnd(nidx) = expt(lidx)
	      bubnd(nidx) = bub(lidx)
	   end if
	end do

* --------------------------------------------------------------------
* compute constant parameters for thermal calculations
* --------------------------------------------------------------------
	do nidx = 1, nsl
	   alpha(nidx) = 0.0
	   beta(nidx)  = 0.0
	   gamma(nidx) = 0.0
	end do

	do nidx = 1, nsl-1
	   alpha(nidx)= dphnd(nidx)+dphnd(nidx+1)
	   beta(nidx) = dphnd(nidx)*dphnd(nidx)+
     &                  dphnd(nidx+1)*dphnd(nidx+1)
	   gamma(nidx)= dphnd(nidx) - dphnd(nidx+1)
	end do

* no flux bottom boundary activated

	alpha(nsl) = 2 * dphnd(nsl)
	beta(nsl)  = 2 * dphnd(nsl)*dphnd(nsl)
	gamma(nsl) = 0.0

	return
	end
