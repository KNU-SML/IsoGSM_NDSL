#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define NSOIL nsoil
      subroutine viclandck(idim,jdim)

*-----------------------------------------------------------------------
*  SUBROUTINE DOCUMENTATION BLOCK
*  SUBPROGRAM viclandck.F
*
*  PRGMMR: Ji Chen   
*  ORG: ECPC/CRD/SIO/UCSD DATE: January 2004
*
*  ABSTRACT: check the land surface with VIC parameters according to 
*            non zero soil density
*-----------------------------------------------------------------------

      implicit none

#include <paramodel.h>
#include <varsfc.h>
#include <comsfcmrg.h>
c
      integer idim, jdim
c
#ifdef VICLSM1

      integer ijdim

      integer nvicland(idim*jdim), landnp(idim,jdim)
      integer landvic(idim,jdim),landsli(idim,jdim)

      integer m, n, np, nf, iv
      integer i, j, ij, k, kk
      integer ii1,jj1,ii21,ii22,jj21,jj22
c
c--- end of variable definition ---
c
#ifdef DBG
      open(1,file='check_vic_Ws.dat',status='unknown')
#endif

      ijdim = idim*jdim
      m = 0                                  ! total number of land grids
      n = 0
      do ij = 1, ijdim
         if(sfcfcs(ij,jsli).eq.1) then       ! land
            m = m + 1
            if(sfcfcs(ij,jsld).le.0.0) then   ! soil density
               n = n + 1
               nvicland(n) = ij             
            end if
         end if
      end do

#ifdef DBG
      close(1)
#endif

      if(n.gt.0) then
         do ij = 1, ijdim
            i = mod(ij,idim)
            if(i.eq.0) i = idim
            j = (ij-i)/idim + 1

            if(sfcfcs(ij,jsld).gt.0.0) then
               landvic(i,j) = 1
            else
               landvic(i,j) = 0
            end if
            if(sfcfcs(ij,jsli).eq.1.) then
               landsli(i,j) = 1
            else
               landsli(i,j) = 0
            end if
         end do
c
c         write(*, *) 'start to fill out VIC grids with sli land'
c
         do i = 1, idim
            do j =1, jdim
               k = 0
               if(landsli(i,j).eq.1.and.landvic(i,j).eq.0)then
                  ii1 = i
                  jj1 = j
                  kk = 1

 55               ii21 = ii1 + kk
                  if(ii21.ge.1.and.ii21.le.idim) then
                     if(landvic(ii21,jj1).eq.1) then
                        k = (jj1 - 1)*idim + ii21
                     end if
                  end if
                  if(k.eq.0) then
                     ii22 = ii1 - kk
                     if(ii22.ge.1.and.ii22.le.idim) then
                        if(landvic(ii22,jj1).eq.1) then
                           k = (jj1 - 1)*idim + ii22
                        end if
                     end if
                  end if
                  if(k.eq.0) then
                     jj21 = jj1 + kk
                     if(jj21.ge.1.and.jj21.le.jdim) then
                        if(landvic(ii1,jj21).eq.1) then
                           k = (jj21 - 1)*idim + ii1
                        end if
                     end if
                  end if
                  if(k.eq.0) then
                     jj22 = jj1 - kk
                     if(jj22.ge.1.and.jj22.le.jdim) then
                        if(landvic(ii1,jj22).eq.1) then
                           k = (jj22 - 1)*idim + ii1
                        end if
                     end if
                  end if

                  if(k.eq.0) then
                     if(ii21.ge.1.and.ii21.le.idim.and.
     &                  jj21.ge.1.and.jj21.le.jdim) then
                        if(landvic(ii21,jj21).eq.1) then
                           k = (jj21 - 1)*idim + ii21
                        end if
                     end if
                  end if
                  if(k.eq.0) then
                     if(ii22.ge.1.and.ii22.le.idim.and.
     &                  jj21.ge.1.and.jj21.le.jdim) then
                        if(landvic(ii22,jj21).eq.1) then
                           k = (jj21 - 1)*idim + ii22
                        end if
                     end if
                  end if
                  if(k.eq.0) then
                     if(ii21.ge.1.and.ii21.le.idim.and.
     &                  jj22.ge.1.and.jj22.le.jdim) then
                        if(landvic(ii21,jj22).eq.1) then
                           k = (jj22 - 1)*idim + ii21
                        end if
                     end if
                  end if
                  if(k.eq.0) then
                     if(ii22.ge.1.and.ii22.le.idim.and.
     &                  jj22.ge.1.and.jj22.le.jdim) then
                        if(landvic(ii22,jj22).eq.1) then
                           k = (jj22 - 1)*idim + ii22
                        end if
                     end if
                  end if
                  if(k.eq.0) then
                     if((ii21.ge.1.and.ii21.le.idim).or.
     &                  (ii22.ge.1.and.ii22.le.idim).or.
     &                  (jj21.ge.1.and.jj21.le.jdim).or.
     &                  (jj22.ge.1.and.jj22.le.jdim)) then
                        kk = kk + 1
                        goto 55
                     else
                        write(*,*) 'no vic land found for sli land'
#ifdef MP
#ifdef RSM
                        call rmpabort
#else
                        call mpabort
#endif
#else
                        call abort
#endif
                     end if
                  end if
               end if
               landnp(i,j) = k
            end do
         end do

#ifdef DBG
         open(77,file='VIC_land.dat',status='unknown')
         open(88,file='SLI_land.dat',status='unknown')
         open(99,file='VIC_SLI_land.dat',status='unknown')
         do j = 1, jdim
            write(77,177) (landvic(i,j),i=1,idim)
            write(88,177) (landsli(i,j),i=1,idim)
            write(99,199) (landnp(i,j),i=1,idim)
 177        format(721i1)
 199        format(721i7)
         end do
         close(77)
         close(88)
         close(99)
         open(66,file='VIC_noland.dat',status='unknown')
         write(66,*) 'the number of land grids without VIC para ',n
#endif
         do ij = 1, n
            i = mod(nvicland(ij),idim)
            if(i.eq.0) i = idim
            j = (nvicland(ij)-i)/idim + 1
#ifdef DBG
            write(66,166) nvicland(ij), landnp(i,j)
 166        format(2i15)
#endif
            np = nvicland(ij)
            nf = landnp(i,j)

            sfcfcs(np,jveg)=sfcfcs(nf,jveg)
            sfcfcs(np,jcpy)=sfcfcs(nf,jcpy)
            sfcfcs(np,jvet)=sfcfcs(nf,jvet)
            sfcfcs(np,jbif)=sfcfcs(nf,jbif)
            sfcfcs(np,jds )=sfcfcs(nf,jds )
            sfcfcs(np,jdsm)=sfcfcs(nf,jdsm)
            sfcfcs(np,jws )=sfcfcs(nf,jws )
            sfcfcs(np,jcef)=sfcfcs(nf,jcef)
            sfcfcs(np,jlai)=sfcfcs(nf,jlai)
            sfcfcs(np,jslz)=sfcfcs(nf,jslz)
            sfcfcs(np,jsnz)=sfcfcs(nf,jsnz)

            do m = 1, lsoil_
               sfcfcs(np,jrot+m-1)=sfcfcs(nf,jrot+m-1)
               sfcfcs(np,jexp+m-1)=sfcfcs(nf,jexp+m-1)
               sfcfcs(np,jkst+m-1)=sfcfcs(nf,jkst+m-1)
               sfcfcs(np,jdph+m-1)=sfcfcs(nf,jdph+m-1)
               sfcfcs(np,jbub+m-1)=sfcfcs(nf,jbub+m-1)
               sfcfcs(np,jqrt+m-1)=sfcfcs(nf,jqrt+m-1)
               sfcfcs(np,jbkd+m-1)=sfcfcs(nf,jbkd+m-1)
               sfcfcs(np,jsld+m-1)=sfcfcs(nf,jsld+m-1)
               sfcfcs(np,jwcr+m-1)=sfcfcs(nf,jwcr+m-1)
               sfcfcs(np,jwpw+m-1)=sfcfcs(nf,jwpw+m-1)
               sfcfcs(np,jsmr+m-1)=sfcfcs(nf,jsmr+m-1)
               sfcfcs(np,jsmx+m-1)=sfcfcs(nf,jsmx+m-1)
            enddo
            do n = 1, nsoil_
               sfcfcs(np,jdpn+n-1)=sfcfcs(nf,jdpn+n-1)
               sfcfcs(np,jsxn+n-1)=sfcfcs(nf,jsxn+n-1)
               sfcfcs(np,jepn+n-1)=sfcfcs(nf,jepn+n-1)
               sfcfcs(np,jbbn+n-1)=sfcfcs(nf,jbbn+n-1)
               sfcfcs(np,japn+n-1)=sfcfcs(nf,japn+n-1)
               sfcfcs(np,jbtn+n-1)=sfcfcs(nf,jbtn+n-1)
               sfcfcs(np,jgmn+n-1)=sfcfcs(nf,jgmn+n-1)
            enddo
         end do
#ifdef DBG
         close(66)
#endif
      end if
#endif
      return
      end
