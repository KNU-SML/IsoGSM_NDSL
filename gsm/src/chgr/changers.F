#define IJCAP ijcap
#define ILATG ilatg
#define ILEVS ilevs
#define ILONF ilonf
#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LTSTP ltstp
#define LVAR lvar
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      program changers
c$$$  main program documentation block
c                .      .    .                                       .
c main program:  chgr        change resolution of mrf sig and sfc files.
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15
c
c abstract: changes the horizontal and vertical resolution
c           of the spectral model sigma and surface files.
c           for the sigma file, first the input data
c           are truncated or extended (with zero coefficients)
c           to the specified output horizontal resolution.
c           then the data are transformed to gridspace.
c           there the surface pressure is interpolated to the
c           specified output orography.  the temperature lapse rate
c           in this calculation is computed from the data
c           except below ground where it is set to -6.5k/km.
c           if necessary, the input moisture profile is extended
c           to the full input model layers using climatology.
c           then the wind, temperature and moisture profiles
c           are interpolated to the pressure locations determined
c           from the output surface pressure and the specified
c           output sigma layer values.  the profiles are never
c           extrapolated but held constant outside the input domain,
c           except below the input surface the temperature lapse rate
c           is again fixed at -6.5k/km and the relative humidity
c           is held fixed instead of the specific humidity.
c           the data are then transformed back to spectral space.
c           for the surface file, the input data are interpolated
c           horizontally to the output grid.  the output sea-land mask
c           is not interpolated but specified.
c
c program history log:
c   91-03-15  mark iredell  docblock written (prehistorical program)
c   01-03-31  hann-ming juang  add w3tag calls for nco implementation
c
c usage:
c   input files:
c     unit05   - namelist namchgr
c                ltrn     = .true. for ibm order files
c                lnew     = .true. for newer moisture climatology
c                lrep     = .true. to use orography from unit12
c     unit11   - input sigma file
c                (if empty, sigma file interpolation is skipped)
c     unit12   - output spectral orography
c     unit13   - output sea-land mask
c     unit14   - input surface file
c                (if empty, surface file interpolation is skipped)
c
c   output files:
c     unit06   - printed output and diagnostics
c     unit51   - output sigma file
c     unit52   - output surface file
c
c   subprograms called:
c     newsig   - get output sigma values
c     glats    - get gaussian latitudes
c     epslon   - compute spectral constants
c     greadi   - read input sigma file
c     chgres   - change horizontal spectral resolution
c     cmpind   - compute indices for transposing
c     transp   - transpose spectral coefficients
c     rmsgt    - print statistics
c     extrap   - transform to grid and interpolate vertically
c     transp   - transpose spectral coefficients
c     gwrite   - write output sigma file
c     gl2gl    - gaussian grid to gaussian grid interpolation
c
c attributes:
c   language: fortran 77.
c   machine:
c
c$$$
c
#include <paramodel.h>
      parameter(idim=lonf_,jdim=latg_,mwave=jcap_,
     1          kdim=levs_,kdimq=levs_*ntotal_,kdims=levs_,
     2          nvadm=lvar_,nptdm=lpnt_,nstdm=ltstp_,
     3          ndiag=ldiag_,ngrid=lgrid_)
c
      parameter(idimt=idim*2,jdimhf=jdim/2,ijdim=idim*jdim,
     1          mdim=(mwave+1)*(mwave+2),
     3          mdimv=(mwave+1)*(mwave+4),
     4          mwvp2=(mwave+1)*2,
     5          mdimhf=mdim/2,mdimvh=mdimv/2,
     6          kdimp=kdim+1,kdimm=kdim-1)
c
#include <parmchgr.h>
c
      parameter(idimti=idimi*2,jdimhfi=jdimi/2,ijdimi=idimi*jdimi,
     1          mdimi=(mwavei+1)*(mwavei+2),
     3          mdimvi=(mwavei+1)*(mwavei+4),
     4          mwvp2i=(mwavei+1)*2,
     5          mdimhfi=mdimi/2,mdimvhi=mdimvi/2,
     6          kdimpi=kdimi+1,kdimmi=kdimi-1)
c
#include <comchgr.h>
c hmhj 2
      common/inpver/ siin(kdimpi),slin(kdimi),delin(kdimi),
     1               ciin(kdimpi),clin(kdimi),rpiin(kdimmi)
     2               ,ak(kdimpi),bk(kdimpi)                                     ! 20070501(12)
c hmhj1
      common/comver/ si  (kdimp ),sl  (kdim ),del  (kdim ),
     1               ci  (kdimp ),cl  (kdim ),rpi  (kdimm),
     $               ak5 (kdimp ),bk5 (kdimp)
c
      common/plncom/eps   (mdimv ),
     1              colrad(jdimhf),wgt(jdimhf ),
     2              wgtcs (jdimhf),rcs2(jdimhf)
c
      common /label/ lab(8)
      character*4 lab
c
      real, allocatable :: ggridi(:,:),ggrido(:,:)
      dimension slmsko(ijdim)
c
      dimension ensemble(2),dummy2(21)
c
      logical ltrn,lnew,lrep
      character*4 sfcftyp
c   
      dimension lev(1000)
      character*16 svar(1000)
c
      namelist/namchgr/ nin,nou,ltrn,lnew,lrep,sfcftyp
c
      data nin,norog,nslmsk,nsfcin/11,12,13,14/
      data nou,nsfcou/51,52/
      data ltrn/.true./
      data lnew/.true./
      data lrep/.true./
c
      read (5,namchgr)
      write(6,namchgr)
c
c  define constants
c
      call glats (jdimhf, colrad, wgt, wgtcs, rcs2)
      call epslon(eps,mwave)
c
c  read in original initial data
c
      call greadi(nin,fhour,idate,gzi,qi,tei,dii,zei,rqi
     &       ,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
     &       ,iret)
      if(iret.ne.0) then
        print *,'hit end of file reading first record of sig file'
        call abort
      endif
c
      print 1101,fhour,idate
 1101 format(1h0,1x,'original initial data read in',/,
     1           3x,'fhour=',f6.1,' ihour=',i2,' month=',i2,' day=',i2,
     2              ' year=',i4)
      print 1102,mwavei,idimi,jdimi,kdimi,kdimqi
 1102 format(1x,'input  file resolution =',5i10)
#ifndef HYBRID      
      call newsig(ci,si,del,sl,cl,rpi)
#else
      call newhyb(ak5,bk5,ci,si,del,sl,cl,rpi)
      do k=1,kdimp
        print 1002, k,ak5(k),bk5(k)
      enddo
 1002 format (1x,'ak bk',i4,2f15.7)
#endif
      print 1103,mwave ,idim ,jdim ,kdim ,kdimq
 1103 format(1x,'output file resolution =',5i10)
c
c  change horizontal resolution by padding 0 or truncation
c
      call chgres(mwavei, qi,mwave, q)
      call chgres(mwavei,gzi,mwave,gz)
      do k=1,kdimi
        call chgres(mwavei,tei(1,k),mwave,te(1,k))
        call chgres(mwavei,dii(1,k),mwave,di(1,k))
        call chgres(mwavei,zei(1,k),mwave,ze(1,k))
      enddo
      do k=1,kdimqi
        call chgres(mwavei,rqi(1,k),mwave,rq(1,k))
      enddo
c
#ifdef CHGR_SIMPLE
         do m=1,mdim
           gzo(m)=gz(m)
           qo (m)=q (m)
         enddo
         if(lrep) then
           read(norog) gzo
         endif
         do k=1,kdimi
           do m=1,mdim
             teo(m,k)=te(m,k)
             dio(m,k)=di(m,k)
             zeo(m,k)=ze(m,k)
             rqo(m,k)=rq(m,k)
           enddo
         enddo
         ltrn=.false.
#else
        call cmpind
        if(ltrn) then
          call transp(gz,1)
          call transp( q,1)
          call transp(te,kdimi)
          call transp(di,kdimi)
          call transp(ze,kdimi)
          call transp(rq,kdimqi)
        endif
        if(lrep) then
          read(norog) gzo
          print *,'gzo of output resolution read in'
          call transp(gzo,1)
        else
          do m=1,mdim
            gzo(m)=gz(m)
          enddo
        endif
        call rmsgti(gz,q,di,te,ze,delin,rq)
        call extrap(lnew)
#endif
c
      call rmsgto(gzo,qo,dio,teo,zeo,del,rqo)
#ifndef HYBRID
      call gwrite(nou,fhour,idate,qo,teo,dio,zeo,rqo,sl,si,gzo,ltrn
#else
      call gwrite(nou,fhour,idate,qo,teo,dio,zeo,rqo,ak5,bk5,gzo,ltrn
#endif
     &       ,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
     &       ,iret)
c
  202 continue
c
c  addition for the change resolution of surface file
c  note that this program uses only bilinear interpol;ation, disregarding
c  whether is it land-ocean mask related or it is index of types (e.g. 
c  vegetaion index).  It is assumed that sfc program is run to overwrite
c  these inproperly interpolated fields from higher resolution climatology.
c  In therms of scripts, SFC0=yes should be specified.
c
      call sfcfld(sfcftyp,1,nrecs,lev,svar,mxlv)
      allocate (ggridi(ijdimi,mxlv),ggrido(ijdim,mxlv))
c
      rewind nsfcin
      read (nsfcin,end=800,err=800) lab
      print *,' end read lab of sfc.'
      write(nsfcou) lab
      read (nsfcin) fhour,idate
      print *,' read fhour,idate of sfc.'
      print *,'fhour,idate=',fhour,idate
      write(nsfcou) fhour,idate
c
      do nfld = 1,nrecs
        read(nsfcin,end=909,err=909) 
     1             ((ggridi(ij,l),ij=1,ijdimi),l=1,lev(nfld))
        print *,' read gaussian grid sfc. rec =',nfld
c
c  special handling of land sea mask
c   interpolate sea ice mask first.  Other areas, use input mask (from mtn).
c
        indslm=11
        if(sfcftyp.eq.'ncp1') indslm=1
        if(nfld.eq.indslm) then
          do i=1,ijdimi
            if(ggridi(i,1).eq.1.) then
              ggridi(i,1)=0.
            endif
          enddo
#ifdef DBG
          print *,'gl2gl (1)'
          call maxmn(ggridi,idimi,jdimi,1,'ggrdi')
#endif
          call gl2gl(ggridi,idimi,jdimi,ggrido,idim,jdim)
          read(nslmsk) slmsko
          do i=1,ijdim
            if(ggrido(i,1).gt.1..and.slmsko(i).eq.0.) then
              ggrido(i,1)=2.
            else
              ggrido(i,1)=slmsko(i)
            endif
          enddo
          call maxmn(ggrido,idim,jdim,1,'ggrdo')
        else
          do l=1,lev(nfld)
            call gl2gl(ggridi(1,l),idimi,jdimi,ggrido(1,l),idim,jdim)
            call maxmn(ggrido(1,l),idim,jdim,1,'ggrido')
          enddo
        endif
        write(nsfcou) ((ggrido(ij,l),ij=1,ijdim),l=1,lev(nfld))
      enddo
      stop
  909 continue
      print *,'hit eof or error reading sfc file body'
      call abort
  800 continue
      print *,'warning warning warning -- surface file empty'
      deallocate (ggridi,ggrido)
      stop
      end
