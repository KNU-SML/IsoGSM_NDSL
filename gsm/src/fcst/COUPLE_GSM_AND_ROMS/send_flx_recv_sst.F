#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1 jgrd1
#define JGRD1P jgrd1p
#define JGRD jgrd
#define JWAV1 jwav1
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLWAVP llwavp
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define PI pi
#define TTP ttp
#include <define.h>
      subroutine send_flx_recv_sst
c
#include <paramodel.h>
#ifdef MP
#undef DBG
#undef DBG_FULL
#ifndef RSM
#include <comio.h>
#endif
c
#ifdef RSM
#include <rscomf.rerun.h>
#include <rscommap.h>
#ifdef MP
#include <comrmpi.h>
#define MPABORT rmpabort
#else
#define MPABORT abort
#endif
#else
#include <comfphys.h>
#include <comfspec.vr.h>
#include <comfver.h>
#include <comfgrid.h>
#include <radiag.h>
#ifdef MP
#include <commpi.h>
#define MPABORT mpabort
#else
#define MPABORT abort
#endif
#endif
c
#ifdef MP
#ifdef RSM
#define MPGP2F rmpgp2f
#define ILEN2S igrd12p_
#define JLEN2S jgrd12p_
#define ILEN igrd1_
#define ILEN2 igrd12_
#define JLEN jgrd1_
#define JLEN2 jgrd12_
#else
#define MPGP2F mpgp2f
#define ILEN2S lonf2p_
#define JLEN2S latg2p_
#define ILEN lonf_
#define ILEN2 lonf2_
#define JLEN latg_
#define JLEN2 latg2_
#endif
#else
#ifdef RSM
#define ILEN2S igrd12_
#define JLEN2S jgrd12_
#define ILEN igrd1_
#define ILEN2 igrd12_
#define JLEN jgrd1_
#define JLEN2 jgrd12_
#else
#define ILEN2S lonf2_
#define JLEN2S latg2_
#define ILEN lonf_
#define ILEN2 lonf2_
#define JLEN latg_
#define JLEN2 latg2_
#endif
#endif
c
c  full arrays
c
      parameter(len=ILEN*JLEN)
      parameter(nradfld=19)
      parameter(lens=ILEN2S*JLEN2S)
c
c  full working arrays
c
      dimension work(len,4)
      dimension slmsep(len)
      real rtimer(nradfld)
c
c     integer, parameter :: send_blocklen=8192 
      integer, parameter :: send_blocklen=524288 
c
      integer, parameter :: send_max_numfld=17
c
      integer status(MPI_STATUS_SIZE)
c
      integer :: send_numfld,send_numblock
      integer, allocatable :: send_reclen(:)
      integer, dimension(5) :: send_header
      real, allocatable :: send_array(:,:)
      character*1, allocatable :: send_var(:,:)
c
      integer :: send_datalen
c
      integer :: recv_numfld,recv_numblock
      integer, allocatable :: recv_reclen(:)
      integer, dimension (5) :: recv_header(5)
      real, allocatable ::  recv_array(:,:)
c
      integer :: ind,ilen,jlen,n,i,j,ij
c
      real, allocatable :: fullon(:,:),fullat(:,:)
c
      real, allocatable :: tseax(:,:)
c
      integer maxnumfld
      parameter( maxnumfld=12)
      character*8,var(maxnumfld)
      data var( 1)/'lon     '/
      data var( 2)/'lat     '/
      data var( 3)/'lsmask  '/
      data var( 4)/'ustr    '/
      data var( 5)/'vstr    '/
      data var( 6)/'shflx   '/
      data var( 7)/'lhflx   '/
      data var( 8)/'dlwfsfc '/
      data var( 9)/'ulwfsfc '/
      data var(10)/'uswfsfc '/
      data var(11)/'dswfsfc '/
      data var(12)/'prate   '/
c
c     parameter( maxnumfld=16)
c     data var(13)/.f.,'u10m    '/
c     data var(14)/.f.,'v10m    '/
c     data var(15)/.f.,'t2m     '/
c     data var(16)/.f.,'q2m     '/
c
      send_numfld=maxnumfld
      allocate (send_array(len,send_numfld))
      allocate (send_var(8,send_numfld))
c
      send_datalen=len*send_numfld
      send_numblock=(send_datalen/send_blocklen)+1
      if(mod(send_datalen,send_blocklen).eq.0) then
        send_numblock=send_numblock-1
      endif
      allocate (send_reclen(send_numblock))
      do n=1,send_numblock
        send_reclen(n)=send_blocklen
        if(send_blocklen*n.gt.send_datalen) then
            send_reclen(n)=send_datalen-
     1                        send_blocklen*(send_numblock-1)
        endif
      enddo
c
      iyr=idate(4)
      imo=idate(2)
      ida=idate(3)
      ihr=idate(1)
      ifhr=nint(fhour)
      ithr=nint(thour)
      if(ithr.eq.0) ither=1
c
      dhour=dtpost/3600.
      if(dtpost.gt.0) then
        rtime=1./dtpost
      else
        rtime=0.
      endif
c
      secswr=max(dhour,dtswav) * 3600.
      seclwr=max(dhour,dtlwav) * 3600.
      if(secswr.gt.0.) then
        rtimsw=1./secswr
      else
        rtimsw=1.
      endif
      if(seclwr.gt.0.) then
        rtimlw=1./seclwr
      else
        rtimlw=1.
      endif
      do n=1,nradfld
      rtimer(n)=rtimsw
      enddo
      rtimer(1)=rtimlw
c
c  get lat/lon
c
      allocate (fullat(ILEN,JLEN),fullon(ILEN,JLEN))
#ifdef RSM
      call MPGP2F(flat,igrd12p_,jgrd12p_,fullat,igrd12_,jgrd12_,1)
      call MPGP2F(flon,igrd12p_,jgrd12p_,fullon,igrd12_,jgrd12_,1)
      rad=180./_pi_
      do j=1,JLEN
        do i=1,ILEN
          fullat(i,j)=fullat(i,j)*rad
          fullon(i,j)=fullon(i,j)*rad
          if(fullon(i,j).gt.180.) fullon(i,j)=fullon(i,j)-360.
        enddo
      enddo
c
c    rproj.eq.0       ... mercater
c    abs(rproj).eq.1  ... polar projection
c
#else
c
c gaussian longitude and latitude
c      
      do j=1,JLEN
        do i=1,ILEN
          fullon(i,j)=float(i-1)*360./float(ILEN)
          if(fullon(i,j).gt.180.) fullon(i,j)=fullon(i,j)-360.
        enddo
      enddo
      do j=1,JLEN/2
        do i=1,ILEN
           fullat(i,j)=90.-colrad(j)*pi_/180.
        enddo
      enddo
      do j=JLEN/2+1,JLEN
        jj=JLEN-j+1
        do i=1,ILEN
          fullat(i,j)=-fullat(i,jj)
        enddo
      enddo
#endif
c
c  get mask
c
#ifdef MP
      call MPGP2F(slmsk,ILEN2S,JLEN2S,slmsep,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#else
      do n=1,len
      slmsep(n)=slmsk(n,1)
      enddo
#endif
#ifdef RSM
      call shalfo(slmsep,1)
#else
      call rowsep(slmsep)
#endif
c
#ifdef MP
      endif
#endif
c
c 1. longitude
c
      ind=1
      do n=1,len
        send_array(n,ind)=fullon(n,1)
      enddo
      do i=1,8
        send_var(i,ind)=var(1)(i:i)
      enddo
      ind=ind+1
c
c 2. latitude
c
      do n=1,len
        send_array(n,ind)=fullat(n,1)
      enddo
      do i=1,8
        send_var(i,ind)=var(2)(i:i)
      enddo
      ind=ind+1
c
c 3  land sea mask
c
#ifdef MP
      if(mype.eq.master) then
#endif
        do n=1,len
          work(n,1)=mod(slmsep(n),2.)
        enddo
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'lsmask')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(3)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 4  surface u-stress
c
#ifdef MP
      call MPGP2F(dusfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
            work(n,1)=work(n,1)*rtime
#else
            work(n,1)=dusfc(n,1)*rtime
#endif
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'uflx')
#endif
        do n=1,len
           send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(4)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 5  surface v-stress
c
#ifdef MP
      call MPGP2F(dvsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
            work(n,1)=work(n,1)*rtime
#else
            work(n,1)=dvsfc(n,1)*rtime
#endif
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'vflx')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(5)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 6  sensible heat flux
c
#ifdef MP
      call MPGP2F(dtsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
          work(n,1)=work(n,1)*rtime
#else
          work(n,1)=dtsfc(n,1)*rtime
#endif
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'tflx')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(6)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 7  latent heat flux
c
#ifdef MP
      call MPGP2F(dqsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
          work(n,1)=work(n,1)*rtime
#else
          work(n,1)=dqsfc(n,1)*rtime
#endif
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'qflx')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(7)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 8 downward longwave radiation flux at the surface
c
#ifdef MP
      call MPGP2F(dlwsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
          work(n,1)=work(n,1)*rtime
#else
          work(n,1)=dlwsfc(n,1)*rtime
#endif
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'dlwsfc')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(8)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 9 upward longwave radiation flux at the surface
c
#ifdef MP
      call MPGP2F(ulwsfc,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
          work(n,1)=work(n,1)*rtime
#else
          work(n,1)=ulwsfc(n,1)*rtime
#endif
      enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'ulwsfc')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(9)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c   10  upward short wave radiation flux at the surface
c   11  downward short wave radiation flux at the surface
c
#ifdef MP
#ifdef RSMVECTORIZE
      call fillarray(fluxr,26)
#endif
#endif
      do k=3,4
        ii=0
#ifdef MP
        call MPGP2F(fluxr(1,1,k),ILEN2S,JLEN2S,work,
     1              ILEN2,JLEN2,1)
        if(mype.eq.master) then
#else
          do j=1,JLEN2
            do i=1,ILEN2
              ii=ii+1
              work(ii,1)=fluxr(i,j,k)
            enddo
          enddo
#endif
          do n=1,len
            work(n,1)=work(n,1)*rtimer(k)
          enddo
#ifdef RSM
          call shalfo(work,1)
#else
          call rowsep(work)
#endif
#ifdef DBG
          call maxmin(work,ILEN*JLEN,1,1,1,'radflxes')
#endif
          do n=1,len
            send_array(n,ind)=work(n,1)
          enddo
          do i=1,8
            send_var(i,ind)=var(7+k)(i:i)
          enddo
          ind=ind+1
#ifdef MP
        endif
#endif
      enddo
c
c 12  precipitation rate
c
#ifdef MP
      call MPGP2F(geshem,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
      if(mype.eq.master) then
#endif
        do n=1,len
#ifdef MP
          work(n,1)=work(n,1)*1.e3*rtime
#else
          work(n,1)=geshem(n,1)*1.e3*rtime
#endif
        enddo
#ifdef RSM
        call shalfo(work,1)
#else
        call rowsep(work)
#endif
#ifdef DBG
        call maxmin(work,ILEN*JLEN,1,1,1,'prate')
#endif
        do n=1,len
          send_array(n,ind)=work(n,1)
        enddo
        do i=1,8
          send_var(i,ind)=var(12)(i:i)
        enddo
        ind=ind+1
#ifdef MP
      endif
#endif
c
c 13  10m wind u
c
c#ifdef MP
c      call MPGP2F(u10m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
c      if(mype.eq.master) then
c#endif
c      do n=1,len
c#ifdef MP
c      work(n,1)=work(n,1)
c#else
c      work(n,1)=u10m(n,1)
c#endif
c      enddo
c#ifdef RSM
c      call shalfo(work,1)
c#else
c      call rowsep(work)
c#endif
c#ifdef DBG
c      call maxmin(work,ILEN*JLEN,1,1,1,'u10m')
c#endif
c      do n=1,len
c        send_array(n,ind)=work(n,1)
c      enddo
c      do i=1,8
c        send_var(i,ind)=var(13)(i:i)
c      enddo
c      ind=ind+1
c#ifdef MP
c      endif
c#endif
cc
cc 14  10m wind v
cc
c#ifdef MP
c      call MPGP2F(v10m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
c      if(mype.eq.master) then
c#endif
c      do n=1,len
c#ifdef MP
c      work(n,1)=work(n,1)
c#else
c      work(n,1)=v10m(n,1)
c#endif
c      enddo
c#ifdef RSM
c      call shalfo(work,1)
c#else
c      call rowsep(work)
c#endif
c#ifdef DBG
c      call maxmin(work,ILEN*JLEN,1,1,1,'v10m')
c#endif
c      do n=1,len
c        send_array(n,ind)=work(n,1)
c      enddo
c      do i=1,8
c        send_var(i,ind)=var(14)(i:i)
c      enddo
c      ind=ind+1
c#ifdef MP
c      endif
c#endif
cc
cc 15  temperature 2m
cc
c#ifdef MP
c      call MPGP2F(t2m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
c      if(mype.eq.master) then
c#endif
c      do n=1,len
c#ifdef MP
c      work(n,1)=work(n,1)
c#else
c      work(n,1)=t2m(n,1)
c#endif
c      enddo
c#ifdef RSM
c      call shalfo(work,1)
c#else
c      call rowsep(work)
c#endif
c      il2k=2
c#ifdef DBG
c      call maxmin(work,ILEN*JLEN,1,1,1,'t2m')
c#endif
c      do n=1,len
c        send_array(n,ind)=work(n,1)
c      enddo
c      do i=1,8
c        send_var(i,ind)=var(15)(i:i)
c      enddo
c      ind=ind+1
c#ifdef MP
c      endif
c#endif
cc
cc 16  specific humidity 2m
cc
c#ifdef MP
c      call MPGP2F(q2m,ILEN2S,JLEN2S,work,ILEN2,JLEN2,1)
c      if(mype.eq.master) then
c#endif
c      do n=1,len
c#ifdef MP
c      work(n,1)=work(n,1)
c#else
c      work(n,1)=q2m(n,1)
c#endif
c      enddo
c#ifdef RSM
c      call shalfo(work,1)
c#else
c      call rowsep(work)
c#endif
c      il2k=2
c#ifdef DBG
c      call maxmin(work,ILEN*JLEN,1,1,1,'q2m')
c#endif
c      do n=1,len
c        send_array(n,ind)=work(n,1)
c      enddo
c      do i=1,8
c        send_var(i,ind)=var(16)(i:i)
c      enddo
c      ind=ind+1
c#ifdef MP
c      endif
c#endif
c
      if(mype.eq.master) then
        if(ind-1.ne.send_numfld) then
          print *,'abnormal termination send_flx_recv_sst.'
          call MPABORT
        endif
      endif
c
c  send flux fields to ocean model
c
      call mpi_barrier(mpi_comm_world,ierr)
      if( ierr.ne.0 ) then
        print *,'PE',mype,':***** Error stop in send_flx_recv_sst **** '
        print *,'PE',mype,':error code from mpi_barrier = ',ierr
        call MPABORT
      endif
c
      if(mype.eq.master) then
        master_ocean=npes_
        call MPI_COMM_RANK(MPI_COMM_WORLD,mype_world,jerr)
#ifdef DBG_FULL
        print *,mype_world,' sending flx to ',master_ocean,
     1   ' in chunks of ',send_blocklen,
     2   ' for ',send_numblock,' times'
#endif
        send_header(1)=send_blocklen
        send_header(2)=send_numblock
        send_header(3)=ILEN
        send_header(4)=JLEN
        send_header(5)=send_numfld
#ifdef DBG_FULL
        print *,'send_header=',send_header
#endif
        call MPI_SSEND(send_header,5,MPI_INTEGER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,ierr)
        call MPI_SSEND(send_reclen,send_numblock,MPI_INTEGER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,ierr)
c
        do n=1,send_numblock
          ind=(n-1)*send_blocklen+1
#ifdef DBG_FULL
          print *,'sending flux block number=',n,
     1            ' length=',send_reclen(n)
#endif
          call MPI_SSEND(send_array(ind,1),send_reclen(n),MPI_REAL8,
     1                  master_ocean,msgtag,MPI_COMM_WORLD,ierr)
        enddo
c
        call MPI_SSEND(send_var,8*send_numfld,MPI_CHARACTER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,ierr)
c
c  receive sst from ocean
c
        call MPI_RECV(recv_header,5,MPI_INTEGER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,status,ierr)
#ifdef DBG_FULL
        print *,'recv_header=',recv_header
#endif
        recv_blocklen=recv_header(1)
        recv_numblock=recv_header(2)
        ilen_recv=recv_header(3)
        jlen_recv=recv_header(4)
        recv_numfld=recv_header(5)
#ifdef DBG
        print *,'recv_header(5)=recv_numfld=',recv_numfld
#endif
        allocate (recv_reclen(recv_numblock))
        allocate (recv_array(len,recv_numfld))
        call MPI_RECV(recv_reclen,recv_numblock,MPI_INTEGER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,status,ierr)
c
        if(ilen_recv.ne.ILEN.or.jlen_recv.ne.JLEN) then
          print *,'ilen_recv.ne.ILEN.or.jlen_recv.ne.JLEN'
          print *,'ilen_recv,jlen_recv,ILEN,JLEN=',
     &             ilen_recv,jlen_recv,ILEN,JLEN
          print *,'interpolation necessary'
          call MPABORT
        endif
c
        do n=1,recv_numblock
          ind=(n-1)*recv_blocklen+1
#ifdef DBG_FULL
          print *,'receiving sst block number=',n,
     1            ' length=',recv_reclen(n)
#endif
          call MPI_RECV(recv_array(ind,1),recv_reclen(n),MPI_REAL8,
     1                 master_ocean,msgtag,MPI_COMM_WORLD,status,ierr)
        enddo
c
        do n=1,recv_numfld
          call rmaxmin(recv_array(1,n),ILEN*JLEN,1,1,1,'tsea')
#ifndef RSM
          call row1ns(recv_array(1,n))
#else
          call shalfi(recv_array(1,n),1)
#endif
        enddo
c
c  change unit of sst from c to k.
c
        do i=1,ILEN*JLEN
          recv_array(i,1)=recv_array(i,1)+_ttp_
        enddo
      endif
c
      call mpi_barrier(mpi_comm_world,ierr)
      if( ierr.ne.0 ) then
        print *,'PE',mype,':***** Error stop in send_flx_recv_sst **** '
        print *,'PE',mype,':error code from mpi_barrier = ',ierr
        call MPABORT
      endif
#ifdef DBG_FULL
c 
      print *,'Atmos: Comm between atmos and ocean now complete'
#endif
c
c  spread to all the rsm nodes
c
#ifdef RSM
      call rmpbcasti(recv_numfld,1)
#else
      call mpbcasti(recv_numfld,1)
#endif
#ifdef DBG_FULL
      print *,'recv_numfld=',recv_numfld
#endif
c     
      if(mype.ne.master) then
        allocate (recv_array(len,recv_numfld))
      endif
      allocate (tseax(ILEN2S,JLEN2S))
c
#ifdef RSM
      call rmpsynall
      call rmpgf2p(recv_array,ILEN2,JLEN2,tseax,ILEN2S,JLEN2S,1)
#else
      call mpsynall
      call mpgf2p(recv_array,ILEN2,JLEN2,tseax,ILEN2S,JLEN2S,1)
#endif
c
c  replace rsm surface temperature only over ocean.
c
      do j=1,JLEN2S
        do i=1,ILEN2S
          if(slmsk(i,j).eq.0.) then
            tsea(i,j)=tseax(i,j)
          endif
        end do
      end do
c
      if(mype.eq.master) then
        deallocate (recv_reclen)
      endif
      deallocate (recv_array)
      deallocate (tseax)
c
      deallocate (send_array)
      deallocate (send_var)
      deallocate (send_reclen)
      deallocate (fullat,fullon)
c
#ifdef DBG
      print *,'mype=',mype,' send_flx_recv_sst finished'
#endif
#endif
      return
      end
