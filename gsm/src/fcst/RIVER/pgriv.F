#define PI pi
#define IO2 io2
#define JO2 jo2
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAP22P lcap22p
#define LCAPP lcapp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNT2P lnt2p
#define LNTP lntp
#define LNUV lnuv
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NCPUS ncpus
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define RD rd
#define RV rv
#define RERTH rerth
#define TWOJ1 twoj1
c
      subroutine pgriv
c
c  river discharge estimation coded by kei yoshimura (07/02/09)
c
c  compute river discharge by TRIP (Oki and Sud, 1998)
c     S(t+1)=S(t)*exp(-dc*dt)+(1-exp(-dc*dt))*I/dc
c     O=(S(t+1)-S(t))/dt+I
c     S:river storage 
c     O:river discharge
c     I:input to the flow (inflow from upstream and total runoff 
c       accumulated in a grid)
c     dc=v/l (v: flow velocity, l: flowpath length)
c
c     input runoff are taken from LSM and output river discharge etc.
c     are written to riv.ftXX in grib format.
c
c  variables:
c     gdriv: river storage [kg/m/m]
c     rflow: river discharge [kg/s]
c     imap : flow direction (0:sea, 9:rivermouth)
c           8  1  2
c           7  *  3
c           6  5  4
c     roff : converted runoff [kg/m/m/s]
c
c  initial and restart file:
c     rivi (if not exsited, all river channel storage is assumed
c           to have 10mm water)
c
c  Note1: 
c     Calculation is in a single processor
c  Note2: 
c     Map(s) of riverflow direction is very important.
c     Currently, 1deg x 1deg version (Oki and Sud, 1998) and
c     0.5deg x 0.5deg version are available.
c     The map should be placed/linked as fort.93
c  Note3: 
c     io_ and jo_ in paramodel.h are used as io2_ and jo2_.
c     t126 (io_=360,jo_=181) and t248 (io_=720,jo_=361) are
c     usable in default, but other resolution io2_ and jo2_ 
c     needs to be selected manually. => in comfriv.h and gau2llc.F
c
c     Grid coordination (1deg) is centered, such as 
c     (0.5E, 89.5N)  - (0.5W, 89.5N) ; i=1,io_
c                    :               ; j=1,jo_-1
c     (0.5E, 89.5S)  - (0.5W, 89.5S) 
c     (        dummy column        ) ; j=jo_
c
c  included h file:
c     comfriv.h (in include directory)
c  subroutine called: 
c     gau2llc: gausian to latlon-centered cordination
c
c  influenced subroutine:
c     gsmstep: call this subroutine (pgriv) and
c              write restart file (rivi)
c     gbphys : calculate total runoff (surface+subsurface)
c     wriflx : writing to riv.ftXX in grib format.
c     idsdef : (in share directory) define output variable IDs.
c              gdriv(199), rflow(200), imap(202), roff(203)
c
**** LET'S GO!! ****
c
#include <paramodel.h>
c
#include <comfver.h>
c
#ifdef MP
#include <commpi.h>
#endif
c
#include <comfriv.h>
c
      parameter(len=lonf_*latg_)
      dimension work(len)
      real runof(io2_,jo2_)
      real otflw(io2_,jo2_)
      real inflw(io2_,jo2_)
      real gdrivo(io2_,jo2_)
      real rdest(io2_,jo2_)
      real area(io2_,jo2_)
      real alon(io2_),dlon(io2_)
      real alat(jo2_),dlat(jo2_)
#ifdef ISOTOPE
      dimension work1(len)
      real runof1(io2_,jo2_)
      real otflw1(io2_,jo2_)
      real inflw1(io2_,jo2_)
      real gdrivo1(io2_,jo2_)
      dimension work2(len)
      real runof2(io2_,jo2_)
      real otflw2(io2_,jo2_)
      real inflw2(io2_,jo2_)
      real gdrivo2(io2_,jo2_)
#endif
c
      integer idest(io2_,jo2_)
      integer jdest(io2_,jo2_)
      integer iofs(0:9)
      integer jofs(0:9)
      integer idx,idy
      integer itmp(io2_)
c                        0  1  2  3  4  5  6  7  8  9
      DATA       IOFS  / 0, 0, 1, 1, 1, 0,-1,-1,-1, 0 /
      DATA       JOFS  / 0,-1,-1, 0, 1, 1, 1, 0,-1, 0 /
c
      character c3*3
      integer ifirst
      data ifirst/0/
      integer nrvmap,nrvini
      data nrvmap,nrvini/93,94/
      real vriver
      data vriver/0.385e0/       !! v=0.5 m/s, meandering ratio=1.3
      save ifirst,idest,jdest,area,rdest,vriver
      real gdrivall,gdrivallo,rflowall,inflwall,rivall,runofall
      real areaall,rbud
c
c  Initial setting
c    
      if (ifirst.eq.0) then
#ifdef MP
      if (mype.eq.master) then
#endif
        ifirst=1
c
c   reading river direction map
c
        close(nrvmap)
        open(nrvmap,file='flw')
        write(c3,'(i3.3)') io2_
        do j=1,jo2_-1
          read(nrvmap,'('//c3//'i1)') (itmp(i),i=1,io2_)
c          print'(180i1)',(itmp(i),i=1,io2_)
          do i=1,io2_
            imap(i,j)=real(itmp(i))
          enddo
        enddo
        close(nrvmap)
c
        do j=1,jo2_-1
          do i=1,io2_
            idest(i,j)=i+iofs(int(imap(i,j)))
            if (idest(i,j).eq.  0  ) idest(i,j)=io2_
            if (idest(i,j).eq.io2_+1) idest(i,j)=1
            jdest(i,j)=j+jofs(int(imap(i,j)))
          enddo
        enddo
cc reading initial
        close(nrvini)
        open(nrvini,file='rivi',status='old',err=99,
     $       form='unformatted')
        print*,'rivi opened'
        read(nrvini,err=99)gdriv
#ifdef ISOTOPE
        read(nrvini,err=99)gdriv1
        read(nrvini,err=99)gdriv2
#endif
        close(nrvini)
        goto 98
 99     print*,'no rivi'
        do j=1,jo2_-1
          do i=1,io2_
            gdriv(i,j)=10.
#ifdef ISOTOPE
            gdriv1(i,j)=10.
            gdriv2(i,j)=10.
#endif
          enddo
        enddo
 98     continue
c 
c   identify river mouth
c
        do j=1,jo2_-1
          do i=1,io2_
            idx=idest(i,j)
            idy=jdest(i,j)
            if (imap(idx,idy).eq.0 .and. (imap(i,j).ne.0)) then
              imap(idx,idy)=9
              print*,'new rivermouth',idx,idy
            endif
          enddo
        enddo
c
c   polarious location
c
        do i=1,io2_
          alon(i)=(real(i)-0.5)/real(io2_)*2.*pi_
          dlon(i)=1./real(io2_)
        enddo
        do j=1,jo2_-1
          alat(j)=(0.5-(real(j)-0.5)/real(jo2_-1))*pi_
          dlat(j)=1./real(jo2_-1)
        enddo
c        
c   calculate rdest
c
        do j=1,jo2_-1
          do i=1,io2_
            rdest(i,j)=0.
            idx=idest(i,j)
            idy=jdest(i,j)
            if (imap(i,j).ge.1 .and. imap(i,j).le.8) then
              cosx=sin(alat(j))*sin(alat(idy))
     $            +cos(alat(j))*cos(alat(idy))
     $            *cos(alon(i)-alon(idx))
              rdest(i,j)=vriver/(acos(cosx)*rerth_)
            elseif (imap(i,j).eq.9) then !! assume flows toward east
              idxx=mod(i,io2_)+1
              cosx=sin(alat(j))*sin(alat(j))
     $            +cos(alat(j))*cos(alat(j))
     $            *cos(alon(i)-alon(idxx))
              rdest(i,j)=vriver/(acos(cosx)*rerth_)
            endif
          enddo
        enddo
c
c  calculate area
        do j=1,jo2_-1
          do i=1,io2_
            area(i,j)=2.*rerth_*rerth_*pi_*pi_*dlon(i)*dlat(j)
     $                *cos(alat(j))
c            if (imap(i,j).ge.1) then
c              print 100,i,j,idest(i,j),jdest(i,j),
c     $              imap(i,j),rdest(i,j),area(i,j)
c            endif
          enddo
        enddo
c
        print*,'end of pgriv first setting'
#ifdef MP
      endif
#endif
      endif !! end of initial setting
 100  format('pgriv first step: ',4i4,f4.0,2e10.2)
c
#ifdef MP
c
c  MPI -> single
c
      call mpgp2f(trunof,lonf2p_,latg2p_,work,lonf2_,latg2_,1)
#ifdef ISOTOPE
      call mpgp2f(trunof1,lonf2p_,latg2p_,work1,lonf2_,latg2_,1)
      call mpgp2f(trunof2,lonf2p_,latg2p_,work2,lonf2_,latg2_,1)
#endif
#else
      do n=1,len
        work(n)=trunof(n,1)
#ifdef ISOTOPE
        work1(n)=trunof1(n,1)
        work2(n)=trunof2(n,1)
#endif
      enddo
#endif
c
c  single
c
#ifdef MP
      if (mype.eq.master) then
#endif
        call rowsep(work)
#ifdef ISOTOPE
        call rowsep(work1)
        call rowsep(work2)
#endif
c
c  interpolate gaussian to center-lonlat
c
        call gau2llc(work,lonf_,latg_,
     &               180./float(io2_),90.-90./float(jo2_-1),
     $               360./float(io2_),180./float(jo2_-1),
     &               runof,io2_,jo2_)
#ifdef ISOTOPE
        call gau2llc(work1,lonf_,latg_,
     &               180./float(io2_),90.-90./float(jo2_-1),
     $               360./float(io2_),180./float(jo2_-1),
     &               runof1,io2_,jo2_)
        call gau2llc(work2,lonf_,latg_,
     &               180./float(io2_),90.-90./float(jo2_-1),
     $               360./float(io2_),180./float(jo2_-1),
     &               runof2,io2_,jo2_)
#endif
c      print*,'check3',work(51,31),runof(51,31)
c
c  reset
c        
        do j=1,jo2_-1
          do i=1,io2_
            gdrivo(i,j)=gdriv(i,j)
            inflw(i,j)=0.
#ifdef ISOTOPE
            gdrivo1(i,j)=gdriv1(i,j)
            inflw1(i,j)=0.
            gdrivo2(i,j)=gdriv2(i,j)
            inflw2(i,j)=0.
#endif
          enddo
        enddo
c      print*,'check4',runof(51,31),gdriv(50,31)
c            
c  trip
c
        do j=1,jo2_-1
          do i=1,io2_
            inflw(i,j)=inflw(i,j)+runof(i,j)
#ifdef ISOTOPE
            inflw1(i,j)=inflw1(i,j)+runof1(i,j)
            inflw2(i,j)=inflw2(i,j)+runof2(i,j)
#endif
            idx=idest(i,j)
            idy=jdest(i,j)
            if (imap(i,j).ge.1) then
              if (rdest(i,j).eq.0.) then
                 print*,'RDEST',i,j,rdest(i,j)
              endif
              gdriv(i,j)=gdriv(i,j)*exp(-rdest(i,j)*deltim)
     $              + (1.e0-exp(-(rdest(i,j)*deltim)))
     $              * runof(i,j)/rdest(i,j)
              otflw(i,j)=((gdrivo(i,j)-gdriv(i,j))/deltim
     $             + inflw(i,j))*area(i,j)
#ifdef ISOTOPE
              gdriv1(i,j)=gdriv1(i,j)*exp(-rdest(i,j)*deltim)
     $              + (1.e0-exp(-(rdest(i,j)*deltim)))
     $              * runof1(i,j)/rdest(i,j)
              otflw1(i,j)=((gdrivo1(i,j)-gdriv1(i,j))/deltim
     $             + inflw1(i,j))*area(i,j)
              gdriv2(i,j)=gdriv2(i,j)*exp(-rdest(i,j)*deltim)
     $              + (1.e0-exp(-(rdest(i,j)*deltim)))
     $              * runof2(i,j)/rdest(i,j)
              otflw2(i,j)=((gdrivo2(i,j)-gdriv2(i,j))/deltim
     $             + inflw2(i,j))*area(i,j)
#endif
              if (imap(i,j).ne.9) then
                gdriv(idx,idy)=gdriv(idx,idy)
     $                + otflw(i,j)/area(idx,idy)*deltim
                inflw(idx,idy)=inflw(idx,idy)
     $               + otflw(i,j)/area(idx,idy)
#ifdef ISOTOPE
                gdriv1(idx,idy)=gdriv1(idx,idy)
     $                + otflw1(i,j)/area(idx,idy)*deltim
                inflw1(idx,idy)=inflw1(idx,idy)
     $               + otflw1(i,j)/area(idx,idy)
                gdriv2(idx,idy)=gdriv2(idx,idy)
     $                + otflw2(i,j)/area(idx,idy)*deltim
                inflw2(idx,idy)=inflw2(idx,idy)
     $               + otflw2(i,j)/area(idx,idy)
#endif
              endif
            else !! at seas or big lakes
              otflw(i,j)=0.
              gdriv(i,j)=0.
#ifdef ISOTOPE
              otflw1(i,j)=0.
              gdriv1(i,j)=0.
              otflw2(i,j)=0.
              gdriv2(i,j)=0.
#endif
            endif
            rflow(i,j)=rflow(i,j)+otflw(i,j)*deltim !! average for output
            roff(i,j)=roff(i,j)+runof(i,j)*deltim   !! averege for output
#ifdef ISOTOPE
            rflow1(i,j)=rflow1(i,j)+otflw1(i,j)*deltim !! average for output
            roff1(i,j)=roff1(i,j)+runof1(i,j)*deltim   !! averege for output
            rflow2(i,j)=rflow2(i,j)+otflw2(i,j)*deltim !! average for output
            roff2(i,j)=roff2(i,j)+runof2(i,j)*deltim   !! averege for output
#endif
          enddo
        enddo
c      print*,'check5',runof(51,31),gdriv(50,31)
c
c  budget check    
c
#ifdef DBG
        gdrivall=0.
        gdrivallo=0.
        rflowall=0.
        inflwall=0.
        rivall=0.
        runofall=0.
        areaall=0.
        do j=1,jo2_-1
          do i=1,io2_
            if (imap(i,j).ge.1) then
              areaall=areaall+area(i,j)
              gdrivall=gdrivall+gdriv(i,j)*area(i,j)
              gdrivallo=gdrivallo+gdrivo(i,j)*area(i,j)
              rflowall=rflowall+otflw(i,j)*deltim
              inflwall=inflwall+inflw(i,j)*area(i,j)*deltim
              runofall=runofall+runof(i,j)*area(i,j)*deltim
              if (imap(i,j).le.8) then
                rivall=rivall+otflw(i,j)*deltim
              endif
              rbud=gdriv(i,j)-gdrivo(i,j)-
     $             (inflw(i,j)-otflw(i,j)/area(i,j))*deltim
              if (abs(rbud).gt.0.1) then
                print'(a15,2i4,f4.0,5e15.7)',
     $                'Riv Inbalance',i,j,imap(i,j),rbud,
     $                gdriv(i,j),gdrivo(i,j),inflw(i,j)*deltim,
     $                otflw(i,j)/area(i,j)*deltim
              endif
            endif
          enddo
        enddo
        print '(a20,2f20.7)','*** pgriv_trip bud ',
     $       (gdrivall-gdrivallo-(inflwall-rflowall))/areaall,
     $       (inflwall-runofall-rivall)/areaall
        print '(a20,e20.7)','*** trip check A   ',areaall
        print '(a20,f20.7)','*** trip check S1  ',gdrivallo/areaall
        print '(a20,f20.7)','*** trip check S2  ',gdrivall/areaall
        print '(a20,f20.7)','*** trip check O   ',rflowall/areaall
        print '(a20,f20.7)','*** trip check I   ',inflwall/areaall
        print '(a20,f20.7)','*** trip check R   ',runofall/areaall
        print '(a20,f20.7)','*** trip check Riv ',rivall/areaall
#endif
c
#ifdef MP
      endif
#endif
c
      return
      end
