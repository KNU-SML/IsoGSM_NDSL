#define LEVP1 levp1
#define LEVS levs
#define MLVARK mlvark
#define SLVARK slvark
      subroutine wriken(npoint,ikfreq,imodk,itnum,svdata,
     1                  lab,fhour,idate,si,sl,nvrken,nptken,nstken,nn)
#include <paramodel.h>
#ifdef KEN
c
       character*8 lab
       dimension svdata(nvrken,nptken,nstken)
       dimension si(levp1_),sl(levs_)
       dimension idate(4)
c
      character*128 fno
c
      if(npoint.le.0) return
c
      call fnam('ken',nchi,fhour,fno,ncho)
#ifdef ASSIGN
      call assign('assign -R')
#endif
      open(unit=nn,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      write(6,*) ' error in opening file ',fno(1:ncho)
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
  901 continue
#ifndef NOPRINT
      write(6,*) ' file ',fno(1:ncho),' opened. unit=',nn
#endif
c
c  note, that in several scenerios, itnum=itnum+1 at the
c  bottom of the 10000 loop, so undo it
c
       if (ikfreq.gt.1) then
         if (imodk.le.0) then
           itnum = itnum - 1
         end if
       else
         itnum = itnum - 1
       end if
#ifndef NOPRINT
       print 1047,itnum,npoint
 1047  format(1h0,i6,' steps of ken(campana) gridpt data saved for ',
     1            i5,' points')
#endif
       do j=1,npoint
        do  k=1,itnum
          if (svdata(44,j,k).le.0.) then
            do i=25,27
              svdata(i,j,k) = svdata(i,j,k-1)
            enddo
            do i=41,46
              svdata(i,j,k) = svdata(i,j,k-1)
            enddo
            do i = 48, 49
              svdata(i,j,k) = svdata(i,j,k-1)
            enddo
            do i = 51, 58
              svdata(i,j,k) = svdata(i,j,k-1)
            enddo
            i=38
            svdata(i,j,k) = svdata(i,j,k-1)
            do lv=1,levs_
              svdata(lv+slvark_+(mlvark_-1)*levs_,j,k)=
     1               svdata(lv+slvark_+(mlvark_-1)*levs_,j,k-1)
            enddo
          endif
          if (svdata(50,j,k).le.0.) then
            i = 47
            svdata(i,j,k) = svdata(i,j,k-1)
            i = 50
            svdata(i,j,k) = svdata(i,j,k-1)
          end if
          if (svdata(70,j,k).le.0.) then
            do i = 70, 74
              svdata(i,j,k) = svdata(i,j,k-1)
            enddo
          end if
        enddo
      enddo
c
       rewind nn
       write(nn) lab
       write(nn) fhour,idate,si,sl
       write(nn) nvrken,nptken,nstken,npoint,itnum
       do j=1,npoint
        write(nn) ((svdata(i,j,k),k=1,itnum),i=1,nvrken)
       enddo
       close(nn)
#endif
      return
      end
