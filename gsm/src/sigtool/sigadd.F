#define JCAP jcap
#define LEVS levs
      program sigadd
#include <paramodel.h>
c     
c     average two sigma-level spectrum files
c     usage: sigadd inputsig1 inputsig2 outputsig
c     
      integer jcap,kdim
      parameter (jcap=jcap_,kdim=levs_)
      real*8 var((jcap+1)*(jcap+2))
      real*8 var2((jcap+1)*(jcap+2))
      character head*32
      character head2*32
      character infile1*80,infile2*80,outfile*80
      real*8 dummy1,dummy2(245)
      real*8 dummy3,dummy4(245)
      integer i,j,k,idum(4),idum2(4)
c
      if (iargc().ne.3) then
         write(6,*) 'Usage: sigadd input1 input2 output'
         stop
      endif
      call getarg(1,infile1)
      call getarg(2,infile2)
      call getarg(3,outfile)
c
      open(1,file=infile1,form='unformatted')
      open(2,file=infile2,form='unformatted')
c      open(10,file=outfile,form='unformatted',
c     $     access='direct',recl=4*(jcap+1)*(jcap+2))
      open(10,file=outfile,form='unformatted')
      read(1) head
      read(1) dummy1,idum,dummy2
      read(2) head2
      read(2) dummy3,idum2,dummy4
c      write(6,'(i4,3e15.7)') 
c     $     (i,dummy2(i),dummy4(i),dummy2(i)-dummy4(i),i=1,245)
      write(10) head
      write(10) dummy3,idum2,dummy4
c
      do k=1,2
         read(1) var
         read(2) var2
c         write(6,*) k,k,var(1),var2(1)
         write(10) ((var(j)+var2(j)),j=1,(jcap+1)*(jcap+2))
c         write(10,rec=k) (real(var(j)),j=1,(jcap+1)*(jcap+2))
      enddo
c
      do k=1,kdim
         read(1) var
         read(2) var2
c         write(6,*) k+3,k+3,var(1),var2(1)
         write(10) ((var(j)+var2(j)),j=1,(jcap+1)*(jcap+2))
c         write(10,rec=k+2) (real(var(j)),j=1,(jcap+1)*(jcap+2))
      enddo
c
      do k=1,kdim*2
         read(1) var
         read(2) var2
c         write(6,*) k+2+kdim,2+kdim+mod(k+1+kdim,2)*kdim+int((k+1)/2),
c     $        var(1),var2(1)
         write(10) ((var(j)+var2(j)),j=1,(jcap+1)*(jcap+2))
c         write(10,rec=2+kdim+mod(k+1+kdim,2)*kdim+int((k+1)/2))
c     $        (real(var(j)),j=1,(jcap+1)*(jcap+2))
      enddo
c
      do k=1,kdim
         read(1) var
         read(2) var2
c         write(6,*) k+2+kdim*3,k+2+kdim*3,var(1),var2(1)
         write(10) ((var(j)+var2(j)),j=1,(jcap+1)*(jcap+2))
c         write(6,*) k+2+kdim*3,k+2+kdim*3,var(1)
c         write(10,rec=k+2+kdim*3) (real(var(j)),j=1,(jcap+1)*(jcap+2))
      enddo
c         
      stop
      end
