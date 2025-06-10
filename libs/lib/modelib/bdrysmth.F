      subroutine bdrysmth(idim,jdim,smthflx)
       implicit none
       integer  idim,jdim
       integer,parameter::nbzgrd=8
       real*8,parameter::gammin=0.01,gammax=1.
       real*8    alpha(nbzgrd)
       real*8    rltb(idim*jdim), smthflx(idim*jdim)
       integer   nitr,itr,i,j,ij

       call nudgecoef(nbzgrd,gammin,gammax,alpha)
       print *,' This is Lehmann nudging coefficients.',
     $         ' zone width=',nbzgrd,' grid points',
     $         ' gammin=',gammin,' gammax=',gammax
       print *,' nudging coeffcients=',alpha

       do ij=1,idim*jdim
        rltb(ij)=1.
       enddo
       do j=1,jdim
        do i=1,nbzgrd
          ij=(j-1)*idim+i
          rltb(ij)=1.-alpha(i)
        enddo
        do i=idim,idim-nbzgrd+1,-1
          ij=(j-1)*idim+i
          rltb(ij)=1.-alpha(idim-i+1)
        enddo
       enddo
       do i=1,idim
        do j=1,nbzgrd
          ij=(j-1)*idim+i
          rltb(ij)=(1.-alpha(j))*rltb(ij)
        enddo
        do j=jdim,jdim-nbzgrd+1,-1
          ij=(j-1)*idim+i
          rltb(ij)=(1.-alpha(jdim-j+1))*rltb(ij)
        enddo
       enddo

       nitr=15
       do itr=1,nitr
        do j=2,jdim-1
         do i=2,idim-1
           ij=(j-1)*idim+i
           if(rltb(ij).ne.1.) then
             smthflx(ij)=(smthflx(i-1+(j-1)*idim)+
     &                    smthflx(i+(j-1-1)*idim)+
     &                    smthflx(i+1+(j-1)*idim)+
     &                    smthflx(i+(j+1-1)*idim)+
     &                    smthflx(i+(j-1)*idim)*4.0)*0.125
           endif
         enddo
        enddo
        do i=1,idim
          smthflx(i              )=smthflx(i+idim)
          smthflx(i+(jdim-1)*idim)=smthflx(i+(jdim-1-1)*idim)
        enddo
        do j=1,jdim
          smthflx(1+(j-1)*idim   )=smthflx(2+(j-1)*idim)
          smthflx(idim+(j-1)*idim)=smthflx(idim-1+(j-1)*idim)
        enddo
       enddo


      end

      subroutine nudgecoef(s,gammin,gammax,alpha)
c
c  Input: s . . . . . . . . . . . . . . width of boundary relaxation zone
c  must be a power of 2.
c  gammin . . . . . . minimal Courant number
c  gammax . . . . . . maximal Courant number
c  (for which "optimal" relaxation coefficients shall be determined)
c  Output: alpha (.) .... weights of externally specified values in the 
c  boundary zone (corresponding to "optimal" relaxation coefficients)
c
      implicit none
      integer s
      real*8 alpha(s)
c
      real*8 gammin,gammax
      real*8, allocatable :: p (:),q (:)
      real*8, allocatable :: pp(:),qq(:)
      real*8 my,kk,kdt2,xxx
      integer n,i,j
c
      allocate (p (0:s),q (0:s))
      allocate (pp(0:s),qq(0:s))
c
c   Computation of P(.), Q(.):
c   p(.) . . . . . . . . . . . . coefficients of polynomial P(b)
c   q(.) . . . . . . . . . . . . coefficients of polynomial Q(b)
c   pp(.), qq(.) .... auxiliary variables for computation of p(.), q(.)
C   n ..... S'
c   "Initialization" for n = 1:
c
      n=1
      p(0)=0.
      p(1)=1.
      q(0)=1.
      q(1)=0.
      my=sqrt(gammax/gammin)
c
c Begin of main loop (Step from n to 2*n):
c
      do while (n.lt.s) 
        my=sqrt((my+1./my)/2.)
        do i=0,n+n
          pp(i)=0.
          qq(i)=0.
        enddo
        do i=0,n
          do j=0,n
            pp(i+j)=pp(i+j)+p(i)*p(j)+q(i)*q(j)
            qq(i+j)=qq(i+j)+2.*my*p(i)*q(j)
          enddo
        enddo
        do i=0,n+n
          p(i)=pp(i)
          q(i)=qq(i)
        enddo
        n=2*n
      enddo
      if(n.ne.s) then
        write (6,*)' ! s is not a power of 2 !'
        call abort
      endif
c
c Computation of K+ and alpha:
c p (.) . . . . . . . . . . . . coefficients of polynomial P_i (b)
c q (.) . . . . . . . . . . . . coefficients of polynomial P_i-1 (b)
c kk . . . . . . . . . . . . . . K+
c kdt2 . . . . . . . . . . . . 2 K dt
c
      do i=n,1,-1
        kk=p(i)/q(i-1)
        do j=i,1,-1
          xxx=q(j)
          q(j)=p(j)-kk*q(j-1)
          p(j)=xxx
        enddo
        xxx=q(0)
        q(0)=p(0)
        p(0)=xxx
        kdt2=kk*sqrt(gammin*gammax)
        alpha(i)=kdt2/(1.+kdt2)
      enddo
c
c REMARK: This alpha corresponds to the leap-frog scheme, 
c whereas kdt2 is independent of the integration scheme.
c
      deallocate (p,q,pp,qq)
      return
      endsubroutine nudgecoef
