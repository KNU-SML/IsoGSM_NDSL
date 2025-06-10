cfpp$ noconcur r
c-----------------------------------------------------------------------
      subroutine tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    tridi2      solves tridiagonal matrix problems.
c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
c
c abstract: this routine solves multiple tridiagonal matrix problems
c   with 2 right-hand-side and solution vectors for every matrix.
c   the solutions are found by eliminating off-diagonal coefficients,
c   marching first foreward then backward along the matrix diagonal.
c   the computations are vectorized around the number of matrices.
c   no checks are made for zeroes on the diagonal or singularity.
c
c program history log:
c   91-05-07  iredell
c
c usage:    call tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
c
c   input argument list:
c     l        - integer number of tridiagonal matrices
c     n        - integer order of the matrices
c     cl       - real (l,2:n) lower diagonal matrix elements
c     cm       - real (l,n) main diagonal matrix elements
c     cu       - real (l,n-1) upper diagonal matrix elements
c                (may be equivalent to au if no longer needed)
c     r1       - real (l,n) 1st right-hand-side vector elements
c                (may be equivalent to a1 if no longer needed)
c     r2       - real (l,n) 2nd right-hand-side vector elements
c                (may be equivalent to a2 if no longer needed)
c
c   output argument list:
c     au       - real (l,n-1) work array
c     a1       - real (l,n) 1st solution vector elements
c     a2       - real (l,n) 2nd solution vector elements
c
c remarks: this routine can be easily modified to solve a different
c   number of right-hand-sides and solutions per matrix besides 2.
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
c
c$$$
      dimension cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),r2(l,n),
     &          au(l,n-1),a1(l,n),a2(l,n)
c-----------------------------------------------------------------------
      do i=1,l
        fk=1./cm(i,1)
        au(i,1)=fk*cu(i,1)
        a1(i,1)=fk*r1(i,1)
        a2(i,1)=fk*r2(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)=fk*cu(i,k)
          a1(i,k)=fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k)=fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        enddo
      enddo
      do i=1,l
        fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n)=fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n)=fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k)=a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k)=a2(i,k)-au(i,k)*a2(i,k+1)
        enddo
      enddo
c-----------------------------------------------------------------------
      return
      end
