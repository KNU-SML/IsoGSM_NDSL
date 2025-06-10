      subroutine rasset(lm, dt, nsphys, fh,  prs
     *,                 rasal, krmin, krmax, kfmax, ncrnd, rannum
     *,                 mct, pctp, pcbt, kctop, dtf, dsfc, len)
!
!     this version of setras assumes same number of levels in the vertical
!     but levels are reverse order
!
!
      parameter (actp=1.7,   facm=1.00)
!
      real dsfc(len), sig(lm+1), prs(lm+1)
     *,    rasal(lm), rannum(200,nsphys)
      real ph(15),    a(15),     ac(16), ad(16)
      integer, allocatable ::  nrnd(:)
!
      common /rasacr/ ac, ad
!
      data ph/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0
     *,       550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/
!
       data a/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677
     *,       0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664
     *,       0.0553, 0.0445, 0.0633/
!
      dimension pctop(mct+1), kctop(mct+1), sgh(mct+1)
!
      data pinc/5.0/
!     data pctp/70.0/, pcbt/900.0/, pinc/5.0/

!
      integer krsize
      logical first
      data iseed/0/, first/.true./, fh0/0.0/
      save iseed, first, fh0, krsize, nrnd
      real seed
!
cfpp$ noconcur r
!
      if (first) then
         first = .false.
!                                   set other parameters
         lmm1 = lm - 1
         tem  = 1.0 / prs(lm+1)
         do l=1,lm+1
            sig(l) = prs(l) * tem
         enddo
!
         tem  = (pcbt - pctp - (mct-1)*mct*pinc/2) / mct
         pctop(1) = pctp * 0.001
         do i=1,mct
           pctop(i+1) = 0.001 * (pctop(i)*1000.0 + tem + pinc*(i-1))
         enddo
         print *,' pctop=',pctop
         do l=1,lm
           if (sig(l) . le. pctop(1)) kctop(1) = l
         enddo
         tem2 = 0.0
         do i=2,mct+1
           if (kctop(i-1) .lt. lm) then
             ii       = kctop(i-1) + 1
             kctop(i) = ii
             do l=ii,lm
               if (sig(l) . le. pctop(i)) kctop(i) = l
             enddo
             if (pctop(i) .gt. sig(kctop(i))) then
               tem1 = sig(kctop(i)) - sig(kctop(i-1))
               if (tem1 .lt. tem2) kctop(i) = min(lm, kctop(i)+1)
               tem2 = sig(kctop(i)) - sig(kctop(i-1))
             endif
           else
             kctop(i) = lm + 1
           endif
         enddo
         do i=1,mct+1
           if (kctop(i) .le. lm) mctm = i
         enddo
!        do i=1,mctm
!          sgh(i) = sig(kctop(i))
!        enddo
!        do i=1,mct
!          sgc(i)  = 0.5 * (sgh(i) + sgh(i+1))
!          dsgh(i) = sgh(i+1) - sgh(i)
!        enddo
         print *,' mct=',mct,' mctm=',mctm
         print *,' kctop=',(kctop(i),i=1,mctm)
!        print *,' sgh=',(sig(kctop(i)),i=1,mctm)
!        print *,' sgc=',sgc
!        print *,' dsgh=',dsgh
         do i=1,mctm-1
           if (kctop(i) .eq. kctop(i+1)) then
             print *,' for cloud-type = ',i, ' kctop same at top'
     *,              ' and bottom - run stopped'
             stop
           endif
         enddo
!
         krmin = 1
         krmax = mctm - 1
         kfmax = mctm - 1
         do l=1,mctm-1
            sgc  = 0.5 * (sig(kctop(l)) + sig(kctop(l+1)))
            if (sgc .le. 0.760) krmax = l
            if (sgc .le. 0.930) kfmax = l
         enddo
!
!        print *,' krmin=',krmin,' krmax=',krmax,' kfmax=',kfmax
!
         call random_seed(size=krsize)
         print *,' krsize=',krsize
         allocate (nrnd(krsize))
!
!     compute ncrnd and the random numbers
!
!        ncrnd   = 42 * (dt/3600) + 0.50001

         ncrnd   = (krmax-krmin+1) * (dtf/1200) + 0.50001
         if (dt .gt. dtf) ncrnd = (5*ncrnd) / 4
!        ncrnd   = (krmax-krmin+1) * (dt/1200) + 0.50001
         ncrnd   = max(ncrnd, 1)

!        print *,' ncrnd=',ncrnd,' dt=',dt,' dtf=',dtf
!    *,          ' rasalf=',rasalf
!
      if (fh .eq. fh0) then
         iseed = iseed + 1
      else
         iseed = fh*3600.0 + 0.0001
         iseed = max(iseed, 1)
         fh0   = fh
      endif
!
!                                                       for c90
!     call ranset(iseed)
!     do nn=1,nsphys
!       do i=1,ncrnd
!         rannum(i,nn) = ranf()
!       enddo
!     enddo
!
!                                                       for ibmsp
!     allocate (nrnd(krsize))
      nrnd = iseed
      call random_seed(put=nrnd)
      call random_number(rannum)
!
!     seed = iseed
!     call srand(seed)     
!     do nn=1,nsphys
!       do i=1,ncrnd
!         rannum(i,nn) = rand()
!       enddo
!     enddo
!
      return
      end
