      subroutine parusr(str,lun,i1,io)

      common /usrstr/ nnod,ncon,nods(20),nodc(10),vals(10),kons(10)

      character*(*) str
      character*80  ust
      character*20  utg(30)
      logical       bump
      real*4 vals

      data maxusr /30/
      data maxnod /20/
      data maxcon /10/

c----------------------------------------------------------------------
c----------------------------------------------------------------------

      if(len(str).gt.80) goto 900
      ust  = str

      ncon = 0
      nnod = 0

c  process string pieces(s) into cond nodes and store nodes
c  --------------------------------------------------------

      call parseq(ust,utg,maxusr,ntot)

      do n=1,ntot
      call parutg(lun,io,utg(n),nod,kon,val,*908)
      if(kon.ne.0) then
         ncon = ncon+1
         if(ncon.gt.maxcon) goto 901
         nodc(ncon) = nod
         kons(ncon) = kon
         vals(ncon) = val
      else
         nnod = nnod+1
         if(nnod.gt.maxnod) goto 902
         nods(nnod) = nod
      endif
      enddo

c  sort cond nodes in jump/link table order
c  ----------------------------------------

      do i=1,ncon
      do j=i+1,ncon
      if(nodc(i).gt.nodc(j)) then
         nod     = nodc(i)
         nodc(i) = nodc(j)
         nodc(j) = nod

         kon     = kons(i)
         kons(i) = kons(j)
         kons(j) = kon

         val     = vals(i)
         vals(i) = vals(j)
         vals(j) = val
      endif
      enddo
      enddo

c  check on special rules for bump nodes
c  -------------------------------------

      bump = .false.

      do n=1,ncon
      if(kons(n).eq.5) then
         if(io.eq.0)   goto 903
         if(n.ne.ncon) goto 904
         bump = .true.
      endif
      enddo

c  check store node count and alignment
c  ------------------------------------

      if(.not.bump .and. nnod.eq.0) goto 905
      if(nnod.gt.i1)                goto 906

      irpc = -1
      do i=1,nnod
      if(nods(i).gt.0) then
         if(irpc.lt.0) irpc = lstrpc(nods(i),lun)
         if(irpc.ne.lstrpc(nods(i),lun)) goto 907
      endif
      enddo

      return
900   call abt('parusr - user string > 80 chars         :'//ust)
901   call abt('parusr - too many cond nodes            :'//ust)
902   call abt('parusr - too many stor nodes            :'//ust)
903   call abt('parusr - bump on input not allowed      :'//ust)
904   call abt('parusr - bump must be on inner node     :'//ust)
905   call abt('parusr - user string has no store nodes :'//ust)
906   call abt('parusr - must be at least i1 store nodes:'//ust)
907   call abt('parusr - store nodes must in one rep grp:'//ust)
908   call abt('parusr - parutg:'                         //ust)
      end
