      subroutine nfdr_getfixf(iunit,condir,f,idim,jdim,ann,ifm,slp)
c
      implicit none
c
      character*128 condir
c
      character*1 mdlltr
      integer idim,jdim,i,j
      integer itmp
      integer nchar
c
      real f(20,20)
      real ann(idim,jdim),slp(idim,jdim)
      integer iunit,ifm(idim,jdim)
c
c  read data from 1988 nfdrfuel model and put it in the f array
c
      do nchar=1,200
        if (condir(nchar:nchar).eq.' ') go to 1
      enddo
    1 continue
      nchar=nchar-1
      if(nchar.le.0) then
        print *,'error in getfixf' 
        call abort
      endif
      open (unit=iunit,file=condir(1:nchar)//'/nfdrgrid25.mdl',
     1      form='formatted',status='old')
      do i=1,20
        read(iunit,10) mdlltr,f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),
     *   f(i,6),f(i,7),f(i,8),f(i,9),f(i,10),f(i,11),f(i,12),f(i,13),
     *   f(i,14),f(i,15),f(i,16),f(i,17),f(i,18),f(i,19),f(i,20)
      enddo
   10 format(a1, 6(f4.2,1x), f5.2,1x, f4.0,2x, f3.0,3x, f2.0,4x,
     * f1.0,1x, 3(f4.0,1x), f3.0,1x, f2.0,1x, f4.2,1x, 3(f2.1,1x),
     * f4.2)
      close(iunit)
c
c  read annual precipitation climatology
c
      open(iunit,file=condir(1:nchar)//'/meanppt25k.fmt',
     *    status='old',form='formatted')
      do j=1,jdim
       do i=1,idim
         read(iunit,'(i3)') itmp
         ann(i,j)=float(itmp)
       enddo
      enddo
      close(iunit)
c
c  read fuelmap
c
      open(iunit,file=condir(1:nchar)//'/nfdrfuel25k.fmt',
     *     status='old',form='formatted')
      do j=1,jdim
       do i=1,idim
         read(iunit,'(i3)') ifm(i,j)
       enddo
      enddo
      close(iunit)
c
c  read slope
c
      open(iunit,file=condir(1:nchar)//'/slope25k.fmt',
     *     status='old',form='formatted')
      do j=1,jdim
       do i=1,idim
         read(iunit,'(i3)') itmp
         slp(i,j)=float(itmp)
       enddo
      enddo
      close(iunit)
c
      return
      end
