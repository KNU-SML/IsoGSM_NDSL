      subroutine makeoa(zavg,zmax,var,glat,oa4,ol,ioa4,
     1           oro,xnsum,xnsum1,xnsum2,xnsum3,xnsum4,
     2           ist,ien,jst,jen,im,jm,imn,jmn,xlat,numi)
      dimension glat(jmn),xlat(jm)
      integer zavg(imn,jmn),zmax(imn,jmn)
      dimension oro(im,jm)
      dimension oa4(im,jm,4),ioa4(im,jm,4)
      dimension ist(im,jm),ien(im,jm),jst(jm),jen(jm)
      dimension xnsum(im,jm),xnsum1(im,jm),xnsum2(im,jm)
      dimension xnsum3(im,jm),xnsum4(im,jm)
      dimension var(im,jm),ol(im,jm,4),numi(jm)
      logical flag
c
c---- global xlat and xlon ( degree )
c
      im1 = im - 1
      jm1 = jm - 1
      delxn = 360./imn      ! mountain data resolution
c
      do j=1,jmn
         glat(j) = -90. + (j-1) * delxn + delxn * 0.5
      enddo
      print *,' im=',im,' jm=',jm,' imn=',imn,' jmn=',jmn
c
c---- find the average of the modes in a grid box
c
      do j=1,jm
      do i=1,numi(j)
         delx  = 360./numi(j)       ! gaussian grid resolution
         faclon  = delx / delxn
c        ist(i,j) = faclon * float(i-1) + faclon * 0.5
c        ist(i,j) = ist(i,j) + 1
c        ien(i,j) = faclon * float(i) + faclon * 0.5
         ist(i,j) = faclon * float(i-1) + 1.0001
         ien(i,j) = faclon * float(i)   + 0.0001
         if (ist(i,j) .le. 0)      ist(i,j) = ist(i,j) + imn
         if (ien(i,j) .lt. ist(i,j)) ien(i,j) = ien(i,j) + imn
cx         print*, ' i j ist ien ',i,j,ist(i,j),ien(i,j)
      enddo
      enddo
      do j=1,jm-1
         flag=.true.
         do j1=1,jmn
            if(flag.and.glat(j1).gt.xlat(j+1)) then
              jst(j) = j1
              flag = .false.
            endif
         enddo
         flag=.true.
         do j1=jst(j),jmn
            if(flag.and.glat(j1).gt.xlat(j)) then
              jen(j) = j1 - 1
              flag = .false.
            endif
         enddo
      enddo
c
      do j=1,jm
        do i=1,numi(j)
          xnsum(i,j) = 0.0
          oro(i,j)   = 0.0
        enddo
      enddo
      do j=1,jm1
        do i=1,numi(j)
          do ii1=ist(i,j),ien(i,j)
             i1 = ii1
             if (ii1 .gt. imn) i1 = i1 - imn
             do j1=jst(j),jen(j)
                xnsum(i,j) = xnsum(i,j) + 1
                height = float(zavg(i1,j1))
                if(height.lt.-990.) height = 0.0
                oro(i,j) = oro(i,j) + height
             enddo
          enddo
        enddo
      enddo
      do j=1,jm1
        do i=1,numi(j)
          if(xnsum(i,j).gt.0.) then
             oro(i,j) = oro(i,j) / xnsum(i,j)
          else
             oro(i,j) = 0.0
          endif
        enddo
      enddo
c
c---- count number of max. peaks higher than the average mode
c     in a grid box
!
      do j=1,jm
        do i=1,numi(j)
          xnsum(i,j) = 0.0
        enddo
      enddo
      do j=1,jm1
        do i=1,numi(j)
           do ii1=ist(i,j),ien(i,j)
              i1 = ii1
              if(ii1.gt.imn) i1 = i1 - imn
             do j1=jst(j),jen(j)
                if(float(zmax(i1,j1)) .gt. oro(i,j))
     1             xnsum(i,j) = xnsum(i,j) + 1
             enddo
           enddo
        enddo
      enddo
c
c---- count number of mode. higher than the hc, critical height
c     in a grid box
      do j=1,jm
        do i=1,numi(j)
          xnsum1(i,j) = 0.0
          xnsum2(i,j) = 0.0
          xnsum3(i,j) = 0.0
          xnsum4(i,j) = 0.0
        enddo
      enddo
      do j=1,jm1
        do i=1,numi(j)
           hc = 1116.2 - 0.878 * var(i,j)
           do ii1=ist(i,j),ien(i,j)
             i1 = ii1
             if(ii1.gt.imn) i1 = i1 - imn
             do j1=jst(j),jen(j)
               if(float(zavg(i1,j1)) .gt. hc)
     1            xnsum1(i,j) = xnsum1(i,j) + 1
               xnsum2(i,j) = xnsum2(i,j) + 1
             enddo
           enddo
c
           inci = nint((ien(i,j)-ist(i,j)) * 0.5)
           isttt = min(max(ist(i,j)-inci,1),imn)
           ieddd = min(max(ien(i,j)-inci,1),imn)
c
           incj = nint((jen(j)-jst(j)) * 0.5)
           jsttt = min(max(jst(j)-incj,1),jmn)
           jeddd = min(max(jen(j)-incj,1),jmn)
c
           do i1=isttt,ieddd
             do j1=jsttt,jeddd
               if(float(zavg(i1,j1)) .gt. hc)
     1            xnsum3(i,j) = xnsum3(i,j) + 1
               xnsum4(i,j) = xnsum4(i,j) + 1
             enddo
           enddo
cx         print*,' i j hc var ',i,j,hc,var(i,j)
cx         print*,'xnsum12 ',xnsum1(i,j),xnsum2(i,j)
cx         print*,'xnsum34 ',xnsum3(i,j),xnsum4(i,j)
        enddo
      enddo
c
c---- calculate the 3d orographic asymmetry for 4 wind directions
c---- and the 3d orographic subgrid orography fraction
c     (kwd = 1  2  3  4)
c     ( wd = w  s sw nw)
c
c
      do kwd = 1, 4
        do j=1,jm
          do i=1,numi(j)
            oa4(i,j,kwd) = 0.0
          enddo
        enddo
      enddo
c
      do j=1,jm-2
        do i=1,numi(j)
        ii = i + 1
        if (ii .gt. numi(j)) ii = ii - numi(j)
          xnpu = xnsum(i,j)    + xnsum(i,j+1)
          xnpd = xnsum(ii,j)   + xnsum(ii,j+1)
          if (xnpd .ne. xnpu) oa4(ii,j+1,1) = 1. - xnpd / max(xnpu , 1.)
          ol(ii,j+1,1) = (xnsum3(i,j+1)+xnsum3(ii,j+1))/
     1                   (xnsum4(i,j+1)+xnsum4(ii,j+1))
        enddo
      enddo
      do j=1,jm-2
        do i=1,numi(j)
        ii = i + 1
        if (ii .gt. numi(j)) ii = ii - numi(j)
          xnpu = xnsum(i,j+1)   + xnsum(ii,j+1)
          xnpd = xnsum(i,j)     + xnsum(ii,j)
          if (xnpd .ne. xnpu) oa4(ii,j+1,2) = 1. - xnpd / max(xnpu , 1.)
          ol(ii,j+1,2) = (xnsum3(ii,j)+xnsum3(ii,j+1))/
     1                   (xnsum4(ii,j)+xnsum4(ii,j+1))
        enddo
      enddo
      do j=1,jm-2
        do i=1,numi(j)
        ii = i + 1
        if (ii .gt. numi(j)) ii = ii - numi(j)
          xnpu = xnsum(i,j+1)  + ( xnsum(i,j) + xnsum(ii,j+1) )*0.5
          xnpd = xnsum(ii,j)   + ( xnsum(i,j) + xnsum(ii,j+1) )*0.5
          if (xnpd .ne. xnpu) oa4(ii,j+1,3) = 1. - xnpd / max(xnpu , 1.)
          ol(ii,j+1,3) = (xnsum1(ii,j)+xnsum1(i,j+1))/
     1                   (xnsum2(ii,j)+xnsum2(i,j+1))
        enddo
      enddo
      do j=1,jm-2
        do i=1,numi(j)
        ii = i + 1
        if (ii .gt. numi(j)) ii = ii - numi(j)
          xnpu = xnsum(i,j)    + ( xnsum(ii,j) + xnsum(i,j+1) )*0.5
          xnpd = xnsum(ii,j+1) + ( xnsum(ii,j) + xnsum(i,j+1) )*0.5
          if (xnpd .ne. xnpu) oa4(ii,j+1,4) = 1. - xnpd / max(xnpu , 1.)
          ol(ii,j+1,4) = (xnsum1(i,j)+xnsum1(ii,j+1))/
     1                   (xnsum2(i,j)+xnsum2(ii,j+1))
        enddo
      enddo
c
      do kwd = 1, 4
        do i=1,numi(j)
          ol(i,1,kwd)  = ol(i,2,kwd)
          ol(i,jm,kwd) = ol(i,jm-1,kwd)
        enddo
      enddo
c
      do kwd=1,4
        do j=1,jm
          do i=1,numi(j)
            t = oa4(i,j,kwd)
            oa4(i,j,kwd) = sign( min( abs(t), 1. ), t )
          enddo
        enddo
      enddo
c
      ns0 = 0
      ns1 = 0
      ns2 = 0
      ns3 = 0
      ns4 = 0
      ns5 = 0
      ns6 = 0
      do kwd=1,4
      do j=1,jm
      do i=1,numi(j)
         t = abs( oa4(i,j,kwd) )
         if(t .eq. 0.) then
            ioa4(i,j,kwd) = 0
            ns0 = ns0 + 1
         else if(t .gt. 0. .and. t .le. 1.) then
            ioa4(i,j,kwd) = 1
            ns1 = ns1 + 1
         else if(t .gt. 1. .and. t .le. 10.) then
            ioa4(i,j,kwd) = 2
            ns2 = ns2 + 1
         else if(t .gt. 10. .and. t .le. 100.) then
            ioa4(i,j,kwd) = 3
            ns3 = ns3 + 1
         else if(t .gt. 100. .and. t .le. 1000.) then
            ioa4(i,j,kwd) = 4
            ns4 = ns4 + 1
         else if(t .gt. 1000. .and. t .le. 10000.) then
            ioa4(i,j,kwd) = 5
            ns5 = ns5 + 1
         else if(t .gt. 10000.) then
            ioa4(i,j,kwd) = 6
            ns6 = ns6 + 1
         endif
      enddo
      enddo
      enddo
c
      write(6,*) "! makeoa exit"
c
      return
      end
