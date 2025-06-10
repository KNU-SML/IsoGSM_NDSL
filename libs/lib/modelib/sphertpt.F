      subroutine sphpt1(jcap,kmx,wave,physv,rlon,ppnm)
c                                                                       
c  this routine performs conversion of spectral coefficients
c  to physical space on arbitrary lat/lon points.
c
c  input:
c     jcap:  spherical coefficient resolution
c      kmx:  number of levels
c     wave:  spherical coefficents
c     rlon:  longituds of the points
c
c  output:
c    physv:  physical space value
c
      dimension wave((jcap+1)*(jcap+2),kmx)
      dimension physv(kmx)
c                                                                       
c  local dimension
c
      real, allocatable :: gfftc(:,:)
      real, allocatable :: sumc(:,:,:)
c
      dimension ppnm((jcap+1)*(jcap+2))
c
      allocate (sumc(2,0:jcap,kmx))
      allocate (gfftc(2,0:jcap))
c
      jcap1=jcap+1
c
      gfftc(1,0)=1.0
      gfftc(2,0)=0.0
      do mm=1,jcap
        waven=float(mm)
        gfftc(1,mm)= cos(waven*rlon)*2.0
        gfftc(2,mm)=-sin(waven*rlon)*2.0
      enddo
c                                                                       
c  do transform                                                 
c                                                                       
      do n=1,2
      do k=1,kmx
      do m=0,jcap
        sumc(n,m,k)=0.0
      enddo
      enddo
      enddo
      do l=0,jcap
        ls=l*((2*jcap+3)-l)
        do i=1,jcap1-l
          ii=(i-1)*2
          do k=1,kmx
            sumc(1,l,k)=sumc(1,l,k)+ppnm(ls+1+ii)*wave(ls+1+ii,k)
            sumc(2,l,k)=sumc(2,l,k)+ppnm(ls+2+ii)*wave(ls+2+ii,k)
          enddo
        enddo
      enddo
      do k=1,kmx
        physv(k)=0.0
      enddo
      do i=1,jcap1
        ii=i-1
        do k=1,kmx
          physv(k)=physv(k)+gfftc(1,ii)*sumc(1,ii,k)
     1                     +gfftc(2,ii)*sumc(2,ii,k)
        enddo
      enddo
c
      deallocate (gfftc)
      deallocate (sumc)
c
      return                                                            
      end                                                               
      subroutine sphpt2(jcap,kmx,wave,physvx,physvy,rlon,ppnm,hhnm)
c
c    1                  rlon,rsnln,rcsln)
c
c  this routine performs conversion of spectral coefficients
c  to space derivative in physical space on arbitrary lat/lon points.
c
c  input:
c     jcap:  spherical coefficient resolution
c      kmx:  number of levels
c     wave:  spherical coefficents
c     rlon:  longituds of the points
c      pnm:  legendre polynomial
c      hnm:  legendre polynomial y-derivative
c not used
c    rsnln:  coordinate rotation sin factor
c    rcsln:  coordinate rotation cos factor
c
c  output:
c    physvx:  d/dx
c    physvy:  d/dy
c   
      dimension wave((jcap+1)*(jcap+2),kmx)
     1         ,physvx(kmx),physvy(kmx)
      dimension ppnm((jcap+1)*(jcap+2)),hhnm((jcap+1)*(jcap+2))
c
c  local variables
c
      dimension qex(0:jcap)
      dimension gfftc(2,0:jcap),gfftx(2,0:jcap)  
      dimension sumc(2,0:jcap,kmx),sumy(2,0:jcap,kmx)
c
      parameter(er_=6.3712e06)
      parameter(eriv_=1./er_)
c
      jcap1=jcap+1
c                                                                       
c  get fer and qimer for spectral transfer
c
      do m=0,jcap
        qex(m)=eriv_*float(m)
      enddo
c
      gfftc(1,0)=1.0
      gfftc(2,0)=0.0
      do mm=1,jcap
        waven=float(mm)
        gfftc(1,mm)= cos(waven*rlon)*2.0
        gfftc(2,mm)=-sin(waven*rlon)*2.0
      enddo
      do mm=0,jcap
        gfftx(2,mm)=-gfftc(1,mm)*qex(mm)
        gfftx(1,mm)= gfftc(2,mm)*qex(mm)
      enddo
c                                                                      
c  do transform                                                 
c                                                                       
      do n=1,2
      do k=1,kmx
      do m=0,jcap
        sumc(n,m,k)=0.0
        sumy(n,m,k)=0.0
      enddo
      enddo
      enddo
      do l=0,jcap
        ls=l*((2*jcap+3)-l)
        do i=1,jcap1-l
          ii=(i-1)*2
          do k=1,kmx
            sumc(1,l,k)=sumc(1,l,k)+ppnm(ls+1+ii)*wave(ls+1+ii,k)
            sumc(2,l,k)=sumc(2,l,k)+ppnm(ls+2+ii)*wave(ls+2+ii,k)
            sumy(1,l,k)=sumy(1,l,k)+hhnm(ls+1+ii)*wave(ls+1+ii,k)
            sumy(2,l,k)=sumy(2,l,k)+hhnm(ls+2+ii)*wave(ls+2+ii,k)
          enddo
        enddo
      enddo
      do k=1,kmx
        physvx(k)=0.0
        physvy(k)=0.0
      enddo
      do i=1,jcap1
        ii=i-1
        do k=1,kmx
          physvx(k)=physvx(k)+gfftx(1,ii)*sumc(1,ii,k)
     1                       +gfftx(2,ii)*sumc(2,ii,k)
          physvy(k)=physvy(k)+gfftc(1,ii)*sumy(1,ii,k)*eriv_
     1                       +gfftc(2,ii)*sumy(2,ii,k)*eriv_
        enddo
      enddo
c
c  rotation if different coordinate
c
c     do k=1,kmx
c       do n=1,lngrdb
c         qy=physvy(n,k)
c         qx=physvx(n,k)
c         physvx(n,k)=-qx*rsnln(n)-qy*rcsln(n)
c         physvy(n,k)= qx*rcsln(n)-qy*rsnln(n)
c       enddo
c     enddo
c
      return                                                            
      end
      subroutine rpln2i(jcap,qlnt,dqlnt,colrad)
c 
c  computes legendre polinomial and its y-derivatives 
c 
      dimension qlnt((jcap+1)*(jcap+2))
      dimension dqlnt((jcap+1)*(jcap+2))
c
c local arrays
c
      dimension qlnv((jcap+1)*(jcap+4))
      dimension x(jcap+1)
      dimension dpln((jcap+1)*(jcap+4))
c
      dimension indxmv((jcap+1)*(jcap+4))
      dimension deps((jcap+1)*(jcap+4)),rdeps((jcap+1)*(jcap+4)),
     1          dx(2*(jcap+1)),y(jcap+1)
      dimension dxab((jcap+1)*(jcap+2),2)
c
      call gpln2i(jcap,indxmv,deps,rdeps,dx,y)
      call ggozrm(jcap,dxab)
c
      jcap1=jcap+1
      jcap2=jcap+2
c
      colr=colrad
      sinlat=cos(colr)
      cos2=1.0-sinlat*sinlat
      prod=1.0
      do ll=1,jcap1
        x(ll)=0.5*prod
        prod=prod*cos2*y(ll)
      enddo
      do ll=1,jcap1
        x(ll)=sqrt(x(ll))
      enddo
      do ll=1,jcap1
        dpln(2*ll-1)=x(ll)
        dpln(2*ll  )=x(ll)
      enddo
      lplus=2*jcap1
      do ll=1,2*jcap1
        dpln(ll+lplus)=dx(ll)*sinlat*dpln(ll)
      enddo
      lp2=0
      lp1=2*jcap1
      lp0=2*2*jcap1
      len=2*jcap1-2
      do n=3,jcap2
        do ll=1,len
          dpln(ll+lp0)=(sinlat*dpln(ll+lp1)
     1             -deps(ll+lp1)*dpln(ll+lp2))*rdeps(ll+lp0)
        enddo
        lp2=lp1
        lp1=lp0
        lp0=lp0+len
        len=len-2
      enddo
c
c  transpose vector dpln array from cra. order to ibm order.
c
      do i=1,jcap1*(jcap+4)
        qlnv(indxmv(i))=dpln(i)
      enddo
c
      lpv=0
      lpt=0
      len=2*jcap1
      do n=1,jcap1
        do ll=1,len
          qlnt(ll+lpt)=qlnv(ll+lpv)
        enddo
        lpv=lpv+len+2
        lpt=lpt+len
        len=len-2
      enddo
      lp0=0
      lp1=2
      len=2*jcap1
      do i=1,jcap1
        do ll=1,len
           dqlnt(ll+lp0)=+qlnv(ll+lp1)*dxab(ll+lp0,2)
        enddo
        lp1=lp1+len+2
        lp0=lp0+len
        len=len-2
      enddo
      lend=jcap1*jcap2-4
      do ll=1,lend
        dqlnt(ll+2)=dqlnt(ll+2)+qlnt(ll)*dxab(ll+2,1)
      enddo
c
      return
      end
      subroutine gpln2i(jcap,indxmv,deps,rdeps,dx,y)
c
      dimension indxmv((jcap+1)*(jcap+4))
      dimension deps((jcap+1)*(jcap+4)),rdeps((jcap+1)*(jcap+4)),
     1          dx(2*(jcap+1)),y(jcap+1)
      dimension x(jcap+1)
c
      jcap1=jcap+1
      jcap2=jcap+2
c
      do ll=1,jcap1
        rdeps(ll) = 0.0
      enddo
c
      lplus=jcap1
      len  =jcap1
      do inde=2,jcap2
        do ll=1,len
          l=ll-1
          n=l+inde-1
          rdeps(ll+lplus)=(n*n-l*l)/(4.0*n*n-1.0)
        enddo
        lplus=lplus+len
        len=len-1
      enddo
      do i=jcap2,jcap1*jcap2/2+jcap1
        rdeps(i)=sqrt(rdeps(i))
      enddo
      do i=1,jcap1*jcap2/2+jcap1
        deps(2*i-1)=rdeps(i)
        deps(2*i  )=rdeps(i)
      enddo
      ibegin=2*jcap1+1
      do i=ibegin,jcap1*(jcap+4)
        rdeps(i)=1.0/deps(i)
      enddo
      do ll=1,jcap1
        x(ll)=ll*2+1
      enddo
      do ll=1,jcap1
        y(ll)=x(ll)/(x(ll)-1.)
      enddo
      do ll=1,jcap1
        x(ll)=sqrt(x(ll))
      enddo
      do ll=1,jcap1
        dx(2*ll-1)=x(ll)
        dx(2*ll  )=x(ll)
      enddo
c
c    set index array for transposing vector array
c   from cray order to ibm order.
c
      l=0
      do nn=1,jcap2
        lln=min0(jcap2-nn+1,jcap1)
        do ll=1,lln
          indx=((jcap+3)*(ll-1)-(ll-1)*ll/2+nn)*2
          l=l+2
          indxmv(l-1)=indx-1
          indxmv(l  )=indx
        enddo
      enddo
c
      return
      end
      subroutine ggozrm(jcap,dxab)
c
      dimension   eps(jcap+2,jcap+1)
      dimension dxint((jcap+1)*(jcap+2))
      dimension    dx((jcap+1)*2,jcap+2)
      dimension  deps((jcap+1)*2,jcap+2)
      dimension  dxab((jcap+1)*(jcap+2),2)
c
      jcap1=jcap+1
      jcap2=jcap+2
      do ll=1,jcap1
        l=ll-1
        do inde=2,jcap2
          n=l+inde-1
          a=(n*n-l*l)/(4.0*n*n-1.0)
          eps(inde,ll)=sqrt(a)
        enddo
      enddo
      do ll=1,jcap1
        eps(1,ll) = 0.0e0
      enddo
c 
      do ll=1,jcap1*2
         dxint(2*ll-1)=ll
         dxint(2*ll  )=ll
      enddo
      lp = 0
      do i=1,jcap2
        do ll=1,jcap1*2
          dx(ll,i)=dxint(ll+lp)
        enddo
        lp=lp+2
      enddo
      do i=1,jcap2
        do ll=1,jcap1
          deps(2*ll-1,i)=eps(i,ll)
          deps(2*ll  ,i)=eps(i,ll)
        enddo
      enddo
      do ll=1,jcap1*2
        dxab(ll,1)=0.0
      enddo
      lp1=jcap1*2
      len=jcap1*2-2
      do i=1,jcap
        do ll=1,len
          dxab(ll+lp1,1)= dx(ll,i+1)*deps(ll,i+1)
          dxab(ll+lp1,2)=-dx(ll,i  )*deps(ll,i+2)
        enddo
        lp1=lp1+len
        len=len-2
      enddo
      do i=1,jcap2
        do ll=1,jcap1*2
          dx(ll,i)=dx(ll,i)-1.e0
        enddo
      enddo
      do ll=1,jcap1*2
        dxab(ll,2)=-dx(ll,1)*deps(ll,2)
      enddo
c
c  transpose scalar arrays dxab from cray order to ibm order.
c
      call transo (dxab,jcap,2)
c
      return
      end
      subroutine transo(a,jcap,kmax)
c
      dimension indxnn((jcap+1)*(jcap+2)),
     1          indxmm((jcap+1)*(jcap+2))
c
      dimension a((jcap+1)*(jcap+2),kmax)
      dimension b((jcap+1)*(jcap+2))
c
      jcap1=jcap+1
      jcap2=jcap+2
c
      call cmpind(jcap,indxnn,indxmm)
c
      do k=1,kmax
        do m=1,jcap1*jcap2
          b(indxmm(m))=a(m,k)
        enddo
        do m=1,jcap1*jcap2
          a(m,k)=b(m)
        enddo
      enddo
c
      return
      end
      subroutine cmpind(jcap,indxnn,indxmm)
c
      dimension indxnn((jcap+1)*(jcap+2)),indxmm((jcap+1)*(jcap+2))
c
c  indxnn :  1-d index of converting input form spher coeff array
c                  to transposed form array
c  indxmm :  1-d index of converting transposed form spher coeff
c                  array to input form spherical coeff array
c
      jcap1=jcap+1
      l=0
      do m=1,jcap1
        nend=jcap1-m+1
        do nn=1,nend
              n=nn+m-1
          l=l+2
          indx=(jcap1*(n-m)-(n-m)*(n-m-1)/2+m)*2-1
          indxnn(l-1)=indx
          indxnn(l  )=indx+1
        enddo
      enddo
c
      l=0
      do nn=1,jcap1
        lln=jcap1-nn+1
        do ll=1,lln
          n=ll+nn-1
          m=ll
          indx=(m*jcap1-(jcap1-n)-(m-1)*m/2)*2-1
          l=l+2
          indxmm(l-1)=indx
          indxmm(l  )=indx+1
        enddo
      enddo
      return
c
      end
