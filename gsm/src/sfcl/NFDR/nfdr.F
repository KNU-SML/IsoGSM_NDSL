#define NPES npes
      subroutine nfdr(idimdl,jdimdl,iyear,month,iday,ihour,fhour,
     *                proj,orient,truth,cotru,
     *                delx,dely,rlat1,rlat2,rlon1,rlon2,
     *                dbgrd,rhgrd,pptgrd,cldgrd,wsgrd,
     *                pptdurgrd,tmaxgrd,tmingrd,rhmaxgrd,rhmingrd,
     *                condir)
c
c  nfdr main driver
c
      implicit none
c
#include <paramodel.h>
c
#ifdef MP
#include <comrmpi.h>
#endif
c
c  fixed fields (file i/o)
c
      integer idim,jdim,ijdim
      parameter(idim=184,jdim=116,ijdim=idim*jdim)
      integer ifm(idim,jdim)
      real ann(idim,jdim),slp(idim,jdim)
      real f(20,20)
      integer jyear,jmonth,jday,jhour
      real prev_fhour
c
      real*4 tmp(idim,jdim)
      integer*4 iyear4,month4,iday4,ihour4
      real*4 fhour4
#ifdef DBG
      real*4 ann4(idim,jdim),slp4(idim,jdim)
      real*4 rfm4(idim,jdim)
#endif
c
c  input (arguments)
c
      integer idimdl,jdimdl,iyear,month,iday,ihour
      real fhour
      real proj,orient,truth,cotru,delx,dely,
     1     rlat1,rlat2,rlon1,rlon2
      character*128  condir
c 
c  input forecast model variables (arguments)
c
      real dbgrd   (idimdl,jdimdl),rhgrd    (idimdl,jdimdl),
     *     pptgrd  (idimdl,jdimdl),cldgrd   (idimdl,jdimdl),
     *     wsgrd   (idimdl,jdimdl),pptdurgrd(idimdl,jdimdl),
     *     tmaxgrd (idimdl,jdimdl),tmingrd  (idimdl,jdimdl),
     *     rhmaxgrd(idimdl,jdimdl),rhmingrd (idimdl,jdimdl)
#ifdef DBG
      real*4 sdbgrd   (idimdl,jdimdl),srhgrd    (idimdl,jdimdl),
     *       spptgrd  (idimdl,jdimdl),scldgrd   (idimdl,jdimdl),
     *       swsgrd   (idimdl,jdimdl),spptdurgrd(idimdl,jdimdl),
     *       stmaxgrd (idimdl,jdimdl),stmingrd  (idimdl,jdimdl),
     *       srhmaxgrd(idimdl,jdimdl),srhmingrd (idimdl,jdimdl)
#endif
c
c  input and output (file i/o)
c
      real cumppt (idim,jdim),avethbndy(idim,jdim),
     1     ythmois(idim,jdim),fm100    (idim,jdim)
      integer kdi(idim,jdim)
c
c  output array (file i/o)
c
      integer isc(idim,jdim),ibi     (idim,jdim),ierc(idim,jdim),
     1         ic(idim,jdim),kd_class(idim,jdim)
c
c  work variables
c 
c  input forecast model variables interopolated to fire danger grid
c
      real  dbpix(idim,jdim),rhpix(idim,jdim),pptpix(idim,jdim),
     1     cldpix(idim,jdim),wspix(idim,jdim),
     2     pptdur(idim,jdim),tmax (idim,jdim),tmin  (idim,jdim),
     3      rhmax(idim,jdim),rhmin(idim,jdim)
#ifdef DBG
      real*4  sdbpix(idim,jdim),srhpix(idim,jdim),spptpix(idim,jdim),
     1        scldpix(idim,jdim),swspix(idim,jdim),
     2     spptdur(idim,jdim),stmax (idim,jdim),stmin  (idim,jdim),
     3     srhmax(idim,jdim),srhmin(idim,jdim)
#endif
      integer i,j,ij,jdate
c
      logical validmdl
      real sgbrt,rire,extmoi,fm1,tf,ros,scm
      real w1d,w1dp,w10d,w100d,w1000d,wwood,
     *     wherb,drfuel,sig1d,sig10d,sg100d,s1000d,
     *     sgwood,sgherb,hd,
     *     depth,wndftr,wfmin,wfmax,
     *     fhprev
      real fm10,bndry,fm1000,fherbc,fmwood,wherbc,phasein
      logical lmask
      real rslmsk(idimdl*jdimdl),slmask(idim*jdim)
c
      integer nfrdint
      parameter(nfrdint=24)
      character rsmprt*80
      integer nrsmprt
c
      integer iunit
      parameter(iunit=33)
c
c  interpolate from model grid to fire danger grid
c
c  land ocean mask not used
c
c  just for safe side, fill rslmsk and slmask with 1.
c
      lmask=.false.
      do ij=1,idimdl*jdimdl
        rslmsk(ij)=1.
      enddo
      do ij=1,idim*jdim
        slmask(ij)=1.
      enddo
      call nfdr_ll2rsm(dbgrd,idimdl,jdimdl,0,dbpix,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(rhgrd,idimdl,jdimdl,0,rhpix,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(pptgrd,idimdl,jdimdl,0,pptpix,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(cldgrd,idimdl,jdimdl,0,cldpix,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(wsgrd ,idimdl,jdimdl,0,wspix ,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(pptdurgrd ,idimdl,jdimdl,0,pptdur,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(tmaxgrd,idimdl,jdimdl,0,tmax,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(tmingrd,idimdl,jdimdl,0,tmin,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(rhmaxgrd,idimdl,jdimdl,0,rhmax,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
      call nfdr_ll2rsm(rhmingrd,idimdl,jdimdl,0,rhmin,
     1             lmask,rslmsk,slmask,
     2             proj,orient,truth,cotru,
     3             delx,dely,rlat1,rlat2,rlon1,rlon2)
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      do ij=1,idimdl*jdimdl
         sdbgrd    (ij,1)=dbgrd(ij,1)
         srhgrd    (ij,1)=rhgrd(ij,1)
         spptgrd   (ij,1)=pptgrd(ij,1)
         scldgrd   (ij,1)=cldgrd(ij,1) 
         swsgrd    (ij,1)=wsgrd(ij,1)
         spptdurgrd(ij,1)=pptdurgrd(ij,1)
         stmaxgrd  (ij,1)=tmaxgrd(ij,1)
         stmingrd  (ij,1)=tmingrd(ij,1) 
         srhmaxgrd (ij,1)=rhmaxgrd(ij,1)
         srhmingrd (ij,1)=rhmingrd(ij,1)
      enddo
      open(unit=66,file='forcing.before_proj',form='unformatted')
      write(66) sdbgrd    
      write(66) srhgrd    
      write(66) spptgrd   
      write(66) scldgrd   
      write(66) swsgrd    
      write(66) spptdurgrd
      write(66) stmaxgrd  
      write(66) stmingrd  
      write(66) srhmaxgrd 
      write(66) srhmingrd 
      close(66)
      do ij=1,idim*jdim
        sdbpix(ij,1)=dbpix(ij,1)
        srhpix(ij,1)=rhpix(ij,1)
        spptpix(ij,1)=pptpix(ij,1)
        scldpix(ij,1)=cldpix(ij,1)
        swspix(ij,1)=wspix(ij,1)
        spptdur(ij,1)=pptdur(ij,1)
        stmax(ij,1)=tmax(ij,1)
        stmin(ij,1)=tmin(ij,1)
        srhmax(ij,1)=rhmax(ij,1)
        srhmin(ij,1)=rhmin(ij,1)
      enddo
      open(unit=66,file='forcing.after_proj',form='unformatted')
      write(66) sdbpix
      write(66) srhpix
      write(66) spptpix
      write(66) scldpix
      write(66) swspix
      write(66) spptdur
      write(66) stmax
      write(66) stmin
      write(66) srhmax
      write(66) srhmin
      close(66)
c
#endif
c
c  change unit
c
c  from m/s to miles/hour
c
      do j=1,jdim
        do i=1,idim
          if(tmax(i,j).gt.0.) then
            wspix(i,j)=wspix(i,j)*3600./1609.3
          endif
        enddo
      enddo
c
c  from m to inch
c
      do j=1,jdim
        do i=1,idim
          if(tmax(i,j).gt.0.) then
            pptpix(i,j)=pptpix(i,j)*1000./25.4
          endif
        enddo
      enddo
c
c  from second to hour
c
      do j=1,jdim
        do i=1,idim
          if(tmax(i,j).gt.0.) then
            pptdur(i,j)=pptdur(i,j)/3600.
          endif
        enddo
      enddo
c
c  from degree K to degree F
c
      do j=1,jdim
        do i=1,idim
          if(tmax(i,j).gt.0.) then
            dbpix(i,j)=(dbpix(i,j)-273.16)*1.80+32.0
            tmax (i,j)=(tmax (i,j)-273.16)*1.80+32.0
            tmin (i,j)=(tmin (i,j)-273.16)*1.80+32.0
          endif
        enddo
      enddo
c
      call nfdr_julian(month,iday,jdate)
      jdate=jdate-1
#ifdef MP
c
      if(mype.eq.master) then
#endif
c
c  read fixed field
c
      call nfdr_getfixf(iunit,condir,f,idim,jdim,ann,ifm,slp)
      print *,'fixed fields for nfdr read in'
#ifdef DBG
      open(unit=88,file='./nfdrfix',form='unformatted')
      do j=1,jdim
        do i=1,idim
          ann4(i,j)=ann(i,j)
          slp4(i,j)=slp(i,j)
          rfm4(i,j)=ifm(i,j)
        enddo
      enddo
      write(88) ann4
      write(88) slp4
      write(88) rfm4
#endif
c
c  read restart file
c
      open (unit=iunit,file='nfdr.restart',
     *      form='unformatted',status='old',err=800)
      read(iunit) jyear,jmonth,jday,jhour,prev_fhour
      print *,'year,month,day,hour,fhour of restart=',
     1         jyear,jmonth,jday,jhour,prev_fhour
      read(iunit) cumppt
      read(iunit) avethbndy
      read(iunit) ythmois
      read(iunit) fm100
      read(iunit) kdi
      print *,' nfdr restartf file nfdr.restart read in'
      go to 801
  800 continue
      print *,'Warning.  nfdr_init called.  This is dead start for nfdr'
      call nfdr_init(idim,jdim,
     1              cumppt,avethbndy,ythmois,kdi,fm100)

  801 continue
c
      do j = 1,jdim
        do i = 1,idim
          ij=(j-1)*idim+i
c
c  missing forcing
c
          if(tmax(i,j).eq.-999.) then
            isc(i,j)=-999
            ibi(i,j)=-999
            ierc(i,j)=-999
            ic(i,j)=-999
            cumppt(i,j)=-999.
            avethbndy(i,j)=-999.
            ythmois(i,j)=-999.
            fm100(i,j)=-999.
            kdi(i,j)=-999
            kd_class(i,j)=-999
            go to 100
          endif
c
c  Determine if this pixel represents a valid fuel model
c  The following ifm values are not valid -- 0(no data),
c  13(ag land), 22(water), 23(barren), 24(marsh).
c
          validmdl = .false.
          if(ifm(i,j).ge. 1.and.ifm(i,j).le.12.or.
     *       ifm(i,j).ge.14.and.ifm(i,j).le.21) validmdl = .true.
#ifdef DBG
c
c  diagnstics
c
          if(ij.eq.8000)then
             write(6,17) ifm(i,j),cumppt(i,j),avethbndy(i,j),
     1                   ythmois(i,j),kdi(i,j),fm100(i,j)
  17         format(' ',i5,3(f6.2,1x),i3,1x,f6.2,1x)
          endif
#endif  
c
c  Set each index to a characteristic value for those pixels that 
c  do not have a valid valid fuel model.
c
          if(.not. validmdl) then
             call nfdr_nofuel(ifm(i,j),isc(i,j),
     1                        ierc(i,j),ibi(i,j),ic(i,j))
          endif
c
c  calc kbdi, but overwrite for water (ifm=22) or outside US (ifm=0)
c
          call nfdr_kbdi(pptpix(i,j),ann(i,j),cumppt(i,j),tmax(i,j),
     1              kdi(i,j))
          kd_class(i,j)=int((kdi(i,j)-1)/50) + 1
          if(ifm(i,j).eq.22) kd_class(i,j)=17
          if(ifm(i,j).eq. 0) kd_class(i,j)=18
c
c  Assign fuel data and do calculations only for valid fuel models.
c
          if(validmdl)then
            call nfdr_setloads(ifm(i,j),f,w1d,w1dp,
     1                         w10d,w100d,w1000d,wwood,
     *                    wherb,drfuel,sig1d,sig10d,sg100d,s1000d,
     *                    sgwood,sgherb,hd,
     *                    scm,extmoi,depth,wndftr,wfmin,wfmax,kdi(i,j))
            call nfdr_m100(jdate,pptdur(i,j),tmax(i,j),tmin(i,j),
     1                rhmax(i,j),rhmin(i,j),ij,bndry,fm100(i,j))
            call nfdr_oneten(dbpix(i,j),rhpix(i,j),cldpix(i,j),
     1                 fm100(i,j),tf,fm1,fm10)
            call nfdr_m1000(bndry,avethbndy(i,j),ythmois(i,j),fm1000)
            call nfdr_curing(fm1,ann(i,j),fm1000,wherb,fherbc,fmwood,
     *                  w1dp,wherbc)
            call nfdr_ercsc(w1dp,w10d,w100d,w1000d,wwood,wherbc,sig1d,
     *                 sig10d,sg100d,s1000d,sgherb,sgwood,depth,fm1,
     *                 fm10,fm100(i,j),fm1000,fmwood,fherbc,extmoi,hd,
     *                 wspix(i,j),wndftr,slp(i,j),isc(i,j),
     *                 ierc(i,j),rire,sgbrt,ros)
            call nfdr_bi_ic(sgbrt,rire,extmoi,fm1,tf,ros,scm,
     *                 ibi(i,j),ic(i,j))
c        
c  Adjust indexes for winter conditions.  Assume winter if max dry bulb
c  is < 40 deg F, transition indexes from 40 to 50 deg F, and leave
c  them as calculated if the max temperature is > 50 deg F.
c
            if(tmax(i,j) .lt. 40.) then
               isc(i,j)=0
               ibi(i,j)=0
               ierc(i,j)=0
               ic(i,j)=0
            endif
c
            if(tmax(i,j) .ge. 40.0 .and. tmax(i,j) .le. 50.0) then
               phasein = (tmax(i,j) - 40.0)/10.0
               isc(i,j) = int(0.0 + float(isc(i,j)) * phasein + 0.5)
               ibi(i,j) = int(0.0 + float(ibi(i,j)) * phasein + 0.5)
               ierc(i,j) = int(0.0 + float(ierc(i,j)) * phasein + 0.5)
               ic(i,j) = int(0.0 + float(ic(i,j)) * phasein + 0.5)
            endif
          endif 
c
#ifdef DBG            
          if(ij .eq. 8000) then
            write(6,18)ifm(i,j),cumppt(i,j),avethbndy(i,j),
     1                 ythmois(i,j),kdi(i,j),fm100(i,j),fm1000
   18       format(' ',i5,3(f6.2,1x),i3,1x,2(f6.2,1x))        
          endif
#endif
  100     continue
        enddo
      enddo
c
c  output to be gribbed?
c
      call fnam('nfdr.ieee',9,fhour,rsmprt,nrsmprt)
      open (unit=iunit,file=rsmprt(1:nrsmprt),
     1      form='unformatted',status='unknown')
      iyear4=iyear
      month4=month
      iday4=iday
      ihour4=ihour
      fhour4=fhour
      write(iunit) iyear4,month4,iday4,ihour4,fhour4
      do j=1,jdim
        do i=1,idim
          tmp(i,j)=isc(i,j)
        enddo
      enddo
      write(iunit) tmp
      do j=1,jdim
        do i=1,idim
          tmp(i,j)=ibi(i,j)
        enddo
      enddo
      write(iunit) tmp
      do j=1,jdim
        do i=1,idim
          tmp(i,j)=ic(i,j)
        enddo
      enddo
      write(iunit) tmp
      do j=1,jdim
        do i=1,idim
          tmp(i,j)=ierc(i,j)
        enddo
      enddo
      write(iunit) tmp
      do j=1,jdim
        do i=1,idim
          tmp(i,j)=kd_class(i,j)
        enddo
      enddo
      write(iunit) tmp
      close(iunit)
      print *,'nfdr file ',rsmprt(1:nrsmprt),' written'
c              
c  restart file for archive              
c              
      call fnam('nfdr.restart',12,fhour,rsmprt,nrsmprt)
      open (unit=iunit,file=rsmprt(1:nrsmprt),
     1      form='unformatted',status='unknown')
      write(iunit) iyear,month,iday,ihour,fhour
      write(iunit) cumppt
      write(iunit) avethbndy
      write(iunit) ythmois
      write(iunit) fm100
      write(iunit) kdi
      close(iunit)
      print *,'nfdr restart file ',rsmprt(1:nrsmprt),' written'
c              
c  also overwrites restart file              
c              
      open (unit=iunit,file='nfdr.restart',
     1      form='unformatted',status='unknown')
      write(iunit) iyear,month,iday,ihour,fhour
      write(iunit) cumppt
      write(iunit) avethbndy
      write(iunit) ythmois
      write(iunit) fm100
      write(iunit) kdi
      close(iunit)
      print *,'nfdr restart file nfdr.restart overwritten'
#ifdef MP
c
c  endif for mype .eq. master
c
      endif
#endif
c
      return
      end
