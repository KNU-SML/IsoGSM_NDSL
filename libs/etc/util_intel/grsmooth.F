      program grsmooth
c
c  grib space smoother program 
c 
c  mtotl ... maxmum number of fields in a given grib file
c  lkpds ... number of pds words to be recorded 
c 
      character*64, allocatable :: label(:,:) 
      integer, allocatable :: lenla(:,:)
      character*8, allocatable :: labbr(:)
      character*8, allocatable :: lvabbr(:)
      real, allocatable :: rlev(:)
      integer, allocatable :: ilindx(:)
      integer, allocatable :: ivindx(:)
      integer, allocatable :: ipindx(:)
      integer, allocatable :: nltotl(:)
      logical llevl,lltyp
c
      real, allocatable :: gridout(:,:)
c
      integer, allocatable :: kpds(:,:),kgds(:,:)
      integer, allocatable :: lskip(:),lgrib(:)
c
      integer kpdsout(25),kgdsout(22)
c
c  fngrib ... grib file name(s). 
c             dimension for history sequence (in forecast time or initial)
c  ldebug  ... true for debugging printout. output goes to fnprint
c
c  namelist variable types
c 
      character*128, allocatable :: fngrib(:)
      character*128 fngrbo
c
      data msk1/32000/,msk2/4000/
c
      character*80 fnkpds5,fnkpds6
c
      data ldebug/.true./
c
c  fort unit numbers used inside the program
c
c  luptr  ... unit unmber for diagnostic/debug print output other than uniplo
c  lupgb  ... unit unmber for input (grib) file
c  lupgo  ... unit unmber for grib output file
c
      data lupgb/21/
      data lupgo/51/
      data luptr/50/
c
      parameter(mbuf=1024*128*64)
      character*1 cbuf(mbuf),dbuf(mbuf)
      character, allocatable :: bbuf(:)
c
      real, allocatable :: grid(:)
      logical*1, allocatable :: lbms(:)
      real*8, allocatable :: fr(:)
      integer, allocatable :: ibm(:)
      real*4 , allocatable :: grid4(:)
      character, allocatable ::  grib(:)
c
      logical ldbg
      common/comdbg/ ldbg
c
      parameter(lenpds=28,lengds=32,mxbit=16)
c
      integer*4 jpds4(25),jgds4(22),jens4(5)
      integer*4 mpds4(25),mgds4(22),mens4(5)
      integer*4 lupgb4,luptr4,msk14,msk24,mnum4,mbuf4
      integer*4 nlen4,nnum4,iret4,ndata4
      integer*4, allocatable :: kpds4(:,:),kgds4(:,:)
      integer*4, allocatable :: lskip4(:),lgrib4(:)
c
c  read main program namelist
c
      read (5,'(a80)') fnkpds5
      read (5,'(a80)') fnkpds6
      write(6,'(a80)') fnkpds5
      write(6,'(a80)') fnkpds6
c
      ldbg=ldebug
c
      read (5,*,end=908) nfiles
      go to 907
  908 continue
      write(6,*) ' eof on data'
      call abort
  907 continue
      write(6,*) ' nfiles=',nfiles
      if(nfiles.eq.0) then
        call abort
      endif 
      print *,'nfiles=',nfiles
      allocate (fngrib(nfiles))
c
      do nfile=1,nfiles
        fngrib(nfile)='        '
      enddo
c
      do nfile=1,nfiles
        read (5,'(a128)') fngrib(nfile)
        write(6,'(8h fngrib=,a80)') fngrib(nfile)
      enddo
c
      call iniparm(fnkpds5)
      call inilevl(fnkpds6)
c
c  loop through grib file
c
      do nfile=1,nfiles
c
        print *,' processing ', fngrib(nfile)(1:80)
c
c  get grib index buffer
c
        call numchar(fngrib(nfile),nch)
        call nainit(lupgb,fngrib(nfile)(1:nch),isize,iret)
        allocate (bbuf(isize))
        call naopen(lupgb,fngrib(nfile)(1:nch),bbuf,isize,iret)
c
c  build variable name table from the first grib file
c
        rewind lupgb
        mnum=0
        lupgb4=lupgb
        msk14=msk1
        msk24=msk2
        mnum4=mnum
        mbuf4=mbuf
c
        call ngetgir(bbuf,isize,msk14,msk24,mnum4,mbuf4,
     1               cbuf,nlen4,nnum4,iret4)
        nlen=nlen4
        nnum=nnum4
        iret=iret4
        if(iret.ne.0) then
          print *,'error.  cbuf length too short in getgir'
          call abort
        endif
        if(nnum.eq.0) then
          print *,'error.  not a grib file found in getgir'
          call abort
        endif
        if(nlen.eq.0) then
          print *,'error.  nlen=0 found in getgir'
          call abort
        endif
c
        mtotl=nnum
        allocate (kpds(25,mtotl),kgds(22,mtotl))
        allocate (lskip(mtotl),lgrib(mtotl))
        allocate (kpds4(25,mtotl),kgds4(22,mtotl))
        allocate (lskip4(mtotl),lgrib4(mtotl))
        allocate (labbr(mtotl))
c
        do n=1,mtotl
          labbr(n)='        '
        enddo
c
        nlen4=nlen
        nnum4=nnum
        call unpindx(cbuf,nlen4,nnum4,kpds4,kgds4,lskip4,lgrib4,iret4)
        print *,'KPDS4(8,1)=',KPDS4(8,1)
        print *,'KPDS4(21,1)=',KPDS4(21,1)
        do n=1,nnum4
c         if(kpds4(8,n).eq.100) then
c           kpds4(21,n)=kpds4(21,n)-1
c         endif
          kpds4(8,n)=kpds4(8,n)+(kpds4(21,n)-1)*100
c         kpds4(8,n)=kpds4(8,n)+1900
          kpds4(25,n)=-1
          do i=1,25
            kpds(i,n)=kpds4(i,n)
          enddo
          do i=1,22
            kgds(i,n)=kgds4(i,n)
          enddo
          lskip(n)=lskip4(n)
          lgrib(n)=lgrib4(n)
        enddo
        iret=iret4
c
        imax=kgds4(2,1)
        jmax=kgds4(3,1)
        ijmax=imax*jmax
        if(ldbg) then
          print *,'imax=',imax,' jmax=',jmax,' ijmax=',ijmax
        endif
c
        mdata=ijmax
        allocate (grid(mdata))
        allocate (gridout(imax,jmax))
        allocate (lbms(mdata))
        allocate (fr(mdata))
        allocate (ibm(mdata))
        allocate (grid4(mdata))
        allocate (grib(30+lenpds+lengds+mdata*(mxbit+1)/8))
c
        mtotl=nnum4
        allocate (label(3,mtotl))
        allocate (lenla(3,mtotl))
        allocate (labbr(mtotl))
        allocate (lvabbr(mtotl))
        allocate (rlev(mtotl))
        allocate (ilindx(mtotl))
        allocate (ivindx(mtotl))
        allocate (ipindx(mtotl))
        allocate (nltotl(mtotl))
c
        ntotl=nnum4
c
        do n=1,ntotl
          call getparm(kpds(1,n),1,label(1,n),lenla(1,n),labbr(n))
          call getlevl(kpds(1,n),1,label(1,n),lenla(1,n),rlev(n),lvabbr)
        enddo
c
c  get rid of unnecessary spaces from label
c
        do n=1,ntotl
          call rmblnk(label(1,n),lenla(1,n),3)
        enddo
c
        call flsrch(ntotl,luptr,
     1      mtotl,label,lenla,labbr,lvabbr,rlev,
     2      ilindx,ivindx,ipindx,
     3      nltotl,nvtotl,nptotl,llevl,lltyp)
c
        if(ldbg) then
          do n=1,ntotl
            print *,n,')',label(1,n)(1:lenla(1,n)),
     1                            ' (',labbr(n),' ) ',
     2                ' at ',label(2,n)(1:lenla(2,n)),' ',
     3                       label(3,n)(1:lenla(3,n))
          enddo
        endif
c
        mnum=0
c
        do ij=1,ijmax
          lbms(ij)=.true.
        enddo
c
        rewind lupgb
        mnum4=0
        lupgb4=lupgb
        msk14=msk1
        msk24=msk2
        mbuf4=mbuf
c
        nn4=0
        do n=1,ntotl
          nn4=0
          do i=1,25
            jpds4(i)=kpds4(i,n)
          enddo
          do i=1,22
            jgds4(i)=kgds4(i,n)
          enddo
          do i=1,5
            jens4(i)=-1
          enddo
          nlen4=nlen
          nnum4=nnum
c
          call getgbss(cbuf,nlen4,nnum4,nn4,jpds4,jgds4,jens4,
     &                 k4,mpds4,mgds4,mens4,mskip4,mgrib4,iret4)
c
          if(ldbg) then
            print *,' debug print.'
            print *,'n=',n,' jpds4=',(jpds4(i),i=1,25)
            print *,'n=',n,' mpds4=',(mpds4(i),i=1,25)
            print *,' '
            print *,'n=',n,' jgds4=',(jgds4(i),i=1,22)
            print *,'n=',n,' mgds4=',(mgds4(i),i=1,22)
          endif
c
          if(mgrib4.eq.0) then
            print *,' error in getgbss.  field not found.'
            call abort
          endif
c
          if(ldbg) then
            print *,'ymdhf=',
     &            mpds4(8),mpds4(9),mpds4(10),mpds4(11),mpds4(14)
          endif
c
          do ij=1,ijmax
            lbms(ij)=.true.
          enddo
c
          lupgb4=lupgb
          luptr4=luptr
c
          call nrdgb(bbuf,isize,mgrib4,mskip4,
     1              mpds4,mgds4,ndata4,lbms,grid4,luptr4)
c
          if(ndata4.ne.ijmax) then
            write(luptr,*) ' error in rdgb. ndata.ne.ijmax'
            write(luptr,*) ' ndata=',ndata4,' ijmax=',ijmax
            write(luptr,*) ' nfile=',nfile
            write(luptr,*) ' fngrib=',fngrib(nfile)
            print *,' error in rdgb.   ndata.ne.ijmax'
            call abort
          endif
c
c  this is the place smoothing should take place
c
          do j=1,jmax
            do i=1,imax
              ijind=imax*(j-1)+i
              gridout(i,j)=grid4(ijind)
            enddo
          enddo
          call bdrysmth(imax,jmax,gridout)
c
c  open output grib file
c
          if(n.eq.1) then
            fngrbo=fngrib(nfile)(1:nch)//'_s'
            open(unit=lupgo,file=fngrbo,form='unformatted',err=968)
            go to 967
  968       continue
            print *,'error in opening file ',fngrbo
            call abort
  967       continue
            print *,' file ',fngrbo(1:128),' opened. unit=',lupgo
          endif
c
          do i=1,25
            kpdsout(i)=kpds(i,n)
          enddo
c         if(kpdsout(8).eq.100) then
c            kpdsout(21)=kpdsout(21)+1  
c         endif
c         kpdsout(8)=mod(kpdsout(8),100)
          do i=1,22
            kgdsout(i)=kgds(i,n)
          enddo
          if(ldbg) then
            print *, 'ijmax=',ijmax
            print *, 'kpdsout=',(kpdsout(i),i=1,25)
            print *, 'kgdsout=',(kgdsout(i),i=1,22)
          endif
          call putgb(lupgo,ijmax,kpdsout,kgdsout,lbms,
     1               gridout,ibm,fr,grib,igrib,luptr,iret)
          if(iret.ne.0) then
            print *, 'error in putgb'
            print *, 'kpdsout=',(kpdsout(i),i=1,25)
            print *, 'kgdsout=',(kgdsout(i),i=1,22)
            call abort
          else
            if(ldbg) then
              print *, 'grib-end=',(grib(i),i=igrib-3,igrib)
            endif
            print *, n,') grib record written. lgrib=',igrib
          endif
        enddo
        close(lupgo)
        deallocate (grid)
        deallocate (gridout)
        deallocate (lbms)
        deallocate (fr)
        deallocate (ibm)
        deallocate (grid4)
        deallocate (grib)
        deallocate (kpds,kgds)
        deallocate (lskip,lgrib)
        deallocate (kpds4,kgds4)
        deallocate (lskip4,lgrib4)
        deallocate (labbr)
        deallocate (bbuf)
        deallocate (label)
        deallocate (lenla)
        deallocate (lvabbr)
        deallocate (rlev)
        deallocate (ilindx)
        deallocate (ivindx)
        deallocate (ipindx)
        deallocate (nltotl)
      enddo
      close(luptr)
c
      stop
      end
      SUBROUTINE INIPARM(FNKPDS5)
C
C     Initialize KPDS(5) table
C
      PARAMETER(LUKP5=90)
C
      CHARACTER*80 FNKPDS5
C
      PARAMETER (NKPDS5=255)
      CHARACTER*64 PINFO(0:NKPDS5)
      CHARACTER*64 PNAME(0:NKPDS5)
      CHARACTER*64 PUNIT(0:NKPDS5)
      DIMENSION LENPI(0:NKPDS5),LENPN(0:NKPDS5),LENPU(0:NKPDS5)
C
      COMMON/COMPRM/ PINFO,PUNIT,PNAME,LENPI,LENPU,LENPN
C
      CHARACTER*80 CPDS5
      CHARACTER*80 STRING
C
C  Open Parameter table (PDS5) file
C
      OPEN(LUKP5,FILE=FNKPDS5,STATUS='OLD',ERR=1)
      GO TO 2
    1 CONTINUE
      print *, 'ERROR IN OPENING FILE ',FNKPDS5
      CALL ABORT
    2 CONTINUE
      print *, 'FILE ',FNKPDS5(1:50),' opened.  Unit=',LUKP5
C
      N=0
  300 CONTINUE
      READ(LUKP5,100,END=200) STRING
  100 FORMAT(A80)
      IS=1
			NDELIM=0
      DO 500 I=1,80
			IF(STRING(I:I).EQ.':') THEN
        NDELIM=NDELIM+1
        NCHAR=I-1-IS+1
        IF(NCHAR.GT.0) THEN
          IF(NDELIM.EQ.1) THEN
            CPDS5(1:NCHAR)=STRING(IS:I-1)
						READ(CPDS5,110) M
  110       FORMAT(I3)
          ELSEIF(NDELIM.EQ.2) THEN
						JJ=0
            DO 505 J=IS+1,I-1
            IF(STRING(J-1:J-1).NE.' '.OR.STRING(J:J).NE.' ') THEN
	  					JJ=JJ+1
              PINFO(M)(JJ:JJ)=STRING(J-1:J-1)
						ENDIF
  505       CONTINUE
            IF(STRING(I-1:I-1).NE.' ') THEN
	  					JJ=JJ+1
              PINFO(M)(JJ:JJ)=STRING(I-1:I-1)
						ENDIF
            LENPI(M)=JJ
          ELSEIF(NDELIM.EQ.3) THEN
						JJ=0
            DO 510 J=IS,I-1
            IF(STRING(J:J).NE.' ') THEN
						  JJ=JJ+1
              PUNIT(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  510       CONTINUE
            LENPU(M)=JJ
						JJ=0
            DO 520 J=I+1,80
            IF(STRING(J:J).NE.' ') THEN
						  JJ=JJ+1
              PNAME(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  520       CONTINUE
            LENPN(M)=JJ
						GO TO 500
          ENDIF
  				IS=I+1
				ENDIF
			ENDIF
  500 CONTINUE
C
			IF(NDELIM.EQ.1) THEN
				PINFO(M)(1:1)='?'
				PUNIT(M)(1:1)='?'
				PNAME(M)(1:1)='?'
        LENPI(M)=1
        LENPU(M)=1
        LENPN(M)=1
			ELSEIF(NDELIM.EQ.2) THEN
				PUNIT(M)(1:1)='?'
				PNAME(M)(1:1)='?'
        LENPU(M)=1
        LENPN(M)=1
			ENDIF
C
      N=N+1
      GO TO 300
  200 CONTINUE
C
      CLOSE(UNIT=LUKP5)
      RETURN
      END
      SUBROUTINE INILEVL(FNKPDS6)
C
      PARAMETER (LUKP6=90)
      PARAMETER (NKPDS6=255)
      CHARACTER*80  PDS6 (0:NKPDS6)
C
      CHARACTER*64 LINFO(0:NKPDS6)
      CHARACTER*64 LUNIT(0:NKPDS6)
      CHARACTER*64 LFACT(0:NKPDS6)
      CHARACTER*64 LBASE(0:NKPDS6)
      CHARACTER*64 LTYPE(0:NKPDS6)
      CHARACTER*64 LEVLR(0:NKPDS6)
      CHARACTER*3  LEVAB(0:NKPDS6)
C
      DIMENSION LENLI(0:NKPDS6),LENLU(0:NKPDS6),LENLF(0:NKPDS6),
     1          LENLB(0:NKPDS6),LENLT(0:NKPDS6),LENLE(0:NKPDS6)
C
      COMMON/COMLVL/ LINFO,LUNIT,LFACT,LBASE,LTYPE,LEVLR,LEVAB,
     1               LENLI,LENLU,LENLF,LENLB,LENLT,LENLE
C
      CHARACTER*80 FNKPDS6
C
      CHARACTER*160 STRING
C
C  Open Level table (PDS6) file
C
      OPEN(LUKP6,FILE=FNKPDS6,STATUS='OLD',ERR=761)
      GO TO 762
  761 CONTINUE
      print *,'ERROR IN OPENING FILE ',FNKPDS6
      CALL ABORT
  762 CONTINUE
      print *, 'FILE ',FNKPDS6(1:50),' opened.  Unit=',LUKP6
C
      N=1
  300 CONTINUE
      READ(LUKP6,100,END=200) STRING
  100 FORMAT(A80)
      IS=1
			NDELIM=0
      DO 500 I=1,80
			IF(STRING(I:I).EQ.':') THEN
        NDELIM=NDELIM+1
        NCHAR=I-1-IS+1
        IF(NCHAR.GT.0) THEN
          IF(NDELIM.EQ.1) THEN
            PDS6(N-1)(1:NCHAR)=STRING(IS:I-1)
            READ(PDS6(N-1),110) M
  110       FORMAT(I3)
          ELSEIF(NDELIM.EQ.2) THEN
						JJ=0
            DO 505 J=IS+1,I-1
            IF(STRING(J-1:J-1).NE.' '.OR.STRING(J:J).NE.' ') THEN
	  					JJ=JJ+1
              LINFO(M)(JJ:JJ)=STRING(J-1:J-1)
						ENDIF
  505       CONTINUE
            IF(STRING(I-1:I-1).NE.' ') THEN
	  					JJ=JJ+1
              LINFO(M)(JJ:JJ)=STRING(I-1:I-1)
						ENDIF
						LENLI(M)=JJ
          ELSEIF(NDELIM.EQ.3) THEN
						JJ=0
            DO 506 J=IS,I-1
            IF(STRING(J-1:J-1).NE.' ') THEN
	  					JJ=JJ+1
              LUNIT(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  506       CONTINUE
						LENLU(M)=JJ
            IF(LUNIT(M)(1:1).EQ.'-') THEN
              LUNIT(M)(1:1)=' '
              LENLU(M)=1
            ENDIF
          ELSEIF(NDELIM.EQ.4) THEN
						JJ=0
            DO 507 J=IS,I-1
            IF(STRING(J-1:J-1).NE.' ') THEN
	  					JJ=JJ+1
              LFACT(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  507       CONTINUE
						LENLF(M)=JJ
          ELSEIF(NDELIM.EQ.5) THEN
						JJ=0
            DO 508 J=IS,I-1
            IF(STRING(J-1:J-1).NE.' ') THEN
	  					JJ=JJ+1
              LBASE(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  508       CONTINUE
						LENLB(M)=JJ
          ELSEIF(NDELIM.EQ.6) THEN
						JJ=0
            DO 510 J=IS,I-1
            IF(STRING(J:J).NE.' ') THEN
						  JJ=JJ+1
              LTYPE(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  510       CONTINUE
						LENLT(M)=JJ
          ELSEIF(NDELIM.EQ.7) THEN
						JJ=0
            DO 520 J=IS,I-1
            IF(STRING(J:J).NE.' ') THEN
						  JJ=JJ+1
              LEVLR(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  520       CONTINUE
						LENLE(M)=JJ
						JJ=0
C           DO 530 J=I+1,80
            DO 530 J=I+1,I+3
            IF(STRING(J:J).NE.' ') THEN
						  JJ=JJ+1
              LEVAB(M)(JJ:JJ)=STRING(J:J)
						ENDIF
  530       CONTINUE
						GO TO 500
          ENDIF
  				IS=I+1
				ENDIF
			ENDIF
  500 CONTINUE
C
			IF(NDELIM.EQ.1) THEN
				LINFO(M)(1:1)='?'
				LUNIT(M)(1:1)='?'
        LFACT(M)(1:1)='?'
        LBASE(M)(1:1)='?'
				LTYPE(M)(1:1)='?'
				LEVLR(M)(1:1)='?'
				LEVAB(M)(1:3)='???'
        LENLI(M)=1
        LENLU(M)=1
        LENLF(M)=1
        LENLB(M)=1
        LENLT(M)=1
        LENLE(M)=1
			ELSEIF(NDELIM.EQ.2) THEN
				LUNIT(M)(1:1)='?'
        LFACT(M)(1:1)='?'
        LBASE(M)(1:1)='?'
				LTYPE(M)(1:1)='?'
				LEVLR(M)(1:1)='?'
				LEVAB(M)(1:3)='???'
        LENLU(M)=1
        LENLF(M)=1
        LENLB(M)=1
        LENLT(M)=1
        LENLE(M)=1
			ELSEIF(NDELIM.EQ.3) THEN
        LFACT(M)(1:1)='?'
        LBASE(M)(1:1)='?'
				LTYPE(M)(1:1)='?'
				LEVLR(M)(1:1)='?'
				LEVAB(M)(1:3)='???'
        LENLF(M)=1
        LENLB(M)=1
        LENLT(M)=1
        LENLE(M)=1
			ELSEIF(NDELIM.EQ.4) THEN
        LBASE(M)(1:1)='?'
				LTYPE(M)(1:1)='?'
				LEVLR(M)(1:1)='?'
				LEVAB(M)(1:3)='???'
        LENLB(M)=1
        LENLT(M)=1
        LENLE(M)=1
			ELSEIF(NDELIM.EQ.5) THEN
				LTYPE(M)(1:1)='?'
				LEVLR(M)(1:1)='?'
				LEVAB(M)(1:3)='???'
        LENLT(M)=1
        LENLE(M)=1
			ELSEIF(NDELIM.EQ.6) THEN
				LEVLR(M)(1:1)='?'
				LEVAB(M)(1:3)='???'
        LENLE(M)=1
			ENDIF
C
      N=N+1
      GO TO 300
  200 CONTINUE
C
      CLOSE(UNIT=LUKP6)
C
      RETURN
      END
      SUBROUTINE GETPARM(KPDS,NGI,LABEL,LENLA,LABBR)
C
      DIMENSION KPDS(22,*)
C
C     Transrates KPDS(5,*) parameter to variable names
C
      CHARACTER*64 LABEL(3,*)
      DIMENSION LENLA(3,*)
      CHARACTER*8 LABBR(*)
C
      PARAMETER (NKPDS5=255)
      CHARACTER*64 PINFO(0:NKPDS5)
      CHARACTER*64 PNAME(0:NKPDS5)
      CHARACTER*64 PUNIT(0:NKPDS5)
      DIMENSION LENPI(0:NKPDS5),LENPN(0:NKPDS5),LENPU(0:NKPDS5)
C
      SAVE NUNDEF
      DATA NUNDEF/0/
C
      COMMON/COMPRM/ PINFO,PUNIT,PNAME,LENPI,LENPU,LENPN
C
      LOGICAL LDEBUG
      COMMON/COMDBG/ LDEBUG
C
      DO N=1,NGI
				K=KPDS(5,N)
C
				IF(K.GT.NKPDS5) THEN
           print *,' FATAL ERROR in INDEX file! ',
     1          ' KPDS(5,N)=',KPDS(5,N),' FOR N=',N
           print *,' KPDS(5,N) CANNOT BE GREATER THAN ',NKPDS5
           PRINT *,' FATAL ERROR in INDEX file! ',
     1          ' KPDS(5,N)=',KPDS(5,N),' FOR N=',N
           PRINT *,' KPDS(5,N) CANNOT BE GREATER THAN ',NKPDS5
           CALL ABORT
        ENDIF
        IF(PNAME(K)(1:1).EQ.'?'.AND.PINFO(K)(1:1).EQ.'?'
     1     .AND.PUNIT(K)(1:1).EQ.'?') THEN
           print *,' SERIOUS PROBLEM in INDEX file',
     1           ' KPDS(5,N)=',KPDS(5,N),' FOR N=',N
           print *,' THE PARAMETER IS NOT DEFINED ',
     1           'IN THE KPDS(5) TABLE FILE '
           PRINT *,' Serious problem in INDEX file',
     1           ' KPDS(5,N)=',KPDS(5,N),' FOR N=',N
           PRINT *,' THE PARAMETER IS NOT DEFINED ',
     1           'IN THE KPDS(5) TABLE FILE '
           NUNDEF=NUNDEF+1
           IF(NUNDEF.GE.100) THEN
             PRINT *,'Too many undefined variables'
             CALL ABORT
           ENDIF
           IF(NUNDEF.LT.10) THEN
             WRITE(PINFO(K),'(6HUNDEF0,I1)') NUNDEF
           ELSEIF(NUNDEF.LT.100) THEN
             WRITE(PINFO(K),'(5HUNDEF,I2)') NUNDEF
           ENDIF
           LENPI(K)=7
           PNAME(K)=PINFO(K)
           LENPN(K)=7
           PUNIT(K)='?'
           LENPU(K)=1
           print *,'NUNDEF=',NUNDEF,' PINFO(K)=',PINFO(K),
     1                    ' PNAME(K)=',PNAME(K),' PUNIT(K)=',PUNIT(K)
		    ENDIF
C
        LABEL(1,N)=PINFO(K)(1:LENPI(K))//' ('//PUNIT(K)(1:LENPU(K))//')'
        LENLA(1,N)=LENPI(K)+2+LENPU(K)+1
C
        LABBR(N)='        '
        LABBR(N)=PNAME(K)(1:LENPN(K))
      ENDDO
C
      RETURN
			END
      SUBROUTINE GETLEVL(KPDS,NGI,LABEL,LENLA,RLEV,LVABBR,LUPTR)
C
C     Transrates KPDS(6,*) parameter to level names
C
      CHARACTER*64 LABEL(3,*)
      DIMENSION LENLA(3,*)
      DIMENSION RLEV(*)
      CHARACTER*8 LVABBR(*)
C
      DIMENSION KPDS(22,*)
C
      CHARACTER*160 FMT
C
      PARAMETER (NKPDS6=255)
      CHARACTER*64 LINFO(0:NKPDS6)
      CHARACTER*64 LUNIT(0:NKPDS6)
      CHARACTER*64 LFACT(0:NKPDS6)
      CHARACTER*64 LBASE(0:NKPDS6)
      CHARACTER*64 LTYPE(0:NKPDS6)
      CHARACTER*64 LEVLR(0:NKPDS6)
      CHARACTER*3  LEVAB(0:NKPDS6)
C
      DIMENSION LENLI(0:NKPDS6),LENLU(0:NKPDS6),LENLF(0:NKPDS6),
     1          LENLB(0:NKPDS6),LENLT(0:NKPDS6),LENLE(0:NKPDS6)
C
      COMMON/COMLVL/ LINFO,LUNIT,LFACT,LBASE,LTYPE,LEVLR,LEVAB,
     1               LENLI,LENLU,LENLF,LENLB,LENLT,LENLE
C
      CHARACTER*64 LTOPL,LBTML
C
      SAVE NUNDEF
      DATA NUNDEF/0/
C
      DO 100 N=1,NGI
				K=KPDS(6,N)
				IF(K.GT.NKPDS6.OR.K.LT.0) THEN
           print *,' FATAL ERROR in INDEX FILE! ',
     1          'KPDS(6,N)=',KPDS(6,N),' FOR N=',N
           print *,' KPDS(6,N) CANNOT BE GREATER THAN ',NKPDS6
           CALL ABORT
        ENDIF
        IF(LINFO(K)(1:1).EQ.'?'.AND.LUNIT(K)(1:1).EQ.'?'
     1     .AND.LEVLR(K)(1:1).EQ.'?') THEN
           print *,' Serious problem in index FILE! ',
     1        'KPDS(6,N)=',KPDS(6,N),' FOR N=',N
           print *,' THE LEVEL IS NOT DEFINED ',
     1        'IN THE FILE KPDS(6) TABLE FILE.'
           NUNDEF=NUNDEF+1
           IF(NUNDEF.GE.100) THEN
           ENDIF
           IF(NUNDEF.LT.10) THEN
             WRITE(LINFO(K),'(6HUNDEF0,I1)') NUNDEF
           ELSEIF(NUNDEF.LT.100) THEN
             WRITE(LINFO(K),'(5HUNDEF,I2)') NUNDEF
           ENDIF
           LENLI(K)=7
           LUNIT(K)='?'
           LENLU(K)=1
           LEVLR(K)='level'
           LTYPE(K)='-'
		    ENDIF
				IF(K.GE.140.AND.K.LE.159) THEN
           print *,' ERROR in utilizing index FILE ! ',
     1        'KPDS(6,N)=',KPDS(6,N),' FOR N=',N
           print *,' Mixed precision level not programed',
     1        'to work in this program.'
           PRINT *,' ERROR in utilizing index FILE ! ',
     1        'KPDS(6,N)=',KPDS(6,N),' FOR N=',N
           PRINT *,' Mixed precision level not programed',
     1        'to work in this program.'
  		     CALL ABORT
        ENDIF
C
        IF(LTYPE(K)(1:1).EQ.'?') THEN
          print *,' FATAL ERROR! KPDS(6,N)=',
     1            KPDS(6,N),' FOR N=',N
          print *,' LVTYP NOT SPECIFIED IN THE LEVEL subroutine'
          PRINT *,' FATAL ERROR! KPDS(6,N)=',
     1           KPDS(6,N),' FOR N=',N
          PRINT *,' LTYPE NOT SPECIFIED IN THE LEVEL subroutine'
          CALL ABORT
        ENDIF
C
        DO I=1,LENLF(K)
          IF(LFACT(K)(I:I).EQ.'.') THEN
            NDEC=LENLF(K)-I
          ENDIF
        ENDDO
        IF(LENLF(K).LT.10) THEN
          IF(NDEC.LT.10) THEN
            WRITE(FMT,'(2H(F,I1,1H.,I1,1H))') LENLF(K),NDEC
          ELSE
            WRITE(FMT,'(2H(F,I1,1H.,I2,1H))') LENLF(K),NDEC
          ENDIF
        ELSEIF(LENLF(K).LT.100) THEN
          IF(NDEC.LT.10) THEN
            WRITE(FMT,'(2H(F,I2,1H.,I1,1H))') LENLF(K),NDEC
          ELSE
            WRITE(FMT,'(2H(F,I2,1H.,I2,1H))') LENLF(K),NDEC
          ENDIF
        ELSE
          IF(NDEC.LT.10) THEN
            WRITE(FMT,'(2H(F,I3,1H.,I1,1H))') LENLF(K),NDEC
          ELSE
            WRITE(FMT,'(2H(F,I3,1H.,I2,1H))') LENLF(K),NDEC
          ENDIF
        ENDIF
C
        READ(LFACT(K),FMT) RFACT
        DO I=1,LENLB(K)
          IF(LBASE(K)(I:I).EQ.'.') THEN
            NDEC=LENLB(K)-I
          ENDIF
        ENDDO
        IF(LENLB(K).LT.10) THEN
          IF(NDEC.LT.10) THEN
            WRITE(FMT,'(2H(F,I1,1H.,I1,1H))') LENLB(K),NDEC
          ELSE
            WRITE(FMT,'(2H(F,I1,1H.,I2,1H))') LENLB(K),NDEC
          ENDIF
        ELSEIF(LENLB(K).LT.100) THEN
          IF(NDEC.LT.10) THEN
            WRITE(FMT,'(2H(F,I2,1H.,I1,1H))') LENLB(K),NDEC
          ELSE
            WRITE(FMT,'(2H(F,I2,1H.,I2,1H))') LENLB(K),NDEC
          ENDIF
        ELSE
          IF(NDEC.LT.10) THEN
            WRITE(FMT,'(2H(F,I3,1H.,I1,1H))') LENLB(K),NDEC
          ELSE
            WRITE(FMT,'(2H(F,I3,1H.,I2,1H))') LENLB(K),NDEC
          ENDIF
        ENDIF
C
        READ(LBASE(K),FMT) RBASE
C
        IF(LEVLR(K)(1:5).EQ.'layer') THEN
          NN=KPDS(7,N)/2**8
          TOPL=NN*RFACT+RBASE
          BTML=(KPDS(7,N)-NN*2**8)*RFACT+RBASE
        ELSEIF(LEVLR(K)(1:5).EQ.'level') THEN
          TOPL=KPDS(7,N)*RFACT+RBASE
          BTML=KPDS(7,N)*RFACT+RBASE
        ELSEIF(LEVLR(K)(1:6).EQ.'slevel') THEN
  		    IF(LTYPE(K)(1:7).EQ.'surface') THEN
            TOPL=0.
            BTML=0.
          ELSE
            TOPL=KPDS(7,N)*RFACT+RBASE
            BTML=KPDS(7,N)*RFACT+RBASE
          ENDIF
        ELSE
           print *,' FATAL ERROR! KPDS(6,N)=',
     1             KPDS(6,N),' FOR N=',N
           print *,' ILLEGAL LEVEL TYPE ENCOUNTERED ',
     1             'IN THE FILE KPDS(6) TABLE FILE'
           print *,' LEVLR(K)=',LEVLR(K)
           PRINT *,' FATAL ERROR! KPDS(6,N)=',
     1             KPDS(6,N),' FOR N=',N
           PRINT *,' ILLEGAL LEVEL TYPE ENCOUNTERED ',
     1             'IN THE FILE KPDS(6) TABLE FILE'
           PRINT *,' LEVLR(K)=',LEVLR(K)
  		     CALL ABORT
		    ENDIF
C
        RLEV(N)=BTML
C
        IF(FLOAT(INT(TOPL)).EQ.TOPL) THEN
          WRITE(LTOPL,'(I5)') INT(TOPL)
          LENT=5
        ELSE
          WRITE(LTOPL,'(F8.2)') TOPL
          LENT=8
        ENDIF
        IF(FLOAT(INT(BTML)).EQ.BTML) THEN
          WRITE(LBTML,'(I5)') INT(BTML)
          LENB=5
        ELSE
          WRITE(LBTML,'(F8.2)') BTML
          LENB=8
        ENDIF
C
        LABEL(2,N)=LINFO(K)(1:LENLI(K))
        LENLA(2,N)=LENLI(K)
C
  		  IF(LEVLR(K)(1:5).EQ.'level'.OR.LEVLR(K)(1:6).EQ.'slevel') THEN
  		    IF(LTYPE(K)(1:7).EQ.'surface'.OR.LTYPE(K)(1:1).EQ.'-') THEN
            LABEL(3,N)=' '
            LENLA(3,N)=1
	  	    ELSE
            LABEL(3,N)=LTOPL(1:LENT)//LUNIT(K)(1:LENLU(K))
            LENLA(3,N)=LENT+LENLU(K)
		  		ENDIF
		 	  ELSEIF(LEVLR(K)(1:5).EQ.'layer') THEN
          LABEL(3,N)=LBTML(1:LENB)//LUNIT(K)(1:LENLU(K))//' and '//
     1               LTOPL(1:LENT)//LUNIT(K)(1:LENLU(K))
          LENLA(3,N)=LENB+LENLU(K)+5+LENT+LENLU(K)
        ENDIF
C
        LVABBR(N)=LEVAB(K)
C
  100 CONTINUE
C
      RETURN
			END
      SUBROUTINE RMBLNK(LABEL,LENLA,N)
C
C  Get rid of two or more consecutive blanks here
C
      CHARACTER*64 LABEL(N)
      DIMENSION LENLA(N)
      CHARACTER*160 LABEX
C
      DO II=1,N
        JJ=0
        DO J=2,LENLA(II)
          IF(LABEL(II)(J-1:J-1).NE.' '.OR.LABEL(II)(J:J).NE.' ') THEN
	          JJ=JJ+1
            LABEX(JJ:JJ)=LABEL(II)(J-1:J-1)
          ENDIF
          IF(JJ.GT.1) THEN
            IF(LABEX(JJ-1:JJ-1).EQ.'='.AND.LABEX(JJ:JJ).EQ.' ') THEN
              JJ=JJ-1
            ENDIF
          ENDIF
        ENDDO
C
        IF(LABEL(II)(LENLA(II):LENLA(II)).NE.' ') THEN
  	    	JJ=JJ+1
          LABEX(JJ:JJ)=LABEL(II)(LENLA(II):LENLA(II))
        ENDIF
C
        IF(JJ.LE.0) THEN
          LENLA(II)=1
          LABEL(II)(1:1)=' '
        ELSE
          LENLA(II)=JJ
          LABEL(II)(1:LENLA(II))=LABEX(1:LENLA(II))
        ENDIF
C
      ENDDO
C
      RETURN
      END
      SUBROUTINE FLSRCH(NTOTL,LUPTR,
     1      MTOTL,LABEL,LENLA,LABBR,LVABBR,RLEV,
     2      ILINDX,IVINDX,IPINDX,
     3      NLTOTL,NVTOTL,NPTOTL,LLEVL,LLTYP)
C
      CHARACTER*64 LABEL
      CHARACTER*8  LABBR
      CHARACTER*8  LVABBR
      LOGICAL LLEVL,LLTYP
      CHARACTER*64 LABTMP
C
      DIMENSION LABEL(3,MTOTL),LENLA(3,MTOTL),LABBR(MTOTL),
     1          LVABBR(MTOTL),RLEV(MTOTL),
     2          ILINDX(MTOTL),IVINDX(MTOTL),IPINDX(MTOTL),
     3          NLTOTL(MTOTL)
C
      LOGICAL, ALLOCATABLE :: LFP(:)
      LOGICAL LPFND
C
      LOGICAL LDEBUG
      COMMON/COMDBG/ LDEBUG
C
C  Modify variable name if same variable name but with differenct level types
C
      ALLOCATE (LFP(MTOTL))
      DO N=1,NTOTL
        LFP(N)=.FALSE.
      ENDDO
      DO N=1,NTOTL
        IF(.NOT.LFP(N)) THEN
          DO M=N+1,NTOTL
            IF(.NOT.LFP(M)) THEN
              IF(LABBR(N).EQ.LABBR(M).AND.
     1           LABEL(2,N)(1:LENLA(2,N)).NE.
     2           LABEL(2,M)(1:LENLA(2,M))) THEN
                J=0
                DO I=1,5
                  IF(LABBR(M)(I:I).NE.' ') THEN
                    J=J+1
                    LABTMP(J:J)=LABBR(M)(I:I)
                  ENDIF
                ENDDO
                LABBR(M)(1:J+3)=LABTMP(1:J)//LVABBR(M)(1:3)
                DO I=J+4,8
                  LABBR(M)(I:I)=' '
                ENDDO
                print *,'Variable "',LABTMP(1:J),'" on ',
     1                        LABEL(2,M)(1:LENLA(2,M)),' modified to "',
     2                        LABBR(M)(1:J+3),'"'
                LFP(M)=.TRUE.
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C  Find number of variables 
C
      DO N=1,NTOTL
        LFP(N)=.FALSE.
      ENDDO
      NVTOTL=0
      NN=1
      NV=1
      DO N=1,NTOTL
        IF(.NOT.LFP(N)) THEN
          NVTOTL=NVTOTL+1
          IVINDX(NVTOTL)=N
          LFP(N)=.TRUE.
          DO M=N+1,NTOTL
            IF(.NOT.LFP(M).AND.
     1          LABBR(N).EQ.LABBR(M)) THEN
              LFP(M)=.TRUE.
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      IF(LDEBUG) THEN
        print *,' NVTOTL=',NVTOTL
        DO K=1,NVTOTL
          print *,LABBR(IVINDX(K))
        ENDDO
      ENDIF
C
C  Find number of levels for each variable
C
      LLEVL=.FALSE.
      DO N=1,NTOTL
        LFP(N)=.FALSE.
      ENDDO
      NNN=0
      DO NV=1,NVTOTL
        NVL=IVINDX(NV)
        NLTOTL(NV)=0
        DO N=1,NTOTL
          IF(.NOT.LFP(N).AND.
     1                       LABEL(1,N)(1:LENLA(1,N)).EQ.
     2                       LABEL(1,NVL)(1:LENLA(1,NVL)).AND.
     3                       LABBR(N).EQ.LABBR(NVL)) THEN
            NNN=NNN+1
            NLTOTL(NV)=NLTOTL(NV)+1
            ILINDX(NNN)=N
            LFP(N)=.TRUE.
            DO M=N+1,NTOTL
              IF(.NOT.LFP(M).AND.
     1                       LABEL(1,M)(1:LENLA(1,M)).EQ.
     1                       LABEL(1,NVL)(1:LENLA(1,NVL)).AND.
     2                       LABBR(M).EQ.LABBR(NVL).AND.
     3                       LABEL(3,M)(1:LENLA(3,M)).EQ.
     4                       LABEL(3,N)(1:LENLA(3,N))) THEN
                  LFP(M)=.TRUE.
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C
      IF(LDEBUG) THEN
        print *,' '
        print *,' NLTOTL(NV)=',(NLTOTL(NV),NV=1,NVTOTL)
        NNN=0
        DO NV=1,NVTOTL
          print *,(LABEL(3,ILINDX(K+NNN))
     1                   (1:LENLA(3,ILINDX(K+NNN))),',',K=1,NLTOTL(NV))
          NNN=NNN+NLTOTL(NV)
        ENDDO
      ENDIF
      LLEVL=.FALSE.
      DO NV=1,NVTOTL
        IF(NLTOTL(NV).GT.1) THEN
          LLEVL=.TRUE.
          GO TO 101
        ENDIF
      ENDDO
  101 CONTINUE
C
C  Find number of multi-level level types
C
      DO N=1,NTOTL
        LFP(N)=.FALSE.
      ENDDO
      NPTOTL=0
      NN=0
      DO N=1,NTOTL
        IF(.NOT.LFP(N).AND.LENLA(3,N).GT.1) THEN
          NPTOTL=NPTOTL+1
          LFP(N)=.TRUE.
          IPINDX(NPTOTL)=N
          DO M=N+1,NTOTL
            IF(.NOT.LFP(M).AND.
     1                         LABEL(2,N)(1:LENLA(2,N)).EQ.
     2                         LABEL(2,M)(1:LENLA(2,M))) THEN
              LFP(M)=.TRUE.
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(NPTOTL.GT.1) THEN
        LLTYP=.TRUE.
      ELSE
        LLTYP=.FALSE.
      ENDIF
C
      IF(LDEBUG) THEN
        print *,' '
        print *,' NPTOTL=',NPTOTL
        DO K=1,NPTOTL
          print *,LABEL(2,IPINDX(K))(1:LENLA(2,IPINDX(K)))
        ENDDO
      ENDIF
C
      print *,' '
      IF(LLEVL) THEN
        print *,' Multi-level fields found in a given file'
        print *,' Number of multi-level fields=',NLTOTL
      ENDIF
      IF(LLTYP) THEN
        print *,' Same variable with different level type found',
     1                 ' in a given file'
        print *,' Number of such fields=',NPTOTL
      ENDIF
C
      DEALLOCATE(LFP)
      RETURN
      END
      SUBROUTINE WRYTE(LU,LC,C)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    WRYTE       WRITE DATA OUT BY BYTES
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: EFFICIENTLY WRITE UNFORMATTED A CHARACETER ARRAY.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL WRYTE(LU,LC,C)
C
C   INPUT ARGUMENT LIST:
C     LU       - INTEGER UNIT TO WHICH TO WRITE
C     LC       - INTEGER NUMBER OF CHARACTERS OR BYTES TO WRITE
C     C        - CHARACETER (LC) DATA TO WRITE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      CHARACTER C(LC)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      WRITE(LU) C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE PUTGB(LUGB,KF,KPDS,KGDS,LB,F,IBM,FR,GRIB,LGRIB,
     1                 LUPTR,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK

C SUBPROGRAM: PUTGB          PACKS AND WRITES A GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 94-04-01
C
C ABSTRACT: PACK AND WRITE A GRIB MESSAGE.
C   THIS SUBPROGRAM IS NEARLY THE INVERSE OF GETGB.
C
C PROGRAM HISTORY LOG:
C   94-04-01  IREDELL
C
C USAGE:    CALL PUTGB(LUGB,KF,KPDS,KGDS,LB,F,IBM,FR,GRIB,IRET)
C   INPUT ARGUMENTS:
C     LUGB         INTEGER UNIT OF THE UNBLOCKED GRIB DATA FILE
C     KF           INTEGER NUMBER OF DATA POINTS
C     KPDS         INTEGER (25) PDS PARAMETERS
C     KGDS         INTEGER (22) GDS PARAMETERS
C     LB           LOGICAL (KF) BITMAP IF PRESENT
C     F            REAL (KF) DATA
C   OUTPUT ARGUMENTS:
C     IRET         INTEGER RETURN CODE
C                    0      ALL OK
C                    OTHER  W3FI72 GRIB PACKER RETURN CODE
C
C SUBPROGRAMS CALLED:
C   R63W72         MAP W3FI63 PARAMETERS ONTO W3FI72 PARAMETERS
C   GTBITS         GET NUMBER OF BITS AND ROUND DATA
C   W3FI72         PACK GRIB
C   WRYTE          WRITE DATA
C
C ATTRIBUTES:
C   LANGUAGE: F77 FORTRAN
C
C$$$
      INTEGER KPDS(25),KGDS(22)
      LOGICAL*1 LB(KF)
      REAL F(KF)
      INTEGER IBM(KF),IPDS(25),IGDS(18),IBDS(9)
      CHARACTER PDS(28)
C
C     CHARACTER GRIB(200+17*KF/8)
      CHARACTER*1 GRIB(*)
C
      REAL*8 FR(KF),FMAX,FMIN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET W3FI72 PARAMETERS
C     PRINT *,'KGDS(1)=',KGDS(1)
      CALL R63W72(KPDS,KGDS,IPDS,IGDS)
C     PRINT *,'IGDS(3)=',IGDS(3)
      DO I=1,9
      IBDS(I)=0
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COUNT VALID DATA
      KBM=KF
      IF(IPDS(7).NE.0) THEN
        KBM=0
        DO I=1,KF
          IF(LB(I)) THEN
            IBM(I)=1
            KBM=KBM+1
          ELSE
            IBM(I)=0
          ENDIF
        ENDDO
        IF(KBM.EQ.KF) IPDS(7)=0
      ENDIF
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET NUMBER OF BITS AND ROUND DATA
      IF(KBM.EQ.0) THEN
        DO I=1,KF
          FR(I)=0.
        ENDDO
        NBIT=0
      ELSE
        CALL GTBITS(IPDS(7),IPDS(25),KF,IBM,F,FR,FMIN,FMAX,NBIT)
        NBIT=MIN(NBIT,16)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PACK AND WRITE GRIB DATA
C
C  FR IS REAL*8
C
      CALL W3FI68(IPDS,PDS)
C
C     CALL W3FI72(0,FR,0,NBIT,0,IPDS,PDS,
      CALL W3FI72(0,FR,0,NBIT,1,IPDS,PDS,
     &            1,255,IGDS,0,0,IBM,KF,IBDS,
     &            KFO,GRIB,LGRIB,IRET)
      IF(IRET.EQ.0) CALL WRYTE(LUGB,LGRIB,GRIB)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE R63W72(KPDS,KGDS,IPDS,IGDS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    R63W72      CONVERT W3FI63 PARMS TO W3FI72 PARMS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: DETERMINES THE INTEGER PDS AND GDS PARAMETERS
C           FOR THE GRIB1 PACKING ROUTINE W3FI72 GIVEN THE PARAMETERS
C           RETURNED FROM THE GRIB1 UNPACKING ROUTINE W3FI63.
C
C PROGRAM HISTORY LOG:
C   91-10-31  MARK IREDELL
C
C USAGE:    CALL R63W72(KPDS,KGDS,IPDS,IGDS)
C
C   INPUT ARGUMENT LIST:
C     KPDS     - INTEGER (25) PDS PARAMETERS FROM W3FI63
C     KGDS     - INTEGER (22) GDS PARAMETERS FROM W3FI63
C
C   OUTPUT ARGUMENT LIST:
C     IPDS     - INTEGER (25) PDS PARAMETERS FOR W3FI72
C     IGDS     - INTEGER (18) GDS PARAMETERS FOR W3FI72
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION KPDS(25),KGDS(22),IPDS(25),IGDS(18)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  DETERMINE PRODUCT DEFINITION SECTION (PDS) PARAMETERS
      IPDS(1)=28			! LENGTH OF PDS
      IPDS(2)=KPDS(19)			! PARAMETER TABLE VERSION
      IPDS(3)=KPDS(1)			! ORIGINATING CENTER
      IPDS(4)=KPDS(2)			! GENERATING MODEL
      IPDS(5)=KPDS(3)			! GRID DEFINITION
      IPDS(6)=MOD(KPDS(4)/128,2)	! GDS FLAG
      IPDS(7)=MOD((KPDS(4)-128)/64,2)		! BMS FLAG
      IPDS(8)=KPDS(5)			! PARAMETER INDICATOR
      IPDS(9)=KPDS(6)			! LEVEL TYPE
      IF(IPDS(9).EQ.101.OR.IPDS(9).EQ.104.OR.IPDS(9).EQ.106.OR.
     &   IPDS(9).EQ.108.OR.IPDS(9).EQ.110.OR.IPDS(9).EQ.112.OR.
     &   IPDS(9).EQ.114.OR.IPDS(9).EQ.121.OR.IPDS(9).EQ.128.OR.
     &   IPDS(9).EQ.141)  THEN
        IPDS(10)=MOD(KPDS(7)/256,256)	! LEVEL VALUE 1
        IPDS(11)=MOD(KPDS(7),256)	! LEVEL VALUE 2
      ELSE
        IPDS(10)=0			! LEVEL VALUE 1
        IPDS(11)=KPDS(7)		! LEVEL VALUE 2
      ENDIF
      IPDS(12)=KPDS(8)		! YEAR OF CENTURY
      IPDS(13)=KPDS(9)			! MONTH
      IPDS(14)=KPDS(10)			! DAY
      IPDS(15)=KPDS(11)			! HOUR
      IPDS(16)=KPDS(12)			! MINUTE
      IPDS(17)=KPDS(13)			! FORECAST TIME UNIT
      IPDS(18)=KPDS(14)			! TIME RANGE 1
      IPDS(19)=KPDS(15)			! TIME RANGE 2
      IPDS(20)=KPDS(16)			! TIME RANGE INDICATOR
      IPDS(21)=KPDS(17)			! NUMBER IN AVERAGE
      IPDS(22)=KPDS(20)			! NUMBER MISSING IN AVERAGE
      IPDS(23)=KPDS(21)			! CENTURY
      IPDS(24)=KPDS(23)			! SUBCENTER
      IPDS(25)=KPDS(22)			! DECIMAL SCALING
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  DETERMINE GRID DEFINITION SECTION (GDS) PARAMETERS
      IGDS(1)=KGDS(19)			! NUMBER OF VERTICAL COORDINATES
      IGDS(2)=KGDS(20)			! VERTICAL COORDINATES
C     PRINT *,'IGDS(3)=',IGDS(3)
      IGDS(3)=KGDS(1)			! DATA REPRESENTATION
      IGDS(4)=KGDS(2)			! (UNIQUE TO REPRESENTATION)
      IGDS(5)=KGDS(3)			! (UNIQUE TO REPRESENTATION)
      IGDS(6)=KGDS(4)			! (UNIQUE TO REPRESENTATION)
      IGDS(7)=KGDS(5)			! (UNIQUE TO REPRESENTATION)
      IGDS(8)=KGDS(6)			! (UNIQUE TO REPRESENTATION)
      IGDS(9)=KGDS(7)			! (UNIQUE TO REPRESENTATION)
      IGDS(10)=KGDS(8)		! (UNIQUE TO REPRESENTATION)
C
C  POSSIBLE ERROR IN W3LIB 
C
      IF(KGDS(1).NE.1) THEN
        IGDS(11)=KGDS(9)			! (UNIQUE TO REPRESENTATION)
        IGDS(12)=KGDS(10)			! (UNIQUE TO REPRESENTATION)
        IGDS(13)=KGDS(11)			! (UNIQUE TO REPRESENTATION)
        IGDS(14)=KGDS(12)			! (UNIQUE TO REPRESENTATION)
      ELSE                   ! Mercator projection
        IGDS(11)=KGDS(13)     ! Y-increments
        IGDS(12)=KGDS(12)			! X-Increments
        IGDS(13)=KGDS(9)      ! Cutting latitude
        IGDS(14)=KGDS(11)     ! Grid Representation type
      ENDIF
C
      IGDS(15)=KGDS(13)			! (UNIQUE TO REPRESENTATION)
      IGDS(16)=KGDS(14)			! (UNIQUE TO REPRESENTATION)
      IGDS(17)=KGDS(15)			! (UNIQUE TO REPRESENTATION)
      IGDS(18)=KGDS(16)			! (UNIQUE TO REPRESENTATION)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
C-----------------------------------------------------------------------
CFPP$ NOCONCUR R
      SUBROUTINE GTBITS(IBM,IDS,LEN,MG,G,GROUND,GMIN,GMAX,NBIT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GTBITS      COMPUTE NUMBER OF BITS AND ROUND FIELD.
C   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
C
C ABSTRACT: THE NUMBER OF BITS REQUIRED TO PACK A GIVEN FIELD
C   AT A PARTICULAR DECIMAL SCALING IS COMPUTED USING THE FIELD RANGE.
C   THE FIELD IS ROUNDED OFF TO THE DECIMAL SCALING FOR PACKING.
C   THE MINIMUM AND MAXIMUM ROUNDED FIELD VALUES ARE ALSO RETURNED.
C   GRIB BITMAP MASKING FOR VALID DATA IS OPTIONALLY USED.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C
C USAGE:    CALL GTBITS(IBM,IDS,LEN,MG,G,GMIN,GMAX,NBIT)
C   INPUT ARGUMENT LIST:
C     IBM      - INTEGER BITMAP FLAG (=0 FOR NO BITMAP)
C     IDS      - INTEGER DECIMAL SCALING
C                (E.G. IDS=3 TO ROUND FIELD TO NEAREST MILLI-VALUE)
C     LEN      - INTEGER LENGTH OF THE FIELD AND BITMAP
C     MG       - INTEGER (LEN) BITMAP IF IBM=1 (0 TO SKIP, 1 TO KEEP)
C     G        - REAL (LEN) FIELD
C
C   OUTPUT ARGUMENT LIST:
C     GROUND   - REAL*8 (LEN) FIELD ROUNDED TO DECIMAL SCALING
C     GMIN     - REAL*8 MINIMUM VALID ROUNDED FIELD VALUE
C     GMAX     - REAL*8 MAXIMUM VALID ROUNDED FIELD VALUE
C     NBIT     - INTEGER NUMBER OF BITS TO PACK
C
C SUBPROGRAMS CALLED:
C   ISRCHNE  - FIND FIRST VALUE IN AN ARRAY NOT EQUAL TO TARGET VALUE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION MG(LEN),G(LEN)
      REAL*8 GROUND(LEN),GMAX,GMIN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ROUND FIELD AND DETERMINE EXTREMES WHERE BITMAP IS ON
c      DS=10.**IDS
c
      XMX=ABS(G(1))
      DO I=2,LEN
        XMX=MAX(XMX,ABS(G(I)))
      ENDDO
      RMAX=2.**31
      DS=10.**IDS
      RXMX=XMX*DS
      DO WHILE(RXMX.GT.RMAX)
        IDS=IDS-1
        DS=10.**IDS
        RXMX=XMX*DS
      ENDDO
c
      IF(IBM.EQ.0) THEN
        GROUND(1)=NINT(G(1)*DS)/DS
        GMAX=GROUND(1)
        GMIN=GROUND(1)
        DO I=2,LEN
          GROUND(I)=NINT(G(I)*DS)/DS
          GMAX=MAX(GMAX,GROUND(I))
          GMIN=MIN(GMIN,GROUND(I))
        ENDDO
      ELSE
        I1=ISRCHNE(LEN,MG,1,0)
        IF(I1.GT.0.AND.I1.LE.LEN) THEN
          GROUND(I1)=NINT(G(I1)*DS)/DS
          GMAX=GROUND(I1)
          GMIN=GROUND(I1)
          DO I=1,I1-1
            GROUND(I)=0.
          ENDDO
          DO I=I1+1,LEN
            IF(MG(I).NE.0) THEN
              GROUND(I)=NINT(G(I)*DS)/DS
              GMAX=MAX(GMAX,GROUND(I))
              GMIN=MIN(GMIN,GROUND(I))
            ELSE
              GROUND(I)=0.
            ENDIF
          ENDDO
        ELSE
          GMAX=0.
          GMIN=0.
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE NUMBER OF BITS
      NBIT=LOG((GMAX-GMIN)*DS+0.9)/LOG(2.)+1.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE INCDTE(IY,IM,ID,IH,JY,JM,JD,JH,INCHR)
C
C  Compute JY,JM,JD,JH of forecast
C
      JY=IY
      JM=IM
      JD=ID
      INCDY=INCHR/24
      JH=IH+MOD(INCHR,24)
      INCDY=INCDY+JH/24
      JH=MOD(JH,24)
      DO INCD=1,INCDY
        JD=JD+1
        IF(JM.EQ.4.OR.JM.EQ.6.OR.JM.EQ.9.OR.JM.EQ.11) THEN
          MONDY=30
        ELSEIF(JM.EQ.2) THEN
          IF(MOD(JY,4).EQ.0) THEN
            MODNY=29
          ELSE
            MONDY=28
          ENDIF
        ELSE
          MONDY=31
        ENDIF
        IF(JD.GT.MONDY) THEN
          JM=JM+1
          JD=1
          IF(JM.GT.12) THEN
            JY=JY+1
            JM=1
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
      SUBROUTINE MAXMIN(F,IDIM,JDIM,IMAX,JMAX,KMAX,LUPTR)
C
      DIMENSION F(IDIM,JDIM,KMAX)
C
      DO 10 K=1,KMAX
C
      FMAX=F(1,1,K)
      FMIN=F(1,1,K)
C
      DO 20 J=1,JMAX
      DO 20 I=1,IMAX
      IF(FMAX.LE.F(I,J,K)) THEN
      FMAX=F(I,J,K)
      IIMAX=I
      JJMAX=J
      ENDIF
      IF(FMIN.GE.F(I,J,K)) THEN
      FMIN=F(I,J,K)
      IIMIN=I
      JJMIN=J
      ENDIF
   20 CONTINUE
C
      WRITE(6,100) K,FMAX,IIMAX,JJMAX,FMIN,IIMIN,JJMIN
  100 FORMAT(2X,'LEVEL=',I2,' MAX=',E10.4,' AT I=',I5,' J=',I5,
     1                      ' MIN=',E10.4,' AT I=',I5,' J=',I5)
C
   10 CONTINUE
C
      RETURN
      END
      FUNCTION ISRCHNE(N,X,INCX,TARGET)
      INTEGER X(*), TARGET
      J=1
      ISRCHNE=0
      IF(N.LE.0) RETURN
      IF(INCX.LT.0) J=1-(N-1)*INCX
      DO I=1,N
        IF(X(J).NE.TARGET) THEN
          ISRCHNE=I
          RETURN
        ENDIF
        J=J+INCX
      ENDDO
      RETURN
      END
      subroutine numchar(char,num)
c
      implicit none
      character*128 char
      integer num,n
c
      num=1
      dowhile(char(num:num).ne.' ')
        num=num+1
      enddo
      num=num-1
      do n=num+1,128
        char(n:n)=' '
      enddo
      return
      end
