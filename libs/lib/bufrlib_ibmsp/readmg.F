      SUBROUTINE READMG(LUNIT,SUBSET,JDATE,IRET)

      COMMON /HRDWRD/ NBYTW,NBITW,NREV,IORD(8)
      COMMON /MSGCWD/ NMSG(10),NSUB(10),MSUB(10),INODE(10),IDATE(10)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)

      CHARACTER*8 SEC0,SUBSET
      CHARACTER*1 TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      IRET = 0

C  CHECK THE FILE STATUS
C  ---------------------

      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      CALL WTSTAT(LUNIT,LUN,IL, 1)

C  READ A MESSAGE INTO A MESSAGE BUFFER
C  ------------------------------------

1     MBIT = 0
      IMSG = 8/NBYTW+1
      READ(LUNIT,ERR=902,END=100) SEC0,(MBAY(I,LUN),I=IMSG,LMSG(SEC0))
C     CALL PKC(SEC0,8,MBAY(1,LUN),MBIT)
      call chrtrn(mbay(1,lun),sec0,8)

C  SECTION 1
C  ---------

      IAD1 = 8
      LEN1 = IUPB(MBAY(1,LUN),IAD1+ 1,24)
      LEN2 = IUPB(MBAY(1,LUN),IAD1+ 8, 1)
      MTYP = IUPB(MBAY(1,LUN),IAD1+ 9, 8)
      MEAR = IUPB(MBAY(1,LUN),IAD1+13, 8)
      MMON = IUPB(MBAY(1,LUN),IAD1+14, 8)
      MDAY = IUPB(MBAY(1,LUN),IAD1+15, 8)
      MOUR = IUPB(MBAY(1,LUN),IAD1+16, 8)

      JDATE = MEAR*10**6 + MMON*10**4 + MDAY*10**2 + MOUR

      IF(MTYP.EQ.11) GOTO 1

C  SECTION 2
C  ---------

      IAD2 = IAD1+LEN1
      LEN2 = IUPB(MBAY(1,LUN),IAD2+1,24) * LEN2

C  SECTION 3
C  ---------

      IAD3 = IAD2+LEN2
      LEN3 = IUPB(MBAY(1,LUN),IAD3+1 ,24)
      JSUB = IUPB(MBAY(1,LUN),IAD3+5 ,16)
      ISUB = IUPB(MBAY(1,LUN),IAD3+10,16)

      CALL NUMTAB(LUN,ISUB,SUBSET,TAB,ITAB)
      CALL NEMTBA(LUN,SUBSET,MSGT,MSTB,INOD)
      IF(MSGT.NE.MTYP) GOTO 903

C  SECTION 4
C  ---------

      IAD4 = IAD3+LEN3
      LEN4 = IUPB(MBAY(1,LUN),IAD4+1,24)

      MBYT(LUN) = IAD4+4

C  NORMAL EXIT
C  -----------

      IF(INOD.NE.INODE(LUN)) CALL STRCLN
      IDATE(LUN) = JDATE
      INODE(LUN) = INOD
      MSUB (LUN) = JSUB
      NSUB (LUN) = 0
      NMSG (LUN) = NMSG(LUN)+1
      IRET = 0
      RETURN

C  EOF ON ATTEMPTED READ
C  ---------------------

100   CALL CLOSBF(LUNIT)
      INODE(LUN) = 0
      IDATE(LUN) = 0
      SUBSET = ' '
      JDATE = 0
      IRET = -1
      RETURN

C  ERROR EXITS
C  -----------

900   CALL ABT('READMG - FILE IS CLOSED                   ')
901   CALL ABT('READMG - FILE IS OPEN FOR OUTPUT          ')
902   CALL ABT('READMG - I/O ERROR READING MESSAGE        ')
903   CALL ABT('READMG - MSGTYPE MISMATCH FOR '//SUBSET    )
      END
