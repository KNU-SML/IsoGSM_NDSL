      SUBROUTINE BFRINI
#include "ufbprm.h"
      integer TAB_,SEQ_,TBA_,TBB_,TBD_
      parameter(INV_=_INV_,TAB_=_TAB_,SEQ_=_SEQ_)
      parameter(TBA_=_TBA_,TBB_=_TBB_,TBD_=_TBD_)

      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(5000),MBYT(10),MBAY(5000,10)
      COMMON /PADESC/ IBCT,IPD1,IPD2,IPD3,IPD4
      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)
      COMMON /SCALEU/ SCLU(-1000:1000)                                        
      COMMON /STBFR / IOLUN(10),IOMSG(10)
      COMMON /TABABD/ NTBA(0:10),NTBB(0:10),NTBD(0:10),MTAB(TBA_,10),
     .                IDNA(TBA_,10,2),IDNB(TBB_,10),IDND(TBD_,10),
     .                TABA(TBA_,10),TABB(TBB_,10),TABD(TBD_,10)
      COMMON /DXTAB / MAXDX,IDXV,NXSTR(10),LDXA(10),LDXB(10),LDXD(10),
     .                LD30(10),DXSTR(10)
      COMMON /TABLES/ MAXTAB,NTAB,TAG(TAB_),TYP(TAB_),KNT(TAB_),
     .                JUMP(TAB_),LINK(TAB_),JMPB(TAB_),
     .                IBT(TAB_),IRF(TAB_),ISC(TAB_),
     .                ITP(TAB_),VALI(TAB_),KNTI(TAB_),
     .                ISEQ(TAB_,2),JSEQ(SEQ_)
      COMMON /BUFRMG/ MSGLEN,MSGTXT(5000)


      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*56  DXSTR
      CHARACTER*10  TAG
      CHARACTER*6   ADSN(5,2),DNDX(25,10)
      CHARACTER*3   TYPX(5,2),TYPS,TYP
      CHARACTER*1   REPX(5,2),REPS
      DIMENSION     NDNDX(10),NLDXA(10),NLDXB(10),NLDXD(10),NLD30(10)
      DIMENSION     LENX(5)

      DATA ADSN   / '101000','360001','360002','360003','360004' ,
     .              '101255','031002','031001','031001','031000' /
      DATA TYPX   /    'REP',   'DRP',   'DRP',   'DRS' ,  'DRB' ,
     .                 'SEQ',   'RPC',   'RPC',   'RPS' ,  'SEQ' /
      DATA REPX   /      '"',     '(',     '{',     '[' ,    '<' ,
     .                   '"',     ')',     '}',     ']' ,    '>' /
      DATA LENX   /       0 ,     16 ,      8 ,      8  ,     1  /

      DATA (DNDX(I,1),I=1,25)/
     .'102000','031001','000001','000002',
     .'110000','031001','000010','000011','000012','000013','000015',
     .                  '000016','000017','000018','000019','000020',
     .'107000','031001','000010','000011','000012','000013','101000',
     .                  '031001','000030'/

      DATA (DNDX(I,2),I=1,15)/
     .'103000','031001','000001','000002','000003',
     .'101000','031001','300004',
     .'105000','031001','300003','205064','101000','031001','000030'/

      DATA NDNDX /  25 ,  15 , 8*0 /
      DATA NLDXA /  35 ,  67 , 8*0 /
      DATA NLDXB /  80 , 112 , 8*0 /
      DATA NLDXD /  38 ,  70 , 8*0 /
      DATA NLD30 /   5 ,   6 , 8*0 /

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  INITIALIZE /BITBUF/
C  -------------------

C     MAXBYT = 9970
      MAXBYT=19950

C  INITIALIZE /PADESC/
C  -------------------

      IBCT = IFXY('063000')
      IPD1 = IFXY('102000')
      IPD2 = IFXY('031001')
      IPD3 = IFXY('206001')
      IPD4 = IFXY('063255')

C  INITIALIZE /STBFR/
C  ------------------

      DO I=1,10
      IOLUN(I) = 0
      IOMSG(I) = 0
      ENDDO

C  INITIALIZE /REPTAB/
C  -------------------

      DO I=1,5
      LENS(I) = LENX(I)
      DO J=1,2
      IDNR(I,J) = IFXY(ADSN(I,J))
      TYPS(I,J) = TYPX(I,J)
      REPS(I,J) = REPX(I,J)
      ENDDO
      ENDDO

C  INITIALIZE /TABABD/
C  -------------------

      NTBA(0) = TBA_
      NTBB(0) = TBB_
      NTBD(0) = TBD_

C  INITIALIZE /DXTAB/
C  ------------------

      MAXDX = MAXBYT
      IDXV  = 1

      DO J=1,10
      LDXA(J)  = NLDXA(J)
      LDXB(J)  = NLDXB(J)
      LDXD(J)  = NLDXD(J)
      LD30(J)  = NLD30(J)
      DXSTR(J) = '      '
      NXSTR(J) = NDNDX(J)*2
      DO I=1,NDNDX(J)
      I1 = I*2-1
      CALL IPKM(DXSTR(J)(I1:I1),2,IFXY(DNDX(I,J)))
      ENDDO
      ENDDO

C  INITIALIZE /TABLES/
C  -------------------

      MAXTAB = TAB_

C  INITIALIZE /BUFRMG/
C  -------------------

      MSGLEN = 0
                                                                        
C  INITIALIZE /SCALEU/                                                  
C  -------------------                                                  
                                                                        
      SCLU(-1) = .1
      SCLU( 0) =  1
      SCLU( 1) = 10
      DO I=2,1000
      SCLU(-I) = SCLU(1-I)*SCLU(-1)
      SCLU( I) = SCLU(I-1)*SCLU( 1)
      ENDDO                            

      RETURN
      END
