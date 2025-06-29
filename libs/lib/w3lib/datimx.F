#include <define.h>
      SUBROUTINE DATIMX(NOW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    DATIMX      PROVIDES DATE/TIME INFORMATION
C   PRGMMR: R.E.JONES        ORG: W/NMC42    DATE: 93-11-10
C
C ABSTRACT: DATIMX PROVIDES YOU WITH DATE AND TIME INFORMATION 
C   THAT CAN BE USED BY YOUR PROGRAMS. IT HAS THE SAME OUTPUT
C   AS IBM SUBROUTINE DATIMX IN IBM VS FORTRAN 2.5
C
C PROGRAM HISTORY LOG:
C   93-11-10  R. E. JONES   SAME AS IBM VS FORTRAN VER 2.5
C
C USAGE:    CALL DATIMX(NOW)
C   OUTPUT ARGUMENT LIST:  
C     NOW       - INTEGER*4 ARRAY WITH 14 WORDS
C     
C     ELEMENT   CONTENTS
C
C     1         THE VALUE OF THE CLOCK EXPRESSED AS A POSTITIVE 
C               INTEGER, SET TO 1, SET TO -1 IF INVALID
C     2         MILLISECONDS (0-990) PRECISE TO THE HUNDREDTHS              
C               POSITION, SET TO ZERO ON SGI
C     3         SECONDS (0-59)
C     4         MINUTES (0-59)
C     5         HOUR USING A 24-HOUR CLOCK (0-23)
C     6         DAY OF MONTH  (1-31)
C     7         MONTH OF YEAR (1-12)
C     8         YEAR (4 DIGIT, FOR EXAMPLE: 1993)
C     9         RESERVED. VALUE RETURNED AS -1.
C               WE ARE GOING TO PUT THE JULIAN DAY NUMBER HERE
C     10        HOUR USING 12-HOUR CLOCK (1-12)
C     11        AM/PM INDICATOR. 1 IS RETURNED FOR AM, 2 FOR PM.
C     12        DAY OF WEEK (1-7) BEGINNING WITH SUNDAY
C     13        DAY OF YEAR (1-366)
C     14        YEAR (2 DIGIT, FOR EXAMPLE: 93)
C
C   SUBPROGRAMS CALLED:
C     LIBRARY:
C      SYSTEM  - IDATE, TIME
C       W3LIB  - IW3JDN
C
C WARNING: THIS SUBROUTINE IS NOT PORTABLE, DATE SUBROUTINES ARE
C          NOT STANDARD FORTRAN 77.
C
C ATTRIBUTES:
C   LANGUAGE: SiliconGraphics 3.5 FORTRAN 77
C   MACHINE:  SiliconGraphics IRIS-4D/25, 35, INDIGO, Indy 
C
C$$$
C
      INTEGER       NOW(14)
C
      CHARACTER * 8 ATIME
C 
C                 YR    MON     DAY
C
      CALL IDATE(NOW(7),NOW(6),NOW(14))
      CALL TIME(ATIME)
C
      NOW(1)  =  1
C 
C     ATIME IS CHARACTER DATA IN HH:MM:SS FORMAT 
C    
      READ(ATIME,10) NOW(5),NOW(4),NOW(3)
10    FORMAT(I2,1X,I2,1X,I2)
C
      IF (NOW(14).LT.55) THEN
        NOW(8) = NOW(14) + 2000
      ELSE
        NOW(8) = NOW(14) + 1900
      END IF
C
      IF (NOW(5).GT.11) THEN
        NOW(11) = 2
        NOW(10) = NOW(5) - 12
      ELSE
        NOW(11) = 1
        NOW(10) = NOW(5)
      END IF
C
C     COMPUTE JULIAN DAY NUMBER  (JAN 1,1960 IS 2,436,935)
C
      JDN    = IW3JDN(NOW(8),NOW(7),NOW(6))
      NOW(9) = JDN
C
C     COMPUTE DAY OF WEEK FROM JULIAN DAY NUMBER
C
      NOW(12) = MOD((JDN + 1),7) + 1
      IYR     = NOW(8)
C
C     COMPUTE DAY OF YEAR (1-366) USING JULIAN DAY NUMBER
C
      NOW(13) = JDN -
     &  (-31739 +1461 * (IYR+4799) / 4 - 3 * ((IYR+4899)/100)/4)
C
C     NO CLOCK FUNCTION TO GET 1000 TH OF SEC, SET TO ZERO ON SGI
C
      SEC1   = SECNDS(0.0)
      SEC2   = AINT(SEC1)
      NOW(2) = (SEC1 - SEC2) * 1000.0
C
      RETURN
      END
