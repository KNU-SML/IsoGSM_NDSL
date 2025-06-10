      subroutine nfdr_julian (month,iday,jdate)
c
      implicit none
c
      integer month,iday,jdate
c
      integer m(12),i
      data (m(i),i=1,12)/0,31,59,90,120,151,181,212,243,273,304,334/
c
      jdate=m(month)+iday
c
      return
      end
