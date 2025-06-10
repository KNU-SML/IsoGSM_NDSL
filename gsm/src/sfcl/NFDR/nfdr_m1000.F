      subroutine nfdr_m1000 (bndry,avethbndy,ythmois,fm1000)
c
      implicit none
c
      real bndry,avethbndy,ythmois,fm1000
c
c  This subroutine calculates the 1000 hrfm using average bndry
c  conditions for the last 7 days and yesterdays 1000 hrfm value.
c  It is assumed that files containing thousand hour and boundary
c  values already exist.  In this routine ythbndy is yesterdays
c  1000 hrfm, and avethbndry is the 7 day average 1000 hr bndry value.
c
c   fr1 equals (1.0 - 0.82*exp(-7.0 * 24.0 / 1000) = 0.3068
c
      real fr1
      parameter (fr1=0.3068)
c
      avethbndy = (avethbndy*6.0 + bndry)/7.0
      fm1000  = ythmois+(avethbndy - ythmois)*fr1
      ythmois = fm1000
c      
      return
      end
