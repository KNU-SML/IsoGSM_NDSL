      subroutine nfdr_setloads(ifm,f,w1d,w1dp,w10d,w100d,w1000d,wwood,
     *                    wherb,drfuel,sig1d,sig10d,sg100d,s1000d,
     *                    sgwood,sgherb,hd,
     *                    scm,extmoi,depth,wndftr,wfmin,wfmax,kdi)
c
      implicit none
c
      real f(20,20)
c
      real w1d,w1dp,w10d,w100d,w1000d,wwood,
     *     wherb,drfuel,sig1d,sig10d,sg100d,s1000d,sgwood,sgherb,hd,
     *     scm,extmoi,depth,wndftr,wfmin,wfmax
      real sumded,ra1,ra10,ra100,ra1000,addslp,wadd
      integer ifm,kdi
c
      real tons_to_lbs
      data tons_to_lbs/0.0459137/
c
c     ifm reduced by 1 when greater than 13 to avoid going outside
c     bounds of fuel array f(20,20) which has data for 20 parameters 
c     in 20 fuel models.
c
           if(ifm.gt.13) ifm=ifm-1
           w1d = f(ifm,1)*tons_to_lbs
           w1dp = w1d
           w10d = f(ifm,2)*tons_to_lbs
           w100d = f(ifm,3)*tons_to_lbs
           w1000d = f(ifm,4)*tons_to_lbs
           wwood = f(ifm,5)*tons_to_lbs
           wherb = f(ifm,6)*tons_to_lbs
           drfuel = f(ifm,7)*tons_to_lbs
           sig1d = f(ifm,8)
           sig10d = f(ifm,9)
           sg100d = f(ifm,10)
           s1000d = f(ifm,11)
           sgwood = f(ifm,12)
           sgherb = f(ifm,13)
           hd = f(ifm,14)
           scm = f(ifm,15)
           extmoi = f(ifm,16)
           depth = f(ifm,17)
           wndftr = f(ifm,18)
           wfmin = f(ifm,19)
           wfmax = f(ifm,20)
           sumded = w1dp+w10d+w100d+w1000d
           ra1 = w1dp/sumded
           ra10 = w10d/sumded
           ra100 = w100d/sumded
           ra1000 = w1000d/sumded
           addslp = drfuel/700.0
           wadd = 0.
           if(kdi.gt.100) then
             wadd = (real(kdi)-100.0)*addslp
             w1dp = w1d+ra1*wadd
             w10d = w10d+ra10*wadd
             w100d = w100d+ra100*wadd
             w1000d = w1000d+ra1000*wadd
           endif
      return
      end
