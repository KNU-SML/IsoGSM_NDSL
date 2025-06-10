      subroutine cintpx(iin,jtwidl,jin,iout,
     1                 inslat,wgtlat,
     2                 cv,cvt,cvb,camt,ctop,cbot,
     3                 xx,wgt,tt,bb,sum,nn,ltwidl,latrd1)
      dimension cv(iin,jtwidl),cvt(iin,jtwidl),cvb(iin,jtwidl)
      dimension camt(iout),ctop(iout),cbot(iout)
      dimension xx(iout,4),wgt(iout,4),tt(iout,4),bb(iout,4),sum(iout,4)
      dimension nn(iout)
c        simpl linear interpolation of cldamt, unless only 1,2 of the
c         surrounding pts has cv. then,if output gridpt not close enuf
c         do not interpolate to it(prevents spreading of cv clds)..
c           for 1 pt convection-intrp wgt ge (.7)**2 ...
c           for 2 pt convection-sum of intrp wgt ge .45...
c              .45 used rather than .5 to give better result for
c              diagonally opposed pts...
c===>    for tops(cvt) and bots(cvb) just take average of surrounding
c         non-zero cv points.....
c         nn will be number of surrounding pts with cld (gt zero)
c---     nhsh = 1,-1 for northern,southern hemisphere
c         here instead of an extrapolation,just do a simple mean....
c
      lonf=iin/2
      if (inslat.lt.0) go to 600
      inth = mod(ltwidl + inslat + jtwidl - latrd1 - 1,jtwidl) + 1
      inth1 = mod(inth,jtwidl) + 1
      if (inslat.eq.jin) go to 105
      do 100 i=1,iout
        ileft=i
        irght=i+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        wgtlon=0.0
        wgtlat=1.0
c----   normalized distance from upper lat to gaussian lat
        xx(i,1) = cv(ileft,inth)
        xx(i,2) = cv(ileft,inth1)
        xx(i,3) = cv(irght,inth)
        xx(i,4) = cv(irght,inth1)
        wgt(i,1) = (1.e 0-wgtlon)*(1.e 0-wgtlat)
        wgt(i,2) = (1.e 0-wgtlon)*wgtlat
        wgt(i,3) = wgtlon*(1.e 0-wgtlat)
        wgt(i,4) = wgtlon*wgtlat
        tt(i,1) = cvt(ileft,inth)
        tt(i,2) = cvt(ileft,inth1)
        tt(i,3) = cvt(irght,inth)
        tt(i,4) = cvt(irght,inth1)
        bb(i,1) = cvb(ileft,inth)
        bb(i,2) = cvb(ileft,inth1)
        bb(i,3) = cvb(irght,inth)
        bb(i,4) = cvb(irght,inth1)
  100 continue
      go to 130
  105 do 110 i=1,iout
        ileft=i
        irght=i+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        wgtlon=0.0
        wgtlat=1.0
c----   normalized distance from upper lat to gaussian lat
        xx(i,1) = cv(ileft,inth)
        xx(i,3) = cv(irght,inth)
        wgt(i,1) = (1.e 0-wgtlon)*(1.e 0-wgtlat)
        wgt(i,2) = (1.e 0-wgtlon)*wgtlat
        wgt(i,3) = wgtlon*(1.e 0-wgtlat)
        wgt(i,4) = wgtlon*wgtlat
        tt(i,1) = cvt(ileft,inth)
        tt(i,3) = cvt(irght,inth)
        bb(i,1) = cvb(ileft,inth)
        bb(i,3) = cvb(irght,inth)
  110 continue
      iout2 = iout / 2
      do 120 i=1,iout2
        ileft=i
        irght=ileft+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        ileft2=i+iout2
        irght2=ileft2+1
        if(mod(ileft2,lonf).eq.0) irght2=irght2-lonf
        xx(i      ,2) = cv(ileft2,inth)
        xx(i+iout2,2) = cv(ileft ,inth)
        xx(i      ,4) = cv(irght2,inth)
        xx(i+iout2,4) = cv(irght ,inth)
        bb(i      ,2) = cvb(ileft2,inth)
        bb(i+iout2,2) = cvb(ileft ,inth)
        bb(i      ,4) = cvb(irght2,inth)
        bb(i+iout2,4) = cvb(irght ,inth)
        tt(i      ,2) = cvt(ileft2,inth)
        tt(i+iout2,2) = cvt(ileft ,inth)
        tt(i      ,4) = cvt(irght2,inth)
        tt(i+iout2,4) = cvt(irght ,inth)
  120 continue
c---      nn will be number of surrounding pts with cld (gt zero)
  130 do 10 i=1,iout
        nn(i) = 0
   10 continue
      do 12 j=1,4
       do 12 i=1,iout
        sum(i,j) = 0.e0
   12 continue
        do 150 kpt=1,4
          do 14 i=1,iout
            if (xx(i,kpt).gt.0.e0) then
              nn(i) = nn(i) + 1
              sum(i,1) = sum(i,1) + wgt(i,kpt)
              sum(i,2) = sum(i,2) + tt(i,kpt)
              sum(i,3) = sum(i,3) + bb(i,kpt)
            endif
   14     continue
          do 15 i=1,iout
            sum(i,4) = sum(i,4) + wgt(i,kpt) * xx(i,kpt)
   15     continue
  150   continue
        do 16 i=1,iout
          if (nn(i).eq.1.and.sum(i,1).gt.0.49e0) go to 17
          if (nn(i).eq.2.and.sum(i,1).ge.0.45e0) go to 17
          if (nn(i).ge.3) go to 17
            ctop(i) = 0.e0
            cbot(i) = 100.e0
            camt(i) = 0.e0
            go to 18
   17     continue
            ltop = sum(i,2)/nn(i) + 0.5e0
            lbot = sum(i,3)/nn(i) + 0.5e0
            lamt = sum(i,4) + 0.5e0
            ctop(i) = ltop
            cbot(i) = lbot
            camt(i) = lamt
   18     continue
   16   continue
      return
c--- polar region-no extrapolation
  600 continue
      ja = iabs(inslat)
      do 200 i=1,iout
        ileft=i
        irght=ileft+1
        if(mod(ileft,lonf).eq.0) irght=irght-lonf
        wgtlon=0.0
c----    get left point on nearest latitude
        xx(i,1) = cv(ileft,ja)
        xx(i,2) = cv(irght,ja)
        wgt(i,1) = 1.e 0-wgtlon
        wgt(i,2) = wgtlon
        tt(i,1) = cvt(ileft,ja)
        tt(i,2) = cvt(irght,ja)
        bb(i,1) = cvb(ileft,ja)
        bb(i,2) = cvb(irght,ja)
  200 continue
c---      nn will be number of surrounding pts with cld (gt zero)
      do 20 i=1,iout
        nn(i) = 0
   20 continue
      do 22 j=1,4
       do 22 i=1,iout
        sum(i,j) = 0.e0
   22 continue
        do 202 kpt=1,2
          do 24 i=1,iout
            if (xx(i,kpt).gt.0.e0) then
              nn(i) = nn(i) + 1
              sum(i,1) = sum(i,1) + wgt(i,kpt)
              sum(i,2) = sum(i,2) + tt(i,kpt)
              sum(i,3) = sum(i,3) + bb(i,kpt)
            endif
   24     continue
          do 25 i=1,iout
            sum(i,4) = sum(i,4) + wgt(i,kpt) * xx(i,kpt)
   25     continue
  202   continue
        do 26 i=1,iout
          if (nn(i).eq.1.and.sum(i,1).gt.0.7e0) go to 27
          if (nn(i).eq.2) go to 27
            ctop(i) = 0.e0
            cbot(i) = 100.e0
            camt(i) = 0.e0
            go to 28
   27     continue
            ltop = sum(i,2)/nn(i) + 0.5e0
            lbot = sum(i,3)/nn(i) + 0.5e0
            lamt = sum(i,4) + 0.5e0
            ctop(i) = ltop
            cbot(i) = lbot
            camt(i) = lamt
   28     continue
   26   continue
      return
      end
