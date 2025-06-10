      subroutine angrad( pres, alm, stla, ctl2, al2
     *,                  pi, tla, prb, wfn, ufn)
      real refp(6), refr(6), tlac(8), plac(8)
      real pres, stla, ctl2
      real alm, al2, tla, tem, tem1
     *,    prb, acr, wfn, ufn
!
      data plac/100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0/

      data tlac/ 25.0,  22.5,  20.0,  17.5,  15.0,   12.5,   10.0, 10.0/

!
      data refp/500.0, 300.0, 250.0, 200.0, 150.0, 100.0/

      data refr/ 0.25,   0.5,   0.75,  1.0,  1.5,  2.0/

!
cfpp$ expand (crtwrk)

      tem = pi / 180.0

      if (tla .lt. 0.0) then

!
          if (pres .le. plac(1)) then
            tla = tlac(1)
          elseif (pres .le. plac(2)) then
            tla = tlac(2) + (pres-plac(2))*((tlac(1)-tlac(2))/
     *                                           (plac(1)-plac(2)))
          elseif (pres .le. plac(3)) then
            tla = tlac(3) + (pres-plac(3))*((tlac(2)-tlac(3))/
     *                                           (plac(2)-plac(3)))
          elseif (pres .le. plac(4)) then
            tla = tlac(4) + (pres-plac(4))*((tlac(3)-tlac(4))/
     *                                           (plac(3)-plac(4)))
          elseif (pres .le. plac(5)) then
            tla = tlac(5) + (pres-plac(5))*((tlac(4)-tlac(5))/
     *                                           (plac(4)-plac(5)))
          elseif (pres .le. plac(6)) then
            tla = tlac(6) + (pres-plac(6))*((tlac(5)-tlac(6))/
     *                                           (plac(5)-plac(6)))
          elseif (pres .le. plac(7)) then
            tla = tlac(7) + (pres-plac(7))*((tlac(6)-tlac(7))/
     *                                           (plac(6)-plac(7)))
          elseif (pres .le. plac(8)) then
            tla = tlac(8) + (pres-plac(8))*((tlac(7)-tlac(8))/
     *                                           (plac(7)-plac(8)))
          else
            tla = tlac(8)
          endif

          stla = sin(tla*tem)
          tem1 = cos(tla*tem)
          ctl2 = tem1 * tem1

      else
          stla = sin(tla*tem)
          tem1 = cos(tla*tem)
          ctl2 = tem1 * tem1

      endif
        if (pres .ge. refp(1)) then
          tem = refr(1)
        elseif (pres .ge. refp(2)) then
          tem = refr(1) + (pres-refp(1)) *
     *                   ((refr(2)-refr(1))/(refp(2)-refp(1)))
        elseif (pres .ge. refp(3)) then
          tem = refr(2) + (pres-refp(2)) *
     *                   ((refr(3)-refr(2))/(refp(3)-refp(2)))
        elseif (pres .ge. refp(4)) then
          tem = refr(3) + (pres-refp(3)) *
     *                   ((refr(4)-refr(3))/(refp(4)-refp(3)))
        elseif (pres .ge. refp(5)) then
          tem = refr(4) + (pres-refp(4)) *
     *                   ((refr(5)-refr(4))/(refp(5)-refp(4)))
        elseif (pres .ge. refp(6)) then
          tem = refr(5) + (pres-refp(5)) *
     *                   ((refr(6)-refr(5))/(refp(6)-refp(5)))
        else
          tem = refr(6)
        endif
        al2 = max(alm, 2.0e-4/tem)

!
      return
      end
