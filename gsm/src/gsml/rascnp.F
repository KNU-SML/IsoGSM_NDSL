      subroutine rascnp(len, lenc, km, lm, nstrp, dt, si, sl
     *,           krmin, krmax, ncrnd, afac, rannum, ufac
     *,           rgas, cp, grav, alhl
     *,           tin, qin, uin, vin, pin
     *,           evapc, kbot, ktop, kuo, lat, cd
     *,           sig, prj, prh, fpk, hpk, sgb, ods, rasal, prns
     *,           cwmin, cmf, detr )
c
      parameter (frac=5.0, crtmsf=0.0)
c
c      input
c
      dimension tin(lenc,km), qin(lenc,km), uin(lenc,km)
     *,         vin(lenc,km), pin(len),     cd(len)
     *,         evapc(len), kbot(len), ktop(len), kuo(len)
     *,         rannum(ncrnd*3)
     *,         si(km+1), sl(km)
      real      cwmin(lenc,km)
      real      cmf(lenc,km), detr(lenc,km)
c
      dimension sig(lm+1),  prj(lm+1), prh(lm),  fpk(lm), hpk(lm)
     *,         sgb(lm),  ods(lm), rasal(lm), prns(nstrp)
c
c     locals
c
      real pib(len),    pik(len),    pki(len),    pcu(len)
      real tns(len),    qns(len),    uos(len),    vos(len)
      real gam(len,lm)
      real uoc(len,lm), voc(len,lm), tn0(len,lm), qn0(len,lm)
      real tcu(len,lm), qcu(len,lm), clt(len,lm)
      real cwm0(len,lm), cmf0(len,lm),detr0(len,lm)

      kmp1 = km + 1
      ir1  = 5 * lm + 1
      lmm1 = lm - 1
      lmp1 = lm + 1
      rkap = rgas / cp
      ia = 1
c
      do i=1,len
         pib(i) = pin(i) * 10.0
         pik(i) = pib(i) ** rkap
         pki(i) = 1.0 / pik(i)
         uos(i) = 0.0
         vos(i) = 0.0
         tns(i) = 0.0
         qns(i) = 0.0
      enddo
      do l=1,lmm1
         ll = km - l + 1
         do i=1,len
            tn0(i,l) = tin(i,ll) * pki(i) * (1.0/prh(l))
            qn0(i,l) = qin(i,ll)
            uoc(i,l) = uin(i,ll)
            voc(i,l) = vin(i,ll)
         enddo
      enddo
      tem2 = 1.0 / (sig(lm+1) - sig(lm))
      do l=1,nstrp
         tem  = (si(l) - si(l+1)) * tem2
         tem1 = tem * (1.0/prns(l))
         do i=1,len
            tns(i) = tns(i) + tin(i,l) * tem1
            qns(i) = qns(i) + qin(i,l) * tem
            uos(i) = uos(i) + uin(i,l) * tem
            vos(i) = vos(i) + vin(i,l) * tem
         enddo
      enddo
      do i=1,len
         tns(i)    = tns(i) * pki(i)
         tn0(i,lm) = tns(i)
         qn0(i,lm) = qns(i)
         uoc(i,lm) = uos(i)
         voc(i,lm) = vos(i)
         pcu(i)    = 0.0
      enddo
      do l=1,lm
         do i=1,len
            tcu(i,l) = 0.0
            qcu(i,l) = 0.0
            clt(i,l) = 0.0
            cwm0(i,l) = 0.0
            cmf0(i,l) = 0.0
            detr0(i,l) = 0.0
         enddo
      enddo
c
      call rasv (len, len, lm,  dt
     *,          ncrnd, krmin, krmax, frac, rasal, .true.
     *,          cp,  alhl, grav, crtmsf, rannum
     *,          sgb,    ods, prh, prj, fpk, hpk, sig(2)
     *,          pib,    pik, pki, gam
     *,          uoc,    voc, tn0, qn0
     *,          pcu,    clt, tcu, qcu
     *,          cd
     *,          cwm0, cmf0, detr0)
c
      tem = dt / 86400000.0
      do i=1,len
        evapc(i) = pcu(i) * tem
      enddo
c
      do l=1,lmm1
         ll = km - l + 1
         do i=1,len
            tin(i,ll) = tn0(i,l) * pik(i) * prh(l)
            qin(i,ll) = qn0(i,l)
            cwmin(i,ll) = cwmin(i,ll) + cwm0(i,l) 
            cmf(i,ll) = cmf(i,ll) + cmf0(i,l)
            detr(i,ll) = detr(i,ll) + detr0(i,l) 
         enddo
      enddo
      do i=1,len
         kuo(i) = 0
         kbot(i) = kmp1
         ktop(i) = 0
         tns(i) = (tn0(i,lm) - tns(i)) * pik(i) * prh(lm)
         qns(i) = qn0(i,lm) - qns(i)
         uos(i) = uoc(i,lm) - uos(i)
         vos(i) = voc(i,lm) - vos(i)
      enddo

      do l=lmm1,1,-1
         do i=1,len
            if (tcu(i,l) .ne. 0.0) then
               kuo(i) = 1
            endif
         enddo
      enddo
c  new test for convective clouds ! added ib 07/25/96
      do l=lmm1,1,-1
         do i=1,len
            if (clt(i,l) .gt. 0) ktop(i) = l
         enddo
      enddo
      do l=1,lmm1
         do i=1,len
            if (clt(i,l) .gt. 0) kbot(i) = l
         enddo
      enddo
      do i=1,len
         ktop(i) = kmp1 - ktop(i)
         kbot(i) = kmp1 - kbot(i)
      enddo
c
      do l=1,nstrp
         do i=1,len
            tin(i,l) = tin(i,l) + tns(i)
            qin(i,l) = qin(i,l) + qns(i)
         enddo
      enddo
      return
      end
