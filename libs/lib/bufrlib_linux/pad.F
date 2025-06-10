      subroutine pad(ibay,ibit,ibyt,ipadb)

      dimension ibay(*)

c----------------------------------------------------------------------
c----------------------------------------------------------------------

c  pad the subset to an ipadb bit boundary
c  ----------------------------------------

      ipad = ipadb - mod(ibit+8,ipadb)
      call pkb(ipad,8,ibay,ibit)
      call pkb(   0,ipad,ibay,ibit)
      ibyt = ibit/8

      if(mod(ibit,ipadb).ne.0) goto 900
      if(mod(ibit,8    ).ne.0) goto 900

      return
900   call abt('pad - bit pad failure              ')
      end
