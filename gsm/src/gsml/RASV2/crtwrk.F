      subroutine crtwrk(pl, acr)
!
      real pl, acr
      integer iwk
      real ac(16), ad(16)
      common /rasacr/ ac, ad
!
      iwk = pl * 0.02 - 0.999999999
      iwk = max(1, min(iwk,16))
      acr = ac(iwk) + pl * ad(iwk)
!
      return
      end
