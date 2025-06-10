      subroutine mvb(ib1,nb1,ib2,nb2,nbm)

      dimension ib1(*),ib2(*)        

c-----------------------------------------------------------------------
cfpp$ expand(upb,pkb)
c-----------------------------------------------------------------------

      jb1 = 8*(nb1-1)
      jb2 = 8*(nb2-1)

      do n=1,nbm
      call upb(nval,8,ib1,jb1)      
      call pkb(nval,8,ib2,jb2)      
      enddo

      return
      end
