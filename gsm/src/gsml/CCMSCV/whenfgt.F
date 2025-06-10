      subroutine whenfgt(n,array,inc,target,index,nval)
c
c $Id: whenfgt.F,v 1.2 2009/10/02 00:35:47 kana Exp $
c $Author: kana $
c
      dimension index(*), array(*)
      ina=1
      nval=0
      if(inc .lt. 0) ina=(-inc)*(n-1)+1
      do 100 i=1,n
        if(array(ina) .gt. target) then
          nval=nval+1
          index(nval)=i
        endif
        ina=ina+inc
 100  continue
      return
      end           
 
