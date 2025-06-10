c
      parameter( ngwri=2*((mwave+1)*(mwave+2)) )
      parameter( ntrns=(mwave+1)*(mwave+2) )
      parameter( nfl22=(mwave+1)*2*1024+
     1                 (mwave+1)*(mwave+2) )
      parameter( npln3=(mwave+1)*(mwave+2)*2+(mwave+1)*(mwave+1)*2+
     1                 (mwave+1)*2*2+
     2                 (((mwave+1)*(mwave+2))/2)*2*2 )
      parameter( nsums=(mwave+1)*2*(mwave+1) )
      parameter( nsyms=idim*2*kdim )
c
      parameter(ndjm1=(nfl22/npln3*nfl22+npln3/nfl22*npln3)/
     1                (nfl22/npln3      +npln3/nfl22      ))
      parameter(nscrch=ndjm1+1)
