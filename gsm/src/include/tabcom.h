c
c     common block tabcom contains quantities precomputed in subroutine
c     table for use in the longwave radiation program:
c          em1     =  e1 function, evaluated over the 0-560 and
c                     1200-2200 cm-1 intervals
c          em1wde  =  e1 function, evaluated over the 160-560 cm-1
c                     interval
c          table1  =  e2 function, evaluated over the 0-560 and
c                     1200-2200 cm-1 intervals
c          table2  =  temperature derivative of table1
c          table3  =  mass derivative of table1
c          em3     =  e3 function, evaluated over the 0-560 and
c                     1200-2200 cm-1 intervals
c          source  =  planck function, evaluated at specified temps. for
c                     bands used in cts calculations
c          dsrce   =  temperature derivative of source
c          ind     =  index, with value ind(i)=i. used in fst88
c          indx2   =  index values used in obtaining "lower triangle"
c                     elements of avephi,etc.,in fst88
c          kmaxv   =  index values used in obtaining "upper triangle"
c                     elements of avephi,etc.,in fst88
c          kmaxvm  =  kmaxv(l),used for do loop indices
c
      common / tabcom / ind(imax),indx2(lp1v),kmaxv(lp1),
     1     kmaxvm,idummy2(imax+lp1v+lp1+1)
      common/tabcom/em1(28,180),em1wde(28,180),table1(28,180),
     1 table2(28,180),table3(28,180),em3(28,180),source(28,nbly),
     2 dsrce(28,nbly)
