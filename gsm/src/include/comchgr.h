c
c  blank commons
c
      common idate(4)
c
c  common for input arrays
c
      common qi(mdimi),tei(mdimi,kdimi),dii(mdimi,kdimi),
     1       zei(mdimi,kdimi),rqi(mdimi,kdimqi),gzi(mdimi)
c
c  common for modified horizontal resolution
c
      common q  (mdim),te (mdim,kdimi),di (mdim,kdimi),ze (mdim,kdimi),
     1       rq (mdim,kdimqi),gz(mdim)
c
      common ps (idimt),tf (idimt,kdimi),rqf (idimt,kdimqi),
     1       uf (idimt,kdimi),vf (idimt,kdimi)
c
c  common for arrays of modified vertical resolution
c
      common qo (mdim),teo(mdim,kdim),dio(mdim,kdim),zeo(mdim,kdim),
     1       rqo(mdim,kdimq),gzo(mdim)
c
      common pso(idimt),tfo(idimt,kdim),rqfo(idimt,kdimq),
     1       ufo(idimt,kdim),vfo(idimt,kdim)
