c
c common block for diagnostics
c
c the number of diagnostic terms + 2
c
#ifdef T
      parameter ( ntermt = 15 )
      parameter ( nsumt = ntermt - 1 )
      parameter ( nbgtt = ntermt - 2 )
      parameter ( levbgtt = levs_ )
#endif
#ifdef V
      parameter ( ntermv = 14 )
      parameter ( nsumv = ntermv - 1 )
      parameter ( nbgtv = ntermv - 2 )
      parameter ( levbgtv = levs_ )
#endif
c
#ifdef DA
      parameter ( nterm2 = 2 )
#endif
#ifdef DT
      parameter ( nglpa_start = 1 )
      parameter ( nglpa_end = 5 )
      parameter ( nglpb_start = 6 )
      parameter ( nglpb_end = 13 )
      parameter ( nimp_start = 3 )
      parameter ( nimp_end = 5 )
      parameter ( ndebug = 4 )
#endif
#ifdef DV
      parameter ( nglpa_start = 1 )
      parameter ( nglpa_end = 9 )
      parameter ( nglpb_start = 10 )
      parameter ( nglpb_end = 12 )
      parameter ( ndebug = 2 )
#endif
c
#ifdef T
      integer levtem (levbgtt)
#endif
c
#ifdef V
      integer levvor (levbgtv)
#endif
c
c common blocks
c
#ifdef T
      common /combgtt/ btermst(lnt22_, levbgtt, ntermt)
#endif
c
#ifdef V
      common /combgtv/ btermsv(lnt22_, levbgtv, ntermv)
#endif
c
#ifdef A
      common /combgtr/
     &         ulnr    (lnt22_,levs_)
     &       , vlnr    (lnt22_,levs_)
#endif
c
#ifdef DT
      common /combgt2/ bterms2(lnt22_, levbgtt, nterm2, ndebug)
#endif
c
#ifdef DV
      common /combgt2/ bterms2(lnt22_, levbgtv, nterm2, ndebug)
#endif
