c
c    common block bandta contains random band parameters for the lw
c    calculations using 10 cm-1 wide bands.the 15 um co2 complex
c    is 2 bands,560-670 and 670-800 cm-1. ozone coefficients are
c    in 3 bands,670-800 (14.1 um),990-1070 and 1070-1200 (9.6 um).
c    the  (nblw) bands now include:
c                56 bands, 10  cm-1 wide    0  -   560  cm-1
c                 2 bands, 15 um complex  560  -   670  cm-1
c                                         670  -   800  cm-1
c                 3 "continuum" bands     800  -   900  cm-1
c                                         900  -   990  cm-1
c                                        1070  -   1200 cm-1
c                 1 band for 9.6 um band  990  -   1070 cm-1
c               100 bands, 10 cm-1 wide  1200  -   2200 cm-1
c                 1 band for 4.3 um src  2270  -   2380 cm-1
c    thus nblw presently equals    163
c    all bands are arranged in order of increasing wavenumber
c
c        arndm   =   random "a" parameter for (nblw) bands
c        brndm   =   random "b" parameter for (nblw) bands
c        betad   =   continuum coefficients for (nblw) bands
c        ap,bp   =   capphi coefficients for (nblw) bands
c        atp,btp =   cappsi coefficients for (nblw) bands
c        bandlo  =   lowest frequency in each of (nblw) freq. bands
c        bandhi  =   highest frequency in each of (nblw) freq. bands
c        ao3rnd  =   random "a" parameter for ozone in (3) ozone
c                    bands
c        bo3rnd  =   random "b" parameter for ozone in (3) ozone
c                    bands
c        ab15    =   the product arndm*brndm for the two bands
c                    representing the 15 um band complex of co2
c     data for arndm,brndm,ap,bp,atp,btp,ao3rnd,bo3rnd are obtained by
c     using the afgl 1982 catalog. continuum coefficients are from
c     roberts (1976).
c
      common / bandta / arndm(nblw),brndm(nblw),betad(nblw),ap(nblw),
     1                  bp(nblw),atp(nblw),btp(nblw),bandlo(nblw),
     2                  bandhi(nblw),ao3rnd(3),bo3rnd(3),ab15(2)
c
c
c    common block bdwide contains random band parameters for specific
c    wide bands. at present,the information consists of 1) random
c    model parameters for the 15 um band,560-800 cm-1; 2) the
c    continuum coefficient for the 800-990,1070-1200 cm-1 band
c        specifically:
c        awide       =   random "a" parameter for  band
c        bwide       =   random "b" parameter for  band
c        betawd      =   continuum coefficients for band
c        apwd,bpwd   =   capphi coefficients for  band
c        atpwd,btpwd =   cappsi coefficients for band
c        bdlowd      =   lowest frequency in each  freq  band
c        bdhiwd      =   highest frequency in each freq  band
c        ab15wd      =   the product arndm*brndm for the one band
c                        representing the 15 um band complex of co2
c        betinw      =   cont.coefficient for a specified wide
c                        freq.band (800-990 and 1070-1200 cm-1).
c        sko2d       =   1./betinw, used in spa88 for cont. coeffs
c        skc1r       =   betawd/betinw, used for cont. coeff. for
c                        15 um band in fst88
c        sko3r       =   ratio of cont. coeff. for 9.9 um band to
c                        betinw, used for 9.6 um cont coeff in fst88
c     data for awide,bwide,apwd,bpwd,atpwd,btpwd,ao3wd,bo3wd are
c     obtained by using the afgl 1982 catalog. continuum coefficients
c     are from roberts (1976).
c
      common / bdwide / awide,bwide,betawd,
     1                  apwd,bpwd,atpwd,btpwd,
     2                  bdlowd,bdhiwd,betinw,
     3                  ab15wd,sko2d,skc1r,sko3r
c
c    common block bdcomb contains random band parameters for the lw
c    calculations using combined wide frequency bands between 160 and
c    1200 cm-1,as well as the 2270-2380 band for source calc.
c        bands 1-8: combined wide frequency bands for 160-560 cm-1
c        bands 9-14: frequency bands,as in bandta (narrow bands)
c                    for 560-1200 cm-1
c        band  15:  frequency band 2270-2380 cm-1,used for source
c                   calculation only
c        thus nbly presently equals   15
c
c        bands are arranged in order of increasing wavenumber
c        acomb       =   random "a" parameter for (nbly) bands
c        bcomb       =   random "b" parameter for (nbly) bands
c        betacm      =   continuum coefficients for (nbly) bands
c        apcm,bpcm   =   capphi coefficients for (nbly) bands
c        atpcm,btpcm =   cappsi coefficients for (nbly) bands
c        bdlocm      =   lowest frequency in each of (nbly) freq. bands
c        bdhicm      =   highest frequency in each of (nbly) freq. bands
c        ao3cm       =   random "a" parameter for ozone in (3) ozone
c                        bands
c        bo3cm       =   random "b" parameter for ozone in (3) ozone
c                        bands
c        ab15cm      =   the product arndm*brndm for the two bands
c                        representing the 15 um band complex of co2
c        betinc      =   cont.coefficient for a specified wide
c                        freq.band (800-990 and 1070-1200 cm-1).
c        iband       =   index no of the 40 wide bands used in
c                        combined wide band calculations. in other
c                        words,index telling which of the 40 wide
c                        bands between 160-560 cm-1 are included in
c                        each of the first 8 combined wide bands
c     data for acomb,bcomb,apcm,bpcm,atpcm,btpcm,ao3cm,bo3cm are
c     obtained by using the afgl 1982 catalog. continuum coefficients
c     are from roberts (1976). iband index values are obtained by
c     experimentation.
c
      common / bdcomb / iband(40),acomb(nbly),bcomb(nbly),
     1                  betacm(nbly),apcm(nbly),bpcm(nbly),atpcm(nbly),
     2                  btpcm(nbly),bdlocm(nbly),bdhicm(nbly),betinc,
     3                  ao3cm(3),bo3cm(3),ab15cm(2)
