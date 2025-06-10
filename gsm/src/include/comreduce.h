#ifdef REDUCE_GRID
c
      integer lcapd,lonfd,lonfdp,lonfds,lcapdp
      common /comreduce/
     *                lcapd (latg_/2) ! regular grid in lat
     *               ,lonfd (latg_/2) ! regular grid in lat
#ifdef MP
     *               ,lcapdp(latg_/2,0:npes_-1) !regular grid in lat
     *               ,lonfdp(latg2p_,0:npes_-1) !symmetry spread grid in lat
     *               ,lonfds(latg2p_,0:npes_-1) !symmetry spread grid in lat
#endif
#endif
