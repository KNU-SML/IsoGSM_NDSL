c
c   begin compspec
c
        common /compspec/
     *    zea    (lln22p_,levsp_),
     *    dia    (lln22p_,levsp_),
     *    tea    (lln22p_,levsp_),
#ifndef NISLQ
     *    rqa    (lln22p_,levhp_),
#endif
     *    ulna   (lln22p_,levsp_),
     *    vlna   (lln22p_,levsp_),
     *    dpdphia(lln22p_),
     *    dpdlama(lln22p_),
     *    qa     (lln22p_),
     *    qlapa  (lln22p_)
       common /compspec/
     *    za     (lln22p_),
     *    uua    (lln22p_,levsp_),
     *    vva    (lln22p_,levsp_),
     *    ya     (lln22p_,levsp_),
#ifndef NISLQ
     *    rta    (lln22p_,levhp_),
#endif
     *    xa     (lln22p_,levsp_),
     *    wa     (lln22p_,levsp_)
     *
#ifdef NISLQ
        common /compspec/ rqa (lln22p_,levhp_)
        common /compspec/ rta (lln22p_,levhp_)
#endif
