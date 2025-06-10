c
#include <constant.h>
#include <define.h>
#ifdef JCAP 
      integer jcap_
      parameter(jcap_=_jcap_)
#endif
#ifdef LONF
      integer lonf_
      parameter(lonf_=_lonf_)
#endif
#ifdef LATG
      integer latg_
      parameter(latg_=_latg_)
#endif
#ifdef LEVS
      integer levs_
      parameter(levs_=_levs_)
#endif

#ifdef LSOIL
      integer lsoil_
      parameter(lsoil_=_lsoil_)
#endif

#ifdef NSOIL
      integer nsoil_
      parameter(nsoil_=_nsoil_)
#endif

#ifdef MSUB
      integer msub_
      parameter(msub_=_msub_)
#endif

#ifdef LALBD
      integer lalbd_
      parameter(lalbd_=_lalbd_)
#endif

#ifdef SAS
#ifdef NCLDTOP
      integer ncldtop_
      parameter(ncldtop_=_ncldtop_)
#endif
#endif

#ifdef NCPUS
      integer ncpus_
      parameter(ncpus_=_ncpus_)
#endif

#ifdef JCAP1 
      integer jcap1_
      parameter(jcap1_=jcap_+1)
#endif
#ifdef JCAP2 
      integer jcap2_
      parameter(jcap2_=jcap_+2)
#endif
#ifdef LNT
      integer lnt_
      parameter(lnt_=(jcap_+1)*(jcap_+2)/2)
#endif
#ifdef LNT2
      integer lnt2_
      parameter(lnt2_=(jcap_+1)*(jcap_+2))
#endif
#ifdef LNT22
      integer lnt22_
      parameter(lnt22_=(jcap_+1)*(jcap_+2)+1)
#endif
#ifdef LNUT22
      integer lnut22_
      parameter(lnut22_=((jcap_+2)*(jcap_+1)/2+(jcap_+1))*2+1)
#endif



#ifdef RSM
#ifdef IGRD
      integer igrd_
      parameter(igrd_=_igrd_)
#endif
#ifdef JGRD
      integer jgrd_
      parameter(jgrd_=_jgrd_)
#endif
#ifdef CIGRD1
      integer cigrd1_
      parameter(cigrd1_=_cigrd1_)
#endif
#ifdef CJGRD1
      integer cjgrd1_
      parameter(cjgrd1_=_cjgrd1_)
#endif
#ifdef BORDER
      integer border_
      parameter(border_=_border_)
#endif
#ifdef BGF
      integer bgf_
      parameter(bgf_=_bgf_)
#endif
#ifdef IGRD1
      integer igrd1_
      parameter(igrd1_=igrd_+1)
#endif
#ifdef JGRD1
      integer jgrd1_
      parameter(jgrd1_=jgrd_+1)
#endif
#ifdef IGRD2
      integer igrd2_
      parameter(igrd2_=igrd_*2)
#endif
#ifdef IGRD12
      integer igrd12_
      parameter(igrd12_=(igrd_+1)*2)
#endif
#ifdef JGRD12
      integer jgrd12_
      parameter(jgrd12_=(jgrd_+1)/2)
#endif
#ifdef IGRDCUT
      integer igrdcut_
      parameter(igrdcut_=igrd_-12)
#endif
#ifdef IWAV
      integer iwav_
      parameter(iwav_=(igrd_-12)/3*2)
#endif
#ifdef IWAV1
      integer iwav1_
      parameter(iwav1_=(igrd_-12)/3*2+1)
#endif
#ifdef JGRDCUT
      integer jgrdcut_
      parameter(jgrdcut_=jgrd_-12)
#endif
#ifdef JWAV
      integer jwav_
      parameter(jwav_=(jgrd_-12)/3*2)
#endif
#ifdef JWAV1
      integer jwav1_
      parameter(jwav1_=(jgrd_-12)/3*2+1)
#endif
#ifdef LEVR
      integer levr_
      parameter(levr_=_levr_)
#endif
#ifdef LNWAV
      integer lnwav_
      parameter(lnwav_=((igrd_-12)/3*2+1)*((jgrd_-12)/3*2+1))
#endif
#ifdef IGRD1I
      integer igrd1i_
      parameter(igrd1i_=igrd_+1)
#endif
#ifdef JGRD1I
      integer jgrd1i_
      parameter(jgrd1i_=jgrd_+1)
#endif
#ifdef IGRD1O
      integer igrd1o_
      parameter(igrd1o_=igrd_+1)
#endif
#ifdef JGRD1O
      integer jgrd1o_
      parameter(jgrd1o_=jgrd_+1)
#endif
#ifdef JUMPR
      integer jumpr_
      parameter(jumpr_=igrd_*2+3)
#endif
#ifdef RGEN
      integer rgen_
      parameter(rgen_=80)
#endif
#ifdef LNGRD
      integer lngrd_
      parameter(lngrd_=(igrd_+1)*(jgrd_+1))
#endif
#ifdef CLNGRD
      integer clngrd_
      parameter(clngrd_=cigrd1_*cjgrd1_)
#endif
#ifdef RLPNT
      integer rlpnt_
      parameter(rlpnt_=50)     
#endif
#ifdef RLTSTP
      integer rltstp_
      parameter(rltstp_=50)     
#endif
#ifdef RSLVARK
      integer rslvark_
      parameter(rslvark_=50)     
#endif
#ifdef RMLVARK
      integer rmlvark_
      parameter(rmlvark_=50)     
#endif
#ifdef RELX
      integer relx_
      parameter(relx_=_relx_)     
#endif
#ifdef DIFUH
      integer difuh_
      parameter(difuh_=_difuh_)     
#endif
#ifdef DIFUM
      integer difum_
      parameter(difum_=_difum_)     
#endif
#ifdef BSMOOTH
      integer bsmooth_
      parameter(bsmooth_=1)     
#endif
#endif

 


#ifdef NTRAC
      integer ntrac_
      parameter(ntrac_=_ntrac_)
#endif
#ifdef NCLDG
      integer ncldg_
      parameter(ncldg_=_ncldg_)
#endif
#ifdef NCLD
      integer ncld_
      parameter(ncld_=_ncld_)
#endif
#ifdef NTOTAL
      integer ntotal_
      parameter(ntotal_=ntrac_+ncldg_)
#endif
#ifdef LEVH
      integer levh_
      parameter(levh_=levs_*(ntrac_+ncldg_))
#endif

#ifdef LEVP1
      integer levp1_
      parameter(levp1_=levs_+1)
#endif
#ifdef LEVM1
      integer levm1_
      parameter(levm1_=levs_-1)
#endif
#ifdef LNUT
      integer lnut_
      parameter(lnut_=(jcap_+1)*(jcap_+2)/2+jcap_+1)
#endif
#ifdef LNUT2
      integer lnut2_
      parameter(lnut2_=2*((jcap_+1)*(jcap_+2)/2+jcap_+1))
#endif
#ifdef LNUV
      integer lnuv_
      parameter(lnuv_=(jcap_+1)*(jcap_+2))
#endif
#ifdef LNU2
      integer lnu2_
      parameter(lnu2_=(jcap_+1)*(jcap_+2)*2)
#endif
#ifdef TWOJ1
      integer twoj1_
      parameter(twoj1_=2*(jcap_+1))
#endif
#ifdef LONF2
      integer lonf2_
      parameter(lonf2_=lonf_*2)
#endif
#ifdef LATG2
      integer latg2_
      parameter(latg2_=latg_/2)
#endif
#ifdef LATG2SSI
      integer latg2_
      parameter(latg2_=latg_/2+1)
#endif
#ifdef LONF22
      integer lonf22_
#ifdef SCM
      parameter(lonf22_=lonf_*2+1)
#else
      parameter(lonf22_=lonf_*2+2)
#endif
#endif
#ifdef LOTANU
      integer lotanu_
      parameter(lotanu_=4*levs_+1)
#endif

#ifdef GSMTNRES
      real mtnres_
      parameter(mtnres_=_gsmtnres_)
#endif

#ifdef RSMTNRES
      real mtnres_
      parameter(mtnres_=_rsmtnres_)
#endif

#ifdef MTNVAR
      integer mtnvar_
      parameter(mtnvar_=_mtnvar_)
#endif
#ifdef IMN
      integer imn_
      parameter(imn_=360*60/mtnres_)
#endif
#ifdef JMN
      integer jmn_
      parameter(jmn_=180*60/mtnres_)
#endif
#ifdef SLVARK
      integer slvark_
      parameter(slvark_=_slvark_)
#endif
#ifdef MLVARK
      integer mlvark_
      parameter(mlvark_=_mlvark_)
#endif
#ifdef LPNT
      integer lpnt_
      parameter(lpnt_=_lpnt_)
#endif
#ifdef LTSTP
      integer ltstp_
      parameter(ltstp_=_ltstp_)
#endif
#ifdef LOTFBL
      integer lotfbl_
      parameter(lotfbl_=4*levs_)
#endif
#ifdef LOTBA
      integer lotba_
      parameter(lotba_=8*levs_)
#endif
#ifdef LVAR
      integer lvar_
      parameter(lvar_=1)
#endif
#ifdef LDIAG
      integer ldiag_
      parameter(ldiag_=1)
#endif
#ifdef LGRID
      integer lgrid_
      parameter(lgrid_=0)
#endif
#ifdef IGEN
      integer igen_
      parameter(igen_=80)
#endif
#ifdef NST
      integer nst_
      parameter(nst_=1)     
#endif
#ifdef NT
      integer nt_
      parameter(nt_=15)     
#endif
#ifdef NQ
      integer nq_
      parameter(nq_=10)     
#endif
#ifdef NU
      integer nu_
      parameter(nu_=11)     
#endif
#ifdef NV
      integer nv_
      parameter(nv_=11)     
#endif
#ifdef NP
      integer np_
      parameter(np_=6)     
#endif

#ifdef LONSSI
      integer lonssi_
      parameter(lonssi_=_lonssi_)
#endif
#ifdef LATSSI
      integer latssi_
      parameter(latssi_=lonssi_/2+1)
#endif
#ifdef MAXLEV
      integer maxlev_
      parameter(maxlev_=_maxlev_)
#endif
#ifdef MAXREP
      integer maxrep_
      parameter(maxrep_=_maxrep_)
#endif
#ifdef MAXCYL
      integer maxcyl_
      parameter(maxcyl_=_maxcyl_)
#endif
#ifdef MAXPAK
      integer maxpak_
      parameter(maxpak_=_maxpak_)
#endif
#ifdef NQC
      integer nqc_
      parameter(nqc_=_nqc_)
#endif
#ifdef NMX
      integer nmx_
      parameter(nmx_=_nmx_)
#endif
#ifdef N
      integer n_
      parameter(n_=_n_)
#endif

#ifdef NTDATA
      integer ntdata_
      parameter(ntdata_=_ntdata_)
#endif
#ifdef NSDATA
      integer nsdata_
      parameter(nsdata_=_nsdata_)
#endif
#ifdef NWDATA
      integer nwdata_
      parameter(nwdata_=_nwdata_)
#endif
#ifdef NPDATA
      integer npdata_
      parameter(npdata_=_npdata_)
#endif
#ifdef NQDATA
      integer nqdata_
      parameter(nqdata_=_nqdata_)
#endif
#ifdef NPWDAT
      integer npwdat_
      parameter(npwdat_=_npwdat_)
#endif
#ifdef NSPROF
      integer nsprof_
      parameter(nsprof_=_nsprof_)
#endif
#ifdef NSIGDIVT
      integer nsigdivt_
      parameter(nsigdivt_=_nsigdivt_)
#endif
#ifdef JCAPDIVT
      integer jcapdivt_
      parameter(jcapdivt_=_jcapdivt_)
#endif
#ifdef JCAPSTAT
      integer jcapstat_
      parameter(jcapstat_=_jcapstat_)
#endif
#ifdef NMDSZH
      integer nmdszh_
      parameter(nmdszh_=_nmdszh_)
#endif
#ifdef NSIGSAT
      integer nsigsat_
      parameter(nsigsat_=_nsigsat_)
#endif

#ifdef ILONF
      integer ilonf_
      parameter(ilonf_=_ilonf_)
#endif
#ifdef ILATG
      integer ilatg_
      parameter(ilatg_=_ilatg_)
#endif
#ifdef IJCAP
      integer ijcap_
      parameter(ijcap_=_ijcap_)
#endif
#ifdef ILEVS
      integer ilevs_
      parameter(ilevs_=_ilevs_)
#endif

#ifdef LEVRI
      integer levri_
      parameter(levri_=_levs_)
#endif
#ifdef LEVRO
      integer levro_
      parameter(levro_=_levs_)
#endif


#ifdef IO
      integer io_
      parameter(io_=_io_)
#endif
#ifdef JO
      integer jo_
      parameter(jo_=_jo_)
#endif
#ifdef KO
      integer ko_
      parameter(ko_=_ko_)
#endif
#ifdef RIVER
#ifdef IO2
      integer io2_
      parameter(io2_=_io2_)
#endif
#ifdef JO2
      integer jo2_
      parameter(jo2_=_jo2_)
#endif
#endif
#ifdef IOSGB
      integer io_
      parameter(io_=_lonf_)
#endif
#ifdef JOSGB
      integer jo_
      parameter(jo_=_latg_)
#endif
#ifdef KT
      integer kt_
      parameter(kt_=_kt_)
#endif
#ifdef KTT
      integer ktt_
      parameter(ktt_=_ktt_)
#endif
#ifdef KZZ
      integer kzz_
      parameter(kzz_=_kzz_)
#endif
#ifdef MOO
      integer moo_
      parameter(moo_=_moo_)
#endif

#ifdef RERTH
      real rerth_
      parameter(rerth_=_rerth_)
#endif
#ifdef RRERTH
      real rrerth_
      parameter(rrerth_=_rrerth_)
#endif
#ifdef G
      real g_
      parameter(g_=_g_)
#endif
#ifdef OMEGA
      real omega_
      parameter(omega_=_omega_)
#endif
#ifdef RD
      real rd_
      parameter(rd_=_rd_)
#endif
#ifdef RV
      real rv_
      parameter(rv_=_rv_)
#endif
#ifdef RVRDM1
      real rvrdm1_
      parameter(rvrdm1_=_rvrdm1_)
#endif
#ifdef CP
      real cp_
      parameter(cp_=_cp_)
#endif
#ifdef CV
      real cv_
      parameter(cv_=_cv_)
#endif
#ifdef CVAP
      real cvap_
      parameter(cvap_=_cvap_)
#endif
#ifdef CLIQ
      real cliq
      parameter(cliq_=_cliq_)
#endif
#ifdef CICE
      real cice_
      parameter(cice_=_cice_)
#endif
#ifdef HVAP
      real hvap_
      parameter(hvap_=_hvap_)
#endif
#ifdef HSUB
      real hsub_
      parameter(hsub_=_hsub_)
#endif
#ifdef HFUS
      real hfus_
      parameter(hfus_=_hfus_)
#endif
#ifdef PSAT
      real psat_
      parameter(psat_=_psat_)
#endif
#ifdef SBC
      real sbc_
      parameter(sbc_=_sbc_)
#endif
#ifdef SOLR
      real solr_
      parameter(solr_=_solr_)
#endif
#ifdef T0C
      real t0c_
      parameter(t0c_=_t0c_)
#endif
#ifdef TTP
      real ttp_
      parameter(ttp_=_ttp_)
#endif
#ifdef CAL
      real cal_
      parameter(cal_=_cal_)
#endif
#ifdef QMIN
      real qmin_
      parameter(qmin_=_qmin_)
#endif
#ifdef RHOH2O
      real rhoh2o_
      parameter(rhoh2o_=_rhoh2o_)
#endif
#ifdef FV
      real fv_
      parameter(fv_=rv_/rd_-1.)
#endif
#ifdef GOCP
      real gocp_
      parameter(gocp_=g_/cp_)
#endif
#ifdef PI
      real pi_
      parameter(pi_=_pi_)
#endif
#ifdef KLOT
      integer klot_
      parameter(klot_=_levs_)
#endif
#ifdef NSTYPE 
      integer nstype_
      parameter(nstype_=_nstype_)
#endif
#ifdef NVTYPE 
      integer nvtype_
      parameter(nvtype_=_nvtype_)
#endif
#ifdef NPES
      integer npes_
      parameter(npes_=_npes_)
#endif


#ifdef MP
  

#ifdef NCOL
      integer ncol_
      parameter(ncol_=_ncol_)
#endif
#ifdef NROW
      integer nrow_
      parameter(nrow_=_nrow_)
#endif
#ifdef LEVSP
      integer levsp_
      parameter(levsp_=_levsp_)
#endif
#ifdef LEVHP
      integer levhp_
      parameter(levhp_=_levsp_*ntotal_)
#endif
#ifdef LEVP1P
      integer levp1p_
      parameter(levp1p_=levsp_+1)
#endif
#ifdef LEVM1P
      integer levm1p_
      parameter(levm1p_=levsp_-1)
#endif
#ifdef LONFP
      integer lonfp_
      parameter(lonfp_=_lonfp_)
#endif
#ifdef LONF2P
      integer lonf2p_
      parameter(lonf2p_=lonfp_*2)
#endif
#ifdef LONF22P
      integer lonf22p_
      parameter(lonf22p_=lonfp_*2+2)
#endif
#ifdef LATGP
      integer latgp_
      parameter(latgp_=_latgp_)
#endif
#ifdef LATG2P
      integer latg2p_
      parameter(latg2p_=latgp_/2)
#endif
#ifdef JCAPP
      integer jcapp_
      parameter(jcapp_=_jcapp_)
#endif
#ifdef JCAP1P
      integer jcap1p_
      parameter(jcap1p_=jcapp_+1)
#endif
#ifdef LCAPP
      integer lcapp_
      parameter(lcapp_=(jcapp_+1)*2)
#endif
#ifdef LCAP22P
      integer lcap22p_
      parameter(lcap22p_=(jcapp_+1)*4+2)
#endif
#ifdef LNTP
      integer lntp_
      parameter(lntp_=_lntp_)
#endif
#ifdef LNT2P
      integer lnt2p_
      parameter(lnt2p_=lntp_*2)
#endif
#ifdef LNT22P
      integer lnt22p_
      parameter(lnt22p_=lntp_*2+1)
#endif
#ifdef LLNP
      integer llnp_
      parameter(llnp_=_llnp_)
#endif
#ifdef LLN2P
      integer lln2p_
      parameter(lln2p_=llnp_*2)
#endif
#ifdef LLN22P
      integer lln22p_
      parameter(lln22p_=llnp_*2+1)
#endif

#ifdef RSM

#ifdef IGRD1P
      integer igrd1p_
      parameter(igrd1p_=_igrd1p_)
#endif
#ifdef IGRD12P
      integer igrd12p_
      parameter(igrd12p_=igrd1p_*2)
#endif
#ifdef JGRD1P
      integer jgrd1p_
      parameter(jgrd1p_=_jgrd1p_)
#endif
#ifdef JGRD12P
      integer jgrd12p_
      parameter(jgrd12p_=jgrd1p_/2)
#endif
#ifdef LNGRDP
      integer lngrdp_
      parameter(lngrdp_=igrd1p_*jgrd1p_)      
#endif
#ifdef IWAV1P
      integer iwav1p_
      parameter(iwav1p_=_iwav1p_)
#endif
#ifdef LEVRP
      integer levrp_
      parameter(levrp_=_levrp_)
#endif
#ifdef LLWAVP
      integer llwavp_
      parameter(llwavp_=_llwavp_)
#endif
#ifdef LNWAVP
      integer lnwavp_
      parameter(lnwavp_=_lnwavp_)
#endif

#endif

#endif
