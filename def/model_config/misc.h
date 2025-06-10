
/* miscellaneous options */
#define STDAMP             /* Spectral tendency damping scheme for rsm */
#define STDAMP3            /* area average log ps correction to be used with stdamp */
#undef  GSMDAMP            /* Spectral nudging scheme for gsm */
#define MRG_POST           /* merge postprocessor into forecast code */
#undef  PGB_PARALLEL       /* pgb parallel io and parallel computation */
#define PGB_SEMIPARALLEL   /* pgb parallel computation but serial io */
#undef  PGB_NOPARALLEL     /* serial pgb */
#undef  SETMTNZERO         /* set orography to zero over ocean */
#undef  NOPRINT            /* suppress print out */
#undef  RESET_FHOUR        /* reset fhour to zero at the end of the year */
#undef  RMPGP2F_GLB        /* use global memory */ 
#undef  BASEOUT            /* rsm output post processed base field */
#undef  ADJSOILW           /* adjustment of soil moisture using observed precipitation */
#undef  NFDR               /* national fire danger rating system */
#undef  RIV1               /* river routing scheme, 1deg */
#undef  RIV05              /* river routing scheme, 0.5deg */
#undef  SPNG_WID	   /* wide (original) sponge zone in rsml/rltbini.F */
#define SPNG_NRW           /* narrow sponge zone in rsml/rltbini.F */
#undef  CO2PROG            /* prognostic co2 */
#define ZERBZDIA           /* zero out precip at the boundary buffer zone */
#define SMTHBZ             /* smoothing in the buffer zone for cosmetics */
#define LEHMANN            /* Lehmann's boundary nudging coeffs */

#define COMR4 true         /* communication with real*4 */
#define ISOTOPE            /* enable isotope prediction (RAS is needed)*/
#define NOFRAC             /* vapor tracer mode (ISOTOPE is needed) */

#define NISLQ
#define NISLQ_MASS
#define NISLQ_PHYS
