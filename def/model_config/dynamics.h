
/* dynamics options */
#undef  REDUCE_GRID  /* reduced grid for gsm.  defined if _jcap_ > 42  in depend.h */
#define RDAMPUX      /* RSM speed dependent diffusion to allow longer timestep */
#define PSPLIT       /* Process splitting time scheme for both gsm and rsm */
#undef  RSMVECTORIZE /* vectorize mpi rsm physics */
#undef  GET          /* use get for mpi communication in rsm */
#undef  VECSUM       /* use vectorized sum in mpi rsm */
#undef  HYBRID       /* use sigma-p hybrid */

