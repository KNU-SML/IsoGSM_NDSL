#define NHM_MP                   /* mpi */
#undef NHM_KF                    /* KF parameterization */
#define NHM_NP                   /* no cloud parameterization */
#undef NHM_HYRDO                 /* Hydrostatic mode */
#define DECOMPOSITION 2          /* MPI decompostion dimension */

#ifdef NHM_KF
#undef NHM_NP
#define MSWSYS22 20
#endif
#ifdef NHM_NP
#define MSWSYS22 0
#endif

#ifdef NHM_HYDRO
#define MSWSYS20 -1
#else
#define MSWSYS20 2
#endif

#if ( DECOMPOSITION == 1 ) 
#define X_NUM 1
#endif
#if ( DECOMPOSITION == 2 ) 
#define X_NUM -1
#endif

