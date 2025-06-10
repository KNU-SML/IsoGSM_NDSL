#define SODA
#undef  LEVITUS
#define R2FRC
#undef  CaRD10V2
#ifdef  SODA
#define ROMS_STARTUP "SODA"
#define ROMS_INI_NLI 40
#elif defined LEVITUS
#define ROMS_STARTUP "LEVITUS"
#define ROMS_INI_NLI 33
#else
#define ROMS_STARTUP "GROMS"
#define ROMS_INI_NLI 20
#endif
#ifdef  SODA
#define ROMS_BNY_NLI 40
#else
#define ROMS_BNY_NLI 20
#endif
#define ROMS_INI_THETA_S 4
#define ROMS_INI_THETA_B 0.4
#define ROMS_INI_TCLINE  75
#define ROMS_INI_TSTART  15
#define ROMS_INI_TEND    15
#define ROMS_INI_SMTHINP 1
#define ROMS_BNY_THETA_S 4
#define ROMS_BNY_THETA_B 0.4
#define ROMS_BNY_TCLINE  75
#define ROMS_BNY_TSTART  15
#define ROMS_BNY_TEND    15
#define ROMS_BNY_SMTHINP 1
#define ROMS_FOR_NLI 20
#define ROMS_FOR_THETA_S 4
#define ROMS_FOR_THETA_B 0.4
#define ROMS_FOR_TCLINE  75
#define ROMS_FOR_TSTART  15
#define ROMS_FOR_TEND    15
#define ROMS_FOR_SMTHINP 1
#define ROMS_THETA_S 4
#define ROMS_THETA_B 0.4
#define ROMS_TCLINE  75
