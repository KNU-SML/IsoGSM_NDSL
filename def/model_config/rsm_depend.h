/* dependency for rsm */

#define IGRDM1 `expr _igrd_ - 1`
#define JGRDM1 `expr _jgrd_ - 1`

#if ( _rdelx_ == 190000 )
#define DELTAT_REG_SPRING 400
#define DELTAT_REG_SUMMER 600
#define DELTAT_REG_FALL 400
#define DELTAT_REG_WINTER 360
#endif
#if ( _rdelx_ == 60000 )
#define DELTAT_REG_SPRING 400
#define DELTAT_REG_SUMMER 600
#define DELTAT_REG_FALL 400
#define DELTAT_REG_WINTER 360
#endif
#if ( _rdelx_ == 50000 )
#define DELTAT_REG_SPRING 400
#define DELTAT_REG_SUMMER 600
#define DELTAT_REG_FALL 400
#define DELTAT_REG_WINTER 360
#endif
#if ( _rdelx_ == 40000 )
#define DELTAT_REG_SPRING 300
#define DELTAT_REG_SUMMER 420
#define DELTAT_REG_FALL 300
#define DELTAT_REG_WINTER 200
#endif
#if ( _rdelx_ == 30000 )
#define DELTAT_REG_SPRING 180
#define DELTAT_REG_SUMMER 200
#define DELTAT_REG_FALL 180
#define DELTAT_REG_WINTER 120
#endif
#if ( _rdelx_ == 25000 )
#define DELTAT_REG_SPRING 150
#define DELTAT_REG_SUMMER 180
#define DELTAT_REG_FALL 150
#define DELTAT_REG_WINTER 100
#endif
#if ( _rdelx_ == 20000 )
#define DELTAT_REG_SPRING 90
#define DELTAT_REG_SUMMER 100
#define DELTAT_REG_FALL 90
#define DELTAT_REG_WINTER 60
#endif
#if ( _rdelx_ == 15000 )
#define DELTAT_REG_SPRING 60
#define DELTAT_REG_SUMMER 75
#define DELTAT_REG_FALL 60
#define DELTAT_REG_WINTER 45
#endif
#if ( _rdelx_ == 10000 )
#define DELTAT_REG_SPRING 40
#define DELTAT_REG_SUMMER 50
#define DELTAT_REG_FALL 40
#define DELTAT_REG_WINTER 30
#endif
