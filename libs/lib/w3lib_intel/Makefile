LIB = w3lib.a
OBJS = aea.o assign.o datimx.o \
gbyte.o gbytes.o getgbss.o iw3jdn.o iw3pds.o pdsens.o pdseup.o \
sbyte.o sbytes.o swap32.o  unpindx.o w3fi01.o w3fi04.o w3fi63.o \
w3fi83.o w3fp11.o w3fs26.o w3pack.o w3tagb.o xmovex.o \
baopen.o baread.o getgir.o ixgb.o rdgb.o skgb.o \
nainit.o naopen.o nbaread.o ngetgir.o nixgb.o nrdgb.o nskgb.o

SHELL = /bin/sh
AR = ar
CP = cp
RM = rm
CD = cd
CHMOD = chmod

LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
F77 = ifort
FORT_FLAGS = -O2 -shared-intel -mcmodel=medium
LOAD_FLAGS = -O2
INCLUDE_DIR = -I.. -I../..

${LIB} : ${OBJS}
	${AR} -r $@ $?
	${CP} -f ${LIB} ${LIB_DIR}; \
	${CD} ${LIB_DIR}; \
	${CHMOD} 755 ${LIB}

.SUFFIXES :

.SUFFIXES : .o .f

.f.o :
	${F77} ${FORT_FLAGS} -c $<

clean :
	- ${RM} -f *.o *.a core
	- ${RM} -rf rii_files

# DO NOT DELETE
