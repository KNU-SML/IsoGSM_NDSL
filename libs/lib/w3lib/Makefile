LIB = w3lib.a
OBJS = aea.o assign.o baread.o datimx.o \
gbyte.o gbytes.o getgbss.o getgir.o iw3jdn.o iw3pds.o ixgb.o pdsens.o pdseup.o \
rdgb.o sbyte.o sbytes.o skgb.o swap32.o unpindx.o w3fi01.o w3fi04.o w3fi63.o \
w3fi83.o w3fp11.o w3fs26.o w3fi59.o w3fi72.o w3fi75.o \
w3fi68.o w3fi73.o w3fi76.o w3fi58.o w3fi71.o w3fi74.o w3fi82.o \
w3tagb.o xmovex.o xstore.o baopen.o mova2i.o

COBJS = cdate.o cdate2.o

SHELL = /bin/sh
AR = ar
CP = cp
CPP = gcc -E
RM = rm
CD = cd
CHMOD = chmod

LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
F77 = ifort
FORT_FLAGS = -O2 -shared-intel -mcmodel=medium
LOAD_FLAGS = -O2
INCLUDE_DIR = -I.. -I../..

${LIB} : ${OBJS} ${COBJS}
	${AR} -r $@ $?
	${CP} -f ${LIB} ${LIB_DIR}; \
	${CD} ${LIB_DIR}; \
	${CHMOD} 755 ${LIB}

.SUFFIXES :

.SUFFIXES : .o .f .F .c

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $<

.c.x :
	${CC} -o $@ ${INCLUDE_DIR} $< -lm

clean :
	- ${RM} -rf *.o *.a core rii_files
	- ${RM} -rf rii_files

# DO NOT DELETE
