EXECUTABLES = grmap.x grmean.x grmean_daily.x grsplit.x grtran.x ibm2ie3.x cray2ie3.x grvar.x \
cosread.x grsmooth.x
OBJS = ${EXECUTABLES:.x=.o}

SHELL = /bin/sh
CP = cp
RM = rm
CD = cd
CHMOD = chmod
CC = cc
CPP = gcc -E

LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
ETC_DIR = /rhome/yeong/test/svrtest/libs/etc
F77 = ifort
FORT_FLAGS = -r8 -O2 -convert big_endian -shared-intel -mcmodel=medium
LOAD_FLAGS = -r8 -O2 -convert big_endian
INCLUDE_DIR = -I.. -I../..

LIB = ${LIB_DIR}/w3lib.a ${LIB_DIR}/libncaru.a ${LIB_DIR}/modelib.a

execs : ${EXECUTABLES}
	${CP} $? ${ETC_DIR}
	${CD} ${ETC_DIR}; \
	${CHMOD} 755 $?

.SUFFIXES : .x .o .f .F 

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.F.x :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f
	${F77} ${LOAD_FLAGS} -o $@ $*.o ${LIB}

.o.x :
	${F77} ${LOAD_FLAGS} -o $@ $< ${LIB}

.f.x :
	${F77} ${FORT_FLAGS} -c $<
	${F77} ${LOAD_FLAGS} -o $@ $*.o ${LIB}

.f.o :
	${F77} ${FORT_FLAGS} -c $<

clean :
	- ${RM} -rf *.f *.o *.x core rii_files
	- ${CD} ${ETC_DIR}; \
	${RM} -f ${EXECUTABLES}

# DO NOT DELETE
