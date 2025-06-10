EXECUTABLES = bufrdte.x jwpl1.x msgjb.x 
CEXECUTABLES = 

SHELL = /bin/sh
CC = cc
CPP = gcc -E
RM = rm
CD = cd
CP = cp
CHMOD = chmod

F77 = ifort
FORT_FLAGS = -r8 -O2 -convert big_endian -shared-intel -mcmodel=medium
LOAD_FLAGS = -r8 -O2 -convert big_endian
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
ETC_DIR = /rhome/yeong/test/svrtest/libs/etc
INCLUDE_DIR = -I.. -I../..

LIB = ${LIB_DIR}/bufr.a

execs : ${EXECUTABLES} ${CEXECUTABLES}
	if [ ! -s ${ETC_DIR} ] ; then \
		mkdir -p ${ETC_DIR}; \
	fi
	${CP} $? ${ETC_DIR}; \
	${CD} ${ETC_DIR}; \
	${CHMOD} 755 $?

cexecs : ${CEXECUTABLES}
	if [ ! -s ${ETC_DIR} ] ; then \
		mkdir -p ${ETC_DIR}; \
	fi
	${CP} $? ${ETC_DIR}; \
	${CD} ${ETC_DIR}; \
	${CHMOD} 755 $?

.SUFFIXES : .x .o .f .F .c

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
	${F77} ${LOAD_FLAGS} -o $@ $< ${LIB}

.o.x :
	${F77} ${LOAD_FLAGS} -o $@ $< ${LIB}

.f.x :
	${F77} ${FORT_FLAGS} -c $<
	${F77} ${LOAD_FLAGS} -o $@ $*.o ${LIB}

.f.o :
	${F77} ${FORT_FLAGS} -c $<

.c :
	${CC} -o $@ $<

dist_clean :
	- ${CD} ${ETC_DIR}; ${RM} -rf *.f *.o *.x core Makefile rii_files
	- ${RM} -f *.o *.x *.a core ${CEXECUTABLES} Makefile

clean :
	- ${RM} -f *.o *.x *.a core ${CEXECUTABLES} Makefile
# DO NOT DELETE
