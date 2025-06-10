#
OBJS = sigadd2.o sigadd.o sigdif.o 
EXECUTABLES = ${OBJS:.o=.x}

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
INSTALL_DIR=/rhome/yeong/test/svrtest/gsm/bin
INCLUDE_DIR = -I../.. -I../include -I/usr/local/mpi/intel18/mvapich2-2.2/include

execs :	 ${EXECUTABLES}
	${CP} $? ${INSTALL_DIR}; \
	${CD} ${INSTALL_DIR}; \
	${CHMOD} 755 $?

.SUFFIXES : .x .o .f .F

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $<

.F.x :
	if [ linux = mac ] ; then \
		ranlib -s ${LIB} ; \
	fi
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f
	${F77} ${LOAD_FLAGS} -o $@ $*.o

.o.x :
	if [ linux = mac ] ; then \
		ranlib -s ${LIB} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ $<

.f.x :
	if [ linux = mac ] ; then \
		ranlib -s ${LIB} ; \
	fi
	${F77} ${FORT_FLAGS} -c $<
	${F77} ${LOAD_FLAGS} -o $@ $*.o

clean :
	- ${RM} -f *.f *.o *.x *.a core ${EXECUTABLES}
	- ${CD} ${INSTALL_DIR}; ${RM} -rf ${EXECUTABLES}

# DO NOT DELETE
