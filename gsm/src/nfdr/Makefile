EXECUTABLE = nfdr_prog.x
OBJS = @OBJS_NFDR_PROG@

MACHINE = linux
SHELL = /bin/sh
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = ifort
FORT_FLAGS = -r8 -O2 -convert big_endian -shared-intel -mcmodel=medium
LOAD_FLAGS = -r8 -O2 -convert big_endian
INCLUDE_DIR = -I../.. -I../include -I/usr/local/mpi/intel18/mvapich2-2.2/include
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
INSTALLDLIB = /rhome/yeong/test/svrtest/gsm/bin
INSTALL_DIR = /rhome/yeong/test/svrtest/gsm/bin
EXTRA_LIBS = 
EXTRA_LIBS_SFC0 = 

LIBS = ${INSTALLDLIB}/sfcl.a ${INSTALLDLIB}/share.a ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a ${EXTRA_LIBS_SFC0}

.DEFAULT :
	ln -fs ../sfcl/$@ $@

${EXECUTABLE} : ${OBJS} ${LIBS}
	if [ linux = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f .F

.F.f :
	sed '/<paramodel.h>/{G;s/$$/#undef MP/;}' $*.F >TEMP.F
	${CPP} -P ${INCLUDE_DIR} TEMP.F >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f ; rm $*.i TEMP.F

.F.o :
	sed '/<paramodel.h>/{G;s/$$/#undef MP/;}' $*.F >TEMP.F
	${CPP} -P ${INCLUDE_DIR} TEMP.F >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f ; rm $*.i TEMP.F
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $*.f

clean :
	- ${RM} -f *.f *.o ${EXECUTABLE} core

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile
	- ${RM} -rf rii_files
