EXECUTABLE = sfc0.x
OBJS = frh2o.o getslc.o ncp1tonoa1.o ncp2tonoa1.o ncp3tonoa1.o noa1toosu1.o noa1toosu2.o noa1tovic1.o noa1tovic2.o nodepara.o osu1tonoa1.o osu1toosu2.o osu1tovic1.o osu1tovic2.o osu2tonoa1.o osu2toosu1.o osu2tovic1.o osu2tovic2.o sfc0.o sfcfcnv.o vic1ini.o vic1tonoa1.o vic1toosu1.o vic1toosu2.o vic1tovic2.o vic2ini.o vic2tonoa1.o vic2toosu1.o vic2toosu2.o vic2tovic1.o 

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

LIBS = ${INSTALLDLIB}/sfcl.a ${INSTALLDLIB}/share.a ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a  ${EXTRA_LIBS_SFC0}

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/sfc/$@ $@

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
