LIB = mpi.a
OBJS = dimset.o equdis.o equdiv.o equlw12.o mainset.o mpabort.o mpbcastc.o mpbcasti.o mpbcastl.o mpbcastr.o mpdimset.o mpfine.o mpgetken.o mpgetlat.o mpgetspd.o mpgf2p.o mpgfpk2fk.o mpgp2f.o mpgpfk2fpk.o mpinit.o mplatall.o mpn2nn.o mpnk2nn.o mpnk2nx.o mpnl2ny.o mpnn2n.o mpnn2nk.o mpnx2nk.o mpnx2x.o mpny2nl.o mpsf2p.o mpsp2f.o mpsynall.o mpsyncol.o mpsynrow.o mpx2nx.o preduceg.o spcshfli.o spcshflo.o 

MACHINE = linux
SHELL = /bin/sh
AR = ar
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = mpif90
FORT_FLAGS = -r8 -O2 -convert big_endian -I/usr/local/mpi/intel18/mvapich2-2.2/include
LOAD_FLAGS = -r8 -O2 -convert big_endian -L/opt/gm/lib
INCLUDE_DIR = -I../.. -I../include -I/usr/local/mpi/intel18/mvapich2-2.2/include
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
INSTALLDLIB = /rhome/yeong/test/svrtest/gsm/bin

#EXP_PARENT = @EXP_PARENT@
#.DEFAULT :
#	ln -fs ${EXP_PARENT}/src/mpi/$@ $@

${LIB} : ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f .F

.F.f :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i

.F.o :
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed '/^ *$$/d' $*.i >$*.f ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f

.f.o :
	${F77} ${FORT_FLAGS} -c $*.f

clean :
	- ${RM} -f *.f *.o *.a core

dist_clean :
	- ${RM} -f *.f *.o *.x core Makefile
	- ${RM} -rf rii_files
