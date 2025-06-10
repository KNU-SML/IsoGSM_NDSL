LIB = sfcl_par.a

OBJS = adjsoilw.o albedo.o checkfn.o checksndph.o fixio.o fixrd_clim.o fixrd.o fixrecio.o fld2msk.o getarea.o getlatlon.o gulfofca.o i2o.o i2oini.o inds2g.o la2ga.o landtyp.o ll2gg.o ll2msk.o ll2rsm.o ll2xyr.o lonlat.o max_unfwat.o mrgfld.o mrgsnow.o nntprt.o numchar.o qcmxmn.o qcsice.o rdgrbsfc.o rdmask.o rdorog.o rgetarea.o rough.o setalbedo.o sfc.o sfcmrg.o subst.o superead.o tsfcor.o vday.o viclandck.o xsum.o xy2llr.o 

SRCS = ${OBJS:o=F}

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

.DEFAULT :
	ln -fs ../sfcl/$@ $@

${LIB} : ${SRCS} ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}

${OBJS} : ${SRCS}

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
	- ${RM} -f *.f *.o *.x *.a core Makefile
	- ${RM} -rf rii_files
