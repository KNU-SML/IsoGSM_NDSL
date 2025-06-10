LIB = bufr.a
OBJS = abt.o adn30.o bfrini.o bufr.o chekstab.o chrtrn.o chrtrna.o closbf.o \
closmg.o conwin.o copybf.o copymg.o cpbfdx.o datebf.o drstpl.o dxinit.o dxmini.o \
elemdx.o getwin.o idn30.o ifbget.o ifxy.o inctab.o invcon.o invtag.o invwin.o \
ipkm.o ireadmg.o irev.o iupb.o iupm.o jstify.o lmsg.o lstjpb.o lstrpc.o lstrps.o \
makestab.o msgini.o msgupd.o msgwrt.o mvb.o nemock.o nemtab.o nemtba.o nemtbb.o \
nemtbd.o nenuck.o newwin.o numbck.o numtab.o nvnwin.o nxtwin.o openbf.o openmb.o \
openmg.o pad.o parseq.o parusr.o parutg.o pkb.o pkc.o pktdd.o posapn.o rcstpl.o \
rdbfdx.o rdmgsb.o rdtree.o rdusdx.o readdx.o readmg.o readsb.o rewnbf.o rjust.o \
rsvfvm.o seqsdx.o status.o strcln.o string.o strnum.o strsuc.o tabent.o tabsub.o \
trybump.o ufbcnt.o ufbcpy.o ufbcup.o ufbdmp.o ufbevn.o ufbget.o ufbint.o ufbinx.o \
ufbqcd.o ufbqcp.o ufbrep.o ufbrew.o ufbrp.o ufbrw.o ufbtab.o upb.o upbb.o upc.o \
uptdd.o usrtpl.o vals.o wrdlen.o writdx.o writsa.o writsb.o wrtree.o wtstat.o

MACHINE = linux
SHELL = /bin/sh
AR = ar
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = cc
CPP = gcc -E
F77 = ifort
FORT_FLAGS = 
LOAD_FLAGS = 
INCLUDE_DIR = -I.. -I../..
LIB_DIR = /rhome/yeong/test/svrtest/libs/lib
EXRTA_LIBS = 

${LIB} : ${OBJS}
	${AR} -r $@ $?
	${CP} $@ ${LIB_DIR} 
	${CD} ${LIB_DIR} ; \
	${CHMOD} 755 ${LIB}

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
	- ${RM} -f *.f *.o *.a core Makefile
