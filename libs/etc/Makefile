#
MACHINE = linux
SHELL = /bin/sh
CP = cp
CD = cd
RM = rm
MAKE = make

ETC_DIR = /rhome/yeong/test/svrtest/libs/etc

all : 
	${CD} util_intel ; ${MAKE} || exit 8
	${CD} utils ; ${MAKE} || exit 8

clean :  
	${CD} util_intel ; ${MAKE} clean
	${CD} utils ; ${MAKE} clean
	- ${RM} -rf rii_files
