#!/bin/sh
#
set -a
CON_DIR=/rhome/yeong/test/svrtest/libs/con
F77=ifort
FTNID=fort.
FORT_OPTS='-r8 -O2 -convert big_endian -shared-intel -mcmodel=medium'
LOAD_OPTS='-r8 -O2 -convert big_endian'
cat >tune.f <<'EOF'
      parameter(nbin=100,im=2,jm=4,km=3,lm=2)
      dimension rhfd(nbin*im*jm*km*lm)
      dimension rfhffd(nbin*im*jm*km*lm)
      dimension kpts(im*jm*km*lm)
      dimension icdays(15),idate(4)
      read(10,100) nbdayi,icdays
      write(50) nbdayi,icdays
      print *,'nbdayi=',nbdayi,' icdays=',icdays
      read(10,200) fhour,idate
      print *,'fhour=',fhour,' idate',idate
      write(50) fhour,idate
      do n=1,nbdayi
C
        print *,'n=',n
        read(10,300) rhfd
C       print *,'rhfd read'
        write(50) rhfd
C
        read(10,300) rfhffd
C       print *,'rfhffd read'
        write(50) rfhffd
C
        read(10,400) kpts
C       print *,'kprs read'
        write(50) kpts
      enddo
100   format(8i10)
200   format(f13.6,4i4)
300   format(6e13.6)
400   format(6i13)
      stop
      end
EOF
#echo "$F77 $FORT_OPTS -o tune.x tune.f"
$F77 $FORT_OPTS -c tune.f || exit 8
$F77 $LOAD_OPTS -o tune.x tune.o || exit 8

for tune in tune1.asc tune.t42l18.amip.ewmrg.asc tune.t62l28.reanl.ewmrg.asc
do
	INPFILE=$tune
	OUTFILE=${CON_DIR}/`echo $tune | sed 's/\.[^\.]*$//g'`
#	echo ' '
#	echo "Converting cloud tuning file: $tune"
	if [ -s $tune.Z ] ; then
		uncompress $tune.Z || exit 8
	fi
	if [ -s $tune.gz ] ; then
		gunzip $tune.gz || exit 8
	fi
	rm -f ${FTNID}[0-9]* 2>/dev/null 
	ln -s $INPFILE		 ${FTNID}10
	ln -s $OUTFILE     ${FTNID}50
	if [ linux = es ] ; then
		dirhere=`pwd`
		cat >~/tune_gb.sh <<EOF
		#!/bin/sh
		cd $dirhere
		./tune.x || exit 8
EOF
		chmod 755 ~/tune_gb.sh
		rexec gbnode000 -l  ./tune_gb.sh 2>tune.out
	else
		./tune.x >tune.out 2>tune.err  || exit 8
	fi
done
#
rm -f tune.o tune.x tune.f 2>/dev/null
rm -f ${FTNID}[0-9]* 2>/dev/null 
exit 0
