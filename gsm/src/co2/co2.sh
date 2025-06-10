#!/bin/sh
#
# @ job_type = serial
# @ output = co2.out.$(jobid).out
# @ error = co2.out.$(jobid).err
# @ queue
#

set -x

CON_DIR=/rhome/yeong/test/svrtest/libs/con
BIN_DIR=/rhome/yeong/test/svrtest/gsm/bin
INSTALL_DIR=/rhome/yeong/test/svrtest/gsm/bin
FTNID=fort.
MODEL_DIR=/rhome/yeong/test/svrtest/gsm

cd $MODEL_DIR/src/co2

co2ppmv=348

rm -f ${FTNID}[0-9]* 2>/dev/null
ln -fs $INSTALL_DIR/co2con  ${FTNID}66
ln -fs co2output  ${FTNID}16
ln -fs $CON_DIR/cns_${co2ppmv}_490850.asc  ${FTNID}61
ln -fs $CON_DIR/cns_${co2ppmv}_490670.asc  ${FTNID}62
ln -fs $CON_DIR/cns_${co2ppmv}_670850.asc  ${FTNID}63

INPUT=co2.parm
echo $co2ppmv. >co2.parm

if [ linux = es ] ; then
dirhere=`pwd`
cat >~/co2_gb.sh <<EOF
#!/bin/sh
cd $dirhere
./co2.x <co2.parm >co2.out 2>&1
EOF
	chmod 755 ~/co2_gb.sh
	rexec gbnode000 -l  ./co2_gb.sh
else
	$BIN_DIR/co2.x <co2.parm >co2.out 2>&1
fi

if [ $? -ne 0 ] ; then
	echo `date` - co2 execution failed
	exit 8
fi
rm -f ${FTNID}[0-9]*
echo "co2con created on `date`" >$INSTALL_DIR/.co2
