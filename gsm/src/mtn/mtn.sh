#!/bin/sh
#
# @ job_type = serial
# @ output = mtn.out.$(jobid).out
# @ error = mtn.out.$(jobid).err
# @ queue
# @ wall_clock_limit = 00:10:00
# @ class = dev
#
set -x
# creates resolution-specific terrain files
# from 10' navy or USGS 8 min or 4 min data
# argument list:
#   1: output sea-and-land mask file
#   2: output gridded orography file
#   3: output mountain variance file
#   4: output spectral orography file
#   5: output unfiltered gridded orography file
#   6: number of longitudes in Gaussian grid
#   7: number of latitudes in Gaussian grid
#   8: spectral truncation
#   9: rhomboidal flag (=0 for triangular truncation)
#
# Unit 52, orogrd is the spectrally truncated grid orography on 
#          Gaussian grid to be used by the sfc codes.
#
FTNID=fort.
CON_DIR=/rhome/yeong/test/svrtest/libs/con
INSTALL_DIR=/rhome/yeong/test/svrtest/gsm/bin
BIN_DIR=/rhome/yeong/test/svrtest/gsm/bin
MODEL_DIR=/rhome/yeong/test/svrtest/gsm
#
jcap=62
levs=28
lonf=192
latg=94
echo $jcap $levs $lonf $latg

cd $MODEL_DIR/src/mtn

mtnres=8

avg=top${mtnres}m_avg.20i4.asc
var=top${mtnres}m_var.20i4.asc
max=top${mtnres}m_max.20i4.asc
slm=top${mtnres}m_slm.80i1.asc

for file in $avg $var $max $slm
do
	ln -fs $CON_DIR/$file $file
done

rm -f ${FTNID}[0-9]* 2>/dev/null
ln -fs $avg                         ${FTNID}11
ln -fs $var                         ${FTNID}12
ln -fs $max                         ${FTNID}13
ln -fs $slm                         ${FTNID}14
ln -fs $INSTALL_DIR/slmsk           ${FTNID}51
ln -fs $INSTALL_DIR/orogrd          ${FTNID}52
ln -fs $INSTALL_DIR/mtnvar          ${FTNID}53
ln -fs $INSTALL_DIR/orospc          ${FTNID}54
ln -fs $INSTALL_DIR/oroggd          ${FTNID}55
ln -fs $INSTALL_DIR/orospc.smth     ${FTNID}61
ln -fs $INSTALL_DIR/orogrd.smth     ${FTNID}62
ln -fs $INSTALL_DIR/oroviw.ieee     ${FTNID}63

echo "$lonf $latg $jcap 0" >mtn.parm

if [ linux = dec ] ; then
	datasize=1048576
	stacksize=32768
elif [ linux = mac ] ; then
	datasize=unlimited
	stacksize=65536
else
	datasize=unlimited
	stacksize=unlimited
fi
if [ linux != nec ] ; then
	ulimit -d $datasize
	ulimit -s $stacksize
fi

$BIN_DIR/mtn.x <mtn.parm >mtn.out 2>&1
if [ $? -ne 0 ] ; then
	echo `date` - mtn.x execution failed
	rm -f ${FTNID}[0-9]*
	exit 8
fi
rm -f ${FTNID}[0-9]*
echo "mtn files created on `date`" >$INSTALL_DIR/.mtn
