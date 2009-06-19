export CURRENT_DIR=`pwd`
cd ~/testx10/fmm/src
find . -name *.x10 -print0  | xargs -r0 x10c -v
cd $CURRENT_DIR


