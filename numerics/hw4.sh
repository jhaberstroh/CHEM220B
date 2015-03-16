#!/bin/bash
set -o errexit
set -o nounset

if [ -z ${CHEM220_IMGDIR+x} ]; then
    CHEM220_IMGDIR=images
fi
if [ ! -e $CHEM220_IMGDIR ]; then
    mkdir $CHEM220_IMGDIR
fi

if [ ! -e $CHEM220_IMGDIR/hw4 ]; then
    mkdir $CHEM220_IMGDIR/hw4
fi
python hw4.py -save $CHEM220_IMGDIR/hw4 -gaw -dmu -gaa > $CHEM220_IMGDIR/hw4/coeff.txt
