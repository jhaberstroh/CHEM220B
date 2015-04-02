#!/bin/bash
set -o errexit
set -o nounset

if [ -z ${CHEM220_IMGDIR+x} ]; then
    CHEM220_IMGDIR=images
fi
if [ ! -e $CHEM220_IMGDIR ]; then
    mkdir $CHEM220_IMGDIR
fi

if [ ! -e $CHEM220_IMGDIR/hw5 ]; then
    mkdir $CHEM220_IMGDIR/hw5
fi
python hw5.py -save $CHEM220_IMGDIR/hw5
