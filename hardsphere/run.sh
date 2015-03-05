#!/bin/sh 
# Run with HW#-SIM or -DATA to run code for the corresponding homework assignment
set -o errexit

if [ -z ${CHEM220_IMGDIR+x} ]; then
    CHEM220_IMGDIR=images
fi
if [ ! -e $CHEM220_IMGDIR ]; then
    mkdir $CHEM220_IMGDIR
fi

if [ -z ${CHEM220_DATDIR+x} ]; then
    CHEM220_DATDIR=data
fi
if [ ! -e $CHEM220_DATDIR ]; then
    mkdir $CHEM220_DATDIR
fi

if [ -z ${1+x} ]; then
    echo "Requires argument HW#-SIM or HW#-DATA to run"
    exit
fi

if [ $1 == "TEST-MC" ]; then
    gcc hardsphere.cc -lstdc++ -o hardsphere-test -DVERBOSE -DXYZOUT -DACCEPTANCE
    if [ $? == 0 ]; then
        ./hardsphere-test -h
        ./hardsphere-test -nsteq 100 -nstmc 100 -step_size .7
    fi
fi

if [ $1 == "TEST-LJ" ]; then
    gcc lj.cc -lstdc++ -o lj-test -DXYZOUT
    if [ $? == 0 ]; then
        ./lj-test -nsteq 10000 -nstmc 100 -dt .1
    fi
fi

# =============================================================================
# HW 1
# =============================================================================
if [ $1 == "HW1-SIM" ]; then
    echo "Submitted one hardsphere-xyz, one hardsphere-sol and three" \
   " hardsphere-lg tasks"
    gcc hardsphere.cc -lstdc++ -o hardsphere-xyz -DXYZOUT
    ./hardsphere-xyz -nsteq 400 -step_size .7 > $CHEM220_DATDIR/atoms.xyz & 
    
    gcc hardsphere.cc -lstdc++ -o hardsphere-sol -DSMALLSPHERE
    ./hardsphere-sol -nsteq 500 -step_size .7 > $CHEM220_DATDIR/small.csv &
    
    gcc hardsphere.cc -lstdc++ -o hardsphere-lg -DLARGESPHERE
    ./hardsphere-lg -nsteq 500 -step_size .7 -probe 1.0 > $CHEM220_DATDIR/large-2.0.csv &
    ./hardsphere-lg -nsteq 500 -step_size .7 -probe 1.5 > $CHEM220_DATDIR/large-3.0.csv &
    ./hardsphere-lg -nsteq 500 -step_size .7 -probe 2.0 > $CHEM220_DATDIR/large-4.0.csv &
fi

if [ $1 == "HW1-DATA" ]; then
    echo "Saving images to $CHEM220_IMGDIR/hw1"
    if [ ! -e $CHEM220_IMGDIR/hw1 ]; then
        mkdir $CHEM220_IMGDIR/hw1
    fi
    python script/plotxyz.py -dir $CHEM220_DATDIR -save $CHEM220_IMGDIR/hw1
    python script/hist.py -dir $CHEM220_DATDIR -save $CHEM220_IMGDIR/hw1
fi

# =============================================================================
# HW 2
# =============================================================================
if [ $1 == "HW2-SIM" ]; then
    echo "Submitted one hardsphere-fourier task"
    gcc hardsphere.cc -lstdc++ -o hardsphere-fourier -DFOURIER
    ./hardsphere-fourier -nsteq 500 -nstmc 10000 -nstfourier 10 -step_size .7 \
        -maxfouriernum 100 > $CHEM220_DATDIR/fourier.csv &
fi

if [ $1 == "HW2-DATA" ]; then
    echo "Saving images to $CHEM220_IMGDIR/hw2"
    if [ ! -e $CHEM220_IMGDIR/hw2 ]; then
        mkdir $CHEM220_IMGDIR/hw2
    fi
    python script/fourier.py -file $CHEM220_DATDIR/fourier.csv -save $CHEM220_IMGDIR/hw2
fi

# =============================================================================
# HW 3
# =============================================================================
if [ $1 == "HW3-SIM" ]; then
    gcc hardsphere.cc -lstdc++ -o hardsphere-gr -DGR
    if [ $? = 0 ]; then
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size  .05 -density  .9 > $CHEM220_DATDIR/gr_9.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size   .1 -density  .8 > $CHEM220_DATDIR/gr_8.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size  .15 -density  .7 > $CHEM220_DATDIR/gr_7.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size   .2 -density  .6 > $CHEM220_DATDIR/gr_6.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size   .4 -density  .5 > $CHEM220_DATDIR/gr_5.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size  1.0 -density  .4 > $CHEM220_DATDIR/gr_4.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size  2.0 -density  .3 > $CHEM220_DATDIR/gr_3.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size  5.0 -density  .2 > $CHEM220_DATDIR/gr_2.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size 10.0 -density  .1 > $CHEM220_DATDIR/gr_1.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc 5000 -step_size 10.0 -density .02 > $CHEM220_DATDIR/gr_02.csv &
    fi
fi

if [ $1 == "HW3-DATA" ]; then
    python script/gr.py $CHEM220_DATDIR/gr_*.csv
fi


# =============================================================================
# HW 5
# =============================================================================
#if [ $1 == "HW5-SIM" ]; then
#    gcc lj.cc -lstdc++ -o lj -DVERBOSE
#    if [ $? = 0 ]; then
#        ./lj
#    fi
#fi
