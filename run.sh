#!/bin/sh 
# Run with HW#-SIM or -DATA to run code for the corresponding homework assignment
set -o errexit
set -o nounset

if [ -z ${CHEM220_IMGDIR+x} ]; then
    CHEM220_IMGDIR=$(pwd)/images
fi
if [ ! -e $CHEM220_IMGDIR ]; then
    mkdir $CHEM220_IMGDIR
fi

if [ -z ${CHEM220_DATDIR+x} ]; then
    CHEM220_DATDIR=$(pwd)/data
fi
if [ ! -e $CHEM220_DATDIR ]; then
    mkdir $CHEM220_DATDIR
fi

if [ -z ${1+x} ]; then
    echo "Requires argument HW#-SIM or HW#-DATA to run"
    exit
fi

if [ $1 == "TEST-MC" ]; then
    cd hardsphere
    gcc hardsphere.cc -lstdc++ -o hardsphere-test -DVERBOSE -DXYZOUT -DACCEPTANCE
    ./hardsphere-test -h
    ./hardsphere-test -nsteq 100 -nstmc 100 -step_size .7
fi

if [ $1 == "TEST-LJ" ]; then
    cd hardsphere
    gcc lj.cc -lstdc++ -o lj-test -DXYZOUT
    ./lj-test -nsteq 1000 -nstmc 100 -nstxyz 10 > $CHEM220_DATDIR/lj.xyz
    # vmd $CHEM220_DATDIR/lj.xyz
fi

# =============================================================================
# HW 1
# =============================================================================
if [ $1 == "HW1-SIM" ]; then
    cd hardsphere
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
    cd hardsphere
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
    cd hardsphere
    echo "Submitted one hardsphere-fourier task"
    gcc hardsphere.cc -lstdc++ -o hardsphere-fourier -DFOURIER
    ./hardsphere-fourier -nsteq 500 -nstmc 10000 -nstfourier 10 -step_size .7 \
        -maxfouriernum 100 > $CHEM220_DATDIR/fourier.csv &
fi

if [ $1 == "HW2-DATA" ]; then
    cd hardsphere
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
    cd hardsphere
    nstmc=5000
    DATA_MB=$(echo "(1000 * 1000 / 2) * ($nstmc / 50) * 10 * 2.3 / 1024 / 1024" | bc)
    echo "This run will generate something like $DATA_MB MB of data"
    read -p "Would you like to proceed (press y to contine)? " -n 1 -r
    echo    # (optional) move to a new line
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit
    fi
    echo "Submitting ten hardsphere-gr jobs"
    gcc hardsphere.cc -lstdc++ -o hardsphere-gr -DGR
    if [ $? = 0 ]; then
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size  .05 -density  .9 > $CHEM220_DATDIR/gr_9.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size   .1 -density  .8 > $CHEM220_DATDIR/gr_8.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size  .15 -density  .7 > $CHEM220_DATDIR/gr_7.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size   .2 -density  .6 > $CHEM220_DATDIR/gr_6.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size   .4 -density  .5 > $CHEM220_DATDIR/gr_5.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size  1.0 -density  .4 > $CHEM220_DATDIR/gr_4.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size  2.0 -density  .3 > $CHEM220_DATDIR/gr_3.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size  5.0 -density  .2 > $CHEM220_DATDIR/gr_2.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size 10.0 -density  .1 > $CHEM220_DATDIR/gr_1.csv  &
        ./hardsphere-gr -nsteq 500 -nstmc $nstmc -step_size 10.0 -density .02 > $CHEM220_DATDIR/gr_02.csv &
    fi
fi

if [ $1 == "HW3-DATA" ]; then
    cd hardsphere
    if [ ! -e $CHEM220_IMGDIR/hw3 ]; then
        mkdir $CHEM220_IMGDIR/hw3
    fi
    python script/gr.py $CHEM220_DATDIR/gr_*.csv -save $CHEM220_IMGDIR/hw3
fi

# =============================================================================
# HW 4
# =============================================================================
if [ $1 == "HW4-DATA" ]; then
    cd numerics
    if [ ! -e $CHEM220_IMGDIR/hw4 ]; then
        mkdir $CHEM220_IMGDIR/hw4
    fi
    python hw4.py -save $CHEM220_IMGDIR/hw4 -gaw -dmu -gaa > $CHEM220_IMGDIR/hw4/coeff.txt
fi

# =============================================================================
# HW 5
# =============================================================================
if [ $1 == "HW4-DATA" ]; then
    cd numerics
    if [ ! -e $CHEM220_IMGDIR/hw5 ]; then
        mkdir $CHEM220_IMGDIR/hw5
    fi
    python hw5.py -save $CHEM220_IMGDIR/hw5
fi

# =============================================================================
# HW 6
# =============================================================================
if [ $1 == "HW6-SIM" ]; then
    cd hardsphere
    gcc lj.cc -lstdc++ -o lj-energy -DENERGY
    gcc lj.cc -lstdc++ -o lj-velocity -DVELOCITY
    gcc lj.cc -lstdc++ -o lj-xyz -DXYZOUT
    gcc lj.cc -lstdc++ -o lj-gr -DGR
    echo "Assuming 8 cores for parallelization."
    echo "If fewer than 8 cores, expect significantly slower runtime."
    echo "Running energy simulations..."
    echo "Running density=.3 simulations..."
    ./lj-energy   -N_linear 10 -density .3             -nsteq 2000 -nstmd 10000 > $CHEM220_DATDIR/lj-e3.txt &
    ./lj-gr       -N_linear 10 -density .3             -nsteq 2000 -nstmd  3000 > $CHEM220_DATDIR/lj-gr_3.csv &
    ./lj-velocity -N_linear 10 -density .3 -seed 90210 -nsteq 2000 -nstmd 10000 > $CHEM220_DATDIR/lj_1-v3.txt &
    ./lj-velocity -N_linear 10 -density .3 -seed 90211 -nsteq 2000 -nstmd 10000 > $CHEM220_DATDIR/lj_2-v3.txt &
    ./lj-velocity -N_linear 10 -density .3 -seed 90212 -nsteq 2000 -nstmd 10000 > $CHEM220_DATDIR/lj_3-v3.txt &
    ./lj-xyz      -N_linear 10 -density .3 -seed 90210 -nsteq 2000 -nstmd 10000 -nstxyz 100 > $CHEM220_DATDIR/lj_1-x3.xyz   &
    ./lj-xyz      -N_linear 10 -density .3 -seed 90211 -nsteq 2000 -nstmd 10000 -nstxyz 100 > $CHEM220_DATDIR/lj_2-x3.xyz   &
    ./lj-xyz      -N_linear 10 -density .3 -seed 90212 -nsteq 2000 -nstmd 10000 -nstxyz 100 > $CHEM220_DATDIR/lj_3-x3.xyz   

    # echo "Running density=.6 simulations..."
    # ./lj-energy   -N_linear 10 -density .6             -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj-e6.txt &
    # ./lj-gr       -N_linear 10 -density .6             -nsteq 1000 -nstmd 3000 > $CHEM220_DATDIR/lj-gr_6.csv &
    # ./lj-velocity -N_linear 10 -density .6 -seed 90210 -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj_1-v6.txt &
    # ./lj-velocity -N_linear 10 -density .6 -seed 90211 -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj_2-v6.txt &
    # ./lj-velocity -N_linear 10 -density .6 -seed 90212 -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj_3-v6.txt &
    # ./lj-xyz      -N_linear 10 -density .6 -seed 90214 -nsteq 1000 -nstmd 5000 -nstxyz 100 > $CHEM220_DATDIR/lj_1-x6.xyz   &
    # ./lj-xyz      -N_linear 10 -density .6 -seed 90215 -nsteq 1000 -nstmd 5000 -nstxyz 100 > $CHEM220_DATDIR/lj_2-x6.xyz   &
    # ./lj-xyz      -N_linear 10 -density .6 -seed 90216 -nsteq 1000 -nstmd 5000 -nstxyz 100 > $CHEM220_DATDIR/lj_3-x6.xyz
    # 
    # echo "Running density=.8 simulations"
    # ./lj-energy   -N_linear 10 -density .8             -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj-e8.txt &
    # ./lj-gr       -N_linear 10 -density .8             -nsteq 1000 -nstmd 3000 > $CHEM220_DATDIR/lj-gr_8.csv     &
    # ./lj-velocity -N_linear 10 -density .8 -seed 90210 -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj_1-v8.txt &
    # ./lj-velocity -N_linear 10 -density .8 -seed 90211 -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj_2-v8.txt &
    # ./lj-velocity -N_linear 10 -density .8 -seed 90212 -nsteq 1000 -nstmd 5000 > $CHEM220_DATDIR/lj_3-v8.txt &
    # ./lj-xyz      -N_linear 10 -density .8 -seed 90210 -nsteq 1000 -nstmd 5000 -nstxyz 100 > $CHEM220_DATDIR/lj_1-x8.xyz   &
    # ./lj-xyz      -N_linear 10 -density .8 -seed 90211 -nsteq 1000 -nstmd 5000 -nstxyz 100 > $CHEM220_DATDIR/lj_2-x8.xyz   &
    # ./lj-xyz      -N_linear 10 -density .8 -seed 90212 -nsteq 1000 -nstmd 5000 -nstxyz 100 > $CHEM220_DATDIR/lj_3-x8.xyz   

    rm lj-gr lj-velocity lj-energy lj-xyz
fi

if [ $1 == "HW6-DATA" ]; then
    cd hardsphere
    if [ ! -e $CHEM220_IMGDIR/hw6 ]; then
        mkdir $CHEM220_IMGDIR/hw6
    fi
    if [ ! -e $CHEM220_IMGDIR/hw6/check ]; then
        mkdir $CHEM220_IMGDIR/hw6/check
    fi
    echo "Running analysis for density=.3..."
    python script/hw6.py -density .3 \
        -velfile $CHEM220_DATDIR/lj_{1..3}-v3.txt \
        -xyzfile $CHEM220_DATDIR/lj_{1..3}-x3.xyz \
        -enerfile $CHEM220_DATDIR/lj-e3.txt \
        -savename d3 -save $CHEM220_IMGDIR/hw6/check \
        > $CHEM220_IMGDIR/hw6/d3_diff.txt
    echo "Running analysis for density=.6..."
    python script/hw6.py -density .6 \
        -velfile $CHEM220_DATDIR/lj_{1..3}-v6.txt \
        -xyzfile $CHEM220_DATDIR/lj_{1..3}-x6.xyz \
        -enerfile $CHEM220_DATDIR/lj-e6.txt \
        -savename d6 -save $CHEM220_IMGDIR/hw6/check \
        > $CHEM220_IMGDIR/hw6/d6_diff.txt
    echo "Running analysis for density=.8..."
    python script/hw6.py -density .8 \
        -velfile $CHEM220_DATDIR/lj_{1..3}-v8.txt \
        -xyzfile $CHEM220_DATDIR/lj_{1..3}-x8.xyz \
        -enerfile $CHEM220_DATDIR/lj-e8.txt \
        -savename d8 -save $CHEM220_IMGDIR/hw6/check \
        > $CHEM220_IMGDIR/hw6/d8_diff.txt
    python script/gr.py $CHEM220_DATDIR/lj-gr_*.csv -save $CHEM220_IMGDIR/hw6 \
        -N 1000 -T 60
fi


# =============================================================================
# HW 7
# =============================================================================

if [ $1 == "HW7-SIM" ]; then
    for i in {1..200}; do
        for str in 3 5 10; do
            echo "Running 200 STR=$str jobs at all field oscillation rates"
            python pymd/ljmd.py -density .5 -N_linear 5 -field 12 $str -nsteq 1000 \
                -nstxyz 10 -nstmd 1000 -seed $i \
                -posfile $CHEM220_DATDIR/field_12_"$str"_$i.xyz -enerfile $CHEM220_DATDIR/field_12_"$str"_$i.ener &
            python pymd/ljmd.py -density .5 -N_linear 5 -field 24 $str -nsteq 1000 \
                -nstxyz 10 -nstmd 1000 -seed $i \
                -posfile $CHEM220_DATDIR/field_24_"$str"_$i.xyz -enerfile $CHEM220_DATDIR/field_24_"$str"_$i.ener &
            python pymd/ljmd.py -density .5 -N_linear 5 -field 36 $str -nsteq 1000 \
                -nstxyz 10 -nstmd 1000 -seed $i \
                -posfile $CHEM220_DATDIR/field_36_"$str"_$i.xyz -enerfile $CHEM220_DATDIR/field_36_"$str"_$i.ener
        done
    done
    for i in {1..3}; do
        python pymd/ljmd.py -density .5 -N_linear 5 -field 0 0 -nsteq 1000 \
            -nstxyz 10 -nstmd 50000 -seed 9021$i \
            -posfile $CHEM220_DATDIR/field_eq_$i.xyz -enerfile $CHEM220_DATDIR/field_eq_$i.ener &
    done
fi
if [ $1 == "HW7-DATA" ]; then
    if [ ! -e "$CHEM220_IMGDIR/hw7" ]; then
        mkdir $CHEM220_IMGDIR/hw7
    fi
    python hardsphere/script/St.py -savedir $CHEM220_IMGDIR/hw7 -savename St \
        -eqfile $CHEM220_DATDIR/field_eq*.xyz
    for str in 3 5 10; do
        python hardsphere/script/St.py -savedir $CHEM220_IMGDIR/hw7 -savename Xt_$str \
            $CHEM220_DATDIR/field_*_"$str"_{1..200}.xyz -invscale $(echo "$str * .6666667" | bc)
    done
    python numerics/hw7.py -save $CHEM220_IMGDIR/hw7
fi
