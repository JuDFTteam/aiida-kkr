#!/usr/bin/env bash

# clone jukkr repository
git clone --depth 1 -b develop --single-branch gitlab@iffgit.fz-juelich.de:kkr/jukkr.git
cd jukkr/

# now codes are build using gfortran and in serial, code executables will be places in the jukkr directory

# build voronoi code
# cannot build with gfortran at the moment, needs fixing on jukkr side
#./install.py --program=voronoi --compiler=ifort --parallelization=serial
./install.py --program=voronoi --compiler=gfortran --parallelization=serial
cd build/ && make -j4 && cp voronoi.exe ../
cd ..

# build kkrhost code
#./install.py --program=kkrhost --compiler=ifort --parallelization=serial
./install.py --program=kkrhost --compiler=gfortran --parallelization=serial -d
cd build/ && make -j4 && cp kkr.x ../
cd ..

# build kkrimp code
# cannot build with gfortran at the moment, needs fixing on jukkr side
#./install.py --program=kkrimp --compiler=ifort --parallelization=serial
./install.py --program=kkrimp --compiler=gfortran --parallelization=serial
cd build/ && make -j4 && cp kkrflex.exe ../
cd ..
