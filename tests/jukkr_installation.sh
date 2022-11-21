#!/usr/bin/env bash

usage(){
  echo "$0 usage:" && grep " .)\ #" $0;
  exit 0;
}

install_jukkr=""
while getopts fh option; do
  case $option in
    f) # Force installation of jukkr codes
       install_jukkr="yes" && echo "Found -f flag: download and install jukkr code" ;;
    h) # Display help
       usage
  esac
done

if [[ -z "$install_jukkr" ]]; then
    # fake installations
    echo "skip compiling code and use caching instead (still create fake files)"
    mkdir -p jukkr/build_new_kkrhost/
    touch jukkr/build_new_kkrhost/kkr.x && chmod +x jukkr/build_new_kkrhost/kkr.x
    touch jukkr/kkr.x && chmod +x jukkr/kkr.x
    touch jukkr/voronoi.exe && chmod +x jukkr/voronoi.exe
    touch jukkr/kkrflex.exe && chmod +x jukkr/kkrflex.exe
else
    # clone jukkr repository
    git clone -b develop gitlab@iffgit.fz-juelich.de:kkr/jukkr.git
    #git clone --depth 1 -b BdG --single-branch gitlab@iffgit.fz-juelich.de:kkr/jukkr.git
    mkdir -p jukkr/
    cd jukkr/

    # now codes are build using gfortran and in serial
    # code executables will be placed inside the jukkr directory
    # build voronoi code
    echo "build voronoi"
    ./install.py --program=voronoi --compiler=gfortran --parallelization=serial
    cd build/ && make -j4 && cp voronoi.exe ../
    cd ..
    echo ""

    # build kkrhost code
    echo "build kkrhost"
    ./install.py --program=kkrhost --compiler=gfortran --parallelization=serial
    cd build/ && make -j4 && cp kkr.x ../
    cd ..
    echo ""

    #  build kkrimp code
    echo "build kkrimp"
    git checkout BdG # use BdG branch to have new ALAT tolerance input
    ./install.py --program=kkrimp --compiler=gfortran --parallelization=serial
    cd build/ && make -j4 && cp kkrflex.exe ../
    cd ../..
    echo ""
fi

echo 'To complete the installation make sure the executables (or fakes) are found in the PATH'
echo 'the following command should work:'
echo 'cd jukkr && export PATH="$PWD:$PATH" && cd ..'
