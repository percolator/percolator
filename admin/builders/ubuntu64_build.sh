#!/bin/bash
# managing input arguments

while getopts “s:b:r:t:” OPTION; do
    case $OPTION in
        s) src_dir=${OPTARG};;
        t) branch=${OPTARG};;
        r) release_dir=${OPTARG};;
        b) build_dir=${OPTARG};;
        d) build_type=${OPTARG};;
        \?) echo "Invalid option: -${OPTARG}" >&2;;
    esac
done

# get current architecture
ARCH=$(uname -m)

sudo apt-get -y install gawk;
if [[ ! -z `echo -e "$(lsb_release -r)" | gawk '($2>="22.04"){print
    $2}'` ]]; then
    sudo apt-get -y install libcurl4-openssl-dev;
fi

if [[ -z ${build_type} ]]; then
    build_type="Release";
fi
if [[ -z ${build_dir} ]]; then
    build_dir="$(mktemp -d --tmpdir ubuntu_build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
    if [[ -n  ${branch} ]]; then
        sudo apt-get install git;
        src_dir="$(mktemp -d --tmpdir ubuntu_build_XXXX)";
        git clone --branch "$1" https://github.com/percolator/percolator.git "${src_dir}/percolator";
    else
        src_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"/../../../
    fi
fi
if [[ -z ${release_dir} ]]; then
    release_dir=${HOME}/release
    if [ ! -d "${release_dir}" ]; then
        mkdir ${release_dir}
    fi
fi

echo "The Builder $0 is building the Percolator packages with src=${src_dir} an\
d build=${build_dir} for the user"
whoami;

#------------------------------------------------------------------------
#------------------------------------------------------------------------
echo "Checking necessary packages for building percolator...";

# Do not apt-upgrade if this is a travis-ci job
sudo apt-get update;
if [ -z "$TRAVIS" ] && [ -z "$CI" ]; then
    # trap 'echo "EXIT (rc: $?)" && exit 1' ERR
    sudo apt-get upgrade;
    sudo apt-get -y install g++ make cmake rpm fakeroot;
fi


# cd ${src_dir}
# read all urls and file names from a centralized kb file
# source percolator/admin/builders/_urls_and_file_names_.sh

mkdir -p $build_dir
cd ${build_dir}

######percolator########

#-----cmake-----
mkdir -p $build_dir/percolator;
cd $build_dir/percolator;
echo "cmake percolator.....";
(set -x;
    cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=${build_type} -DCMAKE_INSTALL_PREFIX=/usr -DGOOGLE_TEST=1  $src_dir/percolator;
)
#-----make------
echo "make percolator (this will take few minutes).....";
make -j 4;
make test;
make -j 4 package;


echo "Finished buildscript execution";
echo "in build directory ${build_dir}";

cp -v $build_dir/percolator/percolator*.deb ${release_dir};
