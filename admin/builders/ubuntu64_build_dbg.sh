#!/bin/bash
# managing input arguments

while getopts “s:b:r:t:” OPTION; do
    case $OPTION in
        s) src_dir=${OPTARG};;
        t) branch=${OPTARG};;
        r) release_dir=${OPTARG};;
        b) build_dir=${OPTARG};;
        \?) echo "Invalid option: -${OPTARG}" >&2;;
    esac
done

sudo apt-get -y install gawk;
if [[ ! -z `echo -e "$(lsb_release -r)" | gawk '($2>="22.04"){print
    $2}'` ]]; then
    sudo apt-get -y install libcurl4-openssl-dev;
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

cd ${src_dir}

# read all urls and file names from a centralized kb file
source percolator/admin/builders/_urls_and_file_names_.sh

mkdir -p $build_dir
cd ${build_dir}

# download and patch xsd
if [[ $(lsb_release -a) == *"14.04"* ]]; then
    if [ ! -d ${ubuntu_xsd} ]; then
        echo "Installing XSD"
        wget --quiet ${ubuntu_xsd_url}
        tar xjf ${ubuntu_xsd}.tar.bz2
        sed -i 's/setg/this->setg/g' ${ubuntu_xsd}/libxsd/xsd/cxx/zc-istream.txx
        sed -i 's/ push_back/ this->push_back/g' ${ubuntu_xsd}/libxsd/xsd/cxx/tree/parsing.txx
        sed -i 's/ push_back/ this->push_back/g' ${ubuntu_xsd}/libxsd/xsd/cxx/tree/stream-extraction.hxx
    fi
else
    sudo apt-get -y install xsdcxx;
fi

# issue with XercesC in Ubuntu 16.04: https://github.com/percolator/percolator/issues/188
if [[ ! -z `echo -e "$(lsb_release -r)" | gawk '($2>="18.04"){print $2}'` ]]; then
    if [[ ! -d ${ubuntu_xerces}/lib ]]; then
        echo "Installing XercesC"
        # download, compile and link xerces
        wget --no-check-certificate --quiet ${ubuntu_xerces_url}
        tar xzf ${ubuntu_xerces}.tar.gz
        cd ${ubuntu_xerces}/
        ./configure --prefix=${build_dir}/${ubuntu_xerces} --disable-netaccessor-curl --disable-transcoder-icu > ../xercesc_config.log 2>&1
        make -j 4 > ../xercesc_make.log 2>&1
        make install > ../xercesc_install.log 2>&1
    fi
else
    sudo apt-get -y install libxerces-c-dev
fi

# end of section to remove
sudo apt-get -y install libboost-dev libboost-filesystem-dev xsdcxx;
sudo apt-get -y install libboost-system-dev libboost-thread-dev libsqlite3-dev libtokyocabinet-dev zlib1g-dev libbz2-dev;

#------------------------------------------------------------------------
mkdir -p $build_dir/percolator-noxml $build_dir/percolator $build_dir/converters;

######percolator########
#-----cmake-----
cd $build_dir/percolator-noxml;
echo "cmake percolator-noxml.....";
(set -x;
    cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=OFF $src_dir/percolator;
)
#-----make------
echo "make percolator (this will take few minutes).....";
make -j 4;
make -j 4 package;

#-----cmake-----
cd $build_dir/percolator;
echo "cmake percolator.....";
(set -x;
    cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/usr -DGOOGLE_TEST=1 -DCMAKE_PREFIX_PATH="${build_dir}/${ubuntu_xerces}/;${build_dir}/${ubuntu_xsd}/" -DXML_SUPPORT=ON $src_dir/percolator;
)
#-----make------
echo "make percolator (this will take few minutes).....";
make -j 4;
make -j 4 package;

#######converters########
cd $build_dir/converters
#-----cmake-----
echo "cmake converters.....";
(set -x;
    cmake -DTARGET_ARCH=amd64 -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH="${build_dir}/${ubuntu_xerces}/;${build_dir}/${ubuntu_xsd}/" -DSERIALIZE="TokyoCabinet" $src_dir/percolator/src/converters;
)
#-----make------
echo "make converters (this will take few minutes).....";

make -j 4;
make -j 4 package;


###########################

echo "Finished buildscript execution";
echo "in build directory ${build_dir}";

cp -v $build_dir/{percolator-noxml,percolator,converters}/*.deb ${release_dir};
