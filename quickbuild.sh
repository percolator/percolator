#!/bin/bash
# managing input arguments

#--------------------------------------------------------

usage()
{
    cat << EOF
This script detects the operating system and executes the corresponding
builder file.

usage: $0 
                  [[-h]]
                  [[-s source_directory]]
                  [[-b build_directory]]
                  [[-r release_directory]]
                  [[-p ubuntu|centos|fedora|osx]]

If no branch and source_directory is provided, the source
code from which the sourcecode is checked out from will be used.
  -h     prints this help page
EOF
}

skip_gui=""
while getopts “hs:b:r:p:g” OPTION; do
    case $OPTION in
        h)  usage;exit 1;;
        s)  src_dir=$OPTARG;;
        b)  build_dir=$OPTARG;;
        r)  release_dir=$OPTARG;;
        g)  skip_gui="-g";;
        p)  case $OPTARG in
               	tarball)
               	    post="tarball64"
               	    ;;
               	ubuntu)
                    post="ubuntu64"
                    ;;
                fedora)	
                    post="fedora64"
                    ;;
                centos) 
                    post="centos64"
                    ;;
                osx) 
                    post="osx64"
                    ;;
                *)
                    if [[ $OPTARG == *,* ]]; then
                      arr=$(echo $OPTARG | tr "," "\n");
                      multi_platform="1";
                    else
                      echo "Platform $OPTARG is undefined."
                      exit 1
                    fi;;
            esac;;
        \?)  echo "Invalid option: -${OPTARG}" >&2;;
    esac
done

if [[ -z $src_dir ]]; then
  src_dir=$(pwd)
fi

if [[ -z $post ]]; then
  if [ -f /etc/redhat-release ]; then
    linux_issue=$(cat /etc/redhat-release)
  elif [ -f /etc/issue ]; then
    linux_issue=$(cat /etc/issue)
  else
    linux_issue=$(uname -s)
  fi
  
  case "$linux_issue" in 
    *Ubuntu*|*Debian*)
      post="ubuntu64"
      ;;
    *Fedora*)
      post="fedora64"
      ;;
    *CentOS*|*Red*)
      post="centos64"
      ;;
    *Darwin*)
      post="osx64"
      ;;
  esac
fi

if [[ -z $post ]]; then
  echo "Could not determine operating system, please set with the -p argument."
  exit 1
fi

builder="${post}_build.sh";

if [[ -z $build_dir ]]; then
  build_dir=$src_dir/../build/${post}
fi
if [[ -z $release_dir ]]; then
  release_dir=$src_dir/../release/${post}
fi

cd ${src_dir}/admin/builders
./${builder} -b ${build_dir} -r ${release_dir} -s ${src_dir}/../
