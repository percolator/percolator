# Centralized place for urls and files for all linux builders ...
# please do not change compression type in urls, since decompression is
# hardcoded in the respective buiding scripts

#
# macOS
#
# -- xsd --
mac_os_xsd='xsd-4.0.0+dep'
mac_os_xsd_url='https://www.codesynthesis.com/download/xsd/4.0/'${mac_os_xsd}'.tar.bz2'
# -- xerces --
mac_os_xerces='xerces-c-3.2.5'
mac_os_xerces_url='https://archive.apache.org/dist/xerces/c/3/sources/'${mac_os_xerces}'.tar.gz'

#
# CentOS
#
# -- xsd --
centos_xsd='xsd-4.0.0-1.x86_64'
centos_xsd_url='http://www.codesynthesis.com/download/xsd/4.0/linux-gnu/x86_64/'${centos_xsd}'.rpm'
# -- boost --
centos_boost='boost_1_67_0'
centos_boost_url='https://dl.bintray.com/boostorg/release/1.67.0/source/'${centos_boost}'.tar.bz2'

#
# Fedora
#
# -- xsd --
fedora_xsd='xsd-4.0.0-1.x86_64'
fedora_xsd_url='http://www.codesynthesis.com/download/xsd/4.0/linux-gnu/x86_64/'${fedora_xsd}'.rpm'

# Ubuntu
#
# -- xerces --
ubuntu_xerces='xerces-c-3.2.3'
ubuntu_xerces_url='https://archive.apache.org/dist/xerces/c/3/sources/'${ubuntu_xerces}'.tar.gz'
# -- xsd --
ubuntu_xsd='xsd-3.3.0-x86_64-linux-gnu'
ubuntu_xsd_url='http://www.codesynthesis.com/download/xsd/3.3/linux-gnu/x86_64/'${ubuntu_xsd}'.tar.bz2'
