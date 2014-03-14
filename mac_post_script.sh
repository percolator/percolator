#!/bin/sh
#
# Author: Jose Fernandez <jc.fernandez.navarro@gmail.com>

if [ "$(id -u)" != "0" ]; then
	echo "Sorry, you are not root."
	exit 1
fi

echo "running $0"

PERCOLATOR_BUNDLE="`echo "$0" | sed -e 's/\/Contents\/MacOS\/Percolator//'`"
PERCOLATOR_RESOURCES="$PERCOLATOR_BUNDLE/src/xml"
PERCOLATOR_EXE="$PERCOLATOR_BUNDLE/percolator"
QVALITY_EXE="$PERCOLATOR_BUNDLE/qvality"
PIN_VERSION_MAJOR="1"
PIN_VERSION_MINOR="3"
POUT_VERSION_MAJOR="1"
POUT_VERSION_MINOR="5"
WRITABLE_DIR="usr/share/percolator/"
PIN_SCHEMA_LOCATION="$WRITABLE_DIR""xml-pin-$PIN_VERSION_MAJOR-$PIN_VERSION_MINOR/"
POUT_SCHEMA_LOCATION="$WRITABLE_DIR""xml-pout-${POUT_VERSION_MAJOR}-${POUT_VERSION_MINOR}/"

sudo cp $PERCOLATOR_EXE ./bin
sudo cp $QVALITY_EXE ./bin
sudo cp $PERCOLATOR_RESOURCES $PIN_SCHEMA_LOCATION
sudo cp $PERCOLATOR_RESOURCES $POUT_SCHEMA_LOCATION



